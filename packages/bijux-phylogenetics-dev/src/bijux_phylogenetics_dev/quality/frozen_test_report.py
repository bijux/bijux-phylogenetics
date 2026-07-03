"""Compare full-suite test results against a frozen repository reference."""

from __future__ import annotations

import argparse
from collections.abc import Iterable
from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path
import shutil
from typing import Any, Literal, Protocol, cast

from defusedxml import ElementTree

from ..trusted_process import run_text

JsonObject = dict[str, Any]
Outcome = Literal["passed", "failed", "error", "skipped"]
OutputChangeKind = Literal["added", "modified", "removed"]
DriftKind = Literal["test", "output"]
DEFAULT_CURRENT_REF = "HEAD"


class XmlElement(Protocol):
    """Minimal XML element contract needed for JUnit normalization."""

    text: str | None

    def find(self, path: str) -> XmlElement | None:
        """Return the first matching descendant element, if one exists."""

    def get(self, key: str, default: str | None = None) -> str | None:
        """Return one attribute value, or the provided default when absent."""

    def iterfind(self, path: str) -> Iterable[XmlElement]:
        """Iterate over every descendant element that matches the path."""


_EXCLUDED_OUTPUT_SEGMENTS = frozenset(
    {
        "__pycache__",
        ".pytest_cache",
        "benchmarks",
        "htmlcov",
        "hypothesis",
        "npm_cache",
        "pycache",
        "runner",
        "tmp",
        "uv_cache",
        "venv",
        "xdg_cache",
    }
)
_EXCLUDED_OUTPUT_FILENAMES = frozenset({".coverage", "coverage.xml", "junit.xml"})


@dataclass(frozen=True, slots=True)
class SuiteTestRecord:
    """Capture one normalized JUnit test record."""

    nodeid: str
    outcome: Outcome
    message: str


@dataclass(frozen=True, slots=True)
class OutputFileRecord:
    """Fingerprint one governed output file."""

    relative_path: str
    sha256: str
    size_bytes: int


@dataclass(frozen=True, slots=True)
class TestDifference:
    """Describe one changed test result between current and frozen runs."""

    nodeid: str
    current_outcome: Outcome | None
    frozen_outcome: Outcome | None
    current_message: str
    frozen_message: str


@dataclass(frozen=True, slots=True)
class OutputDifference:
    """Describe one changed governed output between current and frozen runs."""

    relative_path: str
    change_kind: OutputChangeKind
    current_sha256: str | None
    frozen_sha256: str | None


@dataclass(frozen=True, slots=True)
class ScientificDriftRecord:
    """Describe one observable scientific drift surfaced by frozen validation."""

    kind: DriftKind
    identifier: str
    summary: str


@dataclass(frozen=True, slots=True)
class ReferenceRecord:
    """Describe one Git reference used during frozen validation."""

    input_ref: str
    resolved_ref: str


@dataclass(frozen=True, slots=True)
class SuiteRunRecord:
    """Capture one isolated full-suite run used in frozen validation."""

    repo_root: str
    artifacts_root: str
    junit_xml: str
    command: tuple[str, ...]
    exit_code: int
    stdout_log: str
    stderr_log: str
    test_count: int
    output_count: int
    full_suite_worker_count: int | None


def _normalized_nodeid(testcase: XmlElement) -> str:
    file_name = _attribute_value(testcase, "file")
    test_name = _attribute_value(testcase, "name")
    class_name = _attribute_value(testcase, "classname")
    if file_name and test_name:
        return f"{file_name}::{test_name}"
    if class_name and test_name:
        return f"{class_name}::{test_name}"
    if test_name:
        return test_name
    raise ValueError("junit testcase is missing name metadata")


def _attribute_value(testcase: XmlElement, attribute_name: str) -> str:
    value = testcase.get(attribute_name, "")
    return "" if value is None else value.strip()


def _test_message(testcase: XmlElement) -> str:
    for outcome_tag in ("failure", "error", "skipped"):
        child = testcase.find(outcome_tag)
        if child is None:
            continue
        message = (child.get("message", "") or "").strip()
        body = (child.text or "").strip()
        if message and body:
            return f"{message}\n{body}"
        return message or body
    return ""


def _test_outcome(testcase: XmlElement) -> Outcome:
    if testcase.find("failure") is not None:
        return "failed"
    if testcase.find("error") is not None:
        return "error"
    if testcase.find("skipped") is not None:
        return "skipped"
    return "passed"


def parse_junit_report(junit_xml: Path) -> dict[str, SuiteTestRecord]:
    """Parse a pytest JUnit XML report into normalized test records."""
    root = cast(
        XmlElement, ElementTree.fromstring(junit_xml.read_text(encoding="utf-8"))
    )
    records: dict[str, SuiteTestRecord] = {}
    for testcase in root.iterfind(".//testcase"):
        nodeid = _normalized_nodeid(testcase)
        records[nodeid] = SuiteTestRecord(
            nodeid=nodeid,
            outcome=_test_outcome(testcase),
            message=_test_message(testcase),
        )
    return records


def is_governed_output_file(relative_path: Path) -> bool:
    """Return whether a file should participate in frozen output comparison."""
    if any(segment in _EXCLUDED_OUTPUT_SEGMENTS for segment in relative_path.parts):
        return False
    return relative_path.name not in _EXCLUDED_OUTPUT_FILENAMES


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_output_manifest(artifacts_root: Path) -> dict[str, OutputFileRecord]:
    """Build a manifest for governed output files under one artifacts root."""
    manifest: dict[str, OutputFileRecord] = {}
    if not artifacts_root.exists():
        return manifest
    for path in sorted(
        candidate for candidate in artifacts_root.rglob("*") if candidate.is_file()
    ):
        relative_path = path.relative_to(artifacts_root)
        if not is_governed_output_file(relative_path):
            continue
        manifest[relative_path.as_posix()] = OutputFileRecord(
            relative_path=relative_path.as_posix(),
            sha256=_sha256(path),
            size_bytes=path.stat().st_size,
        )
    return manifest


def diff_test_records(
    current: dict[str, SuiteTestRecord],
    frozen: dict[str, SuiteTestRecord],
) -> list[TestDifference]:
    """Compute test-result drift between current and frozen runs."""
    differences: list[TestDifference] = []
    for nodeid in sorted(set(current) | set(frozen)):
        current_record = current.get(nodeid)
        frozen_record = frozen.get(nodeid)
        if current_record == frozen_record:
            continue
        differences.append(
            TestDifference(
                nodeid=nodeid,
                current_outcome=None
                if current_record is None
                else current_record.outcome,
                frozen_outcome=None if frozen_record is None else frozen_record.outcome,
                current_message=""
                if current_record is None
                else current_record.message,
                frozen_message="" if frozen_record is None else frozen_record.message,
            )
        )
    return differences


def diff_output_manifests(
    current: dict[str, OutputFileRecord],
    frozen: dict[str, OutputFileRecord],
) -> list[OutputDifference]:
    """Compute governed-output drift between current and frozen runs."""
    differences: list[OutputDifference] = []
    for relative_path in sorted(set(current) | set(frozen)):
        current_record = current.get(relative_path)
        frozen_record = frozen.get(relative_path)
        if current_record == frozen_record:
            continue
        if current_record is None:
            change_kind: OutputChangeKind = "removed"
        elif frozen_record is None:
            change_kind = "added"
        else:
            change_kind = "modified"
        differences.append(
            OutputDifference(
                relative_path=relative_path,
                change_kind=change_kind,
                current_sha256=None
                if current_record is None
                else current_record.sha256,
                frozen_sha256=None if frozen_record is None else frozen_record.sha256,
            )
        )
    return differences


def build_scientific_drift_records(
    test_differences: list[TestDifference],
    output_differences: list[OutputDifference],
) -> list[ScientificDriftRecord]:
    """Summarize observable drift from changed tests and governed outputs."""
    drift_records = [
        ScientificDriftRecord(
            kind="test",
            identifier=difference.nodeid,
            summary=(
                "test outcome changed from "
                f"{difference.frozen_outcome or 'missing'} to "
                f"{difference.current_outcome or 'missing'}"
            ),
        )
        for difference in test_differences
    ]
    drift_records.extend(
        ScientificDriftRecord(
            kind="output",
            identifier=difference.relative_path,
            summary=f"governed output was {difference.change_kind}",
        )
        for difference in output_differences
    )
    return drift_records


def serialize_test_records(records: dict[str, SuiteTestRecord]) -> list[JsonObject]:
    """Convert normalized test records into JSON-friendly payloads."""
    return [asdict(record) for record in records.values()]


def serialize_output_manifest(
    manifest: dict[str, OutputFileRecord],
) -> list[JsonObject]:
    """Convert governed output manifest entries into JSON-friendly payloads."""
    return [asdict(record) for record in manifest.values()]


def serialize_test_differences(differences: list[TestDifference]) -> list[JsonObject]:
    """Convert test drift records into JSON-friendly payloads."""
    return [asdict(difference) for difference in differences]


def serialize_output_differences(
    differences: list[OutputDifference],
) -> list[JsonObject]:
    """Convert governed output drift records into JSON-friendly payloads."""
    return [asdict(difference) for difference in differences]


def serialize_scientific_drift(
    drift_records: list[ScientificDriftRecord],
) -> list[JsonObject]:
    """Convert scientific drift records into JSON-friendly payloads."""
    return [asdict(record) for record in drift_records]


def _git_executable() -> str:
    resolved = shutil.which("git")
    if resolved is None:
        raise SystemExit("git executable not found")
    return resolved


def _make_executable() -> str:
    resolved = shutil.which("make")
    if resolved is None:
        raise SystemExit("make executable not found")
    return resolved


def _write_json(path: Path, payload: JsonObject) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _resolve_reference(
    repo_root: Path,
    ref: str,
    *,
    git_executable: str,
) -> ReferenceRecord:
    completed = run_text(
        [git_executable, "-C", str(repo_root), "rev-parse", ref],
        check=True,
        capture_output=True,
    )
    return ReferenceRecord(input_ref=ref, resolved_ref=completed.stdout.strip())


def _prepare_frozen_worktree(
    repo_root: Path,
    worktree_root: Path,
    *,
    frozen_ref: ReferenceRecord,
    git_executable: str,
) -> Path:
    worktree_root.parent.mkdir(parents=True, exist_ok=True)
    run_text(
        [git_executable, "-C", str(repo_root), "worktree", "prune"],
        check=False,
        capture_output=True,
    )
    run_text(
        [
            git_executable,
            "-C",
            str(repo_root),
            "worktree",
            "remove",
            "--force",
            str(worktree_root),
        ],
        check=False,
        capture_output=True,
    )
    shutil.rmtree(worktree_root, ignore_errors=True)
    run_text(
        [
            git_executable,
            "-C",
            str(repo_root),
            "worktree",
            "add",
            "--detach",
            str(worktree_root),
            frozen_ref.resolved_ref,
        ],
        check=True,
        capture_output=True,
    )
    return worktree_root


def _package_make_command(
    repo_root: Path,
    *,
    package: str,
    artifacts_root: Path,
    make_executable: str,
    target: str = "test-all",
    full_suite_worker_count: int | None = None,
) -> list[str]:
    package_root = repo_root / "packages" / package
    profile_path = os.path.relpath(
        repo_root / "makes" / "packages" / f"{package}.mk",
        start=package_root,
    )
    venv_root = artifacts_root / "venv"
    command = [
        make_executable,
        "-C",
        str(package_root),
        "-f",
        profile_path,
        target,
        f"PROJECT_SLUG={package}",
        f"MONOREPO_ROOT={repo_root}",
        f"ROOT_MAKE_DIR={repo_root / 'makes'}",
        f"CONFIG_DIR={repo_root / 'configs'}",
        f"PROJECT_DIR={package_root}",
        f"PROJECT_ARTIFACTS_DIR={artifacts_root}",
        f"API_DIR={repo_root / 'apis' / package}",
        f"VENV={venv_root}",
        f"VENV_PYTHON={venv_root / 'bin' / 'python'}",
        f"PYTHON={venv_root / 'bin' / 'python'}",
        f"ACT={venv_root / 'bin'}",
    ]
    if full_suite_worker_count is not None:
        command.append(f"FULL_SUITE_PYTEST_WORKERS={full_suite_worker_count}")
    return command


def _run_full_suite(
    repo_root: Path,
    *,
    package: str,
    artifacts_root: Path,
    make_executable: str,
    full_suite_worker_count: int | None,
) -> tuple[SuiteRunRecord, dict[str, SuiteTestRecord], dict[str, OutputFileRecord]]:
    shutil.rmtree(artifacts_root, ignore_errors=True)
    artifacts_root.mkdir(parents=True, exist_ok=True)
    command = _package_make_command(
        repo_root,
        package=package,
        artifacts_root=artifacts_root,
        make_executable=make_executable,
        full_suite_worker_count=full_suite_worker_count,
    )
    completed = run_text(
        command,
        check=False,
        capture_output=True,
        cwd=repo_root,
    )
    stdout_log = artifacts_root / "runner" / "stdout.log"
    stderr_log = artifacts_root / "runner" / "stderr.log"
    stdout_log.parent.mkdir(parents=True, exist_ok=True)
    stdout_log.write_text(completed.stdout, encoding="utf-8")
    stderr_log.write_text(completed.stderr, encoding="utf-8")
    junit_xml = artifacts_root / "test" / "junit.xml"
    test_records = parse_junit_report(junit_xml) if junit_xml.exists() else {}
    output_manifest = build_output_manifest(artifacts_root)
    return (
        SuiteRunRecord(
            repo_root=repo_root.as_posix(),
            artifacts_root=artifacts_root.as_posix(),
            junit_xml=junit_xml.as_posix(),
            command=tuple(command),
            exit_code=completed.returncode,
            stdout_log=stdout_log.as_posix(),
            stderr_log=stderr_log.as_posix(),
            test_count=len(test_records),
            output_count=len(output_manifest),
            full_suite_worker_count=full_suite_worker_count,
        ),
        test_records,
        output_manifest,
    )


def build_frozen_test_report(
    repo_root: Path,
    *,
    package: str,
    frozen_ref: str,
    artifacts_root: Path,
    full_suite_worker_count: int | None = None,
    git_executable: str | None = None,
    make_executable: str | None = None,
) -> JsonObject:
    """Run current and frozen full suites and compare their observed outputs."""
    repo_root = repo_root.resolve()
    artifacts_root = artifacts_root.resolve()
    git_bin = git_executable or _git_executable()
    make_bin = make_executable or _make_executable()
    current_reference = _resolve_reference(
        repo_root, DEFAULT_CURRENT_REF, git_executable=git_bin
    )
    frozen_reference = _resolve_reference(repo_root, frozen_ref, git_executable=git_bin)
    frozen_worktree = _prepare_frozen_worktree(
        repo_root,
        artifacts_root / "worktrees" / f"frozen-{frozen_reference.resolved_ref[:12]}",
        frozen_ref=frozen_reference,
        git_executable=git_bin,
    )
    current_run, current_tests, current_outputs = _run_full_suite(
        repo_root,
        package=package,
        artifacts_root=artifacts_root / "current",
        make_executable=make_bin,
        full_suite_worker_count=full_suite_worker_count,
    )
    frozen_run, frozen_tests, frozen_outputs = _run_full_suite(
        frozen_worktree,
        package=package,
        artifacts_root=artifacts_root / "frozen",
        make_executable=make_bin,
        full_suite_worker_count=full_suite_worker_count,
    )
    changed_tests = diff_test_records(current_tests, frozen_tests)
    changed_outputs = diff_output_manifests(current_outputs, frozen_outputs)
    scientific_drift = build_scientific_drift_records(changed_tests, changed_outputs)
    return {
        "schema_version": 1,
        "package": package,
        "current_ref": asdict(current_reference),
        "frozen_ref": asdict(frozen_reference),
        "current_run": asdict(current_run),
        "frozen_run": asdict(frozen_run),
        "changed_test_count": len(changed_tests),
        "changed_tests": serialize_test_differences(changed_tests),
        "changed_output_count": len(changed_outputs),
        "changed_outputs": serialize_output_differences(changed_outputs),
        "scientific_drift_count": len(scientific_drift),
        "scientific_drift": serialize_scientific_drift(scientific_drift),
        "current_run_succeeded": current_run.exit_code == 0,
        "frozen_run_succeeded": frozen_run.exit_code == 0,
    }


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for frozen full-suite comparison."""
    parser = argparse.ArgumentParser(
        description="Compare package full-suite results against a chosen frozen ref."
    )
    parser.add_argument("command", choices=("report", "check"))
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--package", required=True)
    parser.add_argument("--frozen-ref", required=True)
    parser.add_argument("--artifacts-root")
    parser.add_argument("--json-out")
    parser.add_argument("--full-suite-workers", type=int)
    return parser.parse_args()


def main() -> int:
    """Run the frozen full-suite comparison entry point."""
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    artifacts_root = (
        Path(args.artifacts_root).resolve()
        if args.artifacts_root
        else (repo_root / "artifacts" / args.package / "test-frozen")
    )
    payload = build_frozen_test_report(
        repo_root,
        package=args.package,
        frozen_ref=args.frozen_ref,
        artifacts_root=artifacts_root,
        full_suite_worker_count=args.full_suite_workers,
    )
    json_out = (
        Path(args.json_out).resolve()
        if args.json_out
        else artifacts_root / "report.json"
    )
    _write_json(json_out, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    if args.command == "check" and (
        not payload["current_run_succeeded"]
        or not payload["frozen_run_succeeded"]
        or payload["scientific_drift_count"]
    ):
        raise SystemExit("frozen full-suite validation detected drift")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
