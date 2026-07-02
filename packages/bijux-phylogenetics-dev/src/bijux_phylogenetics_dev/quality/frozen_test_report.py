"""Compare full-suite test results against a frozen repository reference."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
from pathlib import Path
from typing import Any, Literal
from xml.etree import ElementTree

JsonObject = dict[str, Any]
Outcome = Literal["passed", "failed", "error", "skipped"]
OutputChangeKind = Literal["added", "modified", "removed"]
DriftKind = Literal["test", "output"]

_EXCLUDED_OUTPUT_SEGMENTS = frozenset(
    {
        "__pycache__",
        ".pytest_cache",
        "benchmarks",
        "htmlcov",
        "hypothesis",
        "npm_cache",
        "pycache",
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


def _normalized_nodeid(testcase: ElementTree.Element) -> str:
    file_name = testcase.get("file", "").strip()
    test_name = testcase.get("name", "").strip()
    class_name = testcase.get("classname", "").strip()
    if file_name and test_name:
        return f"{file_name}::{test_name}"
    if class_name and test_name:
        return f"{class_name}::{test_name}"
    if test_name:
        return test_name
    raise ValueError("junit testcase is missing name metadata")


def _test_message(testcase: ElementTree.Element) -> str:
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


def _test_outcome(testcase: ElementTree.Element) -> Outcome:
    if testcase.find("failure") is not None:
        return "failed"
    if testcase.find("error") is not None:
        return "error"
    if testcase.find("skipped") is not None:
        return "skipped"
    return "passed"


def parse_junit_report(junit_xml: Path) -> dict[str, SuiteTestRecord]:
    """Parse a pytest JUnit XML report into normalized test records."""
    root = ElementTree.fromstring(junit_xml.read_text(encoding="utf-8"))
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
    if relative_path.name in _EXCLUDED_OUTPUT_FILENAMES:
        return False
    return True


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
    for path in sorted(candidate for candidate in artifacts_root.rglob("*") if candidate.is_file()):
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
                current_outcome=None if current_record is None else current_record.outcome,
                frozen_outcome=None if frozen_record is None else frozen_record.outcome,
                current_message="" if current_record is None else current_record.message,
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
                current_sha256=None if current_record is None else current_record.sha256,
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
