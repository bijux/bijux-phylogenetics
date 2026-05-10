"""Run governed evidence reruns inside a detached clean-room worktree."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from uuid import uuid4

DEFAULT_ARTIFACTS_ROOT = Path("artifacts/root/evidence-cleanroom")


@dataclass(frozen=True)
class EvidenceCleanroomReport:
    study_id: str
    selected_evidence_ids: tuple[str, ...]
    cleanroom_path: str
    rerun_command: tuple[str, ...]
    validate_command: tuple[str, ...]
    artifact_command: tuple[str, ...]
    worktree_clean: bool
    updated_path_count: int
    validation_issue_count: int
    artifact_issue_count: int
    status_entries: tuple[str, ...]


def _run(
    command: list[str],
    *,
    cwd: Path,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        check=True,
        cwd=str(cwd),
        env=env,
        text=True,
        capture_output=True,
    )


def _json_payload(output: str) -> dict[str, object]:
    payload = json.loads(output)
    if not isinstance(payload, dict):
        raise ValueError("expected JSON object output")
    return payload


def build_evidence_cleanroom_report(
    repo_root: Path,
    *,
    study_id: str,
    evidence_ids: list[str],
    artifacts_root: Path | None = None,
) -> EvidenceCleanroomReport:
    repo_root = repo_root.resolve()
    selected_evidence_ids = tuple(evidence_ids)
    if not selected_evidence_ids:
        raise ValueError("at least one evidence id is required for a clean-room rerun")
    artifacts_root = (artifacts_root or repo_root / DEFAULT_ARTIFACTS_ROOT).resolve()
    artifacts_root.mkdir(parents=True, exist_ok=True)
    cleanroom_root = artifacts_root / f"{study_id}-{uuid4().hex[:8]}-worktree"

    _run(
        ["git", "-C", str(repo_root), "worktree", "add", "--detach", str(cleanroom_root), "HEAD"],
        cwd=repo_root,
    )
    try:
        pythonpath_parts = [
            str(cleanroom_root / "packages/bijux-phylogenetics/src"),
            str(cleanroom_root / "packages/bijux-phylogenetics-dev/src"),
        ]
        existing_pythonpath = os.environ.get("PYTHONPATH")
        if existing_pythonpath:
            pythonpath_parts.append(existing_pythonpath)
        env = dict(os.environ)
        env["PYTHONPATH"] = os.pathsep.join(pythonpath_parts)

        rerun_command = [
            sys.executable,
            "-m",
            "bijux_phylogenetics",
            "evidence",
            "book",
            "rerun",
            study_id,
        ]
        rerun_command.extend(selected_evidence_ids)
        rerun_command.append("--json")
        rerun_payload = _json_payload(
            _run(rerun_command, cwd=cleanroom_root, env=env).stdout
        )

        validate_command = [
            sys.executable,
            "-m",
            "bijux_phylogenetics",
            "evidence",
            "book",
            "validate",
            "--json",
        ]
        validate_payload = _json_payload(
            _run(validate_command, cwd=cleanroom_root, env=env).stdout
        )

        artifact_command = [
            sys.executable,
            "-m",
            "bijux_phylogenetics_dev.quality.evidence_artifacts",
            "check",
            "--repo-root",
            ".",
        ]
        artifact_payload = _json_payload(
            _run(artifact_command, cwd=cleanroom_root, env=env).stdout
        )

        status_output = _run(
            ["git", "-C", str(cleanroom_root), "status", "--short"],
            cwd=cleanroom_root,
        ).stdout
        status_entries = tuple(
            line.strip() for line in status_output.splitlines() if line.strip()
        )
        validation_metrics = validate_payload.get("metrics", {})
        validation_issue_count = (
            int(validation_metrics.get("issue_count", 0))
            if isinstance(validation_metrics, dict)
            else 0
        )
        artifact_issue_count = int(artifact_payload.get("issue_count", 0))
        rerun_metrics = rerun_payload.get("metrics", {})
        updated_path_count = (
            int(rerun_metrics.get("updated_path_count", 0))
            if isinstance(rerun_metrics, dict)
            else 0
        )
        return EvidenceCleanroomReport(
            study_id=study_id,
            selected_evidence_ids=selected_evidence_ids,
            cleanroom_path=cleanroom_root.as_posix(),
            rerun_command=tuple(rerun_command),
            validate_command=tuple(validate_command),
            artifact_command=tuple(artifact_command),
            worktree_clean=not status_entries,
            updated_path_count=updated_path_count,
            validation_issue_count=validation_issue_count,
            artifact_issue_count=artifact_issue_count,
            status_entries=status_entries,
        )
    finally:
        subprocess.run(
            ["git", "-C", str(repo_root), "worktree", "remove", "--force", str(cleanroom_root)],
            check=True,
            cwd=str(repo_root),
            text=True,
            capture_output=True,
        )
        if cleanroom_root.exists():
            shutil.rmtree(cleanroom_root, ignore_errors=True)


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run an evidence rerun in a detached clean-room worktree."
    )
    parser.add_argument("command", choices=("report", "check"))
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--study-id", required=True)
    parser.add_argument("--evidence-id", dest="evidence_ids", action="append", default=[])
    parser.add_argument("--artifacts-root", default=str(DEFAULT_ARTIFACTS_ROOT))
    parser.add_argument("--json-out", default="")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    artifacts_root = Path(args.artifacts_root)
    if not artifacts_root.is_absolute():
        artifacts_root = repo_root / artifacts_root
    report = build_evidence_cleanroom_report(
        repo_root,
        study_id=args.study_id,
        evidence_ids=args.evidence_ids,
        artifacts_root=artifacts_root,
    )
    payload = asdict(report)
    json_out = args.json_out
    if json_out:
        json_path = Path(json_out)
        if not json_path.is_absolute():
            json_path = repo_root / json_path
    else:
        json_path = artifacts_root / f"{args.study_id}-cleanroom.json"
    _write_json(json_path, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    if args.command == "check":
        if report.validation_issue_count:
            raise SystemExit("evidence clean-room validation failed")
        if report.artifact_issue_count:
            raise SystemExit("evidence clean-room artifact validation failed")
        if not report.worktree_clean:
            raise SystemExit("evidence clean-room rerun produced repository drift")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
