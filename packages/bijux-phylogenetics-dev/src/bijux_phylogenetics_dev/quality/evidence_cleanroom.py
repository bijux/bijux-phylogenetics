"""Run governed evidence reruns inside a detached clean-room worktree."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import shutil
import subprocess  # nosec B404
import sys
import tomllib
from typing import Any
from uuid import uuid4

from ..trusted_process import run_text
from .policies import EXECUTION_SURFACES_POLICY_PATH

DEFAULT_ARTIFACTS_ROOT = Path("artifacts/root/evidence-cleanroom")


@dataclass(frozen=True)
class EvidenceCleanroomReport:
    """Summarize one clean-room evidence rerun and validation pass."""

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


@dataclass(frozen=True)
class SelectedEvidenceCleanroomReport:
    """Summarize the governed set of clean-room evidence reruns."""

    selection_count: int
    reports: tuple[EvidenceCleanroomReport, ...]


def _run(
    command: list[str],
    *,
    cwd: Path,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    return run_text(
        command,
        check=True,
        capture_output=True,
        cwd=cwd,
        env=env,
    )


def _json_payload(output: str) -> dict[str, Any]:
    payload = json.loads(output)
    if not isinstance(payload, dict):
        raise ValueError("expected JSON object output")
    return payload


def _load_selected_cleanroom_runs(repo_root: Path) -> list[tuple[str, list[str]]]:
    with (repo_root / EXECUTION_SURFACES_POLICY_PATH).open("rb") as handle:
        payload = tomllib.load(handle)
    tool = payload.get("tool", {})
    workspace = tool.get("bijux_phylogenetics", {}) if isinstance(tool, dict) else {}
    policy = (
        workspace.get("execution_surfaces", {}) if isinstance(workspace, dict) else {}
    )
    selections = (
        policy.get("cleanroom_selections", []) if isinstance(policy, dict) else []
    )
    results: list[tuple[str, list[str]]] = []
    for selection in selections:
        if not isinstance(selection, dict):
            continue
        study_id = selection.get("study_id")
        evidence_ids = selection.get("evidence_ids", [])
        if not isinstance(study_id, str) or not isinstance(evidence_ids, list):
            continue
        normalized_ids = [
            evidence_id for evidence_id in evidence_ids if isinstance(evidence_id, str)
        ]
        if normalized_ids:
            results.append((study_id, normalized_ids))
    return results


def build_evidence_cleanroom_report(
    repo_root: Path,
    *,
    study_id: str,
    evidence_ids: list[str],
    artifacts_root: Path | None = None,
) -> EvidenceCleanroomReport:
    """Run one evidence selection in a detached worktree and report the outcome."""
    repo_root = repo_root.resolve()
    selected_evidence_ids = tuple(evidence_ids)
    if not selected_evidence_ids:
        raise ValueError("at least one evidence id is required for a clean-room rerun")
    artifacts_root = (artifacts_root or repo_root / DEFAULT_ARTIFACTS_ROOT).resolve()
    artifacts_root.mkdir(parents=True, exist_ok=True)
    cleanroom_root = artifacts_root / f"{study_id}-{uuid4().hex[:8]}-worktree"
    git_bin = shutil.which("git")
    if git_bin is None:
        raise FileNotFoundError("git is required for evidence clean-room reruns")

    _run(
        [
            git_bin,
            "-C",
            str(repo_root),
            "worktree",
            "add",
            "--detach",
            str(cleanroom_root),
            "HEAD",
        ],
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

        sync_command = [
            sys.executable,
            "-m",
            "bijux_phylogenetics_dev.quality.evidence_artifacts",
            "sync",
            "--repo-root",
            ".",
        ]
        _run(sync_command, cwd=cleanroom_root, env=env)

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
            [git_bin, "-C", str(cleanroom_root), "status", "--short"],
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
        run_text(
            [
                git_bin,
                "-C",
                str(repo_root),
                "worktree",
                "remove",
                "--force",
                str(cleanroom_root),
            ],
            check=True,
            cwd=repo_root,
            capture_output=True,
        )
        if cleanroom_root.exists():
            shutil.rmtree(cleanroom_root, ignore_errors=True)


def build_selected_evidence_cleanroom_reports(
    repo_root: Path,
    *,
    artifacts_root: Path | None = None,
) -> SelectedEvidenceCleanroomReport:
    """Run every governed clean-room evidence selection and collect reports."""
    repo_root = repo_root.resolve()
    reports = tuple(
        build_evidence_cleanroom_report(
            repo_root,
            study_id=study_id,
            evidence_ids=evidence_ids,
            artifacts_root=artifacts_root,
        )
        for study_id, evidence_ids in _load_selected_cleanroom_runs(repo_root)
    )
    return SelectedEvidenceCleanroomReport(
        selection_count=len(reports),
        reports=reports,
    )


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the clean-room evidence CLI."""
    parser = argparse.ArgumentParser(
        description="Run an evidence rerun in a detached clean-room worktree."
    )
    parser.add_argument(
        "command", choices=("report", "check", "report-selected", "check-selected")
    )
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--study-id")
    parser.add_argument(
        "--evidence-id", dest="evidence_ids", action="append", default=[]
    )
    parser.add_argument("--artifacts-root", default=str(DEFAULT_ARTIFACTS_ROOT))
    parser.add_argument("--json-out", default="")
    return parser.parse_args()


def main() -> int:
    """Run the clean-room evidence CLI entry point."""
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    artifacts_root = Path(args.artifacts_root)
    if not artifacts_root.is_absolute():
        artifacts_root = repo_root / artifacts_root
    if args.command in {"report-selected", "check-selected"}:
        selected_report = build_selected_evidence_cleanroom_reports(
            repo_root,
            artifacts_root=artifacts_root,
        )
        payload: dict[str, object] = {
            "selection_count": selected_report.selection_count,
            "reports": [asdict(entry) for entry in selected_report.reports],
        }
    else:
        if not args.study_id:
            raise SystemExit(
                "--study-id is required for single-selection clean-room runs"
            )
        single_report = build_evidence_cleanroom_report(
            repo_root,
            study_id=args.study_id,
            evidence_ids=args.evidence_ids,
            artifacts_root=artifacts_root,
        )
        payload = asdict(single_report)
    json_out = args.json_out
    if json_out:
        json_path = Path(json_out)
        if not json_path.is_absolute():
            json_path = repo_root / json_path
    else:
        stem = args.study_id if args.study_id else "selected-evidence"
        json_path = artifacts_root / f"{stem}-cleanroom.json"
    _write_json(json_path, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    if args.command == "check":
        if single_report.validation_issue_count:
            raise SystemExit("evidence clean-room validation failed")
        if single_report.artifact_issue_count:
            raise SystemExit("evidence clean-room artifact validation failed")
        if not single_report.worktree_clean:
            raise SystemExit("evidence clean-room rerun produced repository drift")
    if args.command == "check-selected":
        for entry in selected_report.reports:
            if entry.validation_issue_count:
                raise SystemExit("evidence clean-room validation failed")
            if entry.artifact_issue_count:
                raise SystemExit("evidence clean-room artifact validation failed")
            if not entry.worktree_clean:
                raise SystemExit("evidence clean-room rerun produced repository drift")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
