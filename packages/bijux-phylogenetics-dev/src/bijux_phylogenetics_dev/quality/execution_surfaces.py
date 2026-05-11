"""Audit make, tox, and workflow ownership for repository execution surfaces."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import tomllib
from typing import Any

import yaml

from .policies import EXECUTION_SURFACES_POLICY_PATH

DEFAULT_JSON_OUT = Path("artifacts/root/execution-surfaces.json")
ROOT_MAKEFILE = Path("makes/root.mk")
TOX_FILE = Path("tox.ini")
WORKFLOWS_DIR = Path(".github/workflows")


@dataclass(frozen=True)
class ExecutionSurfaceIssue:
    """Describe one governed execution-surface contract drift."""

    code: str
    path: str
    message: str


def _load_toml(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def _as_dict(value: object) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _as_str_list(value: object) -> list[str]:
    if not isinstance(value, list):
        return []
    return [entry for entry in value if isinstance(entry, str)]


def _load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"expected mapping at {path}")
    normalized: dict[str, Any] = {}
    for key, value in payload.items():
        normalized["on" if key is True else key] = value
    return normalized


def _root_make_targets(text: str) -> dict[str, str]:
    bodies: dict[str, list[str]] = {}
    current_target: str | None = None
    for raw_line in text.splitlines():
        if raw_line.startswith("\t") and current_target is not None:
            bodies[current_target].append(raw_line.strip())
            continue
        stripped = raw_line.strip()
        if not stripped or stripped.startswith((".PHONY:", "#")):
            current_target = None
            continue
        if ":" not in stripped:
            current_target = None
            continue
        current_target = stripped.split(":", 1)[0].strip()
        bodies.setdefault(current_target, [])
    return {target: "\n".join(lines) for target, lines in bodies.items()}


def _parse_envlist(text: str) -> set[str]:
    in_envlist = False
    envs: set[str] = set()
    for line in text.splitlines():
        if line.startswith("envlist"):
            in_envlist = True
            continue
        if not in_envlist:
            continue
        if not line.startswith((" ", "\t")):
            break
        stripped = line.strip()
        if stripped:
            envs.add(stripped)
    return envs


def _tox_commands(text: str) -> dict[str, str]:
    commands_by_env: dict[str, list[str]] = {}
    current_env: str | None = None
    collecting_commands = False
    for raw_line in text.splitlines():
        stripped = raw_line.strip()
        if stripped.startswith("[testenv:") and stripped.endswith("]"):
            current_env = stripped[len("[testenv:") : -1]
            commands_by_env.setdefault(current_env, [])
            collecting_commands = False
            continue
        if stripped.startswith("[") and stripped.endswith("]"):
            current_env = None
            collecting_commands = False
            continue
        if current_env is None:
            continue
        if stripped == "commands =":
            collecting_commands = True
            continue
        if collecting_commands:
            if raw_line.startswith((" ", "\t")) and stripped:
                commands_by_env[current_env].append(stripped)
                continue
            collecting_commands = False
    return {name: "\n".join(lines) for name, lines in commands_by_env.items()}


def _workflow_jobs(path: Path) -> set[str]:
    payload = _load_yaml(path)
    jobs = payload.get("jobs", {})
    return set(jobs) if isinstance(jobs, dict) else set()


def _workflow_run_commands(path: Path) -> str:
    payload = _load_yaml(path)
    jobs = payload.get("jobs", {})
    if not isinstance(jobs, dict):
        return ""
    commands: list[str] = []
    for job in jobs.values():
        if not isinstance(job, dict):
            continue
        steps = job.get("steps", [])
        if not isinstance(steps, list):
            continue
        for step in steps:
            if isinstance(step, dict) and isinstance(step.get("run"), str):
                commands.append(step["run"])
    return "\n".join(commands)


def build_execution_surfaces_report(repo_root: Path) -> dict[str, Any]:
    """Build the execution-surface ownership and workflow report."""
    repo_root = repo_root.resolve()
    payload = _load_toml(repo_root / EXECUTION_SURFACES_POLICY_PATH)
    tool = _as_dict(payload.get("tool"))
    workspace = _as_dict(tool.get("bijux_phylogenetics"))
    policy = _as_dict(workspace.get("execution_surfaces"))
    required_root_make_targets = _as_str_list(policy.get("required_root_make_targets"))
    required_tox_envs = _as_str_list(policy.get("required_tox_envs"))
    tox_commands = {
        env_name: _as_str_list(commands)
        for env_name, commands in _as_dict(policy.get("tox_commands")).items()
    }
    workflow_jobs = {
        workflow_name: _as_str_list(job_names)
        for workflow_name, job_names in _as_dict(policy.get("workflow_jobs")).items()
    }
    cleanroom_selections = [
        entry
        for entry in policy.get("cleanroom_selections", [])
        if isinstance(entry, dict)
    ]

    issues: list[ExecutionSurfaceIssue] = []

    root_make_text = (repo_root / ROOT_MAKEFILE).read_text(encoding="utf-8")
    root_targets = _root_make_targets(root_make_text)
    missing_root_targets = sorted(
        target
        for target in required_root_make_targets
        if target[:-1] not in root_targets
    )
    for target in missing_root_targets:
        issues.append(
            ExecutionSurfaceIssue(
                code="missing-governed-root-target",
                path=ROOT_MAKEFILE.as_posix(),
                message=f"missing governed root make target {target}",
            )
        )

    tox_text = (repo_root / TOX_FILE).read_text(encoding="utf-8")
    envlist = _parse_envlist(tox_text)
    commands_by_env = _tox_commands(tox_text)
    missing_tox_envs = sorted(env for env in required_tox_envs if env not in envlist)
    for env_name in missing_tox_envs:
        issues.append(
            ExecutionSurfaceIssue(
                code="missing-governed-tox-env",
                path=TOX_FILE.as_posix(),
                message=f"missing governed tox environment {env_name}",
            )
        )
    tox_env_reports: list[dict[str, Any]] = []
    for env_name in required_tox_envs:
        actual_commands = commands_by_env.get(env_name, "")
        expected_commands = tox_commands.get(env_name, [])
        missing_commands = sorted(
            command for command in expected_commands if command not in actual_commands
        )
        if missing_commands:
            issues.append(
                ExecutionSurfaceIssue(
                    code="tox-command-drift",
                    path=TOX_FILE.as_posix(),
                    message=(
                        f"tox environment {env_name} is missing governed commands: "
                        + ", ".join(missing_commands)
                    ),
                )
            )
        tox_env_reports.append(
            {
                "env_name": env_name,
                "present": env_name in envlist,
                "expected_commands": expected_commands,
                "missing_commands": missing_commands,
            }
        )

    workflow_reports: list[dict[str, Any]] = []
    for workflow_name, expected_jobs in sorted(workflow_jobs.items()):
        workflow_path = repo_root / WORKFLOWS_DIR / workflow_name
        if not workflow_path.is_file():
            issues.append(
                ExecutionSurfaceIssue(
                    code="missing-governed-workflow",
                    path=(WORKFLOWS_DIR / workflow_name).as_posix(),
                    message="required execution workflow is missing",
                )
            )
            workflow_reports.append(
                {
                    "workflow_name": workflow_name,
                    "exists": False,
                    "missing_jobs": expected_jobs,
                    "has_cleanroom_commands": False,
                }
            )
            continue
        jobs = _workflow_jobs(workflow_path)
        missing_jobs = sorted(
            job_name for job_name in expected_jobs if job_name not in jobs
        )
        if missing_jobs:
            issues.append(
                ExecutionSurfaceIssue(
                    code="workflow-job-drift",
                    path=(WORKFLOWS_DIR / workflow_name).as_posix(),
                    message=(
                        f"workflow is missing governed jobs: {', '.join(missing_jobs)}"
                    ),
                )
            )
        run_commands = _workflow_run_commands(workflow_path)
        has_cleanroom_commands = True
        if workflow_name == "evidence-governance.yml":
            for selection in cleanroom_selections:
                study_id = selection.get("study_id")
                evidence_ids = selection.get("evidence_ids", [])
                if not isinstance(study_id, str) or not isinstance(evidence_ids, list):
                    continue
                ids = " ".join(
                    evidence_id
                    for evidence_id in evidence_ids
                    if isinstance(evidence_id, str)
                )
                expected_snippet = f"make rerun-evidence-cleanroom EVIDENCE_STUDY_ID={study_id} EVIDENCE_IDS={ids}"
                if (
                    expected_snippet not in run_commands
                    and "make rerun-governed-evidence-cleanroom" not in run_commands
                ):
                    has_cleanroom_commands = False
                    issues.append(
                        ExecutionSurfaceIssue(
                            code="workflow-cleanroom-selection-drift",
                            path=(WORKFLOWS_DIR / workflow_name).as_posix(),
                            message=(
                                "evidence workflow is missing governed clean-room rerun "
                                f"selection for {study_id}"
                            ),
                        )
                    )
        workflow_reports.append(
            {
                "workflow_name": workflow_name,
                "exists": True,
                "missing_jobs": missing_jobs,
                "has_cleanroom_commands": has_cleanroom_commands,
            }
        )

    return {
        "schema_version": 1,
        "required_root_make_targets": required_root_make_targets,
        "missing_root_make_targets": missing_root_targets,
        "tox_envs": tox_env_reports,
        "workflows": workflow_reports,
        "cleanroom_selections": cleanroom_selections,
        "issue_count": len(issues),
        "issues": [asdict(issue) for issue in issues],
    }


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def check_execution_surfaces(
    repo_root: Path,
    *,
    json_out: Path | None = None,
) -> dict[str, Any]:
    """Raise when any governed execution-surface requirement is missing."""
    payload = build_execution_surfaces_report(repo_root)
    if json_out is not None:
        _write_json(json_out, payload)
    if payload["issue_count"]:
        raise SystemExit("execution surface audit failed")
    return payload


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the execution-surface audit."""
    parser = argparse.ArgumentParser(
        description="Audit repository make, tox, and workflow execution surfaces."
    )
    parser.add_argument("command", choices=("report", "check"))
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    return parser.parse_args()


def main() -> int:
    """Run the execution-surface CLI entry point."""
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    json_out = Path(args.json_out)
    if not json_out.is_absolute():
        json_out = repo_root / json_out
    if args.command == "check":
        payload = check_execution_surfaces(repo_root, json_out=json_out)
    else:
        payload = build_execution_surfaces_report(repo_root)
        _write_json(json_out, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
