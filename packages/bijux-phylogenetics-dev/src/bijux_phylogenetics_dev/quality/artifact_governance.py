"""Audit repository execution surfaces for governed artifact output paths."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any

DEFAULT_JSON_OUT = Path("artifacts/root/artifact-governance.json")
ROOT_MAKEFILE = Path("makes/root.mk")
TOX_FILE = Path("tox.ini")
REQUIRED_MAKE_TARGETS = {
    "list-evidence-studies": "$(ROOT_ARTIFACTS_DIR)/evidence-studies.json",
    "build-evidence-book": "$(ROOT_ARTIFACTS_DIR)/evidence-book-build.json",
    "build-evidence-study": "$(ROOT_ARTIFACTS_DIR)/evidence-book-build.json",
    "validate-evidence-book": "$(ROOT_ARTIFACTS_DIR)/evidence-book-validation.json",
    "report-evidence-completeness": "$(ROOT_ARTIFACTS_DIR)/evidence-completeness.json",
    "check-evidence-completeness": "$(ROOT_ARTIFACTS_DIR)/evidence-completeness.json",
    "report-evidence-governance": "$(MAKE) report-artifact-governance",
    "check-evidence-governance": "$(MAKE) check-artifact-governance",
    "report-artifact-governance": "$(ROOT_ARTIFACTS_DIR)/artifact-governance.json",
    "check-artifact-governance": "$(ROOT_ARTIFACTS_DIR)/artifact-governance.json",
    "report-execution-surfaces": "$(ROOT_ARTIFACTS_DIR)/execution-surfaces.json",
    "check-execution-surfaces": "$(ROOT_ARTIFACTS_DIR)/execution-surfaces.json",
    "report-package-boundaries": "$(ROOT_ARTIFACTS_DIR)/package-boundaries.json",
    "check-package-boundaries": "$(ROOT_ARTIFACTS_DIR)/package-boundaries.json",
    "report-release-readiness": "$(MAKE) report-publish-readiness",
    "check-release-readiness": "$(MAKE) check-publish-readiness",
}


@dataclass(frozen=True)
class ArtifactGovernanceIssue:
    """Describe one artifact-output contract drift."""

    code: str
    path: str
    message: str


def _target_bodies(root_make_text: str) -> dict[str, str]:
    bodies: dict[str, list[str]] = {}
    current_target: str | None = None
    for raw_line in root_make_text.splitlines():
        if raw_line.startswith("\t") and current_target is not None:
            bodies[current_target].append(raw_line.strip())
            continue
        stripped = raw_line.strip()
        if not stripped or stripped.startswith(".PHONY:"):
            current_target = None
            continue
        if ":" not in stripped:
            current_target = None
            continue
        if stripped.startswith("#"):
            current_target = None
            continue
        target_name = stripped.split(":", 1)[0].strip()
        current_target = target_name
        bodies.setdefault(target_name, [])
    return {target: "\n".join(lines) for target, lines in bodies.items()}


def build_artifact_governance_report(repo_root: Path) -> dict[str, Any]:
    """Build the repository artifact-governance report."""
    repo_root = repo_root.resolve()
    issues: list[ArtifactGovernanceIssue] = []

    tox_text = (repo_root / TOX_FILE).read_text(encoding="utf-8")
    toxworkdir_line = next(
        (
            line.strip()
            for line in tox_text.splitlines()
            if line.strip().startswith("toxworkdir")
        ),
        "",
    )
    expected_toxworkdir = "toxworkdir = {tox_root}/artifacts/root/tox"
    if toxworkdir_line != expected_toxworkdir:
        issues.append(
            ArtifactGovernanceIssue(
                code="tox-artifact-path-drift",
                path=TOX_FILE.as_posix(),
                message="toxworkdir must stay under artifacts/root/tox",
            )
        )

    root_make_text = (repo_root / ROOT_MAKEFILE).read_text(encoding="utf-8")
    target_bodies = _target_bodies(root_make_text)
    make_targets: list[dict[str, Any]] = []
    for target_name, required_snippet in sorted(REQUIRED_MAKE_TARGETS.items()):
        body = target_bodies.get(target_name, "")
        present = bool(body)
        snippet_present = required_snippet in body
        if not present:
            issues.append(
                ArtifactGovernanceIssue(
                    code="missing-root-target",
                    path=ROOT_MAKEFILE.as_posix(),
                    message=f"missing required root make target {target_name}:",
                )
            )
        elif not snippet_present:
            issues.append(
                ArtifactGovernanceIssue(
                    code="make-artifact-path-drift",
                    path=ROOT_MAKEFILE.as_posix(),
                    message=(
                        f"root make target {target_name} must reference {required_snippet}"
                    ),
                )
            )
        make_targets.append(
            {
                "target_name": target_name,
                "present": present,
                "required_snippet": required_snippet,
                "snippet_present": snippet_present,
            }
        )

    return {
        "schema_version": 1,
        "toxworkdir": toxworkdir_line,
        "expected_toxworkdir": expected_toxworkdir,
        "make_targets": make_targets,
        "issue_count": len(issues),
        "issues": [asdict(issue) for issue in issues],
    }


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the artifact-governance audit."""
    parser = argparse.ArgumentParser(
        description="Audit repo execution surfaces for governed artifact output paths."
    )
    parser.add_argument("command", choices=("report", "check"))
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    return parser.parse_args()


def main() -> int:
    """Run the artifact-governance CLI entry point."""
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    payload = build_artifact_governance_report(repo_root)
    json_out = Path(args.json_out)
    if not json_out.is_absolute():
        json_out = repo_root / json_out
    _write_json(json_out, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    if args.command == "check" and payload["issue_count"]:
        raise SystemExit("artifact governance audit failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
