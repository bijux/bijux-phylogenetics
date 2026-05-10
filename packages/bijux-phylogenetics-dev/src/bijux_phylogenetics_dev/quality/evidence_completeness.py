"""Audit that every evidence bundle is structurally complete for CI and review."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any

from .package_bundles import load_publication_readiness_settings

DEFAULT_JSON_OUT = Path("artifacts/root/evidence-completeness.json")


@dataclass(frozen=True)
class EvidenceCompletenessIssue:
    code: str
    path: str
    message: str


def _manifest_payload(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return payload if isinstance(payload, dict) else {}


def _bundle_roots(repo_root: Path) -> list[Path]:
    return sorted(
        path.parent
        for path in (repo_root / "evidence-book" / "studies").glob("*/evidence-*/manifest.json")
    )


def build_evidence_completeness_report(repo_root: Path) -> dict[str, Any]:
    repo_root = repo_root.resolve()
    publication_settings = load_publication_readiness_settings(repo_root)
    required_bundle_artifacts = [
        str(name)
        for name in publication_settings.get("required_evidence_bundle_artifacts", [])
        if isinstance(name, str)
    ]
    required_input_manifest = str(
        publication_settings.get("required_evidence_input_manifest", "inputs.manifest.json")
    )
    issues: list[EvidenceCompletenessIssue] = []
    bundle_reports: list[dict[str, Any]] = []

    for bundle_root in _bundle_roots(repo_root):
        manifest_path = bundle_root / "manifest.json"
        manifest = _manifest_payload(manifest_path)
        missing_paths = [
            name
            for name in [*required_bundle_artifacts, required_input_manifest]
            if not (bundle_root / name).is_file()
        ]
        missing_fields = [
            field
            for field in ("evidence_id", "study_id", "owner_package", "claim_ids", "comparison_mode", "verdict")
            if field not in manifest or manifest[field] in ("", [], {})
        ]
        if missing_paths:
            issues.append(
                EvidenceCompletenessIssue(
                    code="missing-evidence-bundle-surface",
                    path=bundle_root.relative_to(repo_root).as_posix(),
                    message=f"bundle is missing governed files: {', '.join(missing_paths)}",
                )
            )
        if missing_fields:
            issues.append(
                EvidenceCompletenessIssue(
                    code="incomplete-evidence-manifest",
                    path=manifest_path.relative_to(repo_root).as_posix(),
                    message=f"manifest is missing required fields: {', '.join(missing_fields)}",
                )
            )
        bundle_reports.append(
            {
                "bundle_path": bundle_root.relative_to(repo_root).as_posix(),
                "missing_paths": missing_paths,
                "missing_fields": missing_fields,
            }
        )

    return {
        "schema_version": 1,
        "bundle_count": len(bundle_reports),
        "required_bundle_artifacts": required_bundle_artifacts,
        "required_input_manifest": required_input_manifest,
        "bundles": bundle_reports,
        "issue_count": len(issues),
        "issues": [asdict(issue) for issue in issues],
    }


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def check_evidence_completeness(
    repo_root: Path,
    *,
    json_out: Path | None = None,
) -> dict[str, Any]:
    payload = build_evidence_completeness_report(repo_root)
    if json_out is not None:
        _write_json(json_out, payload)
    if payload["issue_count"]:
        raise SystemExit("evidence completeness audit failed")
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit evidence bundle completeness for CI and review."
    )
    parser.add_argument("command", choices=("report", "check"))
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    json_out = Path(args.json_out)
    if not json_out.is_absolute():
        json_out = repo_root / json_out
    if args.command == "check":
        payload = check_evidence_completeness(repo_root, json_out=json_out)
    else:
        payload = build_evidence_completeness_report(repo_root)
        _write_json(json_out, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
