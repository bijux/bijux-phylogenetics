from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics_dev.quality.evidence_completeness import (
    build_evidence_completeness_report,
)

REPO_ROOT = Path(__file__).resolve().parents[3]


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _minimal_repo(tmp_path: Path) -> Path:
    repo_root = tmp_path / "repo"
    _write(
        repo_root / "configs" / "publication_readiness.toml",
        """
[tool.bijux_phylogenetics.publication_readiness]
required_evidence_input_manifest = "inputs.manifest.json"
required_evidence_bundle_artifacts = ["reference.R", "analysis.py", "checks.json", "report.md", "provenance.json"]
""".strip()
        + "\n",
    )
    bundle_root = repo_root / "evidence-book" / "studies" / "demo-study" / "evidence-001"
    _write(
        bundle_root / "manifest.json",
        json.dumps(
            {
                "evidence_id": "evidence-001",
                "study_id": "demo-study",
                "owner_package": "demo-runtime",
                "claim_ids": ["demo-claim"],
                "comparison_mode": "direct_r_parity",
                "verdict": {"status": "matched"},
            },
            indent=2,
        )
        + "\n",
    )
    for name in (
        "reference.R",
        "analysis.py",
        "checks.json",
        "report.md",
        "provenance.json",
        "inputs.manifest.json",
    ):
        _write(bundle_root / name, "{}\n")
    return repo_root


def test_build_evidence_completeness_report_accepts_complete_bundle(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)

    report = build_evidence_completeness_report(repo_root)

    assert report["issue_count"] == 0


def test_build_evidence_completeness_report_flags_missing_bundle_surface(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    (repo_root / "evidence-book" / "studies" / "demo-study" / "evidence-001" / "report.md").unlink()

    report = build_evidence_completeness_report(repo_root)

    issue_codes = {issue["code"] for issue in report["issues"]}
    assert "missing-evidence-bundle-surface" in issue_codes


def test_repository_evidence_completeness_report_is_clean() -> None:
    report = build_evidence_completeness_report(REPO_ROOT)

    assert report["issue_count"] == 0
