from __future__ import annotations

import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
BUNDLE_ROOT = (
    REPO_ROOT
    / "reports"
    / "taxon-trust"
    / "examples"
    / "taxon-identity-and-retention-workflow"
    / "reviewer-evidence"
)


def test_taxon_trust_evidence_bundle_covers_goals_21_through_30() -> None:
    goal_checks = json.loads((BUNDLE_ROOT / "goal_checks.json").read_text(encoding="utf-8"))
    manifest = json.loads((BUNDLE_ROOT / "manifest.json").read_text(encoding="utf-8"))

    assert [row["goal_id"] for row in goal_checks] == list(range(21, 31))
    assert all(row["verdict"] == "verified" for row in goal_checks)
    assert manifest["study_id"] == "taxon-trust"
    assert manifest["all_goals_verified"] is True
    assert sorted(manifest["goals"]) == list(range(21, 31))
    assert len(manifest["input_checksums"]) == 10


def test_taxon_trust_report_manifest_surfaces_workflow_taxon_sections() -> None:
    report_manifest = json.loads(
        (BUNDLE_ROOT / "taxonomy_report_machine_manifest.json").read_text(
            encoding="utf-8"
        )
    )

    assert report_manifest["report_kind"] == "taxonomy"
    assert report_manifest["sections"] == [
        "reviewer-summary",
        "taxon-audit",
        "taxon-identity",
        "taxon-safety",
        "taxon-namespaces",
        "taxon-rank-consistency",
        "taxon-duplicate-identities",
        "taxon-mapping-conflicts",
        "taxon-crosswalk",
        "taxon-exclusions",
        "taxon-loss",
        "taxon-loss-events",
        "taxon-stability",
        "limitations",
    ]
    assert report_manifest["metrics"]["crosswalk_rows"] == 4
    assert report_manifest["metrics"]["excluded_taxa"] == 2
    assert report_manifest["metrics"]["loss_stage_count"] == 3
    assert report_manifest["metrics"]["unstable_taxa"] == 3
