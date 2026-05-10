from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
BUNDLE_ROOT = (
    REPO_ROOT
    / "evidence-book"
    / "studies"
    / "taxon-trust"
    / "evidence-001"
)


def test_taxon_trust_evidence_bundle_records_expected_claims() -> None:
    claim_verdicts = json.loads(
        (BUNDLE_ROOT / "claim_verdicts.json").read_text(encoding="utf-8")
    )
    manifest = json.loads((BUNDLE_ROOT / "manifest.json").read_text(encoding="utf-8"))

    assert [row["claim_id"] for row in claim_verdicts] == [
        "taxon-spelling-variant-audit",
        "taxonomic-synonym-candidate-detection",
        "controlled-synonym-resolution",
        "ambiguous-synonym-rejection",
        "taxon-namespace-classification",
        "mixed-namespace-detection",
        "taxon-crosswalk-table",
        "taxon-exclusion-reasoning",
        "workflow-taxon-loss-report",
        "taxon-stability-report",
    ]
    assert all(row["verdict"] == "verified" for row in claim_verdicts)
    assert manifest["schema_version"] == 1
    assert manifest["study_id"] == "taxon-trust"
    assert manifest["evidence_id"] == "evidence-001"
    assert manifest["verdict"]["status"] == "matched"
    assert sorted(manifest["claim_ids"]) == sorted(
        row["claim_id"] for row in claim_verdicts
    )
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
