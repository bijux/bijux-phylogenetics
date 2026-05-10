from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.evidence.scaffolding import (
    EvidenceBundleTemplateSpec,
    build_evidence_bundle_template,
    write_evidence_bundle_template,
)


def _template_spec() -> EvidenceBundleTemplateSpec:
    return EvidenceBundleTemplateSpec(
        study_id="template-study",
        evidence_id="evidence-001",
        evidence_title="Template Bundle",
        summary="Template scaffold for a new evidence bundle.",
        owner_package="bijux-phylogenetics",
        claim_id="template-claim",
        claim_title="Template claim",
        claim_summary="Template claim summary.",
        analytical_surfaces=("templates", "scaffolding"),
        claim_tags=("template", "scaffolding"),
        source_basis_locators=(
            "packages/bijux-phylogenetics/tests/fixtures/trees/example_tree.nwk",
        ),
        governed_code_paths=(
            "packages/bijux-phylogenetics/src/bijux_phylogenetics/evidence/scaffolding.py",
        ),
        limitation="Template placeholder only.",
        freshness_date="2026-05-10",
    )


def test_build_evidence_bundle_template_uses_governed_portable_locators() -> None:
    template = build_evidence_bundle_template(_template_spec())

    manifest = template["manifest"]
    assert manifest["source_basis"][0]["locator"] == (
        "packages/bijux-phylogenetics/tests/fixtures/trees/example_tree.nwk"
    )
    assert manifest["freshness"]["governed_code_paths"] == [
        "packages/bijux-phylogenetics/src/bijux_phylogenetics/evidence/scaffolding.py"
    ]
    assert manifest["verdict"]["status"] == "not_comparable"
    assert template["report_filename"] == "results/comparison-observation.json"
    assert template["checks"]["required_result_artifacts"] == [
        "results/README.md",
        "results/manifest.json",
    ]


def test_write_evidence_bundle_template_creates_governed_bundle_files(tmp_path: Path) -> None:
    output_root = tmp_path / "bundle"
    written = write_evidence_bundle_template(output_root, _template_spec())

    assert written["manifest"].exists()
    assert written["claims"].exists()
    assert written["readme"].exists()
    assert written["reference_r"].exists()
    assert written["analysis_py"].exists()
    assert written["checks"].exists()
    assert written["human_report"].exists()
    assert written["provenance"].exists()
    assert written["results_readme"].exists()
    assert written["results_manifest"].exists()
    assert written["report"].exists()
    assert json.loads(written["manifest"].read_text(encoding="utf-8"))["study_id"] == (
        "template-study"
    )
    assert written["report"].name == "comparison-observation.json"
