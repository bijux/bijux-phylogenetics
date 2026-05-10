from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.evidence.book import (
    build_evidence_claim_map,
    build_evidence_book_index,
    build_evidence_false_confidence_audit,
    build_evidence_fragile_example_audit,
    build_evidence_mismatch_archive,
    build_evidence_parity_dashboard,
    build_evidence_portability_audit,
    build_evidence_regeneration_contract,
    build_evidence_scientific_debt_register,
    build_evidence_verdict_workflows,
    render_evidence_catalog,
    render_evidence_false_confidence_audit,
    render_evidence_fragile_example_audit,
    render_evidence_mismatch_archive,
    render_evidence_parity_dashboard,
    render_evidence_portability_audit,
    render_evidence_regeneration_contract,
    render_evidence_scientific_debt_register,
    render_evidence_verdict_workflows,
    validate_evidence_book,
    write_evidence_book_index,
)
from bijux_phylogenetics.evidence.bundle_artifacts import build_bundle_governed_artifacts
from bijux_phylogenetics.evidence.bundle_contracts import (
    ARTIFACT_JSON_FILENAMES,
    RESULT_ARTIFACT_JSON_FILENAMES,
)
from bijux_phylogenetics.evidence.closure import (
    build_analytical_surface_coverage,
    build_claim_reaudit,
    build_closure_criteria,
    build_completion_gates,
    build_evidence_maturity_scorecard,
    build_evidence_review_ritual,
    render_analytical_surface_coverage,
    render_claim_reaudit,
    render_closure_criteria,
    render_completion_gates,
    render_evidence_maturity_scorecard,
    render_evidence_review_ritual,
)
from bijux_phylogenetics.evidence.coverage import (
    build_evidence_coverage_gap_report,
    render_evidence_coverage_gap_report,
)
from bijux_phylogenetics.evidence.freshness import (
    build_evidence_freshness_report,
    render_evidence_freshness_report,
)
from bijux_phylogenetics.evidence.integrity import (
    build_evidence_integrity_report,
    render_evidence_integrity_report,
)
from bijux_phylogenetics.evidence.workbench import (
    DOCS_EVIDENCE_OVERVIEW,
    refresh_evidence_book,
)


def _write_book_fixture(root: Path) -> Path:
    book_root = root / "evidence-book"
    study_root = book_root / "studies" / "taxon-trust"
    bundle_root = study_root / "evidence-001"
    index_root = book_root / "index"
    bundle_root.mkdir(parents=True, exist_ok=True)
    index_root.mkdir(parents=True, exist_ok=True)
    (book_root / "README.md").write_text("# Evidence Book\n", encoding="utf-8")
    (study_root / "README.md").write_text("# Taxon Trust\n", encoding="utf-8")
    (study_root / "study.json").write_text(
        json.dumps(
            {
                "study_id": "taxon-trust",
                "study_title": "Taxon Trust",
                "summary": "Fixture-backed taxon evidence.",
                "owner_package": "bijux-phylogenetics",
                "study_categories": ["scientific-validation"],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    (bundle_root / "README.md").write_text("# Evidence 001\n", encoding="utf-8")
    (bundle_root / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "study_id": "taxon-trust",
                "evidence_id": "evidence-001",
                "evidence_title": "Taxon workflow review bundle",
                "summary": "Validates taxon-workflow trust surfaces.",
                "owner_package": "bijux-phylogenetics",
                "claim_ids": ["taxonomy-review"],
                "source_basis": [
                    {
                        "kind": "repository-fixture",
                        "label": "example taxon workflow tree",
                        "locator": "packages/bijux-phylogenetics/tests/fixtures/trees/example_taxon_workflow_tree.nwk",
                    }
                ],
                "freshness": {
                    "last_generated_on": "2026-05-10",
                    "governed_code_paths": [
                        "packages/bijux-phylogenetics/src/bijux_phylogenetics"
                    ],
                    "source_basis_locators": [
                        "packages/bijux-phylogenetics/tests/fixtures/trees/example_taxon_workflow_tree.nwk"
                    ],
                },
                "ownership": {
                    "owner_package": "bijux-phylogenetics",
                    "analytical_surfaces": ["taxonomy"],
                },
                "claim_tags": ["taxonomy", "review"],
                "comparison_mode": "bijux_native_reinterpretation",
                "verdict": {
                    "status": "matched",
                    "summary": "Observed output matches the checked-in fixture expectations.",
                },
                "limitations": ["Covers one workflow family only."],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    (bundle_root / "claims.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "study_id": "taxon-trust",
                "evidence_id": "evidence-001",
                "claim_count": 1,
                "claims": [
                    {
                        "claim_id": "taxonomy-review",
                        "claim_title": "Taxonomy review",
                        "summary": "Fixture-backed taxonomy review contract.",
                        "verdict": "matched",
                        "evidence_ids": ["evidence-001"],
                        "source_fragments": ["fixture-taxonomy-review"],
                    }
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    for filename, payload in build_bundle_governed_artifacts(root, bundle_root).items():
        target = bundle_root / filename
        target.parent.mkdir(parents=True, exist_ok=True)
        if filename in ARTIFACT_JSON_FILENAMES or filename in RESULT_ARTIFACT_JSON_FILENAMES:
            target.write_text(
                json.dumps(payload, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        else:
            if not isinstance(payload, str):
                raise AssertionError(f"expected text artifact for {filename}")
            target.write_text(payload, encoding="utf-8")
    return root


def test_validate_evidence_book_accepts_governed_layout(tmp_path: Path) -> None:
    repo_root = _write_book_fixture(tmp_path)
    refresh_evidence_book(repo_root)

    report = validate_evidence_book(repo_root)

    assert report.valid is True
    assert report.issues == []
    assert [path.name for path in report.bundle_paths] == ["evidence-001"]


def test_validate_evidence_book_rejects_missing_manifest_fields(tmp_path: Path) -> None:
    repo_root = _write_book_fixture(tmp_path)
    manifest_path = (
        repo_root
        / "evidence-book"
        / "studies"
        / "taxon-trust"
        / "evidence-001"
        / "manifest.json"
    )
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    del payload["verdict"]
    manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    report = validate_evidence_book(repo_root, require_index_outputs=False)

    assert report.valid is False
    assert any("missing keys: verdict" in issue.message for issue in report.issues)


def test_validate_evidence_book_rejects_missing_reviewer_summary(
    tmp_path: Path,
) -> None:
    repo_root = _write_book_fixture(tmp_path)
    refresh_evidence_book(repo_root)
    reviewer_summary_path = (
        repo_root
        / "evidence-book"
        / "studies"
        / "taxon-trust"
        / "evidence-001"
        / "reviewer-summary.json"
    )
    reviewer_summary_path.unlink()

    report = validate_evidence_book(repo_root)

    assert report.valid is False
    assert any(
        "reviewer-summary.json" in issue.message
        or "reviewer-summary.json" in issue.path.as_posix()
        for issue in report.issues
    )


def test_validate_evidence_book_rejects_missing_governed_local_artifact(
    tmp_path: Path,
) -> None:
    repo_root = _write_book_fixture(tmp_path)
    reference_path = (
        repo_root
        / "evidence-book"
        / "studies"
        / "taxon-trust"
        / "evidence-001"
        / "reference.R"
    )
    reference_path.unlink()

    report = validate_evidence_book(repo_root, require_index_outputs=False)

    assert report.valid is False
    assert any("reference.R" in issue.path.as_posix() for issue in report.issues)


def test_validate_evidence_book_rejects_missing_freshness_report(
    tmp_path: Path,
) -> None:
    repo_root = _write_book_fixture(tmp_path)
    refresh_evidence_book(repo_root)
    freshness_report_path = (
        repo_root / "evidence-book" / "index" / "freshness-report.json"
    )
    freshness_report_path.unlink()

    report = validate_evidence_book(repo_root)

    assert report.valid is False
    assert any(
        "freshness-report.json" in issue.path.as_posix() for issue in report.issues
    )


def test_validate_evidence_book_rejects_missing_claim_reaudit(tmp_path: Path) -> None:
    repo_root = _write_book_fixture(tmp_path)
    refresh_evidence_book(repo_root)
    claim_reaudit_path = repo_root / "evidence-book" / "index" / "claim-reaudit.json"
    claim_reaudit_path.unlink()

    report = validate_evidence_book(repo_root)

    assert report.valid is False
    assert any("claim-reaudit.json" in issue.path.as_posix() for issue in report.issues)


def test_write_evidence_book_index_renders_catalog_from_index(tmp_path: Path) -> None:
    repo_root = _write_book_fixture(tmp_path)

    refresh_evidence_book(repo_root)
    index_path, catalog_path = write_evidence_book_index(repo_root)
    payload = build_evidence_book_index(repo_root)
    catalog = render_evidence_catalog(payload)
    claim_map = build_evidence_claim_map(repo_root)
    parity_dashboard = build_evidence_parity_dashboard(repo_root)
    mismatch_archive = build_evidence_mismatch_archive(repo_root)
    verdict_workflows = build_evidence_verdict_workflows(repo_root)
    false_confidence_audit = build_evidence_false_confidence_audit(repo_root)
    fragile_example_audit = build_evidence_fragile_example_audit(repo_root)
    scientific_debt_register = build_evidence_scientific_debt_register(repo_root)
    portability_audit = build_evidence_portability_audit(repo_root)
    regeneration_contract = build_evidence_regeneration_contract(repo_root)
    freshness_report = build_evidence_freshness_report(repo_root)
    integrity_report = build_evidence_integrity_report(repo_root)
    coverage_gap_report = build_evidence_coverage_gap_report(repo_root)
    claim_reaudit = build_claim_reaudit(repo_root)
    analytical_surface_coverage = build_analytical_surface_coverage(repo_root)
    closure_criteria = build_closure_criteria(repo_root)
    maturity_scorecard = build_evidence_maturity_scorecard(repo_root)
    review_ritual = build_evidence_review_ritual(repo_root)
    completion_gates = build_completion_gates(repo_root)
    claim_map_path = repo_root / "evidence-book" / "index" / "claim-map.json"
    parity_dashboard_path = (
        repo_root / "evidence-book" / "index" / "parity-dashboard.json"
    )
    parity_summary_path = repo_root / "evidence-book" / "index" / "parity-dashboard.md"
    mismatch_archive_path = (
        repo_root / "evidence-book" / "index" / "mismatch-archive.json"
    )
    mismatch_summary_path = (
        repo_root / "evidence-book" / "index" / "mismatch-archive.md"
    )
    verdict_workflows_path = (
        repo_root / "evidence-book" / "index" / "verdict-workflows.json"
    )
    verdict_workflows_summary_path = (
        repo_root / "evidence-book" / "index" / "verdict-workflows.md"
    )
    false_confidence_audit_path = (
        repo_root / "evidence-book" / "index" / "false-confidence-audit.json"
    )
    false_confidence_summary_path = (
        repo_root / "evidence-book" / "index" / "false-confidence-audit.md"
    )
    scientific_debt_path = (
        repo_root / "evidence-book" / "index" / "scientific-debt-register.json"
    )
    scientific_debt_summary_path = (
        repo_root / "evidence-book" / "index" / "scientific-debt-register.md"
    )
    portability_path = repo_root / "evidence-book" / "index" / "portability-audit.json"
    portability_summary_path = (
        repo_root / "evidence-book" / "index" / "portability-audit.md"
    )
    fragile_example_path = (
        repo_root / "evidence-book" / "index" / "fragile-example-audit.json"
    )
    fragile_example_summary_path = (
        repo_root / "evidence-book" / "index" / "fragile-example-audit.md"
    )
    regeneration_path = (
        repo_root / "evidence-book" / "index" / "regeneration-contract.json"
    )
    regeneration_summary_path = (
        repo_root / "evidence-book" / "index" / "regeneration-contract.md"
    )
    teaching_migration_path = (
        repo_root / "evidence-book" / "index" / "teaching-and-migration.json"
    )
    teaching_migration_summary_path = (
        repo_root / "evidence-book" / "index" / "teaching-and-migration.md"
    )
    freshness_report_path = (
        repo_root / "evidence-book" / "index" / "freshness-report.json"
    )
    freshness_summary_path = (
        repo_root / "evidence-book" / "index" / "freshness-report.md"
    )
    integrity_report_path = (
        repo_root / "evidence-book" / "index" / "integrity-report.json"
    )
    integrity_summary_path = (
        repo_root / "evidence-book" / "index" / "integrity-report.md"
    )
    coverage_gap_path = repo_root / "evidence-book" / "index" / "coverage-gaps.json"
    coverage_gap_summary_path = (
        repo_root / "evidence-book" / "index" / "coverage-gaps.md"
    )
    claim_reaudit_path = repo_root / "evidence-book" / "index" / "claim-reaudit.json"
    claim_reaudit_summary_path = (
        repo_root / "evidence-book" / "index" / "claim-reaudit.md"
    )
    analytical_surface_coverage_path = (
        repo_root / "evidence-book" / "index" / "analytical-surface-coverage.json"
    )
    analytical_surface_coverage_summary_path = (
        repo_root / "evidence-book" / "index" / "analytical-surface-coverage.md"
    )
    closure_criteria_path = (
        repo_root / "evidence-book" / "index" / "closure-criteria.json"
    )
    closure_criteria_summary_path = (
        repo_root / "evidence-book" / "index" / "closure-criteria.md"
    )
    maturity_scorecard_path = (
        repo_root / "evidence-book" / "index" / "evidence-maturity-scorecard.json"
    )
    maturity_scorecard_summary_path = (
        repo_root / "evidence-book" / "index" / "evidence-maturity-scorecard.md"
    )
    review_ritual_path = (
        repo_root / "evidence-book" / "index" / "evidence-review-ritual.json"
    )
    review_ritual_summary_path = (
        repo_root / "evidence-book" / "index" / "evidence-review-ritual.md"
    )
    completion_gates_path = (
        repo_root / "evidence-book" / "index" / "completion-gates.json"
    )
    completion_gates_summary_path = (
        repo_root / "evidence-book" / "index" / "completion-gates.md"
    )
    reviewer_summary_json_path = (
        repo_root
        / "evidence-book"
        / "studies"
        / "taxon-trust"
        / "evidence-001"
        / "reviewer-summary.json"
    )
    reviewer_summary_markdown_path = (
        repo_root
        / "evidence-book"
        / "studies"
        / "taxon-trust"
        / "evidence-001"
        / "reviewer-summary.md"
    )
    docs_evidence_path = repo_root / DOCS_EVIDENCE_OVERVIEW

    assert index_path.exists()
    assert catalog_path.exists()
    assert claim_map_path.exists()
    assert parity_dashboard_path.exists()
    assert parity_summary_path.exists()
    assert mismatch_archive_path.exists()
    assert mismatch_summary_path.exists()
    assert verdict_workflows_path.exists()
    assert verdict_workflows_summary_path.exists()
    assert false_confidence_audit_path.exists()
    assert false_confidence_summary_path.exists()
    assert scientific_debt_path.exists()
    assert scientific_debt_summary_path.exists()
    assert portability_path.exists()
    assert portability_summary_path.exists()
    assert fragile_example_path.exists()
    assert fragile_example_summary_path.exists()
    assert regeneration_path.exists()
    assert regeneration_summary_path.exists()
    assert teaching_migration_path.exists()
    assert teaching_migration_summary_path.exists()
    assert freshness_report_path.exists()
    assert freshness_summary_path.exists()
    assert integrity_report_path.exists()
    assert integrity_summary_path.exists()
    assert coverage_gap_path.exists()
    assert coverage_gap_summary_path.exists()
    assert claim_reaudit_path.exists()
    assert claim_reaudit_summary_path.exists()
    assert analytical_surface_coverage_path.exists()
    assert analytical_surface_coverage_summary_path.exists()
    assert closure_criteria_path.exists()
    assert closure_criteria_summary_path.exists()
    assert maturity_scorecard_path.exists()
    assert maturity_scorecard_summary_path.exists()
    assert review_ritual_path.exists()
    assert review_ritual_summary_path.exists()
    assert completion_gates_path.exists()
    assert completion_gates_summary_path.exists()
    assert reviewer_summary_json_path.exists()
    assert reviewer_summary_markdown_path.exists()
    assert docs_evidence_path.exists()
    assert payload["study_count"] == 1
    assert payload["evidence_count"] == 1
    assert "Taxon Trust" in catalog
    assert catalog_path.read_text(encoding="utf-8") == catalog
    assert json.loads(claim_map_path.read_text(encoding="utf-8")) == claim_map
    assert (
        json.loads(parity_dashboard_path.read_text(encoding="utf-8"))
        == parity_dashboard
    )
    assert parity_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_parity_dashboard(parity_dashboard)
    assert (
        json.loads(mismatch_archive_path.read_text(encoding="utf-8"))
        == mismatch_archive
    )
    assert mismatch_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_mismatch_archive(mismatch_archive)
    assert (
        json.loads(verdict_workflows_path.read_text(encoding="utf-8"))
        == verdict_workflows
    )
    assert verdict_workflows_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_verdict_workflows(verdict_workflows)
    assert (
        json.loads(false_confidence_audit_path.read_text(encoding="utf-8"))
        == false_confidence_audit
    )
    assert false_confidence_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_false_confidence_audit(false_confidence_audit)
    assert (
        json.loads(scientific_debt_path.read_text(encoding="utf-8"))
        == scientific_debt_register
    )
    assert scientific_debt_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_scientific_debt_register(scientific_debt_register)
    assert json.loads(portability_path.read_text(encoding="utf-8")) == portability_audit
    assert portability_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_portability_audit(portability_audit)
    assert (
        json.loads(fragile_example_path.read_text(encoding="utf-8"))
        == fragile_example_audit
    )
    assert fragile_example_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_fragile_example_audit(fragile_example_audit)
    assert (
        json.loads(regeneration_path.read_text(encoding="utf-8"))
        == regeneration_contract
    )
    assert regeneration_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_regeneration_contract(regeneration_contract)
    assert (
        json.loads(freshness_report_path.read_text(encoding="utf-8"))
        == freshness_report
    )
    assert freshness_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_freshness_report(freshness_report)
    assert (
        json.loads(integrity_report_path.read_text(encoding="utf-8"))
        == integrity_report
    )
    assert integrity_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_integrity_report(integrity_report)
    assert (
        json.loads(coverage_gap_path.read_text(encoding="utf-8")) == coverage_gap_report
    )
    assert coverage_gap_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_coverage_gap_report(coverage_gap_report)
    assert json.loads(claim_reaudit_path.read_text(encoding="utf-8")) == claim_reaudit
    assert claim_reaudit_summary_path.read_text(
        encoding="utf-8"
    ) == render_claim_reaudit(claim_reaudit)
    assert (
        json.loads(analytical_surface_coverage_path.read_text(encoding="utf-8"))
        == analytical_surface_coverage
    )
    assert analytical_surface_coverage_summary_path.read_text(
        encoding="utf-8"
    ) == render_analytical_surface_coverage(analytical_surface_coverage)
    assert (
        json.loads(closure_criteria_path.read_text(encoding="utf-8"))
        == closure_criteria
    )
    assert closure_criteria_summary_path.read_text(
        encoding="utf-8"
    ) == render_closure_criteria(closure_criteria)
    assert (
        json.loads(maturity_scorecard_path.read_text(encoding="utf-8"))
        == maturity_scorecard
    )
    assert maturity_scorecard_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_maturity_scorecard(maturity_scorecard)
    assert json.loads(review_ritual_path.read_text(encoding="utf-8")) == review_ritual
    assert review_ritual_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_review_ritual(review_ritual)
    assert (
        json.loads(completion_gates_path.read_text(encoding="utf-8"))
        == completion_gates
    )
    assert completion_gates_summary_path.read_text(
        encoding="utf-8"
    ) == render_completion_gates(completion_gates)
