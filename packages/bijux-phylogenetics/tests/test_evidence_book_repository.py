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
from bijux_phylogenetics.evidence.teaching import (
    build_migration_guide,
    build_student_safe_reproducibility_contract,
    build_teaching_and_migration_index,
    build_teaching_guide,
    render_migration_guide_markdown,
    render_student_safe_reproducibility_markdown,
    render_teaching_and_migration_index_markdown,
    render_teaching_guide_markdown,
)


REPO_ROOT = Path(__file__).resolve().parents[3]


def test_repository_evidence_book_passes_validation() -> None:
    report = validate_evidence_book(REPO_ROOT)

    assert report.valid is True, [
        f"{issue.path.as_posix()}: {issue.message}" for issue in report.issues
    ]


def test_repository_evidence_book_index_matches_generated_payload() -> None:
    index_path = REPO_ROOT / "evidence-book" / "index" / "evidence-index.json"
    catalog_path = REPO_ROOT / "evidence-book" / "index" / "catalog.md"
    claim_map_path = REPO_ROOT / "evidence-book" / "index" / "claim-map.json"
    parity_dashboard_path = (
        REPO_ROOT / "evidence-book" / "index" / "parity-dashboard.json"
    )
    parity_summary_path = REPO_ROOT / "evidence-book" / "index" / "parity-dashboard.md"
    mismatch_archive_path = (
        REPO_ROOT / "evidence-book" / "index" / "mismatch-archive.json"
    )
    mismatch_summary_path = (
        REPO_ROOT / "evidence-book" / "index" / "mismatch-archive.md"
    )
    verdict_workflows_path = (
        REPO_ROOT / "evidence-book" / "index" / "verdict-workflows.json"
    )
    verdict_workflows_summary_path = (
        REPO_ROOT / "evidence-book" / "index" / "verdict-workflows.md"
    )
    false_confidence_audit_path = (
        REPO_ROOT / "evidence-book" / "index" / "false-confidence-audit.json"
    )
    false_confidence_summary_path = (
        REPO_ROOT / "evidence-book" / "index" / "false-confidence-audit.md"
    )
    scientific_debt_path = (
        REPO_ROOT / "evidence-book" / "index" / "scientific-debt-register.json"
    )
    scientific_debt_summary_path = (
        REPO_ROOT / "evidence-book" / "index" / "scientific-debt-register.md"
    )
    portability_path = REPO_ROOT / "evidence-book" / "index" / "portability-audit.json"
    portability_summary_path = (
        REPO_ROOT / "evidence-book" / "index" / "portability-audit.md"
    )
    fragile_example_path = (
        REPO_ROOT / "evidence-book" / "index" / "fragile-example-audit.json"
    )
    fragile_example_summary_path = (
        REPO_ROOT / "evidence-book" / "index" / "fragile-example-audit.md"
    )
    regeneration_path = (
        REPO_ROOT / "evidence-book" / "index" / "regeneration-contract.json"
    )
    regeneration_summary_path = (
        REPO_ROOT / "evidence-book" / "index" / "regeneration-contract.md"
    )
    freshness_report_path = (
        REPO_ROOT / "evidence-book" / "index" / "freshness-report.json"
    )
    freshness_summary_path = (
        REPO_ROOT / "evidence-book" / "index" / "freshness-report.md"
    )
    integrity_report_path = (
        REPO_ROOT / "evidence-book" / "index" / "integrity-report.json"
    )
    integrity_summary_path = (
        REPO_ROOT / "evidence-book" / "index" / "integrity-report.md"
    )
    coverage_gap_path = REPO_ROOT / "evidence-book" / "index" / "coverage-gaps.json"
    coverage_gap_summary_path = (
        REPO_ROOT / "evidence-book" / "index" / "coverage-gaps.md"
    )

    payload = build_evidence_book_index(REPO_ROOT)
    catalog = render_evidence_catalog(payload)
    claim_map = build_evidence_claim_map(REPO_ROOT)
    parity_dashboard = build_evidence_parity_dashboard(REPO_ROOT)
    mismatch_archive = build_evidence_mismatch_archive(REPO_ROOT)
    verdict_workflows = build_evidence_verdict_workflows(REPO_ROOT)
    false_confidence_audit = build_evidence_false_confidence_audit(REPO_ROOT)
    fragile_example_audit = build_evidence_fragile_example_audit(REPO_ROOT)
    scientific_debt_register = build_evidence_scientific_debt_register(REPO_ROOT)
    portability_audit = build_evidence_portability_audit(REPO_ROOT)
    regeneration_contract = build_evidence_regeneration_contract(REPO_ROOT)
    freshness_report = build_evidence_freshness_report(REPO_ROOT)
    integrity_report = build_evidence_integrity_report(REPO_ROOT)
    coverage_gap_report = build_evidence_coverage_gap_report(REPO_ROOT)

    assert json.loads(index_path.read_text(encoding="utf-8")) == payload
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


def test_repository_teaching_and_migration_surfaces_match_generated_payloads() -> None:
    study_ids = [
        "primate-longevity-signal",
        "primate-pgls-and-signal",
    ]
    teaching_guides = []
    migration_guides = []
    reproducibility_contracts = []

    for study_id in study_ids:
        study_root = REPO_ROOT / "evidence-book" / "studies" / study_id
        study_manifest = json.loads(
            (study_root / "study.json").read_text(encoding="utf-8")
        )
        family_index = json.loads(
            (study_root / "family-index.json").read_text(encoding="utf-8")
        )
        source_fragment_map = json.loads(
            (study_root / "source-fragment-map.json").read_text(encoding="utf-8")
        )
        bundle_manifests = [
            json.loads(path.read_text(encoding="utf-8"))
            for path in sorted(study_root.glob("evidence-*/manifest.json"))
        ]

        teaching_payload = build_teaching_guide(
            study_manifest,
            family_index,
            source_fragment_map,
        )
        migration_payload = build_migration_guide(
            study_manifest,
            source_fragment_map,
            bundle_manifests,
        )
        reproducibility_payload = build_student_safe_reproducibility_contract(
            study_manifest
        )

        assert (
            json.loads((study_root / "teaching-guide.json").read_text(encoding="utf-8"))
            == teaching_payload
        )
        assert (study_root / "teaching-guide.md").read_text(
            encoding="utf-8"
        ) == render_teaching_guide_markdown(teaching_payload)
        assert (
            json.loads(
                (study_root / "migration-guide.json").read_text(encoding="utf-8")
            )
            == migration_payload
        )
        assert (study_root / "migration-guide.md").read_text(
            encoding="utf-8"
        ) == render_migration_guide_markdown(migration_payload)
        assert (
            json.loads(
                (study_root / "student-safe-reproducibility.json").read_text(
                    encoding="utf-8"
                )
            )
            == reproducibility_payload
        )
        assert (study_root / "student-safe-reproducibility.md").read_text(
            encoding="utf-8"
        ) == render_student_safe_reproducibility_markdown(reproducibility_payload)

        teaching_guides.append(teaching_payload)
        migration_guides.append(migration_payload)
        reproducibility_contracts.append(reproducibility_payload)

    landing_payload = build_teaching_and_migration_index(
        teaching_guides,
        migration_guides,
        reproducibility_contracts,
    )
    landing_root = REPO_ROOT / "evidence-book" / "index"

    assert (
        json.loads(
            (landing_root / "teaching-and-migration.json").read_text(encoding="utf-8")
        )
        == landing_payload
    )
    assert (landing_root / "teaching-and-migration.md").read_text(
        encoding="utf-8"
    ) == render_teaching_and_migration_index_markdown(landing_payload)
