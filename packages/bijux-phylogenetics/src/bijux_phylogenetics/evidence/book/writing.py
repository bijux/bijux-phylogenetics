from __future__ import annotations

import json
from pathlib import Path

from ..teaching import (
    TEACHING_AND_MIGRATION_INDEX_FILENAME,
    TEACHING_AND_MIGRATION_SUMMARY_FILENAME,
    build_teaching_and_migration_index,
    render_teaching_and_migration_index_markdown,
)
from .aggregation import (
    build_evidence_book_index,
    build_evidence_claim_map,
    build_evidence_parity_dashboard,
)
from .audits import (
    build_evidence_false_confidence_audit,
    build_evidence_fragile_example_audit,
    build_evidence_mismatch_archive,
    build_evidence_portability_audit,
    build_evidence_regeneration_contract,
    build_evidence_scientific_debt_register,
    build_evidence_verdict_workflows,
)
from .layout import (
    EVIDENCE_CATALOG,
    EVIDENCE_CLAIM_MAP,
    EVIDENCE_FALSE_CONFIDENCE_AUDIT,
    EVIDENCE_FALSE_CONFIDENCE_SUMMARY,
    EVIDENCE_FRAGILE_EXAMPLE_AUDIT,
    EVIDENCE_FRAGILE_EXAMPLE_SUMMARY,
    EVIDENCE_INDEX,
    EVIDENCE_INDEX_DIRNAME,
    EVIDENCE_MISMATCH_ARCHIVE,
    EVIDENCE_MISMATCH_SUMMARY,
    EVIDENCE_PARITY_DASHBOARD,
    EVIDENCE_PARITY_SUMMARY,
    EVIDENCE_PORTABILITY_AUDIT,
    EVIDENCE_PORTABILITY_SUMMARY,
    EVIDENCE_REGENERATION_CONTRACT,
    EVIDENCE_REGENERATION_SUMMARY,
    EVIDENCE_SCIENTIFIC_DEBT_REGISTER,
    EVIDENCE_SCIENTIFIC_DEBT_SUMMARY,
    EVIDENCE_VERDICT_WORKFLOWS,
    EVIDENCE_VERDICT_WORKFLOWS_SUMMARY,
    evidence_book_root,
)
from .rendering import (
    render_evidence_catalog,
    render_evidence_false_confidence_audit,
    render_evidence_fragile_example_audit,
    render_evidence_mismatch_archive,
    render_evidence_parity_dashboard,
    render_evidence_portability_audit,
    render_evidence_regeneration_contract,
    render_evidence_scientific_debt_register,
    render_evidence_verdict_workflows,
)


def write_evidence_book_index(repo_root: Path) -> tuple[Path, Path]:
    root = evidence_book_root(repo_root)
    index_root = root / EVIDENCE_INDEX_DIRNAME
    index_root.mkdir(parents=True, exist_ok=True)
    teaching_migration_path = index_root / TEACHING_AND_MIGRATION_INDEX_FILENAME
    teaching_migration_summary_path = (
        index_root / TEACHING_AND_MIGRATION_SUMMARY_FILENAME
    )
    teaching_migration_payload = build_teaching_and_migration_index(
        [],
        [],
        [],
    )
    teaching_migration_path.write_text(
        json.dumps(teaching_migration_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    teaching_migration_summary_path.write_text(
        render_teaching_and_migration_index_markdown(teaching_migration_payload),
        encoding="utf-8",
    )
    payload = build_evidence_book_index(repo_root)
    claim_map_payload = build_evidence_claim_map(repo_root)
    parity_dashboard_payload = build_evidence_parity_dashboard(repo_root)
    mismatch_archive_payload = build_evidence_mismatch_archive(repo_root)
    verdict_workflows_payload = build_evidence_verdict_workflows(repo_root)
    false_confidence_payload = build_evidence_false_confidence_audit(repo_root)
    scientific_debt_payload = build_evidence_scientific_debt_register(repo_root)
    fragile_examples_payload = build_evidence_fragile_example_audit(repo_root)
    regeneration_payload = build_evidence_regeneration_contract(repo_root)
    index_path = index_root / EVIDENCE_INDEX
    catalog_path = index_root / EVIDENCE_CATALOG
    claim_map_path = index_root / EVIDENCE_CLAIM_MAP
    parity_dashboard_path = index_root / EVIDENCE_PARITY_DASHBOARD
    parity_summary_path = index_root / EVIDENCE_PARITY_SUMMARY
    mismatch_archive_path = index_root / EVIDENCE_MISMATCH_ARCHIVE
    mismatch_summary_path = index_root / EVIDENCE_MISMATCH_SUMMARY
    verdict_workflows_path = index_root / EVIDENCE_VERDICT_WORKFLOWS
    verdict_workflows_summary_path = index_root / EVIDENCE_VERDICT_WORKFLOWS_SUMMARY
    false_confidence_audit_path = index_root / EVIDENCE_FALSE_CONFIDENCE_AUDIT
    false_confidence_summary_path = index_root / EVIDENCE_FALSE_CONFIDENCE_SUMMARY
    scientific_debt_path = index_root / EVIDENCE_SCIENTIFIC_DEBT_REGISTER
    scientific_debt_summary_path = index_root / EVIDENCE_SCIENTIFIC_DEBT_SUMMARY
    portability_path = index_root / EVIDENCE_PORTABILITY_AUDIT
    portability_summary_path = index_root / EVIDENCE_PORTABILITY_SUMMARY
    fragile_examples_path = index_root / EVIDENCE_FRAGILE_EXAMPLE_AUDIT
    fragile_examples_summary_path = index_root / EVIDENCE_FRAGILE_EXAMPLE_SUMMARY
    regeneration_path = index_root / EVIDENCE_REGENERATION_CONTRACT
    regeneration_summary_path = index_root / EVIDENCE_REGENERATION_SUMMARY
    index_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    catalog_path.write_text(
        render_evidence_catalog(payload),
        encoding="utf-8",
    )
    claim_map_path.write_text(
        json.dumps(claim_map_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    parity_dashboard_path.write_text(
        json.dumps(parity_dashboard_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    parity_summary_path.write_text(
        render_evidence_parity_dashboard(parity_dashboard_payload),
        encoding="utf-8",
    )
    mismatch_archive_path.write_text(
        json.dumps(mismatch_archive_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    mismatch_summary_path.write_text(
        render_evidence_mismatch_archive(mismatch_archive_payload),
        encoding="utf-8",
    )
    verdict_workflows_path.write_text(
        json.dumps(verdict_workflows_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    verdict_workflows_summary_path.write_text(
        render_evidence_verdict_workflows(verdict_workflows_payload),
        encoding="utf-8",
    )
    false_confidence_audit_path.write_text(
        json.dumps(false_confidence_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    false_confidence_summary_path.write_text(
        render_evidence_false_confidence_audit(false_confidence_payload),
        encoding="utf-8",
    )
    scientific_debt_path.write_text(
        json.dumps(scientific_debt_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    scientific_debt_summary_path.write_text(
        render_evidence_scientific_debt_register(scientific_debt_payload),
        encoding="utf-8",
    )
    portability_payload = build_evidence_portability_audit(repo_root)
    portability_path.write_text(
        json.dumps(portability_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    portability_summary_path.write_text(
        render_evidence_portability_audit(portability_payload),
        encoding="utf-8",
    )
    fragile_examples_path.write_text(
        json.dumps(fragile_examples_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    fragile_examples_summary_path.write_text(
        render_evidence_fragile_example_audit(fragile_examples_payload),
        encoding="utf-8",
    )
    regeneration_path.write_text(
        json.dumps(regeneration_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    regeneration_summary_path.write_text(
        render_evidence_regeneration_contract(regeneration_payload),
        encoding="utf-8",
    )
    return index_path, catalog_path
