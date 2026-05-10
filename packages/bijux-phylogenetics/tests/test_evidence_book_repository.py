from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.evidence.book import (
    build_evidence_claim_map,
    build_evidence_book_index,
    build_evidence_false_confidence_audit,
    build_evidence_mismatch_archive,
    build_evidence_parity_dashboard,
    build_evidence_verdict_workflows,
    render_evidence_catalog,
    render_evidence_false_confidence_audit,
    render_evidence_mismatch_archive,
    render_evidence_parity_dashboard,
    render_evidence_verdict_workflows,
    validate_evidence_book,
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
    parity_dashboard_path = REPO_ROOT / "evidence-book" / "index" / "parity-dashboard.json"
    parity_summary_path = REPO_ROOT / "evidence-book" / "index" / "parity-dashboard.md"
    mismatch_archive_path = REPO_ROOT / "evidence-book" / "index" / "mismatch-archive.json"
    mismatch_summary_path = REPO_ROOT / "evidence-book" / "index" / "mismatch-archive.md"
    verdict_workflows_path = REPO_ROOT / "evidence-book" / "index" / "verdict-workflows.json"
    verdict_workflows_summary_path = REPO_ROOT / "evidence-book" / "index" / "verdict-workflows.md"
    false_confidence_audit_path = REPO_ROOT / "evidence-book" / "index" / "false-confidence-audit.json"
    false_confidence_summary_path = REPO_ROOT / "evidence-book" / "index" / "false-confidence-audit.md"

    payload = build_evidence_book_index(REPO_ROOT)
    catalog = render_evidence_catalog(payload)
    claim_map = build_evidence_claim_map(REPO_ROOT)
    parity_dashboard = build_evidence_parity_dashboard(REPO_ROOT)
    mismatch_archive = build_evidence_mismatch_archive(REPO_ROOT)
    verdict_workflows = build_evidence_verdict_workflows(REPO_ROOT)
    false_confidence_audit = build_evidence_false_confidence_audit(REPO_ROOT)

    assert json.loads(index_path.read_text(encoding="utf-8")) == payload
    assert catalog_path.read_text(encoding="utf-8") == catalog
    assert json.loads(claim_map_path.read_text(encoding="utf-8")) == claim_map
    assert json.loads(parity_dashboard_path.read_text(encoding="utf-8")) == parity_dashboard
    assert parity_summary_path.read_text(encoding="utf-8") == render_evidence_parity_dashboard(
        parity_dashboard
    )
    assert json.loads(mismatch_archive_path.read_text(encoding="utf-8")) == mismatch_archive
    assert mismatch_summary_path.read_text(encoding="utf-8") == render_evidence_mismatch_archive(
        mismatch_archive
    )
    assert json.loads(verdict_workflows_path.read_text(encoding="utf-8")) == verdict_workflows
    assert verdict_workflows_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_verdict_workflows(verdict_workflows)
    assert json.loads(
        false_confidence_audit_path.read_text(encoding="utf-8")
    ) == false_confidence_audit
    assert false_confidence_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_false_confidence_audit(false_confidence_audit)
