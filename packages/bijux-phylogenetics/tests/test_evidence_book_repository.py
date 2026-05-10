from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.evidence.book import (
    build_evidence_claim_map,
    build_evidence_book_index,
    render_evidence_catalog,
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

    payload = build_evidence_book_index(REPO_ROOT)
    catalog = render_evidence_catalog(payload)
    claim_map = build_evidence_claim_map(REPO_ROOT)

    assert json.loads(index_path.read_text(encoding="utf-8")) == payload
    assert catalog_path.read_text(encoding="utf-8") == catalog
    assert json.loads(claim_map_path.read_text(encoding="utf-8")) == claim_map
