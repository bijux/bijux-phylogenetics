from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.evidence.book import (
    build_evidence_book_index,
    build_evidence_claim_map,
    build_evidence_parity_dashboard,
)


REPO_ROOT = Path(__file__).resolve().parents[3]


def test_repository_index_outputs_match_generated_payloads() -> None:
    index_root = REPO_ROOT / "evidence-book" / "index"

    assert json.loads((index_root / "evidence-index.json").read_text(encoding="utf-8")) == build_evidence_book_index(REPO_ROOT)
    assert json.loads((index_root / "claim-map.json").read_text(encoding="utf-8")) == build_evidence_claim_map(REPO_ROOT)
    assert json.loads((index_root / "parity-dashboard.json").read_text(encoding="utf-8")) == build_evidence_parity_dashboard(REPO_ROOT)
