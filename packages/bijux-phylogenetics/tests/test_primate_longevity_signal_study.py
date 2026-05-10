from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.evidence.studies.primate_longevity_signal import (
    build_primate_claim_registry,
    build_primate_family_index,
    build_primate_parity_policy,
    build_primate_scalar_parity_table,
    build_primate_source_fragment_map,
    render_primate_scalar_parity_table_markdown,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
STUDY_ROOT = (
    REPO_ROOT / "evidence-book" / "studies" / "primate-longevity-signal"
)
BUNDLE_ROOT = STUDY_ROOT / "evidence-001"


def test_primate_study_fragment_outputs_match_generated_payloads() -> None:
    source_fragment_map = json.loads(
        (STUDY_ROOT / "source-fragment-map.json").read_text(encoding="utf-8")
    )
    family_index = json.loads(
        (STUDY_ROOT / "family-index.json").read_text(encoding="utf-8")
    )
    claim_registry = json.loads(
        (BUNDLE_ROOT / "claims.json").read_text(encoding="utf-8")
    )

    assert source_fragment_map == build_primate_source_fragment_map(REPO_ROOT)
    assert family_index == build_primate_family_index(REPO_ROOT)
    assert claim_registry == build_primate_claim_registry(REPO_ROOT)


def test_primate_study_parity_outputs_match_generated_payloads() -> None:
    parity_policy = json.loads(
        (STUDY_ROOT / "parity-policy.json").read_text(encoding="utf-8")
    )
    scalar_parity_table = json.loads(
        (BUNDLE_ROOT / "scalar-parity-table.json").read_text(encoding="utf-8")
    )
    scalar_parity_markdown = (BUNDLE_ROOT / "scalar-parity-table.md").read_text(
        encoding="utf-8"
    )

    generated_table = build_primate_scalar_parity_table(REPO_ROOT)
    assert parity_policy == build_primate_parity_policy(REPO_ROOT)
    assert scalar_parity_table == generated_table
    assert scalar_parity_markdown == render_primate_scalar_parity_table_markdown(
        generated_table
    )
