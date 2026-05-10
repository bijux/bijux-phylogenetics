from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.evidence.studies.primate_pgls_and_signal import (
    build_primate_pgls_signal_bundles,
    build_primate_pgls_signal_claim_registry,
    build_primate_pgls_signal_external_sources,
    build_primate_pgls_signal_family_index,
    build_primate_pgls_signal_parity_policy,
    build_primate_pgls_signal_scalar_parity_table,
    build_primate_pgls_signal_source_fragment_map,
    render_primate_pgls_signal_scalar_parity_table_markdown,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
STUDY_ROOT = REPO_ROOT / "evidence-book" / "studies" / "primate-pgls-and-signal"


def test_primate_pgls_and_signal_study_indexes_match_generated_payloads() -> None:
    external_sources = json.loads(
        (STUDY_ROOT / "provenance" / "lund-course-sources.json").read_text(
            encoding="utf-8"
        )
    )
    source_fragment_map = json.loads(
        (STUDY_ROOT / "source-fragment-map.json").read_text(encoding="utf-8")
    )
    family_index = json.loads(
        (STUDY_ROOT / "family-index.json").read_text(encoding="utf-8")
    )
    parity_policy = json.loads(
        (STUDY_ROOT / "parity-policy.json").read_text(encoding="utf-8")
    )
    claim_registry = json.loads(
        (STUDY_ROOT / "claim-registry.json").read_text(encoding="utf-8")
    )

    assert external_sources == build_primate_pgls_signal_external_sources()
    assert source_fragment_map == build_primate_pgls_signal_source_fragment_map()
    assert family_index == build_primate_pgls_signal_family_index(REPO_ROOT)
    assert parity_policy == build_primate_pgls_signal_parity_policy()
    assert claim_registry == build_primate_pgls_signal_claim_registry(REPO_ROOT)


def test_primate_pgls_and_signal_bundles_match_generated_payloads() -> None:
    generated_bundles = build_primate_pgls_signal_bundles(REPO_ROOT)

    for evidence_id, generated in generated_bundles.items():
        bundle_root = STUDY_ROOT / evidence_id
        assert json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8")) == generated["manifest"]
        assert json.loads((bundle_root / "claims.json").read_text(encoding="utf-8")) == generated["claims"]
        assert json.loads(
            (bundle_root / generated["report_filename"]).read_text(encoding="utf-8")
        ) == generated["report_payload"]
        assert (bundle_root / "README.md").read_text(encoding="utf-8") == generated["readme"]


def test_primate_pgls_and_signal_scalar_table_matches_generated_payload() -> None:
    table = json.loads(
        (
            STUDY_ROOT
            / "evidence-001"
            / "scalar-parity-table.json"
        ).read_text(encoding="utf-8")
    )
    markdown = (
        STUDY_ROOT / "evidence-001" / "scalar-parity-table.md"
    ).read_text(encoding="utf-8")

    generated = build_primate_pgls_signal_scalar_parity_table(REPO_ROOT)

    assert table == generated
    assert markdown == render_primate_pgls_signal_scalar_parity_table_markdown(
        generated
    )
