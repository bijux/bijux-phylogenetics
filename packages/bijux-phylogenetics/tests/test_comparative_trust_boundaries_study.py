from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.evidence.studies.comparative_trust_boundaries import (
    build_comparative_trust_boundaries_bundles,
    build_comparative_trust_boundaries_claim_registry,
    build_comparative_trust_boundaries_family_index,
    build_comparative_trust_boundaries_provenance,
    build_comparative_trust_boundaries_source_fragment_map,
    render_comparative_trust_boundaries_study_manifest,
    render_comparative_trust_boundaries_study_readme,
    render_weak_signal_traits_tsv,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
STUDY_ROOT = REPO_ROOT / "evidence-book" / "studies" / "comparative-trust-boundaries"


def test_comparative_trust_boundaries_study_indexes_match_generated_payloads() -> None:
    study_manifest = json.loads((STUDY_ROOT / "study.json").read_text(encoding="utf-8"))
    readme = (STUDY_ROOT / "README.md").read_text(encoding="utf-8")
    provenance = json.loads(
        (STUDY_ROOT / "provenance" / "runtime-sources.json").read_text(encoding="utf-8")
    )
    fragment_map = json.loads(
        (STUDY_ROOT / "source-fragment-map.json").read_text(encoding="utf-8")
    )
    family_index = json.loads(
        (STUDY_ROOT / "family-index.json").read_text(encoding="utf-8")
    )
    claim_registry = json.loads(
        (STUDY_ROOT / "claim-registry.json").read_text(encoding="utf-8")
    )

    assert study_manifest == render_comparative_trust_boundaries_study_manifest()
    assert readme == render_comparative_trust_boundaries_study_readme()
    assert provenance == build_comparative_trust_boundaries_provenance()
    assert fragment_map == build_comparative_trust_boundaries_source_fragment_map()
    assert family_index == build_comparative_trust_boundaries_family_index(REPO_ROOT)
    assert claim_registry == build_comparative_trust_boundaries_claim_registry(REPO_ROOT)


def test_comparative_trust_boundaries_bundles_match_generated_payloads() -> None:
    generated_bundles = build_comparative_trust_boundaries_bundles(REPO_ROOT)

    for evidence_id, generated in generated_bundles.items():
        bundle_root = STUDY_ROOT / evidence_id
        assert json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8")) == generated["manifest"]
        assert json.loads((bundle_root / "claims.json").read_text(encoding="utf-8")) == generated["claims"]
        assert json.loads(
            (bundle_root / generated["report_filename"]).read_text(encoding="utf-8")
        ) == generated["report_payload"]
        assert (bundle_root / "README.md").read_text(encoding="utf-8") == generated["readme"]


def test_comparative_trust_boundaries_weak_signal_case_stays_threshold_sensitive() -> None:
    traits_table = (
        STUDY_ROOT / "evidence-002" / "weak_signal_traits.tsv"
    ).read_text(encoding="utf-8")
    instability_audit = json.loads(
        (
            STUDY_ROOT / "evidence-002" / "result-instability-audit.json"
        ).read_text(encoding="utf-8")
    )

    assert traits_table == render_weak_signal_traits_tsv()
    assert instability_audit["crosses_alpha_boundary"] is True
    assert instability_audit["below_threshold_count"] == 7
    assert instability_audit["at_or_above_threshold_count"] == 13
