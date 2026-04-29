from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines.validation import audit_alignment_inference_readiness


FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_audit_alignment_inference_readiness_prefers_ml_for_aligned_variable_data() -> None:
    report = audit_alignment_inference_readiness(fixture("example_alignment.fasta"))
    assert report.overall_decision == "ready"
    assert report.recommended_workflow == "maximum_likelihood"
    assert any(decision.workflow == "bayesian" and decision.ready for decision in report.decisions)


def test_audit_alignment_inference_readiness_blocks_unaligned_raw_sequences() -> None:
    report = audit_alignment_inference_readiness(fixture("example_sequences_raw.fasta"))
    assert report.overall_decision == "blocked"
    assert report.recommended_workflow == "unsuitable"
    assert any("not yet aligned" in blocker for decision in report.decisions for blocker in decision.blockers)
