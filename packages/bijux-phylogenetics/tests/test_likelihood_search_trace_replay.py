from __future__ import annotations

from copy import deepcopy
from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_nni_from_alignment,
)
from bijux_phylogenetics.phylo.likelihood.search_replay import (
    replay_nucleotide_likelihood_nni_search_trace,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_likelihood_nni_search_trace_replay_restores_stored_final_topology() -> None:
    report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
    )

    replay_report = replay_nucleotide_likelihood_nni_search_trace(report)

    assert replay_report.source_search_algorithm == report.algorithm
    assert replay_report.accepted_trace_event_count == report.accepted_move_count == 2
    assert replay_report.replayed_step_count == 2
    assert replay_report.replay_failed is False
    assert replay_report.failure_reason is None
    assert replay_report.final_topology_matches is True
    assert (
        replay_report.replayed_final_topology_fingerprint
        == replay_report.stored_final_topology_fingerprint
    )
    assert all(step.step_replayed for step in replay_report.step_rows)

def test_likelihood_search_trace_replay_rejects_tampered_accepted_move_hash() -> None:
    report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
    )
    tampered_report = deepcopy(report)
    tampered_report.trace_rows[1].candidate_topology_fingerprint = "not-a-real-topology"

    replay_report = replay_nucleotide_likelihood_nni_search_trace(tampered_report)

    assert replay_report.accepted_trace_event_count == 1
    assert replay_report.replayed_step_count == 0
    assert replay_report.replay_failed is True
    assert replay_report.final_topology_matches is False
    assert (
        replay_report.failure_reason
        == "accepted trace row did not resolve to exactly one rooted NNI move"
    )
    assert replay_report.step_rows[0].matched_candidate_count == 0
    assert replay_report.step_rows[0].step_replayed is False
