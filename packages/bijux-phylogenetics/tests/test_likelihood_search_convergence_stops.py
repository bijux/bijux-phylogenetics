from __future__ import annotations

from pathlib import Path

import bijux_phylogenetics.phylo.likelihood.nni_search as nni_search_module
from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_nni_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_likelihood_nni_search_stops_when_best_improvement_is_within_tolerance() -> None:
    report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
        search_improvement_tolerance=20.0,
    )

    assert report.accepted_move_count == 0
    assert report.evaluated_neighbor_count == 4
    assert report.stopping_reason == "improvement-within-tolerance"
    assert report.trace_rows[-1].stopping_reason == "improvement-within-tolerance"
    assert report.final_tree_newick == report.start_tree_newick
    assert any(row.improving_move for row in report.candidate_rows)
    assert not any(row.selected_best_move for row in report.candidate_rows)


def test_likelihood_nni_search_stops_cleanly_on_mid_search_failure(
    monkeypatch,
) -> None:
    original = nni_search_module.reoptimize_nucleotide_topology_tree
    call_count = 0

    def fail_after_start(*args, **kwargs):
        nonlocal call_count
        call_count += 1
        if call_count == 1:
            return original(*args, **kwargs)
        raise RuntimeError("synthetic likelihood failure")

    monkeypatch.setattr(
        nni_search_module,
        "reoptimize_nucleotide_topology_tree",
        fail_after_start,
    )

    report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
    )

    assert report.accepted_move_count == 0
    assert report.stopping_reason == "search-failure"
    assert report.trace_rows[-1].stopping_reason == "search-failure"
    assert report.final_tree_newick == report.start_tree_newick
