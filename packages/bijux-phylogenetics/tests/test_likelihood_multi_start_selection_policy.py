from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_multi_start_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_likelihood_multi_start_best_selection_policy_searches_best_scored_pool_start() -> (
    None
):
    report = search_nucleotide_likelihood_multi_start_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        start_tree_count=4,
        starting_tree_selection_policy="best",
    )

    assert report.start_tree_source_policy == "scored-starting-tree-pool"
    assert report.starting_tree_selection_policy == "best"
    assert report.available_start_tree_count == 4
    assert report.start_tree_count == 1
    assert report.best_run_source_label == "likelihood-stepwise-addition-tree"
    assert [row.start_tree_source_label for row in report.run_summaries] == [
        "likelihood-stepwise-addition-tree"
    ]


def test_likelihood_multi_start_strategy_priority_policy_searches_one_start_per_strategy() -> (
    None
):
    report = search_nucleotide_likelihood_multi_start_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        start_tree_count=4,
        starting_tree_selection_policy="strategy-priority",
    )

    assert report.starting_tree_selection_policy == "strategy-priority"
    assert report.available_start_tree_count == 4
    assert report.start_tree_count == 3
    assert [row.start_tree_source_kind for row in report.run_summaries] == [
        "likelihood-stepwise-addition-tree",
        "input-tree",
        "random-tree",
    ]


def test_likelihood_multi_start_random_k_selection_policy_is_deterministic() -> None:
    first_report = search_nucleotide_likelihood_multi_start_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        start_tree_count=4,
        starting_tree_selection_policy="random-k",
        selected_start_tree_count=2,
        starting_tree_selection_seed=17,
    )
    second_report = search_nucleotide_likelihood_multi_start_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        start_tree_count=4,
        starting_tree_selection_policy="random-k",
        selected_start_tree_count=2,
        starting_tree_selection_seed=17,
    )

    assert first_report.starting_tree_selection_policy == "random-k"
    assert first_report.available_start_tree_count == 4
    assert first_report.start_tree_count == 2
    assert [row.start_tree_source_label for row in first_report.run_summaries] == [
        row.start_tree_source_label for row in second_report.run_summaries
    ]
