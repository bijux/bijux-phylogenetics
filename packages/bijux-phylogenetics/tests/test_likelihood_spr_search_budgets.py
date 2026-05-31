from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from bijux_phylogenetics.io.newick import dumps_newick, loads_newick
from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_spr_from_alignment,
    validate_nucleotide_likelihood_spr_search_budget,
)
from bijux_phylogenetics.phylo.topology import (
    apply_rooted_spr_move,
    iter_rooted_spr_move_candidates,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_validate_nucleotide_likelihood_spr_search_budget_rejects_conflicting_candidate_aliases() -> (
    None
):
    try:
        validate_nucleotide_likelihood_spr_search_budget(
            evaluation_budget=2,
            max_candidate_count=3,
        )
    except ValueError as error:
        assert str(error) == (
            "evaluation_budget must match max_candidate_count when both are provided"
        )
    else:
        raise AssertionError("conflicting candidate budget aliases must fail")


def test_validate_nucleotide_likelihood_spr_search_budget_rejects_nonpositive_elapsed_time() -> (
    None
):
    try:
        validate_nucleotide_likelihood_spr_search_budget(max_elapsed_seconds=0.0)
    except ValueError as error:
        assert str(error) == (
            "max_elapsed_seconds must be greater than zero when provided"
        )
    else:
        raise AssertionError("nonpositive elapsed-time budgets must fail")


def test_likelihood_spr_search_stops_for_iteration_budget_and_reports_remaining_neighbors() -> (
    None
):
    report = search_nucleotide_likelihood_spr_from_alignment(
        fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_spr_alignment_5_taxa.fasta"),
        model_name="jc69",
        max_iteration_count=1,
        upper_branch_length_bound=1.0,
    )

    assert report.accepted_move_count == 1
    assert report.iteration_count == 1
    assert report.stopping_reason == "iteration-budget-exhausted"
    assert report.search_budget.max_iteration_count == 1
    assert report.unsearched_candidate_count == count_unique_spr_neighbors(
        report.final_tree_newick
    )
    assert report.trace_rows[-1].stopping_reason == "iteration-budget-exhausted"
    assert (
        report.trace_rows[-1].unsearched_candidate_count
        == report.unsearched_candidate_count
    )


def test_likelihood_spr_search_stops_for_accepted_move_budget() -> None:
    report = search_nucleotide_likelihood_spr_from_alignment(
        fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_spr_alignment_5_taxa.fasta"),
        model_name="jc69",
        max_accepted_move_count=1,
        upper_branch_length_bound=1.0,
    )

    assert report.accepted_move_count == 1
    assert report.iteration_count == 1
    assert report.stopping_reason == "accepted-move-budget-exhausted"
    assert report.search_budget.max_accepted_move_count == 1
    assert report.unsearched_candidate_count == count_unique_spr_neighbors(
        report.final_tree_newick
    )


def test_likelihood_spr_search_stops_for_time_budget_before_second_candidate() -> None:
    clock_values = iter([0.0, 0.0, 0.0, 0.2])

    with patch(
        "bijux_phylogenetics.phylo.likelihood.spr_search.time.monotonic",
        side_effect=lambda: next(clock_values),
    ):
        report = search_nucleotide_likelihood_spr_from_alignment(
            fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk"),
            fixture("alignments", "jc69_likelihood_spr_alignment_5_taxa.fasta"),
            model_name="jc69",
            max_elapsed_seconds=0.1,
            upper_branch_length_bound=1.0,
        )

    assert report.accepted_move_count == 0
    assert report.iteration_count == 1
    assert report.evaluated_neighbor_count == 1
    assert report.stopping_reason == "time-budget-exhausted"
    assert report.search_budget.max_elapsed_seconds == 0.1
    assert report.unsearched_candidate_count > 0
    assert (
        report.trace_rows[-1].unsearched_candidate_count
        == report.unsearched_candidate_count
    )


def count_unique_spr_neighbors(tree_newick: str) -> int:
    tree = loads_newick(tree_newick)
    seen_neighbor_newicks: set[str] = set()
    for candidate in iter_rooted_spr_move_candidates(tree):
        seen_neighbor_newicks.add(dumps_neighbor_newick(tree, candidate))
    return len(seen_neighbor_newicks)


def dumps_neighbor_newick(tree, candidate) -> str:
    return dumps_newick(apply_rooted_spr_move(tree, candidate))
