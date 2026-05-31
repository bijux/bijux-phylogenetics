from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    build_nucleotide_likelihood_starting_tree_pool_from_alignment,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_selection import (
    select_nucleotide_likelihood_starting_tree_pool,
    validate_nucleotide_likelihood_starting_tree_selection_count,
    validate_nucleotide_likelihood_starting_tree_selection_policy,
    validate_nucleotide_likelihood_starting_tree_strategy_priority,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_likelihood_starting_tree_selection_policies_pick_distinct_subsets() -> None:
    report = build_nucleotide_likelihood_starting_tree_pool_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        random_start_tree_count=2,
        random_start_tree_seed=17,
    )

    all_rows = select_nucleotide_likelihood_starting_tree_pool(
        report,
        starting_tree_selection_policy="all",
    )
    best_rows = select_nucleotide_likelihood_starting_tree_pool(
        report,
        starting_tree_selection_policy="best",
    )
    random_rows = select_nucleotide_likelihood_starting_tree_pool(
        report,
        starting_tree_selection_policy="random-k",
        selected_start_tree_count=2,
        selection_seed=17,
    )
    repeated_random_rows = select_nucleotide_likelihood_starting_tree_pool(
        report,
        starting_tree_selection_policy="random-k",
        selected_start_tree_count=2,
        selection_seed=17,
    )
    strategy_rows = select_nucleotide_likelihood_starting_tree_pool(
        report,
        starting_tree_selection_policy="strategy-priority",
    )

    assert [row.tree_id for row in all_rows] == [
        "input-tree",
        "likelihood-stepwise-addition-tree",
        "random-tree-seed-17",
        "random-tree-seed-18",
    ]
    assert [row.tree_id for row in best_rows] == [
        "likelihood-stepwise-addition-tree"
    ]
    assert len(random_rows) == 2
    assert [row.tree_id for row in random_rows] == [
        row.tree_id for row in repeated_random_rows
    ]
    assert [row.source_strategy for row in strategy_rows] == [
        "likelihood-stepwise-addition-tree",
        "input-tree",
        "random-tree",
    ]
    assert len({tuple(row.tree_id for row in rows) for rows in [
        all_rows,
        best_rows,
        random_rows,
        strategy_rows,
    ]}) == 4


def test_likelihood_starting_tree_selection_validators_reject_invalid_requests() -> None:
    assert validate_nucleotide_likelihood_starting_tree_selection_policy("BEST") == "best"
    assert validate_nucleotide_likelihood_starting_tree_selection_count("all", None) is None
    assert validate_nucleotide_likelihood_starting_tree_strategy_priority(
        ["random-tree", "input-tree"]
    ) == ("random-tree", "input-tree")

    try:
        validate_nucleotide_likelihood_starting_tree_selection_policy("first")
    except ValueError as error:
        assert "starting_tree_selection_policy must be one of" in str(error)
    else:
        raise AssertionError("unsupported starting-tree selection policies must fail")

    try:
        validate_nucleotide_likelihood_starting_tree_selection_count("random-k", None)
    except ValueError as error:
        assert "selected_start_tree_count is required" in str(error)
    else:
        raise AssertionError("random-k selection must require one explicit count")

    try:
        validate_nucleotide_likelihood_starting_tree_strategy_priority(
            ["random-tree", "random-tree"]
        )
    except ValueError as error:
        assert str(error) == "strategy_priority must not repeat strategies"
    else:
        raise AssertionError("duplicate strategy-priority declarations must fail")
