from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import build_parsimony_stepwise_addition_tree
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_stepwise_addition_surface() -> None:
    assert (
        parsimony_api.build_parsimony_stepwise_addition_tree
        is build_parsimony_stepwise_addition_tree
    )


def test_build_parsimony_stepwise_addition_tree_scores_each_insertion_under_fitch() -> (
    None
):
    tree, report = build_parsimony_stepwise_addition_tree(
        fixture("nni_search_matrix.tsv"),
    )

    assert report.algorithm == "greedy-stepwise-addition-tree"
    assert report.objective_name == "parsimony-fitch"
    assert report.objective_direction == "minimize"
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert report.all_requested_taxa_present_once is True
    assert report.final_score == 2.0
    assert len(report.trace_rows) == 2
    assert report.trace_rows[0].taxon == "C"
    assert report.trace_rows[0].best_score == 2.0
    assert len(report.trace_rows[0].tested_edge_rows) == 3
    assert sorted({row.score for row in report.trace_rows[1].tested_edge_rows}) == [
        2.0,
        4.0,
    ]
    assert tree.tip_names == ["A", "B", "C", "D"]


def test_build_parsimony_stepwise_addition_tree_rejects_insertion_order_taxa_mismatch() -> (
    None
):
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        build_parsimony_stepwise_addition_tree(
            fixture("nni_search_matrix.tsv"),
            insertion_order=["A", "B", "C", "X"],
        )

    assert (
        error_info.value.code
        == "parsimony_stepwise_addition_insertion_order_taxa_mismatch"
    )
    assert error_info.value.details["missing_taxa"] == ["D"]
    assert error_info.value.details["unexpected_taxa"] == ["X"]
