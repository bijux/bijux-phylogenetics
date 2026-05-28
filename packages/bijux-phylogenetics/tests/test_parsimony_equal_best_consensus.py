from __future__ import annotations

from pathlib import Path

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    ParsimonyConsensusSummary,
    ParsimonyEqualBestConsensusReport,
    ParsimonyEqualBestTree,
    summarize_equal_best_parsimony_trees,
)

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_equal_best_consensus_contracts() -> None:
    assert parsimony_api.ParsimonyEqualBestTree is ParsimonyEqualBestTree
    assert (
        parsimony_api.ParsimonyEqualBestConsensusReport
        is ParsimonyEqualBestConsensusReport
    )
    assert parsimony_api.ParsimonyConsensusSummary is ParsimonyConsensusSummary
    assert (
        parsimony_api.summarize_equal_best_parsimony_trees
        is summarize_equal_best_parsimony_trees
    )


def test_equal_best_parsimony_consensus_reports_exact_best_tree_set() -> None:
    report = summarize_equal_best_parsimony_trees(
        fixture("bootstrap_matrix.tsv"),
        method="fitch",
    )

    assert report.algorithm == "parsimony-equal-best-consensus"
    assert report.method == "fitch"
    assert report.candidate_tree_count == 15
    assert report.best_score == 5.0
    assert report.equal_best_tree_count == 5
    assert report.retained_equal_best_tree_count == 5
    assert report.retained_all_equal_best_trees is True
    assert [row.tree_newick for row in report.equal_best_tree_rows] == [
        "(((A,B),C),D);",
        "(((A,B),D),C);",
        "((A,(C,D)),B);",
        "((A,B),(C,D));",
        "(A,(B,(C,D)));",
    ]
    assert report.strict_consensus is not None
    assert report.majority_consensus is not None
    assert report.strict_consensus.consensus_newick == "(A,B,C,D);"
    assert report.strict_consensus.included_clade_count == 0
    assert report.majority_consensus.consensus_newick == "((A,B)60,(C,D)60);"
    assert report.majority_consensus.included_clade_count == 2
