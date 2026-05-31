from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    ParsimonyConsensusSummary,
    ParsimonyEqualBestConsensusReport,
    ParsimonyEqualBestTree,
    summarize_equal_best_parsimony_trees,
    write_parsimony_equal_best_consensus_artifacts,
)
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

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
    assert (
        parsimony_api.write_parsimony_equal_best_consensus_artifacts
        is write_parsimony_equal_best_consensus_artifacts
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


def test_write_equal_best_parsimony_consensus_artifacts_materializes_outputs(
    tmp_path: Path,
) -> None:
    report = summarize_equal_best_parsimony_trees(
        fixture("bootstrap_matrix.tsv"),
        method="fitch",
    )

    outputs = write_parsimony_equal_best_consensus_artifacts(
        tmp_path / "equal-best-consensus",
        report,
    )

    assert set(outputs) == {
        "equal_best_trees_path",
        "equal_best_scores_path",
        "strict_consensus_tree_path",
        "majority_consensus_tree_path",
        "clade_frequencies_path",
        "run_json_path",
    }
    assert (
        outputs["equal_best_scores_path"]
        .read_text(encoding="utf-8")
        .startswith("tree_index\ttotal_score\ttree_newick\n")
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "parsimony-equal-best-consensus"
    assert payload["equal_best_tree_count"] == 5
    assert payload["retained_all_equal_best_trees"] is True
    assert payload["strict_consensus"]["consensus_newick"] == "(A,B,C,D);"
    assert payload["majority_consensus"]["consensus_newick"] == "((A,B)60,(C,D)60);"


def test_equal_best_parsimony_consensus_suppresses_consensus_when_cap_truncates() -> (
    None
):
    report = summarize_equal_best_parsimony_trees(
        fixture("bootstrap_matrix.tsv"),
        method="fitch",
        max_retained_equal_best_trees=3,
    )

    assert report.equal_best_tree_count == 5
    assert report.retained_equal_best_tree_count == 3
    assert report.retained_all_equal_best_trees is False
    assert report.strict_consensus is None
    assert report.majority_consensus is None
    assert [row.tree_newick for row in report.equal_best_tree_rows] == [
        "(((A,B),C),D);",
        "(((A,B),D),C);",
        "((A,(C,D)),B);",
    ]


def test_write_equal_best_parsimony_consensus_artifacts_omits_partial_consensus_files(
    tmp_path: Path,
) -> None:
    report = summarize_equal_best_parsimony_trees(
        fixture("bootstrap_matrix.tsv"),
        method="fitch",
        max_retained_equal_best_trees=3,
    )

    outputs = write_parsimony_equal_best_consensus_artifacts(
        tmp_path / "equal-best-consensus-truncated",
        report,
    )

    assert set(outputs) == {
        "equal_best_trees_path",
        "equal_best_scores_path",
        "run_json_path",
    }
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["retained_all_equal_best_trees"] is False
    assert payload["strict_consensus"] is None
    assert payload["majority_consensus"] is None


def test_equal_best_parsimony_consensus_rejects_nonpositive_tree_cap() -> None:
    with pytest.raises(
        ParsimonyAnalysisError,
        match="retained-tree cap of at least one",
    ):
        summarize_equal_best_parsimony_trees(
            fixture("bootstrap_matrix.tsv"),
            method="fitch",
            max_retained_equal_best_trees=0,
        )
