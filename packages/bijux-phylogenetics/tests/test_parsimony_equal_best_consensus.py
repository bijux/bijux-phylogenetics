from __future__ import annotations

from pathlib import Path
import json

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    ParsimonyConsensusSummary,
    ParsimonyEqualBestConsensusReport,
    ParsimonyEqualBestTree,
    summarize_equal_best_parsimony_trees,
    write_parsimony_equal_best_consensus_artifacts,
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
    assert outputs["equal_best_scores_path"].read_text(encoding="utf-8").startswith(
        "tree_index\ttotal_score\ttree_newick\n"
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "parsimony-equal-best-consensus"
    assert payload["equal_best_tree_count"] == 5
    assert payload["retained_all_equal_best_trees"] is True
    assert payload["strict_consensus"]["consensus_newick"] == "(A,B,C,D);"
    assert payload["majority_consensus"]["consensus_newick"] == "((A,B)60,(C,D)60);"
