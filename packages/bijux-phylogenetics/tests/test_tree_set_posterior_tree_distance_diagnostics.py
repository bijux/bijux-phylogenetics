from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
from bijux_phylogenetics.trees import (
    compute_posterior_tree_distance_diagnostics,
    write_posterior_tree_distance_artifacts,
    write_posterior_tree_distance_diagnostic_table,
    write_posterior_tree_distance_distribution_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_compute_posterior_tree_distance_diagnostics_ranks_outliers() -> None:
    report = compute_posterior_tree_distance_diagnostics(
        fixture("maximum_clade_credibility_tree_set.nwk")
    )
    rows = {row.source_tree_index: row for row in report.rows}

    assert report.tree_count == 5
    assert report.shared_taxa == ["A", "B", "C", "D", "E", "F"]
    assert report.maximum_clade_credibility_tree_index == 2
    assert (
        report.maximum_clade_credibility_newick
        == "(((A:1,B:1):1,(E:1,F:1):1):1,(C:1,D:1):2);"
    )
    assert report.consensus_method == "majority-rule"
    assert report.consensus_newick == "(A:1,B:1,(C:1,D:1)60:1.66666666666667,E:1,F:1);"
    assert report.row_count == 5
    assert report.distribution_row_count == 19

    assert rows[4].mcc_outlier_rank == 1
    assert rows[5].mcc_outlier_rank == 2
    assert rows[3].consensus_outlier_rank == 5

    assert rows[2].mcc_robinson_foulds_distance == 0
    assert rows[2].mcc_normalized_robinson_foulds == 0.0
    assert rows[2].mcc_branch_score_distance == 0.0

    assert rows[4].mcc_robinson_foulds_distance == 8
    assert rows[4].mcc_normalized_robinson_foulds == 1.0
    assert rows[4].mcc_branch_score_distance == pytest.approx(4.69041575982343)
    assert rows[4].consensus_robinson_foulds_distance == 5
    assert rows[4].consensus_normalized_robinson_foulds == 1.0

    assert rows[1].mcc_robinson_foulds_distance == 2
    assert rows[1].mcc_normalized_robinson_foulds == 0.25
    assert rows[1].mcc_branch_score_distance == pytest.approx(2.8284271247461903)


def test_write_posterior_tree_distance_tables_write_expected_rows(
    tmp_path: Path,
) -> None:
    report = compute_posterior_tree_distance_diagnostics(
        fixture("maximum_clade_credibility_tree_set.nwk")
    )

    diagnostic_path = tmp_path / "posterior-tree-distance-diagnostics.tsv"
    distribution_path = tmp_path / "posterior-tree-distance-distribution.tsv"
    write_posterior_tree_distance_diagnostic_table(diagnostic_path, report)
    write_posterior_tree_distance_distribution_table(distribution_path, report)

    assert diagnostic_path.read_text(encoding="utf-8").splitlines()[:3] == [
        "source_tree_index\trooted_topology_id\tmcc_outlier_rank\tconsensus_outlier_rank\tmcc_robinson_foulds_distance\tmcc_normalized_robinson_foulds\tmcc_branch_score_distance\tconsensus_robinson_foulds_distance\tconsensus_normalized_robinson_foulds\tconsensus_branch_score_distance",
        "1\tA|B||A|B|C|D||C|D||E|F\t4\t3\t2\t0.25\t2.82842712474619\t3\t0.6\t3.23178657161089",
        "2\tA|B||A|B|E|F||C|D||E|F\t5\t4\t0\t0\t0\t3\t0.6\t1.9436506316151",
    ]
    assert distribution_path.read_text(encoding="utf-8").splitlines()[:5] == [
        "reference_tree_kind\tdistance_metric\tobserved_value\ttree_count\tfrequency",
        "consensus\tbranch-score\t1.9436506316151\t2\t0.4",
        "consensus\tbranch-score\t3.23178657161089\t1\t0.2",
        "consensus\tbranch-score\t3.71184290855335\t2\t0.4",
        "consensus\tnormalized-robinson-foulds\t0.6\t3\t0.6",
    ]


def test_write_posterior_tree_distance_artifacts_writes_reference_trees_and_tables(
    tmp_path: Path,
) -> None:
    report = compute_posterior_tree_distance_diagnostics(
        fixture("maximum_clade_credibility_tree_set.nwk")
    )

    from bijux_phylogenetics.trees import (
        compute_consensus_tree,
        compute_maximum_clade_credibility_tree,
    )

    mcc_tree, _mcc_report = compute_maximum_clade_credibility_tree(
        fixture("maximum_clade_credibility_tree_set.nwk")
    )
    consensus_tree, _consensus_report = compute_consensus_tree(
        fixture("maximum_clade_credibility_tree_set.nwk")
    )
    paths = write_posterior_tree_distance_artifacts(
        tmp_path,
        mcc_tree,
        consensus_tree,
        report,
    )

    assert sorted(paths) == [
        "consensus_tree_path",
        "maximum_clade_credibility_tree_path",
        "posterior_tree_distance_diagnostic_table_path",
        "posterior_tree_distance_distribution_table_path",
    ]
    assert (
        paths["maximum_clade_credibility_tree_path"].read_text(encoding="utf-8").strip()
        == "(((A:1,B:1):1,(E:1,F:1):1):1,(C:1,D:1):2);"
    )
    assert paths["consensus_tree_path"].read_text(encoding="utf-8").strip() == (
        "(A:1,B:1,(C:1,D:1)60:1.66666666666667,E:1,F:1);"
    )


def test_compute_posterior_tree_distance_diagnostics_requires_exact_taxa() -> None:
    with pytest.raises(InvalidAlignmentError, match="exact same taxon set"):
        compute_posterior_tree_distance_diagnostics(
            fixture("example_tree_set_mismatched.nwk")
        )
