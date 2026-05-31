from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
from bijux_phylogenetics.trees import (
    compute_maximum_clade_credibility_tree,
    write_maximum_clade_credibility_artifacts,
    write_maximum_clade_credibility_score_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_compute_maximum_clade_credibility_tree_selects_known_best_candidate() -> None:
    tree, report = compute_maximum_clade_credibility_tree(
        fixture("maximum_clade_credibility_tree_set.nwk")
    )

    assert dumps_newick(tree) == "(((A:1,B:1):1,(E:1,F:1):1):1,(C:1,D:1):2);"
    assert report.tree_count == 5
    assert report.shared_taxa == ["A", "B", "C", "D", "E", "F"]
    assert report.rooted_topology_count == 4
    assert report.selected_tree_index == 2
    assert report.maximum_clade_credibility_newick == dumps_newick(tree)
    assert report.rows[0].score_rank == 1
    assert report.rows[0].source_tree_index == 2
    assert report.rows[0].raw_tree_count == 1
    assert report.rows[0].clade_credibility_score == pytest.approx(-3.259697819388456)
    assert [row.score_rank for row in report.rows] == [1, 2, 3, 4, 5]
    assert max(report.rows, key=lambda row: row.raw_tree_count).source_tree_index in {
        4,
        5,
    }
    assert max(row.raw_tree_count for row in report.rows) == 2
    assert report.rows[-1].source_tree_index == 3


def test_write_maximum_clade_credibility_score_table_writes_expected_columns(
    tmp_path: Path,
) -> None:
    _tree, report = compute_maximum_clade_credibility_tree(
        fixture("maximum_clade_credibility_tree_set.nwk")
    )

    output_path = tmp_path / "candidate-score-table.tsv"
    write_maximum_clade_credibility_score_table(output_path, report)

    lines = output_path.read_text(encoding="utf-8").splitlines()

    assert lines[0] == (
        "score_rank\tsource_tree_index\trooted_topology_id\traw_tree_count\t"
        "raw_tree_frequency\tclade_credibility_score\tcandidate_newick"
    )
    assert lines[1].startswith("1\t2\t")
    assert "\t1\t0.2\t-3.25969781938846\t" in lines[1]


def test_write_maximum_clade_credibility_artifacts_writes_tree_and_score_table(
    tmp_path: Path,
) -> None:
    tree, report = compute_maximum_clade_credibility_tree(
        fixture("maximum_clade_credibility_tree_set.nwk")
    )

    paths = write_maximum_clade_credibility_artifacts(tmp_path, tree, report)

    assert sorted(paths) == [
        "candidate_score_table_path",
        "maximum_clade_credibility_tree_path",
    ]
    assert (
        paths["maximum_clade_credibility_tree_path"].read_text(encoding="utf-8").strip()
        == "(((A:1,B:1):1,(E:1,F:1):1):1,(C:1,D:1):2);"
    )
    assert (
        paths["candidate_score_table_path"]
        .read_text(encoding="utf-8")
        .startswith("score_rank\tsource_tree_index\trooted_topology_id")
    )


def test_compute_maximum_clade_credibility_tree_requires_exact_taxa() -> None:
    with pytest.raises(InvalidAlignmentError, match="exact same taxon set"):
        compute_maximum_clade_credibility_tree(
            fixture("example_tree_set_mismatched.nwk")
        )
