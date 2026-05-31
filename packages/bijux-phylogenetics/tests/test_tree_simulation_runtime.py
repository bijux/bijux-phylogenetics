from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.simulation import (
    simulate_coalescent_tree,
    simulate_coalescent_trees,
    simulate_random_tree,
    simulate_random_trees,
    write_coalescent_waiting_time_table,
    write_tree_simulation_envelope_table,
    write_tree_simulation_record_table,
)


def test_simulate_random_trees_reports_uniform_branch_length_envelope() -> None:
    trees, report = simulate_random_trees(tree_count=3, tip_count=5, seed=11)

    assert [tree.tip_count for tree in trees] == [5, 5, 5]
    assert report.model == "random-tree"
    assert report.tree_count == 3
    assert report.tip_count == 5
    assert report.seed == 11
    assert report.branch_length_model == "uniform"
    assert report.population_size is None
    assert report.rooted is True
    assert report.binary is True
    assert report.pooled_branch_count == 24
    assert len(report.records) == 3
    assert {row.metric for row in report.envelope_metrics} == {
        "tree_height_branch_length",
        "total_branch_length",
        "branch_length",
        "cherry_count",
        "sackin_imbalance_index",
        "normalized_colless_imbalance",
    }
    for tree in trees:
        for node in tree.iter_nodes():
            if node is tree.root:
                assert node.branch_length is None
            elif node.branch_length is not None:
                assert 0.0 <= node.branch_length <= 1.0


def test_simulate_coalescent_trees_reports_population_size_and_envelope() -> None:
    _trees, report = simulate_coalescent_trees(
        tree_count=2,
        tip_count=4,
        population_size=2.5,
        seed=13,
    )

    assert report.model == "coalescent"
    assert report.population_size == 2.5
    assert report.branch_length_model == "coalescent-waiting-times"
    assert report.pooled_branch_count == 12
    assert [row.observation_count for row in report.envelope_metrics[:2]] == [2, 2]
    assert report.coalescent_waiting_time_tolerance == 0.2
    assert [row.lineage_count for row in report.coalescent_waiting_time_rows] == [
        4,
        3,
        2,
    ]


def test_simulate_coalescent_trees_match_kingman_waiting_time_expectations() -> None:
    _trees, report = simulate_coalescent_trees(
        tree_count=256,
        tip_count=6,
        population_size=3.0,
        waiting_time_tolerance=0.15,
        seed=17,
    )

    assert report.coalescent_waiting_time_tolerance == 0.15
    assert [row.lineage_count for row in report.coalescent_waiting_time_rows] == [
        6,
        5,
        4,
        3,
        2,
    ]
    assert [row.coalescent_rate for row in report.coalescent_waiting_time_rows] == [
        5.0,
        3.333333333333333,
        2.0,
        1.0,
        0.333333333333333,
    ]
    assert [
        row.expected_waiting_time for row in report.coalescent_waiting_time_rows
    ] == [
        0.2,
        0.3,
        0.5,
        1.0,
        3.0,
    ]
    assert all(
        row.observation_count == 256 for row in report.coalescent_waiting_time_rows
    )
    assert all(row.within_tolerance for row in report.coalescent_waiting_time_rows)
    assert max(row.relative_error for row in report.coalescent_waiting_time_rows) < 0.15


def test_simulate_coalescent_trees_reject_negative_waiting_time_tolerance() -> None:
    with pytest.raises(
        ValueError,
        match="waiting_time_tolerance must be nonnegative",
    ):
        simulate_coalescent_trees(
            tree_count=4,
            tip_count=5,
            population_size=2.0,
            waiting_time_tolerance=-0.01,
        )


def test_simulate_random_tree_matches_single_tree_batch_surface() -> None:
    tree, report = simulate_random_tree(tip_count=5, seed=11)
    trees, batch_report = simulate_random_trees(tree_count=1, tip_count=5, seed=11)

    assert tree.to_newick() == trees[0].to_newick()
    assert report.records == batch_report.records
    assert report.envelope_metrics == batch_report.envelope_metrics


def test_simulate_coalescent_tree_matches_single_tree_batch_surface() -> None:
    tree, report = simulate_coalescent_tree(tip_count=4, population_size=2.5, seed=13)
    trees, batch_report = simulate_coalescent_trees(
        tree_count=1,
        tip_count=4,
        population_size=2.5,
        seed=13,
    )

    assert tree.to_newick() == trees[0].to_newick()
    assert report.records == batch_report.records
    assert report.envelope_metrics == batch_report.envelope_metrics
    assert report.coalescent_waiting_time_tolerance == (
        batch_report.coalescent_waiting_time_tolerance
    )
    assert (
        report.coalescent_waiting_time_rows == batch_report.coalescent_waiting_time_rows
    )


def test_write_tree_simulation_tables_emit_reviewable_ledgers(tmp_path: Path) -> None:
    _trees, report = simulate_random_trees(tree_count=2, tip_count=4, seed=7)

    record_path = write_tree_simulation_record_table(
        tmp_path / "simulation-records.tsv", report
    )
    envelope_path = write_tree_simulation_envelope_table(
        tmp_path / "simulation-envelope.tsv",
        report,
    )

    record_text = record_path.read_text(encoding="utf-8")
    envelope_text = envelope_path.read_text(encoding="utf-8")

    assert record_path.is_file()
    assert envelope_path.is_file()
    assert "tree_height_branch_length" in record_text
    assert "normalized_colless_imbalance" in record_text
    assert "branch_length\tedge\t12\t" in envelope_text
    assert "cherry_count\ttree\t2\t" in envelope_text


def test_write_coalescent_waiting_time_table_emits_lineage_ledgers(
    tmp_path: Path,
) -> None:
    _trees, report = simulate_coalescent_trees(
        tree_count=64,
        tip_count=5,
        population_size=2.5,
        waiting_time_tolerance=0.2,
        seed=19,
    )

    waiting_time_path = write_coalescent_waiting_time_table(
        tmp_path / "coalescent-waiting-times.tsv",
        report,
    )
    waiting_time_text = waiting_time_path.read_text(encoding="utf-8")
    waiting_time_rows = [
        line.split("\t") for line in waiting_time_text.strip().splitlines()[1:]
    ]

    assert waiting_time_path.is_file()
    assert waiting_time_text.startswith(
        "lineage_count\tcoalescent_rate\texpected_waiting_time\tobservation_count\t"
    )
    assert "\twithin_tolerance\twaiting_time_tolerance\n" in waiting_time_text
    assert len(waiting_time_rows) == 4
    assert waiting_time_rows[0][:4] == ["5", "4", "0.25", "64"]
    assert all(row[-1] == "0.2" for row in waiting_time_rows)
