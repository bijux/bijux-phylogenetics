from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.simulation import (
    simulate_coalescent_trees,
    simulate_random_trees,
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


def test_write_tree_simulation_tables_emit_reviewable_ledgers(tmp_path: Path) -> None:
    _trees, report = simulate_random_trees(tree_count=2, tip_count=4, seed=7)

    record_path = write_tree_simulation_record_table(tmp_path / "simulation-records.tsv", report)
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
