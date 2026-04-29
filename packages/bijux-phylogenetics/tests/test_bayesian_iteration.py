from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.bayesian.posterior import (
    compare_bayesian_tree_sets,
    summarize_posterior_node_ages,
    summarize_maximum_clade_credibility_tree,
    thin_posterior_tree_set,
)


FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def test_summarize_maximum_clade_credibility_tree_selects_best_supported_posterior_tree() -> None:
    tree, report = summarize_maximum_clade_credibility_tree(
        fixture("trees/example_tree_set_left.nwk"),
        burnin_fraction=0.0,
    )

    assert tree.tip_count == 4
    assert report.total_tree_count == 3
    assert report.selected_tree_index == 1
    assert report.rooted_topology_count == 2
    assert report.filtered_tree_set_path.parent != fixture("trees").resolve()
    assert "A:0.1" in report.mcc_newick


def test_thin_posterior_tree_set_retains_every_nth_tree_after_burnin(tmp_path: Path) -> None:
    output_path = tmp_path / "thinned.nwk"
    report = thin_posterior_tree_set(
        fixture("trees/example_tree_set_left.nwk"),
        output_path,
        thinning_interval=2,
        burnin_fraction=0.0,
    )

    lines = [line for line in output_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    assert report.retained_tree_count == 2
    assert report.retained_indices == [1, 3]
    assert len(lines) == 2


def test_compare_bayesian_tree_sets_reports_mcc_and_posterior_differences() -> None:
    report = compare_bayesian_tree_sets(
        fixture("trees/example_tree_set_left.nwk"),
        fixture("trees/example_tree_set_right.nwk"),
        burnin_fraction=0.0,
    )

    assert report.left_mcc.selected_tree_index == 1
    assert report.right_mcc.selected_tree_index == 2
    assert report.tree_set_comparison.shared_rooted_topology_count == 1
    assert report.mcc_topology.topology_equal is False


def test_summarize_posterior_node_ages_reports_clade_heights() -> None:
    report = summarize_posterior_node_ages(
        fixture("trees/example_tree_set_left.nwk"),
        burnin_fraction=0.0,
    )

    rows = {row.clade: row for row in report.rows}
    assert report.kept_tree_count == 3
    assert rows["A|B"].mean_height == 0.2
    assert rows["C|D"].maximum_height == 0.2
