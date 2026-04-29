from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.bayesian.comparison import (
    compare_independent_bayesian_runs,
    compare_posterior_tree_sets_by_clock,
    compare_posterior_tree_sets_by_prior,
)
from bijux_phylogenetics.bayesian.posterior import (
    compare_bayesian_tree_sets,
    summarize_posterior_node_ages,
    summarize_maximum_clade_credibility_tree,
    thin_posterior_tree_set,
)


FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def _write_mrbayes_trace(path: Path, rows: list[tuple[int, float, float, float]]) -> Path:
    path.write_text(
        "Gen\tLnL\tTL\talpha\n"
        + "".join(f"{generation}\t{lnl}\t{tl}\t{alpha}\n" for generation, lnl, tl, alpha in rows),
        encoding="utf-8",
    )
    return path


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


def test_compare_independent_bayesian_runs_reports_parameter_shifts_and_topology_conflict(tmp_path: Path) -> None:
    left_trace = _write_mrbayes_trace(
        tmp_path / "left.run1.p",
        [
            (0, -110.0, 0.40, 0.90),
            (100, -109.7, 0.41, 0.92),
            (200, -109.6, 0.42, 0.91),
            (300, -109.5, 0.41, 0.93),
            (400, -109.4, 0.42, 0.92),
            (500, -109.3, 0.41, 0.91),
        ],
    )
    right_trace = _write_mrbayes_trace(
        tmp_path / "right.run1.p",
        [
            (0, -108.0, 0.55, 1.50),
            (100, -107.8, 0.54, 1.52),
            (200, -107.7, 0.56, 1.49),
            (300, -107.6, 0.55, 1.51),
            (400, -107.5, 0.54, 1.50),
            (500, -107.4, 0.55, 1.48),
        ],
    )

    report = compare_independent_bayesian_runs(
        fixture("trees/example_tree_set_left.nwk"),
        fixture("trees/example_tree_set_right.nwk"),
        left_trace_path=left_trace,
        right_trace_path=right_trace,
        burnin_fraction=0.0,
        ess_threshold=2.0,
        mean_shift_threshold=1.0,
    )

    assert report.trace_kind == "mrbayes"
    assert report.tree_comparison.mcc_topology.topology_equal is False
    assert [row.parameter for row in report.parameter_differences] == ["LnL", "TL", "alpha"]
    assert report.parameter_differences[-1].mean_delta > 0.5
    assert "independent runs select different maximum clade credibility topologies" in report.warnings
    assert "one or more posterior parameter means differ materially across independent runs" in report.warnings


def test_compare_posterior_tree_sets_by_prior_reports_age_shifts() -> None:
    report = compare_posterior_tree_sets_by_prior(
        fixture("trees/example_tree_set_left.nwk"),
        fixture("trees/example_tree_set_right.nwk"),
        left_label="strict-prior",
        right_label="broad-prior",
        burnin_fraction=0.0,
    )

    assert report.comparison_axis == "prior"
    assert report.left_label == "strict-prior"
    assert report.right_label == "broad-prior"
    assert report.tree_comparison.mcc_topology.topology_equal is False
    assert report.age_differences[0].clade in {"A|B", "A|C", "B|D", "C|D"}
    assert "prior comparison changes the maximum clade credibility topology" in report.warnings


def test_compare_posterior_tree_sets_by_clock_reports_clock_labels() -> None:
    report = compare_posterior_tree_sets_by_clock(
        fixture("trees/example_tree_set_left.nwk"),
        fixture("trees/example_tree_set_right.nwk"),
        left_label="strict-clock",
        right_label="relaxed-clock",
        burnin_fraction=0.0,
    )

    assert report.comparison_axis == "clock"
    assert report.left_label == "strict-clock"
    assert report.right_label == "relaxed-clock"
    assert report.tree_comparison.tree_set_comparison.shared_rooted_topology_count == 1
    assert report.age_differences
