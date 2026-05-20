from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.beast import (
    assess_calibration_dominance,
    assess_time_tree_readiness,
)
from bijux_phylogenetics.bayesian.posterior_sets.comparison import (
    compare_independent_bayesian_runs,
    compare_ml_tree_to_bayesian_posterior,
    compare_posterior_tree_sets_by_clock,
    compare_posterior_tree_sets_by_prior,
)
from bijux_phylogenetics.bayesian.posterior_sets.tree_sets import (
    compare_bayesian_tree_sets,
    subsample_posterior_tree_set,
    summarize_maximum_clade_credibility_tree,
    summarize_posterior_node_ages,
    thin_posterior_tree_set,
    write_posterior_tree_subsample_table,
)
from bijux_phylogenetics.bayesian.presentation.html_reports import (
    render_bayesian_run_comparison_report,
    render_ml_vs_bayesian_tree_report,
    render_time_tree_readiness_report,
)
from bijux_phylogenetics.bayesian.presentation.posterior_uncertainty import (
    build_posterior_uncertainty_figure_package,
    write_bayesian_limitations_text,
    write_bayesian_methods_summary_text,
    write_supplementary_bayesian_diagnostics_table,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def _write_mrbayes_trace(
    path: Path, rows: list[tuple[int, float, float, float]]
) -> Path:
    path.write_text(
        "Gen\tLnL\tTL\talpha\n"
        + "".join(
            f"{generation}\t{lnl}\t{tl}\t{alpha}\n"
            for generation, lnl, tl, alpha in rows
        ),
        encoding="utf-8",
    )
    return path


def test_summarize_maximum_clade_credibility_tree_selects_best_supported_posterior_tree() -> (
    None
):
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


def test_thin_posterior_tree_set_retains_every_nth_tree_after_burnin(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "thinned.nwk"
    report = thin_posterior_tree_set(
        fixture("trees/example_tree_set_left.nwk"),
        output_path,
        thinning_interval=2,
        burnin_fraction=0.0,
    )

    lines = [
        line
        for line in output_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    assert report.retained_tree_count == 2
    assert report.retained_indices == [1, 3]
    assert len(lines) == 2


def test_subsample_posterior_tree_set_reports_seeded_random_selection(
    tmp_path: Path,
) -> None:
    report = subsample_posterior_tree_set(
        fixture("trees/example_tree_set_left.nwk"),
        method="random",
        sample_count=2,
        burnin_fraction=0.0,
        random_seed=11,
    )
    table_path = tmp_path / "posterior-subsample.tsv"
    write_posterior_tree_subsample_table(table_path, report)

    rows = table_path.read_text(encoding="utf-8").splitlines()
    assert report.selection_method == "random"
    assert report.requested_tree_count == 2
    assert report.random_seed == 11
    assert report.retained_tree_count == 2
    assert report.retained_source_indices == [2, 3]
    assert report.trees[0].retained_order == 1
    assert report.trees[0].source_tree_index == 2
    assert rows[0].startswith("retained_order\tsource_tree_index")
    assert "\trandom\t" in rows[1]


def test_subsample_posterior_tree_set_rejects_ambiguous_method_arguments() -> None:
    with pytest.raises(ValueError, match="sample_count is not supported"):
        subsample_posterior_tree_set(
            fixture("trees/example_tree_set_left.nwk"),
            method="evenly-spaced",
            thinning_interval=2,
            sample_count=2,
        )


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
    assert rows["A|B"].median_height == 0.2
    assert rows["A|B"].lower_95_credible_interval == 0.2
    assert rows["A|B"].upper_95_credible_interval == 0.2
    assert rows["C|D"].maximum_height == 0.2


def test_compare_ml_tree_to_bayesian_posterior_reports_topology_and_branch_length_differences() -> (
    None
):
    report = compare_ml_tree_to_bayesian_posterior(
        fixture("trees/example_tree.nwk"),
        fixture("trees/example_tree_set_right.nwk"),
        burnin_fraction=0.0,
    )

    assert report.mcc_tree_index == 2
    assert report.topology.topology_equal is False
    assert any(
        "maximum-likelihood and Bayesian" in warning for warning in report.warnings
    )


def test_assess_calibration_dominance_flags_single_dominant_calibration(
    tmp_path: Path,
) -> None:
    calibration_path = tmp_path / "dominant-calibration.tsv"
    calibration_path.write_text(
        "calibration_id\tclade_name\tminimum_age\tmaximum_age\tdistribution\n"
        "cal-1\tMammals\t0.0\t0.3\tuniform\n",
        encoding="utf-8",
    )

    report = assess_calibration_dominance(
        fixture("trees/example_tree_named_clades.nwk"),
        calibration_path,
    )

    assert report.valid_calibration_count == 1
    assert report.dominant_calibration_ids == ["cal-1"]
    assert report.warnings


def test_assess_time_tree_readiness_blocks_invalid_tip_dates() -> None:
    report = assess_time_tree_readiness(
        fixture("trees/example_tree.nwk"),
        tip_dates_path=fixture("metadata/example_tip_dates_invalid.tsv"),
    )

    assert report.decision == "blocked"
    assert (
        "tip-date table contains missing, invalid, or mismatched dated taxa"
        in report.blockers
    )


def test_render_time_tree_readiness_report_writes_sections(tmp_path: Path) -> None:
    output_path = tmp_path / "time-tree-readiness.html"
    report = render_time_tree_readiness_report(
        tree_path=fixture("trees/example_tree.nwk"),
        calibration_path=fixture("metadata/example_calibrations.tsv"),
        tip_dates_path=fixture("metadata/example_tip_dates.tsv"),
        out_path=output_path,
    )

    html = output_path.read_text(encoding="utf-8")
    assert report.output_path == output_path
    assert report.method_tier.tier == "parser-only"
    assert "method-tier" in html
    assert "readiness" in html
    assert "calibration-dominance" in html
    assert "limitations" in html


def test_compare_independent_bayesian_runs_reports_parameter_shifts_and_topology_conflict(
    tmp_path: Path,
) -> None:
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
    assert [row.parameter for row in report.parameter_differences] == [
        "LnL",
        "TL",
        "alpha",
    ]
    assert report.parameter_differences[-1].mean_delta > 0.5
    assert (
        "independent runs select different maximum clade credibility topologies"
        in report.warnings
    )
    assert (
        "one or more posterior parameter means differ materially across independent runs"
        in report.warnings
    )


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
    assert (
        "prior comparison changes the maximum clade credibility topology"
        in report.warnings
    )


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


def test_render_bayesian_run_comparison_report_writes_tree_and_trace_sections(
    tmp_path: Path,
) -> None:
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
    output_path = tmp_path / "bayesian-run-comparison.html"

    report = render_bayesian_run_comparison_report(
        left_tree_set_path=fixture("trees/example_tree_set_left.nwk"),
        right_tree_set_path=fixture("trees/example_tree_set_right.nwk"),
        left_trace_path=left_trace,
        right_trace_path=right_trace,
        out_path=output_path,
        burnin_fraction=0.0,
        ess_threshold=2.0,
        mean_shift_threshold=1.0,
    )

    html = output_path.read_text(encoding="utf-8")
    assert report.output_path == output_path
    assert report.warning_count >= 1
    assert report.method_tier.tier == "parser-only"
    assert "method-tier" in html
    assert "run-comparison" in html
    assert "tree-comparison" in html
    assert "parameter-differences" in html
    assert "limitations" in html
    assert "limitations" in report.machine_manifest["sections"]


def test_render_ml_vs_bayesian_tree_report_includes_limitations(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "ml-vs-bayesian-report.html"

    report = render_ml_vs_bayesian_tree_report(
        ml_tree_path=fixture("trees/example_tree.nwk"),
        posterior_tree_path=fixture("trees/example_tree_set_right.nwk"),
        out_path=output_path,
        burnin_fraction=0.0,
    )

    html = output_path.read_text(encoding="utf-8")
    assert report.output_path == output_path
    assert report.warning_count >= 1
    assert "ml-versus-bayesian-summary" in html
    assert "limitations" in html
    assert "limitations" in report.machine_manifest["sections"]


def test_build_posterior_uncertainty_figure_package_writes_consensus_plot_and_tables(
    tmp_path: Path,
) -> None:
    result = build_posterior_uncertainty_figure_package(
        fixture("trees/example_tree_set_left.nwk"),
        out_dir=tmp_path / "posterior-uncertainty-package",
    )

    assert result.consensus_tree_path.exists()
    assert result.consensus_figure_path.exists()
    assert result.clade_support_plot_path.exists()
    assert result.unstable_taxa_plot_path.exists()
    assert result.topology_clusters_plot_path.exists()
    assert result.unstable_taxa_table_path.exists()
    assert result.topology_clusters_table_path.exists()
    assert result.uncertainty_conclusions_table_path.exists()
    assert result.conclusion_summary_path.exists()
    assert result.review_path.exists()
    summary = result.conclusion_summary_path.read_text(encoding="utf-8")
    assert "Tree-Set Uncertainty Summary" in summary
    assert "Conflict-prone clades" in summary


def test_write_supplementary_bayesian_diagnostics_table_writes_burnin_and_chain_rows(
    tmp_path: Path,
) -> None:
    second_chain = tmp_path / "chain-2.log"
    second_chain.write_text(
        "# BEAST fixture log\n"
        "state\tposterior\tlikelihood\tclockRate\ttreeHeight\n"
        "0\t-501.0\t-481.0\t0.0010\t13.0\n"
        "1000\t-500.8\t-480.8\t0.0011\t13.1\n"
        "2000\t-500.6\t-480.6\t0.0012\t13.1\n"
        "3000\t-500.5\t-480.5\t0.0011\t13.2\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "bayesian-diagnostics-table.tsv"

    report = write_supplementary_bayesian_diagnostics_table(
        output_path,
        posterior_tree_path=fixture("trees/example_tree_set_left.nwk"),
        primary_log_path=fixture("metadata/example_beast.log"),
        additional_log_paths=[second_chain],
        burnin_fractions=(0.0, 0.25),
        ess_threshold=2.0,
        mean_shift_threshold=1.0,
        cross_chain_mean_shift_threshold=5.0,
    )

    text = output_path.read_text(encoding="utf-8")
    assert report.output_path == output_path
    assert report.chain_count == 2
    assert "row_kind\tchain\tparameter" in text
    assert "burnin-summary\tprimary" in text
    assert "chain-parameter\tchain_1\tposterior" in text


def test_write_bayesian_methods_summary_text_describes_clock_prior_and_diagnostics(
    tmp_path: Path,
) -> None:
    second_chain = tmp_path / "chain-2.log"
    second_chain.write_text(
        "# BEAST fixture log\n"
        "state\tposterior\tlikelihood\tclockRate\ttreeHeight\n"
        "0\t-501.0\t-481.0\t0.0010\t13.0\n"
        "1000\t-500.8\t-480.8\t0.0011\t13.1\n"
        "2000\t-500.6\t-480.6\t0.0012\t13.1\n"
        "3000\t-500.5\t-480.5\t0.0011\t13.2\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "bayesian-methods-summary.md"

    report = write_bayesian_methods_summary_text(
        output_path,
        posterior_tree_path=fixture("trees/example_tree_set_left.nwk"),
        primary_log_path=fixture("metadata/example_beast.log"),
        additional_log_paths=[second_chain],
        analysis_xml_path=fixture("metadata/beast2_strict_yule_posterior.xml"),
        tree_prior="birth-death",
        clock_model="relaxed-lognormal",
        calibration_path=fixture("metadata/example_calibrations.tsv"),
        tip_dates_path=fixture("metadata/example_tip_dates.tsv"),
        burnin_fractions=(0.0, 0.25),
        ess_threshold=2.0,
        mean_shift_threshold=1.0,
        cross_chain_mean_shift_threshold=5.0,
    )

    text = output_path.read_text(encoding="utf-8")
    assert report.output_path == output_path
    assert "strict" in text
    assert "yule" in text
    assert "example_calibrations.tsv" in text
    assert "effective sample size" in text
    assert "beast2_strict_yule_posterior.xml" in text
    assert "chain length" in text.lower()
    assert "did not execute BEAST itself" in text
    assert "only summarized the prepared XML" in text


def test_write_bayesian_limitations_text_includes_diagnostics_and_dating_risks(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "bayesian-limitations.md"
    report = write_bayesian_limitations_text(
        output_path,
        posterior_tree_path=fixture("trees/example_tree_set_left.nwk"),
        primary_log_path=fixture("metadata/example_beast.log"),
        tree_path=fixture("trees/example_tree.nwk"),
        tip_dates_path=fixture("metadata/example_tip_dates_invalid.tsv"),
    )

    text = output_path.read_text(encoding="utf-8")
    assert report.output_path == output_path
    assert "Bayesian Analysis Limitations" in text
    assert "tip-date" in text.lower()
