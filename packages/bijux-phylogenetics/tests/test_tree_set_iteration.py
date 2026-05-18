from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.reports.service import render_tree_uncertainty_report
from bijux_phylogenetics.trees import (
    assess_tree_set_maturity,
    assess_tree_set_storage_risk,
    assess_tree_set_thinning_sensitivity,
    benchmark_tree_set_uncertainty,
    compare_consensus_thresholds,
    detect_unstable_clades,
    summarize_posterior_topology_diversity,
    write_unstable_clade_table,
)

FIXTURES = Path(__file__).parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_benchmark_tree_set_uncertainty_reports_rows() -> None:
    report = benchmark_tree_set_uncertainty(
        tree_counts=[3], taxon_counts=[4], replicates=1, seed=7
    )

    assert report.tree_counts == [3]
    assert report.taxon_counts == [4]
    assert len(report.rows) == 1
    assert report.rows[0].tree_count == 3
    assert report.rows[0].taxon_count == 4
    assert report.rows[0].peak_memory_bytes >= 0


def test_assess_tree_set_storage_risk_reports_low_risk_for_small_fixture() -> None:
    report = assess_tree_set_storage_risk(fixture("example_tree_set_left.nwk"))

    assert report.tree_count == 3
    assert report.risk_level == "low"


def test_summarize_posterior_topology_diversity_reports_dispersion_metrics(
    tmp_path: Path,
) -> None:
    report = summarize_posterior_topology_diversity(
        fixture("example_tree_set_left.nwk")
    )
    unstable_path = tmp_path / "unstable-clades.tsv"
    write_unstable_clade_table(
        unstable_path,
        detect_unstable_clades(fixture("example_tree_set_left.nwk")),
    )

    assert report.tree_count == 3
    assert report.rooted_topology_count == 2
    assert report.dominant_topology_frequency == pytest.approx(2 / 3)
    assert report.pair_count == 3
    assert report.mean_normalized_robinson_foulds_distance > 0.0
    assert report.maximum_robinson_foulds_distance >= 0
    assert report.unstable_clade_count >= 1
    assert "support_classification" in unstable_path.read_text(encoding="utf-8")


def test_assess_tree_set_thinning_sensitivity_reports_retained_counts() -> None:
    report = assess_tree_set_thinning_sensitivity(
        fixture("example_tree_set_left.nwk"),
        thinning_intervals=[2],
    )

    assert report.original_tree_count == 3
    assert report.rows[0].retained_tree_count == 2
    assert report.rows[0].thinning_interval == 2


def test_compare_consensus_thresholds_reports_threshold_rows() -> None:
    report = compare_consensus_thresholds(
        fixture("example_tree_set_left.nwk"),
        thresholds=[0.5, 0.9],
    )

    assert [row.threshold for row in report.rows] == [0.5, 0.9]
    assert report.rows[0].consensus_newick


def test_assess_tree_set_maturity_reports_decision_and_checks() -> None:
    report = assess_tree_set_maturity(
        fixture("example_tree_set_left.nwk"),
        thinning_intervals=[2],
        consensus_thresholds=[0.5, 0.9],
    )

    assert report.decision in {"experimental", "usable", "production_capable"}
    assert len(report.checks) == 5


def test_render_tree_uncertainty_report_includes_new_iteration_sections(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "tree-uncertainty.html"
    report = render_tree_uncertainty_report(
        tree_set_path=fixture("example_tree_set_left.nwk"),
        out_path=output_path,
    )

    html = output_path.read_text(encoding="utf-8")
    assert report.output_path == output_path
    assert "storage-risk" in html
    assert "thinning-sensitivity" in html
    assert "consensus-threshold-sensitivity" in html
    assert "maturity-gate" in html
    assert "limitations" in html
    assert report.machine_manifest["limitations"]
