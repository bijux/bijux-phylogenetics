from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.benchmark as benchmark_api
from bijux_phylogenetics.benchmark.model_fitting import (
    LargeTreeModelFittingBenchmarkBundle,
    LargeTreeModelFittingThreshold,
    _case_definitions_for_tier,
    _evaluate_threshold,
    benchmark_large_tree_model_fitting,
    write_large_tree_model_fitting_bundle,
)
from tests.support.scientific_output_assertions import (
    assert_selected_scientific_outputs_equivalent,
)


@pytest.mark.slow
def test_benchmark_large_tree_model_fitting_small_tier_reports_governed_cases() -> None:
    report = benchmark_large_tree_model_fitting(tier="small")

    assert report.tier == "small"
    assert report.case_count == 2
    assert report.geiger_match_case_count == 2
    assert report.threshold_pass_case_count == 2
    assert report.too_slow_case_count == 0
    assert report.unstable_case_count == 0
    observations = {row.case_id: row for row in report.observations}
    assert set(observations) == {
        "fitcontinuous-pagel-lambda-100-taxa",
        "fitdiscrete-er-binary-100-taxa",
    }

    continuous = observations["fitcontinuous-pagel-lambda-100-taxa"]
    assert continuous.trait_kind == "continuous"
    assert continuous.taxon_count == 100
    assert continuous.runtime_seconds is not None and continuous.runtime_seconds > 0.0
    assert (
        continuous.peak_memory_bytes is not None and continuous.peak_memory_bytes >= 0
    )
    assert continuous.optimizer_iteration_count == 19
    assert continuous.performance_threshold_passed is True
    assert continuous.matches_geiger_reference is True
    assert continuous.parameter_delta is not None and continuous.parameter_delta <= 0.4

    discrete = observations["fitdiscrete-er-binary-100-taxa"]
    assert discrete.trait_kind == "discrete"
    assert discrete.taxon_count == 100
    assert discrete.runtime_seconds is not None and discrete.runtime_seconds > 0.0
    assert discrete.optimizer_iteration_count is not None
    assert discrete.performance_threshold_passed is True
    assert discrete.matches_geiger_reference is True
    assert discrete.rate_delta is not None and discrete.rate_delta <= 0.5


def test_case_definitions_for_heavy_tier_include_512_taxon_continuous_review() -> None:
    cases = _case_definitions_for_tier("heavy")
    case_ids = {case.case_id for case in cases}
    assert "fitcontinuous-brownian-512-taxa" in case_ids
    assert "fitcontinuous-pagel-lambda-100-taxa" in case_ids
    assert "fitdiscrete-er-binary-100-taxa" in case_ids


def test_evaluate_threshold_flags_runtime_regression() -> None:
    threshold = LargeTreeModelFittingThreshold(
        max_runtime_seconds=5.0,
        max_peak_memory_bytes=1024,
        max_optimizer_step_count=10,
    )
    (
        runtime_within_threshold,
        peak_memory_within_threshold,
        optimizer_step_within_threshold,
        performance_threshold_passed,
    ) = _evaluate_threshold(
        threshold=threshold,
        runtime_seconds=6.0,
        peak_memory_bytes=512,
        optimizer_step_count=8,
    )
    assert runtime_within_threshold is False
    assert peak_memory_within_threshold is True
    assert optimizer_step_within_threshold is True
    assert performance_threshold_passed is False


def test_public_runtime_exports_include_large_tree_model_fitting_benchmark() -> None:
    assert (
        benchmark_api.benchmark_large_tree_model_fitting
        is benchmark_large_tree_model_fitting
    )
    assert (
        benchmark_api.write_large_tree_model_fitting_bundle
        is write_large_tree_model_fitting_bundle
    )


@pytest.mark.slow
def test_write_large_tree_model_fitting_bundle_matches_expected_outputs(
    tmp_path: Path,
) -> None:
    bundle = write_large_tree_model_fitting_bundle(tmp_path / "large-tree-benchmark")

    assert isinstance(bundle, LargeTreeModelFittingBenchmarkBundle)
    expected_root = (
        Path(__file__).resolve().parents[1]
        / "src"
        / "bijux_phylogenetics"
        / "resources"
        / "benchmarks"
        / "large_tree_model_fitting"
        / "expected"
    )
    generated = {
        Path("summary.tsv"): bundle.summary_path,
        Path("observations.tsv"): bundle.observation_table_path,
    }
    assert_selected_scientific_outputs_equivalent(expected_root, generated)
