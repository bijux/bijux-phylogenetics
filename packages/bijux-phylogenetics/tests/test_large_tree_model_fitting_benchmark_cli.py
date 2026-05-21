from __future__ import annotations

import json

import pytest

from bijux_phylogenetics.benchmark.model_fitting import (
    LargeTreeModelFittingBenchmarkReport,
    LargeTreeModelFittingObservation,
    LargeTreeModelFittingThreshold,
)
import bijux_phylogenetics.command_line as command_line_api
from bijux_phylogenetics.command_line import main


def test_cli_benchmark_large_tree_model_fitting_reports_metrics(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    observation = LargeTreeModelFittingObservation(
        case_id="fitcontinuous-pagel-lambda-100-taxa",
        tier="small",
        trait_kind="continuous",
        fit_surface="fitcontinuous-pagel-lambda",
        taxon_count=100,
        status="ok",
        runtime_seconds=12.5,
        peak_memory_bytes=4096,
        memory_observation_kind="python-tracemalloc",
        optimizer_name="bounded-profile",
        optimizer_iteration_count=10,
        optimizer_function_evaluation_count=10,
        converged=True,
        hit_lower_parameter_boundary=False,
        hit_upper_parameter_boundary=False,
        unstable_review=False,
        too_slow_review=False,
        stable_conclusion_supported=True,
        threshold=LargeTreeModelFittingThreshold(
            max_runtime_seconds=35.0,
            max_peak_memory_bytes=1024 * 1024,
            max_optimizer_step_count=12,
        ),
        runtime_within_threshold=True,
        peak_memory_within_threshold=True,
        optimizer_step_within_threshold=True,
        performance_threshold_passed=True,
        geiger_reference_available=True,
        geiger_runtime_seconds=0.34,
        geiger_optimizer_step_count=100,
        geiger_parameter_name="lambda",
        geiger_parameter_value=0.9,
        geiger_rate=0.8,
        geiger_log_likelihood=-77.3,
        geiger_aic=160.6,
        geiger_aicc=160.8,
        parameter_delta=0.03,
        rate_delta=0.05,
        log_likelihood_delta=0.05,
        aic_delta=0.11,
        geiger_match_tolerance=0.4,
        matches_geiger_reference=True,
        notes=["governed test fixture"],
    )
    report = LargeTreeModelFittingBenchmarkReport(
        tier="small",
        observations=[observation],
        case_count=1,
        geiger_match_case_count=1,
        threshold_pass_case_count=1,
        too_slow_case_count=0,
        unstable_case_count=0,
        limitations=["governed test fixture"],
    )
    monkeypatch.setattr(
        command_line_api,
        "benchmark_large_tree_model_fitting",
        lambda tier="small": report,
    )

    exit_code = main(["benchmark", "large-tree-model-fitting", "--json"])
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["command"] == "benchmark"
    assert payload["metrics"]["tier"] == "small"
    assert payload["metrics"]["observation_count"] == 1
    assert payload["metrics"]["case_count"] == 1
    assert payload["metrics"]["geiger_match_case_count"] == 1
    assert payload["metrics"]["threshold_pass_case_count"] == 1
    assert payload["metrics"]["max_taxon_count"] == 100
    assert payload["data"]["tier"] == "small"
