from __future__ import annotations

import json

import pytest

from bijux_phylogenetics.benchmark import (
    MaximumLikelihoodWrapperBenchmarkObservation,
    NativeMaximumLikelihoodAccuracyBenchmarkObservation,
    NativeMaximumLikelihoodAccuracyBenchmarkReport,
)
import bijux_phylogenetics.command_line as command_line_api
from bijux_phylogenetics.command_line import main


def test_cli_benchmark_native_maximum_likelihood_accuracy_reports_metrics_and_overrides(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    captured_call: dict[str, object] = {}
    report = NativeMaximumLikelihoodAccuracyBenchmarkReport(
        observations=[
            NativeMaximumLikelihoodAccuracyBenchmarkObservation(
                case_id="k80-four-taxon-transition-bias",
                model_name="k80",
                search_method="nni",
                taxon_count=4,
                site_count=16,
                simulation_seed=203,
                selected_model_name="K80",
                true_topology_fingerprint="A|B",
                inferred_topology_fingerprint="A|B",
                topology_recovered=True,
                topology_recovery_status="topology-recovered",
                robinson_foulds_distance=0,
                normalized_robinson_foulds=0.0,
                parameter_truth={"kappa": 4.0},
                estimated_parameter_values={"kappa": 2.9174759792656726},
                parameter_absolute_deltas={"kappa": 1.082524020734327},
                parameter_recovery_status="parameter-recovered",
                true_topology_log_likelihood=-56.071600988319304,
                inferred_log_likelihood=-56.179489112246024,
                likelihood_shortfall=0.10788812392672,
                likelihood_quality_status="no-likelihood-shortfall",
            )
        ],
        wrapper_observations=[
            MaximumLikelihoodWrapperBenchmarkObservation(
                case_id="iqtree-jc69-maximum-likelihood",
                wrapper_engine="IQ-TREE-style",
                native_surface="infer_nucleotide_likelihood_tree_from_alignment",
                wrapper_surface="cached IQ-TREE JC maximum-likelihood reference",
                comparison_policy="topology-and-likelihood-tolerance",
                status="tolerance-match",
                supported=True,
                blocking=False,
                native_runtime_seconds=0.75,
                native_peak_memory_bytes=4096,
                native_memory_observation_kind="python-tracemalloc",
                wrapper_runtime_seconds=0.018,
                wrapper_peak_memory_bytes=8110080,
                wrapper_memory_observation_kind="macos-time-maximum-resident-set-size",
                runtime_delta_seconds=0.732,
                peak_memory_delta_bytes=-8105984,
                native_log_likelihood=-33.27106468396422,
                wrapper_log_likelihood=-33.271,
                log_likelihood_delta=-6.468396421865918e-05,
                native_tree_newick="((A:0.0,B:0.0):1.0,(C:0.0,D:0.0):1.0);",
                wrapper_tree_newick="(A:0.0,B:0.0,(C:0.0,D:0.0):1.0);",
                native_topology_signature=("A|B",),
                wrapper_topology_signature=("A|B",),
                robinson_foulds_distance=0,
                normalized_robinson_foulds=0.0,
                topology_equal=True,
                support_comparisons=[],
                notes=["governed test fixture"],
            )
        ],
        model_names=["k80"],
        search_methods=["nni"],
        topology_recovery_case_count=1,
        parameter_recovery_case_count=1,
        parameter_not_applicable_case_count=0,
        likelihood_quality_case_count=1,
        wrapper_supported_case_count=1,
        wrapper_clear_case_count=1,
        all_truth_cases_clear=True,
        all_supported_wrapper_cases_clear=True,
        limitations=["governed test fixture"],
    )

    def fake_benchmark_native_maximum_likelihood_accuracy(**kwargs):
        captured_call.update(kwargs)
        return report

    monkeypatch.setattr(
        command_line_api,
        "benchmark_native_maximum_likelihood_accuracy",
        fake_benchmark_native_maximum_likelihood_accuracy,
    )

    exit_code = main(
        [
            "benchmark",
            "native-maximum-likelihood-accuracy",
            "--model",
            "k80",
            "--search-method",
            "nni",
            "--start-tree-count",
            "1",
            "--start-tree-seed",
            "17",
            "--upper-branch-length-bound",
            "1.0",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert captured_call == {
        "model_names": ["k80"],
        "search_methods": ["nni"],
        "start_tree_count": 1,
        "start_tree_seed": 17,
        "upper_branch_length_bound": 1.0,
    }
    assert payload["command"] == "benchmark"
    assert payload["metrics"]["observation_count"] == 1
    assert payload["metrics"]["topology_recovery_case_count"] == 1
    assert payload["metrics"]["parameter_recovery_case_count"] == 1
    assert payload["metrics"]["likelihood_quality_case_count"] == 1
    assert payload["metrics"]["wrapper_clear_case_count"] == 1
    assert payload["metrics"]["model_count"] == 1
    assert payload["metrics"]["search_method_count"] == 1
    assert payload["data"]["model_names"] == ["k80"]
    assert payload["data"]["search_methods"] == ["nni"]
