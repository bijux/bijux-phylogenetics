from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics.benchmark as benchmark_api
from bijux_phylogenetics.benchmark.likelihood_tree_inference import (
    NativeMaximumLikelihoodAccuracyBenchmarkObservation,
    NativeMaximumLikelihoodAccuracyBenchmarkReport,
    benchmark_native_maximum_likelihood_accuracy,
    write_native_maximum_likelihood_accuracy_observation_table,
    write_native_maximum_likelihood_accuracy_run_json,
    write_native_maximum_likelihood_accuracy_summary_table,
    write_native_maximum_likelihood_accuracy_wrapper_table,
)
from bijux_phylogenetics.benchmark.likelihood_tree_inference.contracts import (
    MaximumLikelihoodWrapperBenchmarkObservation,
)

pytestmark = pytest.mark.slow


def test_benchmark_native_maximum_likelihood_accuracy_reports_truth_and_wrapper_lanes() -> (
    None
):
    report = benchmark_native_maximum_likelihood_accuracy()

    assert report.model_names == ["k80"]
    assert report.search_methods == ["nni"]
    assert len(report.observations) == 1
    assert len(report.wrapper_observations) == 4
    assert report.topology_recovery_case_count == 1
    assert report.parameter_recovery_case_count == 1
    assert report.parameter_not_applicable_case_count == 0
    assert report.likelihood_quality_case_count == 1
    assert report.wrapper_supported_case_count == 3
    assert report.wrapper_clear_case_count == 3
    assert report.all_truth_cases_clear is True
    assert report.all_supported_wrapper_cases_clear is True
    assert report.limitations

    row = report.observations[0]
    assert row.case_id == "k80-four-taxon-transition-bias"
    assert row.model_name == "k80"
    assert row.search_method == "nni"
    assert row.selected_model_name == "K80"
    assert row.topology_recovered is True
    assert row.topology_recovery_status == "topology-recovered"
    assert row.robinson_foulds_distance == 0
    assert row.parameter_recovery_status == "parameter-recovered"
    assert row.parameter_absolute_deltas["kappa"] <= 2.0
    assert row.likelihood_quality_status == "no-likelihood-shortfall"
    assert row.likelihood_shortfall <= 0.25


def test_native_maximum_likelihood_accuracy_rejects_models_without_truth_cases() -> (
    None
):
    with pytest.raises(ValueError, match="supports simulated truth cases"):
        benchmark_native_maximum_likelihood_accuracy(model_names=["jc69"])


def test_native_maximum_likelihood_accuracy_writers_emit_governed_outputs(
    tmp_path: Path,
) -> None:
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

    summary_path = write_native_maximum_likelihood_accuracy_summary_table(
        tmp_path / "summary.tsv",
        report,
    )
    observation_path = write_native_maximum_likelihood_accuracy_observation_table(
        tmp_path / "observations.tsv",
        report,
    )
    wrapper_path = write_native_maximum_likelihood_accuracy_wrapper_table(
        tmp_path / "wrapper.tsv",
        report,
    )
    run_json_path = write_native_maximum_likelihood_accuracy_run_json(
        tmp_path / "run.json",
        report,
    )

    assert summary_path.read_text(encoding="utf-8").startswith(
        "observation_count\twrapper_observation_count\tmodel_names\tsearch_methods\ttopology_recovery_case_count\tparameter_recovery_case_count\tparameter_not_applicable_case_count\tlikelihood_quality_case_count\twrapper_supported_case_count\twrapper_clear_case_count\tall_truth_cases_clear\tall_supported_wrapper_cases_clear\tlimitations\n"
    )
    assert observation_path.read_text(encoding="utf-8").startswith(
        "case_id\tmodel_name\tsearch_method\ttaxon_count\tsite_count\tsimulation_seed\tselected_model_name\ttopology_recovered\ttopology_recovery_status\trobinson_foulds_distance\tnormalized_robinson_foulds\tparameter_recovery_status\tparameter_truth\testimated_parameter_values\tparameter_absolute_deltas\ttrue_topology_log_likelihood\tinferred_log_likelihood\tlikelihood_shortfall\tlikelihood_quality_status\n"
    )
    assert wrapper_path.read_text(encoding="utf-8").startswith(
        "case_id\twrapper_engine\tstatus\tsupported\tblocking\ttopology_equal\tnative_log_likelihood\twrapper_log_likelihood\tlog_likelihood_delta\tsupport_split_count\tnotes\n"
    )
    payload = json.loads(run_json_path.read_text(encoding="utf-8"))
    assert payload["model_names"] == ["k80"]
    assert payload["search_methods"] == ["nni"]
    assert len(payload["observations"]) == 1
    assert len(payload["wrapper_observations"]) == 1


def test_public_runtime_exports_include_native_maximum_likelihood_accuracy_benchmark() -> (
    None
):
    assert (
        benchmark_api.benchmark_native_maximum_likelihood_accuracy
        is benchmark_native_maximum_likelihood_accuracy
    )
    assert (
        benchmark_api.write_native_maximum_likelihood_accuracy_summary_table
        is write_native_maximum_likelihood_accuracy_summary_table
    )
    assert (
        benchmark_api.write_native_maximum_likelihood_accuracy_observation_table
        is write_native_maximum_likelihood_accuracy_observation_table
    )
    assert (
        benchmark_api.write_native_maximum_likelihood_accuracy_wrapper_table
        is write_native_maximum_likelihood_accuracy_wrapper_table
    )
    assert (
        benchmark_api.write_native_maximum_likelihood_accuracy_run_json
        is write_native_maximum_likelihood_accuracy_run_json
    )
