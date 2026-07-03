from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics.benchmark as benchmark_api
from bijux_phylogenetics.benchmark.likelihood_tree_inference import (
    MaximumLikelihoodWrapperBenchmarkObservation,
    MaximumLikelihoodWrapperBenchmarkReport,
    MaximumLikelihoodWrapperBenchmarkSummaryRow,
    NativeMaximumLikelihoodAccuracyBenchmarkObservation,
    NativeMaximumLikelihoodAccuracyBenchmarkReport,
    NativeMaximumLikelihoodMemoryBenchmarkReport,
    NativeMaximumLikelihoodSuiteObservation,
    NativeMaximumLikelihoodSuiteReport,
    NativeMaximumLikelihoodSuiteSummaryRow,
    NativeMaximumLikelihoodSpeedBenchmarkReport,
    benchmark_native_maximum_likelihood_suite,
    write_native_maximum_likelihood_suite_observation_table,
    write_native_maximum_likelihood_suite_run_json,
    write_native_maximum_likelihood_suite_summary_table,
)


pytestmark = pytest.mark.slow


def test_benchmark_native_maximum_likelihood_suite_reports_all_required_dimensions() -> (
    None
):
    report = benchmark_native_maximum_likelihood_suite(
        replicates=1,
        size_classes=[("tiny", 4, 12), ("medium", 4, 16)],
        scaling_model_names=["jc69"],
        scaling_search_methods=["nni"],
        accuracy_model_names=["k80"],
        accuracy_search_methods=["nni"],
    )

    assert report.speed_report.model_names == ["jc69"]
    assert report.memory_report.model_names == ["jc69"]
    assert report.accuracy_report.model_names == ["k80"]
    assert report.wrapper_report.case_count == 4
    assert report.limitations
    assert report.all_supported_cases_clear is True
    assert {
        row.dimension for row in report.observations
    } == {"speed", "memory", "likelihood", "topology", "support", "accuracy"}

    observations = {(row.case_id, row.dimension): row for row in report.observations}
    assert (
        observations[
            ("iqtree-jc69-maximum-likelihood", "likelihood")
        ].status
        == "native-advantage"
    )
    assert (
        observations[
            ("iqtree-jc69-ultrafast-bootstrap-support", "support")
        ].status
        == "native-advantage"
    )
    assert (
        observations[
            ("fasttree-jc69-approximate-maximum-likelihood", "likelihood")
        ].status
        == "expected-assumption-difference"
    )
    assert (
        observations[
            ("raxml-small-maximum-likelihood-reference", "likelihood")
        ].status
        == "unsupported-case"
    )
    assert (
        observations[("k80-four-taxon-transition-bias", "accuracy")].status
        == "native-advantage"
    )
    assert report.native_advantage_case_count > 0
    assert report.expected_assumption_difference_case_count > 0
    assert report.unsupported_case_count == 1
    assert report.native_bug_case_count == 0


def test_public_runtime_exports_include_native_maximum_likelihood_suite_benchmark() -> (
    None
):
    assert (
        benchmark_api.benchmark_native_maximum_likelihood_suite
        is benchmark_native_maximum_likelihood_suite
    )
    assert (
        benchmark_api.write_native_maximum_likelihood_suite_summary_table
        is write_native_maximum_likelihood_suite_summary_table
    )
    assert (
        benchmark_api.write_native_maximum_likelihood_suite_observation_table
        is write_native_maximum_likelihood_suite_observation_table
    )
    assert (
        benchmark_api.write_native_maximum_likelihood_suite_run_json
        is write_native_maximum_likelihood_suite_run_json
    )


def test_native_maximum_likelihood_suite_writers_emit_governed_outputs(
    tmp_path: Path,
) -> None:
    report = NativeMaximumLikelihoodSuiteReport(
        speed_report=NativeMaximumLikelihoodSpeedBenchmarkReport(
            replicates=1,
            observations=[],
            trend_rows=[],
            size_class_labels=["tiny", "medium"],
            model_names=["jc69"],
            search_methods=["nni"],
            limitations=["governed test fixture"],
        ),
        memory_report=NativeMaximumLikelihoodMemoryBenchmarkReport(
            replicates=1,
            observations=[],
            trend_rows=[],
            size_class_labels=["tiny", "medium"],
            model_names=["jc69"],
            search_methods=["nni"],
            stable_memory_case_count=1,
            elevated_memory_case_count=0,
            memory_blowup_case_count=0,
            limitations=["governed test fixture"],
        ),
        accuracy_report=NativeMaximumLikelihoodAccuracyBenchmarkReport(
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
                    estimated_parameter_values={"kappa": 2.9},
                    parameter_absolute_deltas={"kappa": 1.1},
                    parameter_recovery_status="parameter-recovered",
                    true_topology_log_likelihood=-56.07,
                    inferred_log_likelihood=-56.18,
                    likelihood_shortfall=0.11,
                    likelihood_quality_status="no-likelihood-shortfall",
                )
            ],
            wrapper_observations=[],
            model_names=["k80"],
            search_methods=["nni"],
            topology_recovery_case_count=1,
            parameter_recovery_case_count=1,
            parameter_not_applicable_case_count=0,
            likelihood_quality_case_count=1,
            wrapper_supported_case_count=0,
            wrapper_clear_case_count=0,
            all_truth_cases_clear=True,
            all_supported_wrapper_cases_clear=True,
            limitations=["governed test fixture"],
        ),
        wrapper_report=MaximumLikelihoodWrapperBenchmarkReport(
            observations=[
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
            summary_rows=[
                MaximumLikelihoodWrapperBenchmarkSummaryRow(
                    status="tolerance-match",
                    case_count=1,
                    blocking_case_count=0,
                    case_ids=["iqtree-jc69-maximum-likelihood"],
                )
            ],
            case_count=1,
            supported_case_count=1,
            tolerance_match_case_count=1,
            expected_model_assumption_difference_case_count=0,
            unsupported_case_count=0,
            native_bug_case_count=0,
            blocking_case_count=0,
            all_supported_cases_clear=True,
            limitations=["governed test fixture"],
        ),
        observations=[
            NativeMaximumLikelihoodSuiteObservation(
                case_id="iqtree-jc69-maximum-likelihood",
                dimension="likelihood",
                status="native-advantage",
                native_surface="infer_nucleotide_likelihood_tree_from_alignment",
                comparison_surface="cached IQ-TREE JC maximum-likelihood reference",
                measurement_unit="log-likelihood",
                native_value=-33.27106468396422,
                comparison_value=-33.271,
                delta_value=-6.468396421865918e-05,
                notes=["governed test fixture"],
            ),
            NativeMaximumLikelihoodSuiteObservation(
                case_id="k80-four-taxon-transition-bias",
                dimension="accuracy",
                status="native-advantage",
                native_surface="benchmark_native_maximum_likelihood_accuracy",
                comparison_surface="governed simulated truth",
                measurement_unit="normalized-robinson-foulds",
                native_value=0.0,
                comparison_value=0.0,
                delta_value=0.11,
                notes=["parameter status: parameter-recovered"],
            ),
        ],
        summary_rows=[
            NativeMaximumLikelihoodSuiteSummaryRow(
                status="native-advantage",
                observation_count=2,
                dimensions=["accuracy", "likelihood"],
                case_ids=[
                    "iqtree-jc69-maximum-likelihood",
                    "k80-four-taxon-transition-bias",
                ],
            )
        ],
        native_advantage_case_count=2,
        native_bug_case_count=0,
        expected_assumption_difference_case_count=0,
        unsupported_case_count=0,
        all_supported_cases_clear=True,
        limitations=["governed test fixture"],
    )

    summary_path = write_native_maximum_likelihood_suite_summary_table(
        tmp_path / "summary.tsv",
        report,
    )
    observation_path = write_native_maximum_likelihood_suite_observation_table(
        tmp_path / "observations.tsv",
        report,
    )
    run_json_path = write_native_maximum_likelihood_suite_run_json(
        tmp_path / "run.json",
        report,
    )

    assert summary_path.read_text(encoding="utf-8").startswith(
        "observation_count\tsummary_row_count\tnative_advantage_case_count\tnative_bug_case_count\texpected_assumption_difference_case_count\tunsupported_case_count\tall_supported_cases_clear\tlimitations\n"
    )
    assert observation_path.read_text(encoding="utf-8").startswith(
        "case_id\tdimension\tstatus\tnative_surface\tcomparison_surface\tmeasurement_unit\tnative_value\tcomparison_value\tdelta_value\tnotes\n"
    )
    payload = json.loads(run_json_path.read_text(encoding="utf-8"))
    assert payload["native_advantage_case_count"] == 2
    assert len(payload["observations"]) == 2
