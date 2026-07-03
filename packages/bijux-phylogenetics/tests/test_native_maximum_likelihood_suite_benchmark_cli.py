from __future__ import annotations

import json

import pytest

from bijux_phylogenetics.benchmark import (
    MaximumLikelihoodWrapperBenchmarkReport,
    NativeMaximumLikelihoodAccuracyBenchmarkReport,
    NativeMaximumLikelihoodMemoryBenchmarkReport,
    NativeMaximumLikelihoodSuiteObservation,
    NativeMaximumLikelihoodSuiteReport,
    NativeMaximumLikelihoodSuiteSummaryRow,
    NativeMaximumLikelihoodSpeedBenchmarkReport,
)
import bijux_phylogenetics.command_line as command_line_api
from bijux_phylogenetics.command_line import main


def test_cli_native_maximum_likelihood_suite_reports_metrics(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
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
            observations=[],
            wrapper_observations=[],
            model_names=["k80"],
            search_methods=["nni"],
            topology_recovery_case_count=0,
            parameter_recovery_case_count=0,
            parameter_not_applicable_case_count=0,
            likelihood_quality_case_count=0,
            wrapper_supported_case_count=0,
            wrapper_clear_case_count=0,
            all_truth_cases_clear=True,
            all_supported_wrapper_cases_clear=True,
            limitations=["governed test fixture"],
        ),
        wrapper_report=MaximumLikelihoodWrapperBenchmarkReport(
            observations=[],
            summary_rows=[],
            case_count=0,
            supported_case_count=0,
            tolerance_match_case_count=0,
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
            )
        ],
        summary_rows=[
            NativeMaximumLikelihoodSuiteSummaryRow(
                status="native-advantage",
                observation_count=1,
                dimensions=["likelihood"],
                case_ids=["iqtree-jc69-maximum-likelihood"],
            )
        ],
        native_advantage_case_count=1,
        native_bug_case_count=0,
        expected_assumption_difference_case_count=0,
        unsupported_case_count=0,
        all_supported_cases_clear=True,
        limitations=["governed test fixture"],
    )
    monkeypatch.setattr(
        command_line_api,
        "benchmark_native_maximum_likelihood_suite",
        lambda **_: report,
    )

    exit_code = main(["benchmark", "native-maximum-likelihood-suite", "--json"])
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["command"] == "benchmark"
    assert payload["metrics"]["observation_count"] == 1
    assert payload["metrics"]["summary_row_count"] == 1
    assert payload["metrics"]["native_advantage_case_count"] == 1
    assert payload["data"]["native_advantage_case_count"] == 1
    assert payload["data"]["observations"][0]["dimension"] == "likelihood"
