from __future__ import annotations

import json

import pytest

from bijux_phylogenetics.benchmark import (
    MaximumLikelihoodWrapperBenchmarkObservation,
    MaximumLikelihoodWrapperBenchmarkReport,
    MaximumLikelihoodWrapperBenchmarkSummaryRow,
    MaximumLikelihoodWrapperSupportComparison,
)
import bijux_phylogenetics.command_line as command_line_api
from bijux_phylogenetics.command_line import main


def test_cli_benchmark_maximum_likelihood_wrapper_comparison_reports_metrics(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    report = MaximumLikelihoodWrapperBenchmarkReport(
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
                support_comparisons=[
                    MaximumLikelihoodWrapperSupportComparison(
                        split_taxa=["A", "B"],
                        native_support_percent=100.0,
                        wrapper_support_percent=100.0,
                        support_percent_delta=0.0,
                    )
                ],
                notes=["governed test fixture"],
            )
        ],
        summary_rows=[
            MaximumLikelihoodWrapperBenchmarkSummaryRow(
                status="tolerance-match",
                case_count=1,
                blocking_case_count=0,
                case_ids=["iqtree-jc69-maximum-likelihood"],
            ),
            MaximumLikelihoodWrapperBenchmarkSummaryRow(
                status="expected-model-assumption-difference",
                case_count=0,
                blocking_case_count=0,
                case_ids=[],
            ),
            MaximumLikelihoodWrapperBenchmarkSummaryRow(
                status="unsupported-case",
                case_count=0,
                blocking_case_count=0,
                case_ids=[],
            ),
            MaximumLikelihoodWrapperBenchmarkSummaryRow(
                status="native-bug",
                case_count=0,
                blocking_case_count=0,
                case_ids=[],
            ),
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
    )
    monkeypatch.setattr(
        command_line_api,
        "benchmark_maximum_likelihood_wrapper_comparison",
        lambda: report,
    )

    exit_code = main(["benchmark", "maximum-likelihood-wrapper-comparison", "--json"])
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["command"] == "benchmark"
    assert payload["metrics"]["case_count"] == 1
    assert payload["metrics"]["observation_count"] == 1
    assert payload["metrics"]["summary_row_count"] == 4
    assert payload["data"]["case_count"] == 1
    assert payload["data"]["tolerance_match_case_count"] == 1
