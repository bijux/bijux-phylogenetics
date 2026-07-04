from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics.benchmark as benchmark_api
from bijux_phylogenetics.benchmark.likelihood_tree_inference import (
    benchmark_maximum_likelihood_wrapper_comparison,
    write_maximum_likelihood_wrapper_comparison_observation_table,
    write_maximum_likelihood_wrapper_comparison_run_json,
    write_maximum_likelihood_wrapper_comparison_summary_table,
)

pytestmark = pytest.mark.slow


def test_benchmark_maximum_likelihood_wrapper_comparison_reports_cached_ml_lanes() -> (
    None
):
    report = benchmark_maximum_likelihood_wrapper_comparison()

    assert report.case_count == 4
    assert report.supported_case_count == 3
    assert report.tolerance_match_case_count == 2
    assert report.expected_model_assumption_difference_case_count == 1
    assert report.unsupported_case_count == 1
    assert report.native_bug_case_count == 0
    assert report.blocking_case_count == 0
    assert report.all_supported_cases_clear is True

    observations = {row.case_id: row for row in report.observations}
    iqtree_ml = observations["iqtree-jc69-maximum-likelihood"]
    assert iqtree_ml.status == "tolerance-match"
    assert iqtree_ml.topology_equal is True
    assert abs(iqtree_ml.log_likelihood_delta) <= 1e-3

    iqtree_support = observations["iqtree-jc69-ultrafast-bootstrap-support"]
    assert iqtree_support.status == "tolerance-match"
    assert len(iqtree_support.support_comparisons) == 1
    support_row = iqtree_support.support_comparisons[0]
    assert support_row.split_taxa == ["A", "B"]
    assert support_row.native_support_percent == 100.0
    assert support_row.wrapper_support_percent == 100.0
    assert support_row.support_percent_delta == 0.0

    fasttree = observations["fasttree-jc69-approximate-maximum-likelihood"]
    assert fasttree.status == "expected-model-assumption-difference"
    assert fasttree.topology_equal is True
    assert abs(fasttree.log_likelihood_delta) < 0.01

    unsupported = observations["raxml-small-maximum-likelihood-reference"]
    assert unsupported.status == "unsupported-case"
    assert unsupported.supported is False
    assert unsupported.wrapper_runtime_seconds is None


def test_maximum_likelihood_wrapper_comparison_writers_emit_governed_outputs(
    tmp_path: Path,
) -> None:
    report = benchmark_maximum_likelihood_wrapper_comparison()

    summary_path = write_maximum_likelihood_wrapper_comparison_summary_table(
        tmp_path / "summary.tsv",
        report,
    )
    observation_path = write_maximum_likelihood_wrapper_comparison_observation_table(
        tmp_path / "observations.tsv",
        report,
    )
    run_json_path = write_maximum_likelihood_wrapper_comparison_run_json(
        tmp_path / "run.json",
        report,
    )

    assert summary_path.read_text(encoding="utf-8").startswith(
        "case_count\tsupported_case_count\ttolerance_match_case_count\texpected_model_assumption_difference_case_count\tunsupported_case_count\tnative_bug_case_count\tblocking_case_count\tlimitations\n"
    )
    assert observation_path.read_text(encoding="utf-8").startswith(
        "case_id\twrapper_engine\tnative_surface\twrapper_surface\tcomparison_policy\tstatus\tsupported\tblocking\tnative_runtime_seconds\twrapper_runtime_seconds\truntime_delta_seconds\tnative_peak_memory_bytes\twrapper_peak_memory_bytes\tpeak_memory_delta_bytes\tnative_log_likelihood\twrapper_log_likelihood\tlog_likelihood_delta\ttopology_equal\trobinson_foulds_distance\tsupport_split_count\tmax_support_percent_delta\tnotes\n"
    )
    payload = json.loads(run_json_path.read_text(encoding="utf-8"))
    assert payload["case_count"] == 4
    assert payload["native_bug_case_count"] == 0
    assert payload["summary_rows"][0]["status"] == "tolerance-match"
    assert {row["case_id"] for row in payload["observations"]} == {
        "iqtree-jc69-maximum-likelihood",
        "iqtree-jc69-ultrafast-bootstrap-support",
        "fasttree-jc69-approximate-maximum-likelihood",
        "raxml-small-maximum-likelihood-reference",
    }


def test_public_runtime_exports_include_maximum_likelihood_wrapper_benchmark() -> None:
    assert (
        benchmark_api.benchmark_maximum_likelihood_wrapper_comparison
        is benchmark_maximum_likelihood_wrapper_comparison
    )
    assert (
        benchmark_api.write_maximum_likelihood_wrapper_comparison_summary_table
        is write_maximum_likelihood_wrapper_comparison_summary_table
    )
    assert (
        benchmark_api.write_maximum_likelihood_wrapper_comparison_observation_table
        is write_maximum_likelihood_wrapper_comparison_observation_table
    )
    assert (
        benchmark_api.write_maximum_likelihood_wrapper_comparison_run_json
        is write_maximum_likelihood_wrapper_comparison_run_json
    )
