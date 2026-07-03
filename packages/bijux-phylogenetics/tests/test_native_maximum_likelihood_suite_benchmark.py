from __future__ import annotations

import pytest

import bijux_phylogenetics.benchmark as benchmark_api
from bijux_phylogenetics.benchmark.likelihood_tree_inference import (
    benchmark_native_maximum_likelihood_suite,
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
