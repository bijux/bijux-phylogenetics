from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics.benchmark as benchmark_api
from bijux_phylogenetics.benchmark.likelihood_tree_inference import (
    benchmark_native_maximum_likelihood_speed,
    write_native_maximum_likelihood_speed_observation_table,
    write_native_maximum_likelihood_speed_run_json,
    write_native_maximum_likelihood_speed_summary_table,
    write_native_maximum_likelihood_speed_trend_table,
)


pytestmark = pytest.mark.slow


def test_benchmark_native_maximum_likelihood_speed_reports_lanes_and_trends() -> None:
    report = benchmark_native_maximum_likelihood_speed(
        replicates=1,
        size_classes=[("tiny", 4, 12), ("medium", 4, 16)],
        model_names=["jc69"],
        search_methods=["nni", "tbr"],
    )

    assert report.replicates == 1
    assert report.size_class_labels == ["tiny", "medium"]
    assert report.model_names == ["jc69"]
    assert report.search_methods == ["nni", "tbr"]
    assert len(report.observations) == 4
    assert len(report.trend_rows) == 2
    assert report.limitations

    observations = {
        (row.size_class, row.search_method): row for row in report.observations
    }
    assert set(observations) == {
        ("tiny", "nni"),
        ("tiny", "tbr"),
        ("medium", "nni"),
        ("medium", "tbr"),
    }
    assert all(row.model_name == "jc69" for row in report.observations)
    assert all(row.taxon_count == 4 for row in report.observations)
    assert {row.site_count for row in report.observations} == {12, 16}
    assert all(row.runtime_seconds >= 0.0 for row in report.observations)
    assert all(row.peak_memory_bytes >= 0 for row in report.observations)
    assert all(row.selected_model_name == "JC69" for row in report.observations)
    assert observations[("tiny", "nni")].search_algorithm == (
        "nucleotide-likelihood-nni-search"
    )
    assert observations[("tiny", "tbr")].search_algorithm == (
        "nucleotide-likelihood-tbr-search"
    )
    assert {row.search_method for row in report.trend_rows} == {"nni", "tbr"}
    assert all(row.baseline_size_class == "tiny" for row in report.trend_rows)
    assert all(row.comparison_size_class == "medium" for row in report.trend_rows)
    assert all(row.runtime_ratio > 0.0 for row in report.trend_rows)


def test_native_maximum_likelihood_speed_rejects_unsupported_fixed_model() -> None:
    with pytest.raises(ValueError, match="supports fixed models"):
        benchmark_native_maximum_likelihood_speed(
            replicates=1,
            size_classes=[("tiny", 4, 12), ("medium", 4, 16)],
            model_names=["gtr"],
            search_methods=["nni"],
        )


def test_native_maximum_likelihood_speed_writers_emit_governed_outputs(
    tmp_path: Path,
) -> None:
    report = benchmark_native_maximum_likelihood_speed(
        replicates=1,
        size_classes=[("tiny", 4, 12), ("medium", 4, 16)],
        model_names=["jc69"],
        search_methods=["nni"],
    )

    summary_path = write_native_maximum_likelihood_speed_summary_table(
        tmp_path / "summary.tsv",
        report,
    )
    observation_path = write_native_maximum_likelihood_speed_observation_table(
        tmp_path / "observations.tsv",
        report,
    )
    trend_path = write_native_maximum_likelihood_speed_trend_table(
        tmp_path / "trends.tsv",
        report,
    )
    run_json_path = write_native_maximum_likelihood_speed_run_json(
        tmp_path / "run.json",
        report,
    )

    assert summary_path.read_text(encoding="utf-8").startswith(
        "replicates\tobservation_count\ttrend_row_count\tsize_class_labels\tmodel_names\tsearch_methods\tlimitations\n"
    )
    assert observation_path.read_text(encoding="utf-8").startswith(
        "label\tsize_class\tmodel_name\tsearch_method\ttaxon_count\tsite_count\truntime_seconds\tpeak_memory_bytes\tselected_model_name\tbest_run_source_label\tsearch_algorithm\tfinal_log_likelihood\n"
    )
    assert trend_path.read_text(encoding="utf-8").startswith(
        "model_name\tsearch_method\tbaseline_size_class\tcomparison_size_class\tbaseline_taxon_count\tcomparison_taxon_count\tbaseline_site_count\tcomparison_site_count\tbaseline_runtime_seconds\tcomparison_runtime_seconds\truntime_delta_seconds\truntime_ratio\n"
    )
    payload = json.loads(run_json_path.read_text(encoding="utf-8"))
    assert payload["replicates"] == 1
    assert payload["size_class_labels"] == ["tiny", "medium"]
    assert payload["model_names"] == ["jc69"]
    assert payload["search_methods"] == ["nni"]
    assert len(payload["observations"]) == 2
    assert len(payload["trend_rows"]) == 1


def test_public_runtime_exports_include_native_maximum_likelihood_speed_benchmark() -> (
    None
):
    assert (
        benchmark_api.benchmark_native_maximum_likelihood_speed
        is benchmark_native_maximum_likelihood_speed
    )
    assert (
        benchmark_api.write_native_maximum_likelihood_speed_summary_table
        is write_native_maximum_likelihood_speed_summary_table
    )
    assert (
        benchmark_api.write_native_maximum_likelihood_speed_observation_table
        is write_native_maximum_likelihood_speed_observation_table
    )
    assert (
        benchmark_api.write_native_maximum_likelihood_speed_trend_table
        is write_native_maximum_likelihood_speed_trend_table
    )
    assert (
        benchmark_api.write_native_maximum_likelihood_speed_run_json
        is write_native_maximum_likelihood_speed_run_json
    )
