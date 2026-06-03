from __future__ import annotations

import json

import pytest

import bijux_phylogenetics.benchmark as benchmark_api
from bijux_phylogenetics.benchmark import benchmark_large_dataset_stress_suite
from bijux_phylogenetics.command_line import main


def _observations_by_workload(report):
    return {row.workload: row for row in report.observations}


@pytest.mark.evaluation
@pytest.mark.stress_small
@pytest.mark.slow
def test_benchmark_large_dataset_stress_suite_small_tier_reports_all_workloads() -> (
    None
):
    report = benchmark_large_dataset_stress_suite(tier="small")

    assert report.tier == "small"
    assert len(report.observations) == 5
    observations = _observations_by_workload(report)
    assert set(observations) == {
        "large-alignment-inference",
        "multi-locus-supermatrix",
        "posterior-tree-set-consensus",
        "comparative-trait-contrasts",
        "tree-annotation-tables",
    }
    alignment = observations["large-alignment-inference"]
    assert alignment.sequence_count == 256
    assert alignment.alignment_length == 512
    assert alignment.output_row_count > 0
    assert alignment.peak_memory_bytes >= 0
    assert alignment.runtime_seconds >= 0.0
    assert alignment.memory_observation_kind in {
        "python-tracemalloc",
        "sampled-process-rss",
        "mixed",
    }

    supermatrix = observations["multi-locus-supermatrix"]
    assert supermatrix.locus_count == 3
    assert supermatrix.taxon_count == 256
    assert supermatrix.output_row_count > supermatrix.taxon_count

    tree_set = observations["posterior-tree-set-consensus"]
    assert tree_set.tree_count == 256
    assert tree_set.taxon_count == 64
    assert tree_set.output_row_count > 0

    comparative = observations["comparative-trait-contrasts"]
    assert comparative.taxon_count == 256
    assert comparative.output_row_count > 0

    tables = observations["tree-annotation-tables"]
    assert tables.taxon_count == 256
    assert tables.output_row_count > tables.taxon_count

    assert all(row.timeout_seconds == 30.0 for row in report.observations)
    assert report.limitations


@pytest.mark.evaluation
@pytest.mark.stress_small
@pytest.mark.slow
def test_cli_benchmark_stress_suite_reports_tier_and_observation_count(
    capsys: pytest.CaptureFixture[str],
) -> None:
    exit_code = main(["benchmark", "stress-suite", "--tier", "small", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["tier"] == "small"
    assert payload["metrics"]["observation_count"] == 5
    assert payload["data"]["tier"] == "small"


@pytest.mark.evaluation
@pytest.mark.stress_heavy
@pytest.mark.slow
def test_benchmark_large_dataset_stress_suite_heavy_tier_meets_large_input_thresholds() -> (
    None
):
    report = benchmark_large_dataset_stress_suite(tier="heavy")
    observations = _observations_by_workload(report)

    alignment = observations["large-alignment-inference"]
    assert alignment.sequence_count is not None and alignment.sequence_count >= 1000
    assert alignment.alignment_length is not None and alignment.alignment_length >= 1000

    tree_set = observations["posterior-tree-set-consensus"]
    assert tree_set.tree_count is not None and tree_set.tree_count >= 1000

    comparative = observations["comparative-trait-contrasts"]
    assert comparative.taxon_count is not None and comparative.taxon_count >= 500

    supermatrix = observations["multi-locus-supermatrix"]
    assert supermatrix.locus_count is not None and supermatrix.locus_count >= 4

    tables = observations["tree-annotation-tables"]
    assert tables.taxon_count is not None and tables.taxon_count >= 1000
    assert all(row.output_row_count > 0 for row in report.observations)


def test_public_runtime_exports_stress_suite_surface() -> None:
    assert (
        benchmark_api.benchmark_large_dataset_stress_suite
        is benchmark_large_dataset_stress_suite
    )
