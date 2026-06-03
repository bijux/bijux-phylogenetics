from __future__ import annotations

import json

import pytest

import bijux_phylogenetics.benchmark as benchmark_api
from bijux_phylogenetics.benchmark import benchmark_large_alignment_scaling
from bijux_phylogenetics.command_line import main


def _workflows_by_name(report):
    return {workflow.workflow: workflow for workflow in report.workflows}


@pytest.mark.slow
def test_benchmark_large_alignment_scaling_reports_all_review_workflows() -> None:
    report = benchmark_large_alignment_scaling(
        replicates=1,
        size_classes=[
            ("sequences-4-sites-16", 4, 16),
            ("sequences-6-sites-24", 6, 24),
        ],
    )

    assert report.replicates == 1
    assert report.sequence_counts == [4, 6]
    assert report.alignment_lengths == [16, 24]
    assert report.limitations
    workflows = _workflows_by_name(report)
    assert set(workflows) == {
        "alignment-diagnostics",
        "alignment-trimming",
        "distance-analysis",
        "alignment-readiness",
    }
    assert all(
        workflow.scaling_axis == "aligned_sites" for workflow in report.workflows
    )
    for workflow in report.workflows:
        assert [row.sequence_count for row in workflow.observations] == [4, 6]
        assert [row.alignment_length for row in workflow.observations] == [16, 24]
        assert all(row.aligned_site_count > 0 for row in workflow.observations)
        assert all(row.runtime_seconds >= 0.0 for row in workflow.observations)
        assert all(row.peak_memory_bytes >= 0 for row in workflow.observations)
        assert workflow.notes


def test_benchmark_large_alignment_scaling_rejects_invalid_size_classes() -> None:
    with pytest.raises(ValueError, match="at least two sequences and two sites"):
        benchmark_large_alignment_scaling(
            replicates=1,
            size_classes=[("too-small", 1, 8)],
        )


@pytest.mark.slow
def test_cli_benchmark_large_alignment_scaling_reports_workflow_and_size_metrics(
    capsys: pytest.CaptureFixture[str],
) -> None:
    exit_code = main(
        [
            "benchmark",
            "large-alignment-scaling",
            "--replicates",
            "1",
            "--sequence-count",
            "4",
            "--alignment-length",
            "16",
            "--sequence-count",
            "6",
            "--alignment-length",
            "24",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["replicates"] == 1
    assert payload["metrics"]["workflow_count"] == 4
    assert payload["metrics"]["max_sequence_count"] == 6
    assert payload["metrics"]["max_alignment_length"] == 24
    assert payload["metrics"]["observation_count"] == 8
    assert payload["data"]["sequence_counts"] == [4, 6]
    assert payload["data"]["alignment_lengths"] == [16, 24]


def test_public_runtime_exports_large_alignment_scaling_surface() -> None:
    assert (
        benchmark_api.benchmark_large_alignment_scaling
        is benchmark_large_alignment_scaling
    )
