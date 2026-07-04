from __future__ import annotations

import json

import pytest

import bijux_phylogenetics.benchmark as benchmark_api
from bijux_phylogenetics.benchmark import benchmark_large_tree_set_scaling
from bijux_phylogenetics.command_line import main


def _workflows_by_name(report):
    return {workflow.workflow: workflow for workflow in report.workflows}


def test_benchmark_large_tree_set_scaling_reports_all_review_workflows() -> None:
    report = benchmark_large_tree_set_scaling(
        replicates=1,
        size_classes=[
            ("trees-8-taxa-6", 8, 6),
            ("trees-12-taxa-8", 12, 8),
        ],
    )

    assert report.replicates == 1
    assert report.tree_counts == [8, 12]
    assert report.tip_counts == [6, 8]
    assert report.limitations
    workflows = _workflows_by_name(report)
    assert set(workflows) == {
        "tree-set-consensus",
        "pairwise-rf-diversity",
        "topology-clustering",
        "uncertainty-summaries",
    }
    assert all(
        workflow.scaling_axis == "posterior_samples" for workflow in report.workflows
    )
    for workflow in report.workflows:
        assert [row.tree_count for row in workflow.observations] == [8, 12]
        assert [row.tip_count for row in workflow.observations] == [6, 8]
        assert all(row.pair_count > 0 for row in workflow.observations)
        assert all(row.runtime_seconds >= 0.0 for row in workflow.observations)
        assert all(row.peak_memory_bytes >= 0 for row in workflow.observations)
        assert workflow.notes


def test_benchmark_large_tree_set_scaling_rejects_invalid_size_classes() -> None:
    with pytest.raises(ValueError, match="at least two trees and two taxa"):
        benchmark_large_tree_set_scaling(
            replicates=1,
            size_classes=[("too-small", 1, 8)],
        )


def test_cli_benchmark_large_tree_set_scaling_reports_workflow_and_size_metrics(
    capsys: pytest.CaptureFixture[str],
) -> None:
    exit_code = main(
        [
            "benchmark",
            "large-tree-set-scaling",
            "--replicates",
            "1",
            "--tree-count",
            "8",
            "--tip-count",
            "6",
            "--tree-count",
            "12",
            "--tip-count",
            "8",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["replicates"] == 1
    assert payload["metrics"]["workflow_count"] == 4
    assert payload["metrics"]["max_tree_count"] == 12
    assert payload["metrics"]["max_tip_count"] == 8
    assert payload["metrics"]["observation_count"] == 8
    assert payload["data"]["tree_counts"] == [8, 12]
    assert payload["data"]["tip_counts"] == [6, 8]


def test_public_runtime_exports_large_tree_set_scaling_surface() -> None:
    assert (
        benchmark_api.benchmark_large_tree_set_scaling
        is benchmark_large_tree_set_scaling
    )
