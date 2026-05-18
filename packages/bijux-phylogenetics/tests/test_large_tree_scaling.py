from __future__ import annotations

import json

import pytest

import bijux_phylogenetics.benchmark as benchmark_api
from bijux_phylogenetics.benchmark import benchmark_large_tree_scaling
from bijux_phylogenetics.command_line import main


def _workflows_by_name(report):
    return {workflow.workflow: workflow for workflow in report.workflows}


def test_benchmark_large_tree_scaling_reports_all_review_workflows() -> None:
    report = benchmark_large_tree_scaling(replicates=1, tip_counts=[8, 16])

    assert report.replicates == 1
    assert report.tip_counts == [8, 16]
    assert report.limitations
    workflows = _workflows_by_name(report)
    assert set(workflows) == {
        "tree-validation",
        "tree-comparison",
        "tree-rendering",
        "tree-reporting",
    }
    assert all(workflow.scaling_axis == "taxa" for workflow in report.workflows)
    for workflow in report.workflows:
        assert [row.item_count for row in workflow.observations] == [8, 16]
        assert all(row.runtime_seconds >= 0.0 for row in workflow.observations)
        assert all(row.peak_memory_bytes >= 0 for row in workflow.observations)
        assert workflow.notes


def test_benchmark_large_tree_scaling_rejects_invalid_tip_counts() -> None:
    with pytest.raises(ValueError, match="at least 2"):
        benchmark_large_tree_scaling(replicates=1, tip_counts=[1])


def test_cli_benchmark_large_tree_scaling_reports_workflow_and_tip_metrics(
    capsys: pytest.CaptureFixture[str],
) -> None:
    exit_code = main(
        [
            "benchmark",
            "large-tree-scaling",
            "--replicates",
            "1",
            "--tip-count",
            "8",
            "--tip-count",
            "16",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["replicates"] == 1
    assert payload["metrics"]["workflow_count"] == 4
    assert payload["metrics"]["max_tip_count"] == 16
    assert payload["metrics"]["observation_count"] == 8
    assert payload["data"]["tip_counts"] == [8, 16]


def test_public_runtime_exports_large_tree_scaling_surface() -> None:
    assert benchmark_api.benchmark_large_tree_scaling is benchmark_large_tree_scaling
