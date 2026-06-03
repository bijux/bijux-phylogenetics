from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics.benchmark as benchmark_api
from bijux_phylogenetics.benchmark import benchmark_large_tree_scaling
from bijux_phylogenetics.benchmark._fixtures import (
    build_balanced_tree,
    interleaved_taxa,
    write_named_balanced_tree,
)
from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.compare.topology import compare_tree_paths
from bijux_phylogenetics.io.newick import write_newick


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


@pytest.mark.slow
def test_large_tree_scaling_comparison_shape_supports_2048_taxa(
    tmp_path: Path,
) -> None:
    left_path = write_newick(
        tmp_path / "large-tree-balanced-2048.nwk",
        build_balanced_tree(2048, prefix="LargeTaxon"),
    )
    right_path = write_named_balanced_tree(
        tmp_path / "large-tree-permuted-balanced-2048.nwk",
        interleaved_taxa(2048, prefix="LargeTaxon"),
    )

    report = compare_tree_paths(left_path, right_path)

    assert len(report.shared_taxa) == 2048
    assert report.robinson_foulds_distance > 0
