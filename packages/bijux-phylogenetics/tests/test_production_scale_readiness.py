from __future__ import annotations

import pytest

import bijux_phylogenetics.validation as validation_api
from bijux_phylogenetics.validation import (
    ProductionScaleThreshold,
    build_production_scale_readiness_report,
)


def _entries_by_workflow(report):
    return {entry.workflow: entry for entry in report.entries}


@pytest.mark.slow
def test_build_production_scale_readiness_report_classifies_workflow_specific_scales() -> (
    None
):
    report = build_production_scale_readiness_report(
        replicates=1,
        tree_tip_counts=[8, 16],
        alignment_size_classes=[
            ("sequences-4-sites-16", 4, 16),
            ("sequences-6-sites-24", 6, 24),
        ],
        tree_set_size_classes=[
            ("trees-8-taxa-6", 8, 6),
            ("trees-12-taxa-8", 12, 8),
        ],
        stress_tiers=["small"],
    )

    assert report.goal_id == 225
    assert [threshold.scale for threshold in report.scale_definitions] == [
        "small",
        "medium",
        "large",
        "hpc",
    ]
    entries = _entries_by_workflow(report)

    tree_validation = entries["tree-validation"]
    assert tree_validation.scale_dimensions == ["taxa"]
    assert tree_validation.highest_ready_scale == "small"

    alignment = entries["alignment-diagnostics"]
    assert alignment.scale_dimensions == ["taxa", "sites"]
    assert alignment.highest_ready_scale == "below-small"

    posterior_tree_set = entries["posterior-tree-set-consensus"]
    assert posterior_tree_set.scale_dimensions == [
        "taxa",
        "tree_count",
        "posterior_size",
    ]
    assert posterior_tree_set.highest_ready_scale == "small"

    tree_annotation = entries["tree-annotation-tables"]
    assert tree_annotation.scale_dimensions == ["taxa"]
    assert tree_annotation.highest_ready_scale == "medium"

    assert report.limitations


def test_build_production_scale_readiness_report_rejects_empty_scale_definitions() -> (
    None
):
    with pytest.raises(
        ValueError, match="scale_definitions must contain at least one scale threshold"
    ):
        build_production_scale_readiness_report(
            replicates=1,
            tree_tip_counts=[8, 16],
            alignment_size_classes=[("sequences-4-sites-16", 4, 16)],
            tree_set_size_classes=[("trees-8-taxa-6", 8, 6)],
            stress_tiers=["small"],
            scale_definitions=[],
        )


def test_public_runtime_exports_production_scale_readiness_surface() -> None:
    assert (
        validation_api.build_production_scale_readiness_report
        is build_production_scale_readiness_report
    )
    assert validation_api.ProductionScaleThreshold is ProductionScaleThreshold
