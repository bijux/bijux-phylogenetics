from __future__ import annotations

from dataclasses import dataclass

from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .average_linkage import (
    AverageLinkageBuildResult,
    build_average_linkage_tree,
)
from .average_linkage import (
    AverageLinkageClusterHeightRow as WPGMAClusterHeightRow,
)
from .average_linkage import (
    AverageLinkageMergeRow as WPGMAMergeRow,
)


@dataclass(slots=True)
class WPGMABuildReport:
    """Deterministic merge and assumption report for one WPGMA tree build."""

    merge_history: list[WPGMAMergeRow]
    cluster_heights: list[WPGMAClusterHeightRow]
    ultrametric_compatible: bool
    assumption_warnings: list[str]


def _build_report(result: AverageLinkageBuildResult) -> WPGMABuildReport:
    return WPGMABuildReport(
        merge_history=result.merge_history,
        cluster_heights=result.cluster_heights,
        ultrametric_compatible=result.ultrametric_compatible,
        assumption_warnings=result.assumption_warnings,
    )


def build_wpgma_tree(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> tuple[PhyloTree, WPGMABuildReport]:
    """Build one deterministic rooted ultrametric WPGMA tree from a distance lookup."""
    result = build_average_linkage_tree(
        identifiers,
        distance_lookup,
        method_name="wpgma",
        cluster_weighting="equal",
    )
    return result.tree, _build_report(result)
