from __future__ import annotations

from dataclasses import dataclass

from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .average_linkage import (
    AverageLinkageBuildResult,
    build_average_linkage_tree,
)
from .average_linkage import (
    AverageLinkageClusterHeightRow as UPGMAClusterHeightRow,
)
from .average_linkage import (
    AverageLinkageMergeRow as UPGMAMergeRow,
)


@dataclass(slots=True)
class UPGMABuildReport:
    """Deterministic merge and assumption report for one UPGMA tree build."""

    merge_history: list[UPGMAMergeRow]
    cluster_heights: list[UPGMAClusterHeightRow]
    ultrametric_compatible: bool
    assumption_warnings: list[str]


def _build_report(result: AverageLinkageBuildResult) -> UPGMABuildReport:
    return UPGMABuildReport(
        merge_history=result.merge_history,
        cluster_heights=result.cluster_heights,
        ultrametric_compatible=result.ultrametric_compatible,
        assumption_warnings=result.assumption_warnings,
    )


def build_upgma_tree(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> tuple[PhyloTree, UPGMABuildReport]:
    """Build one deterministic rooted ultrametric UPGMA tree from a distance lookup."""
    result = build_average_linkage_tree(
        identifiers,
        distance_lookup,
        method_name="upgma",
        cluster_weighting="taxon-count",
    )
    return result.tree, _build_report(result)
