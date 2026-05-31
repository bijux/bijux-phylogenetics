from __future__ import annotations

from dataclasses import dataclass

from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .extremal_linkage import build_extremal_linkage_tree
from .linkage_common import (
    LinkageClusterHeightRow as CompleteLinkageClusterHeightRow,
)
from .linkage_common import (
    LinkageMergeRow as CompleteLinkageMergeRow,
)


@dataclass(slots=True)
class CompleteLinkageBuildReport:
    """Deterministic merge report for one complete-linkage tree build."""

    merge_history: list[CompleteLinkageMergeRow]
    cluster_heights: list[CompleteLinkageClusterHeightRow]


def build_complete_linkage_tree(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> tuple[PhyloTree, CompleteLinkageBuildReport]:
    """Build one deterministic rooted complete-linkage clustering tree from a distance lookup."""
    result = build_extremal_linkage_tree(
        identifiers,
        distance_lookup,
        method_name="complete-linkage",
        linkage_rule="maximum",
    )
    return result.tree, CompleteLinkageBuildReport(
        merge_history=result.merge_history,
        cluster_heights=result.cluster_heights,
    )
