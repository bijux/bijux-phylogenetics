from __future__ import annotations

from dataclasses import dataclass

from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .extremal_linkage import build_extremal_linkage_tree
from .linkage_common import (
    LinkageClusterHeightRow as SingleLinkageClusterHeightRow,
)
from .linkage_common import (
    LinkageMergeRow as SingleLinkageMergeRow,
)


@dataclass(slots=True)
class SingleLinkageBuildReport:
    """Deterministic merge report for one single-linkage tree build."""

    merge_history: list[SingleLinkageMergeRow]
    cluster_heights: list[SingleLinkageClusterHeightRow]


def build_single_linkage_tree(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> tuple[PhyloTree, SingleLinkageBuildReport]:
    """Build one deterministic rooted single-linkage clustering tree from a distance lookup."""
    result = build_extremal_linkage_tree(
        identifiers,
        distance_lookup,
        method_name="single-linkage",
        linkage_rule="minimum",
    )
    return result.tree, SingleLinkageBuildReport(
        merge_history=result.merge_history,
        cluster_heights=result.cluster_heights,
    )
