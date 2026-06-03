from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .agglomerative import build_agglomerative_clustering_tree
from .linkage_common import (
    LinkageClusterHeightRow,
    LinkageMergeRow,
)

ExtremalLinkageRule = Literal["minimum", "maximum"]


@dataclass(slots=True)
class ExtremalLinkageBuildResult:
    """Shared deterministic result for one min/max linkage clustering tree build."""

    tree: PhyloTree
    merge_history: list[LinkageMergeRow]
    cluster_heights: list[LinkageClusterHeightRow]


def build_extremal_linkage_tree(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
    *,
    method_name: str,
    linkage_rule: ExtremalLinkageRule,
) -> ExtremalLinkageBuildResult:
    """Build one rooted linkage tree whose updates use min or max intercluster distance."""
    if linkage_rule not in {"minimum", "maximum"}:
        raise ValueError(f"unsupported extremal linkage rule: {linkage_rule}")
    result = build_agglomerative_clustering_tree(
        identifiers,
        distance_lookup,
        method_name=method_name,
        update_rule=linkage_rule,
        assess_ultrametric_assumption=False,
    )
    return ExtremalLinkageBuildResult(
        tree=result.tree,
        merge_history=result.merge_history,
        cluster_heights=result.cluster_heights,
    )
