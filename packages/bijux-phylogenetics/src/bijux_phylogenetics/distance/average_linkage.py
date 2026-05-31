from __future__ import annotations

from dataclasses import dataclass

from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .agglomerative import (
    build_agglomerative_clustering_tree,
)
from .linkage_common import (
    LinkageClusterHeightRow as AverageLinkageClusterHeightRow,
)
from .linkage_common import (
    LinkageMergeRow as AverageLinkageMergeRow,
)
from .linkage_common import (
    validate_linkage_distance_lookup,
)


@dataclass(slots=True)
class AverageLinkageBuildResult:
    """Shared deterministic result for one rooted average-linkage tree build."""

    tree: PhyloTree
    merge_history: list[AverageLinkageMergeRow]
    cluster_heights: list[AverageLinkageClusterHeightRow]
    ultrametric_compatible: bool
    assumption_warnings: list[str]


def validate_average_linkage_distance_lookup(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> None:
    """Validate one full symmetric distance lookup before average-linkage clustering."""
    validate_linkage_distance_lookup(identifiers, distance_lookup)


def build_average_linkage_tree(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
    *,
    method_name: str,
    cluster_weighting: str,
) -> AverageLinkageBuildResult:
    """Build one rooted average-linkage tree with explicit cluster-weight policy."""
    if cluster_weighting not in {"taxon-count", "equal"}:
        raise ValueError(f"unsupported average-linkage weighting: {cluster_weighting}")
    result = build_agglomerative_clustering_tree(
        identifiers,
        distance_lookup,
        method_name=method_name,
        update_rule=cluster_weighting,
        assess_ultrametric_assumption=True,
    )
    return AverageLinkageBuildResult(
        tree=result.tree,
        merge_history=result.merge_history,
        cluster_heights=result.cluster_heights,
        ultrametric_compatible=result.ultrametric_compatible,
        assumption_warnings=result.assumption_warnings,
    )
