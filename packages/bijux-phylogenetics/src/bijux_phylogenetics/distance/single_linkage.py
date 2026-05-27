from __future__ import annotations

from dataclasses import dataclass

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .linkage_common import (
    ActiveCluster,
    ClusterKey,
    LinkageClusterHeightRow as SingleLinkageClusterHeightRow,
    LinkageMergeRow as SingleLinkageMergeRow,
    choose_linkage_merge_pair,
    cluster_label,
    cluster_pair_key,
    normalized_branch_length,
    validate_linkage_distance_lookup,
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
    validate_linkage_distance_lookup(identifiers, distance_lookup)
    active_clusters = {
        (identifier,): ActiveCluster(
            key=(identifier,),
            node=TreeNode(name=identifier),
            height=0.0,
            size=1,
        )
        for identifier in identifiers
    }
    active_distances: dict[tuple[ClusterKey, ClusterKey], float] = {}
    for left_index, left_identifier in enumerate(identifiers):
        for right_identifier in identifiers[left_index + 1 :]:
            active_distances[
                cluster_pair_key((left_identifier,), (right_identifier,))
            ] = distance_lookup[(left_identifier, right_identifier)]

    merge_history: list[SingleLinkageMergeRow] = []
    cluster_heights: list[SingleLinkageClusterHeightRow] = []
    merge_index = 1
    while len(active_clusters) > 1:
        active_keys = sorted(active_clusters)
        left_key, right_key = choose_linkage_merge_pair(
            active_distances,
            active_keys,
            method_name="single-linkage",
        )
        left_cluster = active_clusters[left_key]
        right_cluster = active_clusters[right_key]
        pair_distance = active_distances[cluster_pair_key(left_key, right_key)]
        merge_height = pair_distance / 2.0
        left_cluster.node.branch_length = normalized_branch_length(
            merge_height - left_cluster.height
        )
        right_cluster.node.branch_length = normalized_branch_length(
            merge_height - right_cluster.height
        )
        merged_key = tuple(sorted(left_key + right_key))
        merged_node = TreeNode(
            name=f"Inner{merge_index}",
            children=[left_cluster.node, right_cluster.node],
        )
        merged_cluster = ActiveCluster(
            key=merged_key,
            node=merged_node,
            height=merge_height,
            size=left_cluster.size + right_cluster.size,
        )
        merge_history.append(
            SingleLinkageMergeRow(
                merge_index=merge_index,
                left_cluster=cluster_label(left_key),
                right_cluster=cluster_label(right_key),
                left_cluster_size=left_cluster.size,
                right_cluster_size=right_cluster.size,
                pair_distance=pair_distance,
                merge_height=merge_height,
                resulting_cluster=cluster_label(merged_key),
                resulting_cluster_size=merged_cluster.size,
            )
        )
        cluster_heights.append(
            SingleLinkageClusterHeightRow(
                merge_index=merge_index,
                cluster=cluster_label(merged_key),
                height=merge_height,
            )
        )
        updated_distances = {
            pair_key: value
            for pair_key, value in active_distances.items()
            if left_key not in pair_key and right_key not in pair_key
        }
        for other_key in active_keys:
            if other_key in {left_key, right_key}:
                continue
            updated_distances[cluster_pair_key(merged_key, other_key)] = min(
                active_distances[cluster_pair_key(left_key, other_key)],
                active_distances[cluster_pair_key(right_key, other_key)],
            )
        del active_clusters[left_key]
        del active_clusters[right_key]
        active_clusters[merged_key] = merged_cluster
        active_distances.clear()
        active_distances.update(updated_distances)
        merge_index += 1

    root = next(iter(active_clusters.values())).node
    return PhyloTree(root=root, source_format="newick", rooted=True), SingleLinkageBuildReport(
        merge_history=merge_history,
        cluster_heights=cluster_heights,
    )
