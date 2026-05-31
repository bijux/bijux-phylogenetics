from __future__ import annotations

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .distance_joining import (
    TIE_TOLERANCE,
    ActiveDistanceJoinCluster,
    ClusterKey,
    build_three_taxon_join_tree,
    build_two_taxon_join_tree,
    choose_join_pair,
    cluster_pair_key,
    distance_sums,
    normalized_branch_length,
    validate_distance_lookup,
)


def _normalized_reduced_distance(value: float) -> float:
    if abs(value) <= TIE_TOLERANCE:
        return 0.0
    return value


def _join_clusters(
    active_clusters: dict[ClusterKey, ActiveDistanceJoinCluster],
    active_distances: dict[tuple[ClusterKey, ClusterKey], float],
    active_variances: dict[tuple[ClusterKey, ClusterKey], float],
    *,
    join_counter: int,
) -> int:
    active_keys = sorted(active_clusters)
    left_key, right_key = choose_join_pair(
        active_distances,
        active_keys,
        method_name="bionj",
    )
    row_sums = distance_sums(active_distances, active_keys)
    active_count = len(active_keys)
    pair_key = cluster_pair_key(left_key, right_key)
    pair_distance = active_distances[pair_key]
    left_length = 0.5 * pair_distance + (
        (row_sums[left_key] - row_sums[right_key]) / (2 * (active_count - 2))
    )
    right_length = pair_distance - left_length
    active_clusters[left_key].node.branch_length = normalized_branch_length(left_length)
    active_clusters[right_key].node.branch_length = normalized_branch_length(
        right_length
    )
    pair_variance = active_variances[pair_key]
    if abs(pair_variance) <= TIE_TOLERANCE:
        lambda_weight = 0.5
    else:
        variance_delta = sum(
            active_variances[cluster_pair_key(right_key, other_key)]
            - active_variances[cluster_pair_key(left_key, other_key)]
            for other_key in active_keys
            if other_key not in {left_key, right_key}
        )
        lambda_weight = 0.5 + (
            variance_delta / (2 * (active_count - 2) * pair_variance)
        )
        lambda_weight = max(0.0, min(1.0, lambda_weight))

    merged_key = tuple(sorted(left_key + right_key))
    merged_node = TreeNode(
        name=f"Inner{join_counter}",
        children=[
            active_clusters[left_key].node,
            active_clusters[right_key].node,
        ],
    )
    updated_distances = {
        key: value
        for key, value in active_distances.items()
        if left_key not in key and right_key not in key
    }
    updated_variances = {
        key: value
        for key, value in active_variances.items()
        if left_key not in key and right_key not in key
    }
    for other_key in active_keys:
        if other_key in {left_key, right_key}:
            continue
        left_distance = active_distances[cluster_pair_key(left_key, other_key)]
        right_distance = active_distances[cluster_pair_key(right_key, other_key)]
        left_variance = active_variances[cluster_pair_key(left_key, other_key)]
        right_variance = active_variances[cluster_pair_key(right_key, other_key)]
        updated_distances[cluster_pair_key(merged_key, other_key)] = (
            _normalized_reduced_distance(
                lambda_weight * (left_distance - left_length)
                + ((1.0 - lambda_weight) * (right_distance - right_length))
            )
        )
        updated_variances[cluster_pair_key(merged_key, other_key)] = (
            lambda_weight * left_variance
            + ((1.0 - lambda_weight) * right_variance)
            - (lambda_weight * (1.0 - lambda_weight) * pair_variance)
        )
    del active_clusters[left_key]
    del active_clusters[right_key]
    active_clusters[merged_key] = ActiveDistanceJoinCluster(
        key=merged_key,
        node=merged_node,
    )
    active_distances.clear()
    active_distances.update(updated_distances)
    active_variances.clear()
    active_variances.update(updated_variances)
    return join_counter + 1


def build_bionj_tree(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> PhyloTree:
    """Build one deterministic BIONJ tree from a full distance lookup."""
    validate_distance_lookup(identifiers, distance_lookup)
    active_clusters = {
        (identifier,): ActiveDistanceJoinCluster(
            key=(identifier,),
            node=TreeNode(name=identifier),
        )
        for identifier in identifiers
    }
    active_distances: dict[tuple[ClusterKey, ClusterKey], float] = {}
    active_variances: dict[tuple[ClusterKey, ClusterKey], float] = {}
    for left_index, left_identifier in enumerate(identifiers):
        for right_identifier in identifiers[left_index + 1 :]:
            pair_key = cluster_pair_key((left_identifier,), (right_identifier,))
            pair_distance = distance_lookup[(left_identifier, right_identifier)]
            active_distances[pair_key] = pair_distance
            active_variances[pair_key] = pair_distance

    join_counter = 1
    while len(active_clusters) > 3:
        join_counter = _join_clusters(
            active_clusters,
            active_distances,
            active_variances,
            join_counter=join_counter,
        )
    if len(active_clusters) == 3:
        return build_three_taxon_join_tree(
            active_clusters,
            active_distances,
            join_counter=join_counter,
        )
    return build_two_taxon_join_tree(
        active_clusters,
        active_distances,
        join_counter=join_counter,
    )
