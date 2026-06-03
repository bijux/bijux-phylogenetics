from __future__ import annotations

from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_clade_id,
    robinson_foulds_metrics,
    split_sort_key,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .missing_distance_policy import apply_missing_distance_policy
from .models import ImportedDistanceEntry, MissingDistancePolicy
from .patristic_residuals import compute_patristic_residual_diagnostics
from .shared import _build_distance_tree_from_lookup, _pair_key

ROUND_DIGITS = 12


def defined_distance_lookup_from_entries(
    entries: list[ImportedDistanceEntry],
) -> dict[tuple[str, str], float]:
    """Reduce imported matrix entries to one defined undirected distance lookup."""
    defined_lookup: dict[tuple[str, str], float] = {}
    for entry in entries:
        if entry.left_identifier == entry.right_identifier:
            continue
        defined_lookup.setdefault(
            _pair_key(entry.left_identifier, entry.right_identifier),
            float(entry.distance),
        )
    return defined_lookup


def subset_defined_lookup(
    identifiers: list[str],
    defined_lookup: dict[tuple[str, str], float],
) -> dict[tuple[str, str], float]:
    """Keep only distance pairs whose taxa remain in the requested subset."""
    retained = set(identifiers)
    return {
        pair: distance
        for pair, distance in defined_lookup.items()
        if pair[0] in retained and pair[1] in retained
    }


def prune_tree_to_identifiers_or_raise(
    tree_path,
    identifiers: list[str],
    *,
    missing_taxa_message_prefix: str,
) -> PhyloTree:
    """Prune one tree to the requested taxa and fail if any requested taxon is absent."""
    pruned_tree, pruning_report = prune_tree_to_requested_taxa(
        tree_path,
        identifiers,
    )
    if pruning_report.absent_requested_taxa:
        missing = ", ".join(pruning_report.absent_requested_taxa)
        raise ValueError(missing_taxa_message_prefix + missing)
    return pruned_tree


def rooted_rf_metrics(
    tree: PhyloTree,
    reference_tree: PhyloTree,
) -> tuple[int, float]:
    """Compute rooted RF distance over the exact shared taxon scope."""
    shared_taxa = set(tree.tip_names) & set(reference_tree.tip_names)
    metrics = robinson_foulds_metrics(
        tree,
        reference_tree,
        shared_taxa,
        rf_mode="rooted",
    )
    return metrics.distance, round(metrics.normalized_distance, ROUND_DIGITS)


def rooted_rf_clade_differences(
    tree: PhyloTree,
    reference_tree: PhyloTree,
) -> tuple[list[str], list[str], list[str], int, float]:
    """Report rooted RF clade differences using stable clade identifiers."""
    shared_taxa = set(tree.tip_names) & set(reference_tree.tip_names)
    metrics = robinson_foulds_metrics(
        tree,
        reference_tree,
        shared_taxa,
        rf_mode="rooted",
    )
    reference_only_clades = [
        canonical_clade_id(signature)
        for signature in sorted(metrics.right_only_signatures, key=split_sort_key)
    ]
    rebuilt_only_clades = [
        canonical_clade_id(signature)
        for signature in sorted(metrics.left_only_signatures, key=split_sort_key)
    ]
    affected_clades = sorted(
        {*reference_only_clades, *rebuilt_only_clades},
        key=lambda clade_id: (len(clade_id.split("|")), clade_id),
    )
    return (
        reference_only_clades,
        rebuilt_only_clades,
        affected_clades,
        metrics.distance,
        round(metrics.normalized_distance, ROUND_DIGITS),
    )


def resolve_tree_and_residuals(
    identifiers: list[str],
    defined_lookup: dict[tuple[str, str], float],
    *,
    method: str,
    missing_distance_policy: MissingDistancePolicy,
) -> tuple[PhyloTree, float]:
    """Resolve a distance tree and its patristic residual RSS over one taxon set."""
    resolved_lookup, _policy_report = apply_missing_distance_policy(
        identifiers,
        defined_lookup,
        policy=missing_distance_policy,
    )
    tree = _build_distance_tree_from_lookup(
        identifiers,
        resolved_lookup,
        method=method,
    )
    residuals = compute_patristic_residual_diagnostics(
        tree,
        identifiers,
        resolved_lookup,
    )
    return tree, residuals.residual_sum_squares


def resolve_distance_lookup(
    identifiers: list[str],
    defined_lookup: dict[tuple[str, str], float],
    *,
    missing_distance_policy: MissingDistancePolicy,
) -> dict[tuple[str, str], float]:
    """Resolve one possibly incomplete imported lookup under the chosen missing-distance policy."""
    resolved_lookup, _policy_report = apply_missing_distance_policy(
        identifiers,
        defined_lookup,
        policy=missing_distance_policy,
    )
    return resolved_lookup


def compute_tree_residual_sum_squares(
    tree: PhyloTree,
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> float:
    """Compute one tree's patristic residual RSS on an already-resolved distance lookup."""
    residuals = compute_patristic_residual_diagnostics(
        tree,
        identifiers,
        distance_lookup,
    )
    return residuals.residual_sum_squares
