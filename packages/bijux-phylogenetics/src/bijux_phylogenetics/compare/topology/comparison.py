from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.diagnostics.validation.structure import _load_tree
from bijux_phylogenetics.phylo.topology.clades import robinson_foulds_metrics
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .models import (
    InMemoryTopologyComparison,
    RobinsonFouldsComparisonReport,
    RobinsonFouldsMode,
    TaxonOverlapPolicy,
    TreeComparisonReport,
)


def _validate_rf_mode(rf_mode: RobinsonFouldsMode) -> None:
    if rf_mode not in {"rooted", "unrooted"}:
        raise ValueError(
            f"rf_mode must be one of {{'rooted', 'unrooted'}}, got {rf_mode!r}"
        )


def _validate_taxon_overlap_policy(taxon_overlap_policy: TaxonOverlapPolicy) -> None:
    if taxon_overlap_policy not in {"prune-to-shared", "require-identical"}:
        raise ValueError(
            "taxon_overlap_policy must be one of "
            "{'prune-to-shared', 'require-identical'}, "
            f"got {taxon_overlap_policy!r}"
        )


def _resolve_shared_taxa(
    left_taxa: set[str],
    right_taxa: set[str],
    *,
    taxon_overlap_policy: TaxonOverlapPolicy,
) -> tuple[set[str], list[str], list[str]]:
    _validate_taxon_overlap_policy(taxon_overlap_policy)
    shared_taxa = left_taxa & right_taxa
    left_only_taxa = sorted(left_taxa - right_taxa)
    right_only_taxa = sorted(right_taxa - left_taxa)
    if len(shared_taxa) < 2:
        raise ValueError("tree comparison requires at least two shared taxa")
    if taxon_overlap_policy == "require-identical" and (
        left_only_taxa or right_only_taxa
    ):
        raise ValueError(
            "tree comparison requires identical taxon sets when "
            "taxon_overlap_policy='require-identical'"
        )
    return shared_taxa, left_only_taxa, right_only_taxa


def _robinson_foulds_metrics(
    left: PhyloTree,
    right: PhyloTree,
    shared_taxa: set[str],
    *,
    rf_mode: RobinsonFouldsMode,
) -> tuple[int, int, int, float]:
    _validate_rf_mode(rf_mode)
    metrics = robinson_foulds_metrics(
        left,
        right,
        shared_taxa,
        rf_mode=rf_mode,
    )
    return (
        metrics.left_count,
        metrics.right_count,
        metrics.distance,
        metrics.normalized_distance,
    )


def _compare_tree_objects(
    left: PhyloTree,
    right: PhyloTree,
    *,
    rf_mode: RobinsonFouldsMode = "rooted",
    taxon_overlap_policy: TaxonOverlapPolicy = "prune-to-shared",
) -> InMemoryTopologyComparison:
    left_taxa = set(left.tip_names)
    right_taxa = set(right.tip_names)
    shared_taxa, left_only_taxa, right_only_taxa = _resolve_shared_taxa(
        left_taxa,
        right_taxa,
        taxon_overlap_policy=taxon_overlap_policy,
    )
    (
        left_rooted_split_count,
        right_rooted_split_count,
        rooted_distance,
        rooted_normalized,
    ) = _robinson_foulds_metrics(left, right, shared_taxa, rf_mode="rooted")
    (
        left_unrooted_split_count,
        right_unrooted_split_count,
        unrooted_distance,
        unrooted_normalized,
    ) = _robinson_foulds_metrics(left, right, shared_taxa, rf_mode="unrooted")
    topology_equal = rooted_distance == 0
    same_unrooted_topology = unrooted_distance == 0
    selected_distance = rooted_distance if rf_mode == "rooted" else unrooted_distance
    selected_normalized = (
        rooted_normalized if rf_mode == "rooted" else unrooted_normalized
    )
    same_taxa_different_rooting = (
        left_taxa == right_taxa and same_unrooted_topology and not topology_equal
    )
    return InMemoryTopologyComparison(
        shared_taxa=sorted(shared_taxa),
        left_only_taxa=left_only_taxa,
        right_only_taxa=right_only_taxa,
        taxon_overlap_policy=taxon_overlap_policy,
        rf_mode=rf_mode,
        left_informative_clades=left_rooted_split_count,
        right_informative_clades=right_rooted_split_count,
        left_unrooted_splits=left_unrooted_split_count,
        right_unrooted_splits=right_unrooted_split_count,
        robinson_foulds_distance=selected_distance,
        normalized_robinson_foulds=selected_normalized,
        rooted_robinson_foulds_distance=rooted_distance,
        rooted_normalized_robinson_foulds=rooted_normalized,
        unrooted_robinson_foulds_distance=unrooted_distance,
        unrooted_normalized_robinson_foulds=unrooted_normalized,
        topology_equal=topology_equal,
        same_unrooted_topology=same_unrooted_topology,
        same_taxa_different_rooting=same_taxa_different_rooting,
    )


def compare_robinson_foulds(
    left_path: Path,
    right_path: Path,
    *,
    rf_mode: RobinsonFouldsMode = "rooted",
    taxon_overlap_policy: TaxonOverlapPolicy = "prune-to-shared",
) -> RobinsonFouldsComparisonReport:
    """Compare two trees using rooted or unrooted Robinson-Foulds distance."""
    left = _load_tree(left_path)
    right = _load_tree(right_path)
    comparison = _compare_tree_objects(
        left,
        right,
        rf_mode=rf_mode,
        taxon_overlap_policy=taxon_overlap_policy,
    )
    return RobinsonFouldsComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=comparison.shared_taxa,
        left_only_taxa=comparison.left_only_taxa,
        right_only_taxa=comparison.right_only_taxa,
        taxon_overlap_policy=comparison.taxon_overlap_policy,
        rf_mode=comparison.rf_mode,
        left_split_count=(
            comparison.left_informative_clades
            if rf_mode == "rooted"
            else comparison.left_unrooted_splits
        ),
        right_split_count=(
            comparison.right_informative_clades
            if rf_mode == "rooted"
            else comparison.right_unrooted_splits
        ),
        robinson_foulds_distance=comparison.robinson_foulds_distance,
        normalized_robinson_foulds=comparison.normalized_robinson_foulds,
        topology_equal=(
            comparison.topology_equal
            if rf_mode == "rooted"
            else comparison.same_unrooted_topology
        ),
    )


def compare_tree_paths(
    left_path: Path,
    right_path: Path,
    *,
    rf_mode: RobinsonFouldsMode = "rooted",
    taxon_overlap_policy: TaxonOverlapPolicy = "prune-to-shared",
) -> TreeComparisonReport:
    """Compare two trees over their shared taxa."""
    left = _load_tree(left_path)
    right = _load_tree(right_path)
    return _build_tree_comparison_report(
        left_path,
        right_path,
        left,
        right,
        rf_mode=rf_mode,
        taxon_overlap_policy=taxon_overlap_policy,
    )


def _build_tree_comparison_report(
    left_path: Path,
    right_path: Path,
    left: PhyloTree,
    right: PhyloTree,
    *,
    rf_mode: RobinsonFouldsMode,
    taxon_overlap_policy: TaxonOverlapPolicy,
) -> TreeComparisonReport:
    from . import _compare_branch_lengths_for_trees

    comparison = _compare_tree_objects(
        left,
        right,
        rf_mode=rf_mode,
        taxon_overlap_policy=taxon_overlap_policy,
    )
    branch_report = _compare_branch_lengths_for_trees(
        left_path,
        right_path,
        left,
        right,
        taxon_overlap_policy=taxon_overlap_policy,
    )
    same_topology_different_branch_lengths = comparison.topology_equal and any(
        row.left_length != row.right_length for row in branch_report.shared_splits
    )
    return TreeComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=comparison.shared_taxa,
        left_only_taxa=comparison.left_only_taxa,
        right_only_taxa=comparison.right_only_taxa,
        taxon_overlap_policy=comparison.taxon_overlap_policy,
        rf_mode=comparison.rf_mode,
        left_informative_clades=comparison.left_informative_clades,
        right_informative_clades=comparison.right_informative_clades,
        left_unrooted_splits=comparison.left_unrooted_splits,
        right_unrooted_splits=comparison.right_unrooted_splits,
        robinson_foulds_distance=comparison.robinson_foulds_distance,
        normalized_robinson_foulds=comparison.normalized_robinson_foulds,
        rooted_robinson_foulds_distance=comparison.rooted_robinson_foulds_distance,
        rooted_normalized_robinson_foulds=(
            comparison.rooted_normalized_robinson_foulds
        ),
        unrooted_robinson_foulds_distance=(
            comparison.unrooted_robinson_foulds_distance
        ),
        unrooted_normalized_robinson_foulds=(
            comparison.unrooted_normalized_robinson_foulds
        ),
        topology_equal=comparison.topology_equal,
        same_unrooted_topology=comparison.same_unrooted_topology,
        same_taxa_different_rooting=comparison.same_taxa_different_rooting,
        same_topology_different_branch_lengths=same_topology_different_branch_lengths,
    )
