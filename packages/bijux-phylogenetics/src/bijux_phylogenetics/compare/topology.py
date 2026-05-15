from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.pruning import (
    RequestedTaxaPruningReport,
    _prune_tree_against_taxa,
    prune_tree_to_requested_taxa,
)
from bijux_phylogenetics.core.clade_sets import (
    canonical_bipartition,
    canonical_clade_id,
    informative_rooted_clade_nodes,
    informative_rooted_clades,
    informative_unrooted_splits,
    node_support_value,
    robinson_foulds_metrics,
    split_sort_key,
)
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.diagnostics.validation import _load_tree
from bijux_phylogenetics.io.iqtree_support import support_fraction

RobinsonFouldsMode = str
TaxonOverlapPolicy = str
BranchScoreStatus = str

_STRONG_SUPPORT_THRESHOLD = 0.9
_WEAK_SUPPORT_THRESHOLD = 0.7
_SUPPORT_DISAGREEMENT_THRESHOLD = 0.15


@dataclass(slots=True)
class RobinsonFouldsComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    taxon_overlap_policy: str
    rf_mode: str
    left_split_count: int
    right_split_count: int
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    topology_equal: bool


@dataclass(slots=True)
class TreeComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    taxon_overlap_policy: str
    rf_mode: str
    left_informative_clades: int
    right_informative_clades: int
    left_unrooted_splits: int
    right_unrooted_splits: int
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    rooted_robinson_foulds_distance: int
    rooted_normalized_robinson_foulds: float
    unrooted_robinson_foulds_distance: int
    unrooted_normalized_robinson_foulds: float
    topology_equal: bool
    same_unrooted_topology: bool
    same_taxa_different_rooting: bool
    same_topology_different_branch_lengths: bool


@dataclass(slots=True)
class SharedTaxaPruningReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    left_pruning: RequestedTaxaPruningReport
    right_pruning: RequestedTaxaPruningReport
    post_pruning_comparison: TreeComparisonReport


@dataclass(slots=True)
class CladeSupportPair:
    split_id: str
    left_support: float | None
    right_support: float | None
    left_support_fraction: float | None
    right_support_fraction: float | None
    support_fraction_delta: float | None
    support_disagreement: bool


@dataclass(slots=True)
class SupportConflictRow:
    split_id: str
    comparison_status: str
    left_present: bool
    right_present: bool
    left_support: float | None
    right_support: float | None
    left_support_fraction: float | None
    right_support_fraction: float | None
    strongest_support_fraction: float | None
    support_strength: str
    conflict_classification: str
    detail: str


@dataclass(slots=True)
class SupportComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    strong_support_threshold: float
    weak_support_threshold: float
    support_disagreement_threshold: float
    shared_clades: list[CladeSupportPair]
    conflicting_clades: list[SupportConflictRow]


@dataclass(slots=True)
class BranchLengthPair:
    split_id: str
    left_length: float | None
    right_length: float | None
    delta: float | None
    ratio: float | None


@dataclass(slots=True)
class BranchScoreSplit:
    split_id: str
    comparison_status: str
    left_length: float | None
    right_length: float | None
    branch_score_difference: float | None
    squared_difference: float | None


@dataclass(slots=True)
class BranchScoreComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    taxon_overlap_policy: str
    same_taxon_set: bool
    branch_score_distance: float | None
    split_count: int
    shared_split_count: int
    left_only_split_count: int
    right_only_split_count: int
    missing_length_split_count: int
    splits: list[BranchScoreSplit]


@dataclass(slots=True)
class BranchLengthComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    taxon_overlap_policy: str
    same_taxon_set: bool
    shared_splits: list[BranchLengthPair]
    branch_score: BranchScoreComparisonReport


@dataclass(slots=True)
class CladeSetComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    shared_clades: list[str]
    left_only_clades: list[str]
    right_only_clades: list[str]


@dataclass(slots=True)
class CladeOverlapObservation:
    tree_path: Path
    present: bool
    support: float | None


@dataclass(slots=True)
class CladeOverlapRow:
    clade_id: str
    present_in_all_trees: bool
    present_tree_count: int
    absent_tree_count: int
    observations: list[CladeOverlapObservation]


@dataclass(slots=True)
class TreeCladeOverlapSummary:
    tree_path: Path
    clade_count: int
    support_clade_count: int
    unique_clades: list[str]
    excluded_taxa: list[str]


@dataclass(slots=True)
class CladeOverlapComparisonReport:
    tree_paths: list[Path]
    shared_taxa: list[str]
    shared_clades: list[str]
    conflicting_clades: list[str]
    tree_summaries: list[TreeCladeOverlapSummary]
    clade_rows: list[CladeOverlapRow]


@dataclass(slots=True)
class CladeChangeReport:
    left_path: Path
    right_path: Path
    lost_clades: list[str]
    gained_clades: list[str]


@dataclass(slots=True)
class InMemoryTopologyComparison:
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    taxon_overlap_policy: str
    rf_mode: str
    left_informative_clades: int
    right_informative_clades: int
    left_unrooted_splits: int
    right_unrooted_splits: int
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    rooted_robinson_foulds_distance: int
    rooted_normalized_robinson_foulds: float
    unrooted_robinson_foulds_distance: int
    unrooted_normalized_robinson_foulds: float
    topology_equal: bool
    same_unrooted_topology: bool
    same_taxa_different_rooting: bool


@dataclass(slots=True)
class InMemoryBranchLengthComparison:
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    taxon_overlap_policy: str
    shared_splits: list[BranchLengthPair]
    branch_score_distance: float | None
    branch_score_splits: list[BranchScoreSplit]
    missing_length_split_count: int

def _format_clade_set(clades: set[frozenset[str]]) -> list[str]:
    return sorted(canonical_clade_id(clade) for clade in clades)


def _split_id(signature: frozenset[str]) -> str:
    return canonical_clade_id(signature)


def _support_strength(
    value: float | None,
    *,
    strong_support_threshold: float,
    weak_support_threshold: float,
) -> str:
    fraction = support_fraction(value)
    if fraction is None:
        return "unavailable"
    if fraction >= strong_support_threshold:
        return "strong"
    if fraction >= weak_support_threshold:
        return "moderate"
    return "low"


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


def _sum_branch_lengths(left: float | None, right: float | None) -> float | None:
    if left is None or right is None:
        return None
    return left + right


def _unrooted_branch_score_lengths(
    tree: PhyloTree,
) -> dict[frozenset[str], float | None]:
    taxa = set(tree.tip_names)
    if len(taxa) < 2:
        raise ValueError("branch-score comparison requires at least two retained taxa")

    split_lengths: dict[frozenset[str], float | None] = {}
    suppress_binary_root = len(tree.root.children) == 2

    def visit(node: TreeNode, *, parent_is_root: bool) -> set[str]:
        if node.is_leaf():
            if node.name is None or node.name not in taxa:
                return set()
            descendant_taxa = {node.name}
        else:
            descendant_taxa: set[str] = set()
            for child in node.children:
                descendant_taxa.update(
                    visit(
                        child, parent_is_root=node is tree.root and suppress_binary_root
                    )
                )
        if not descendant_taxa:
            return set()
        if node is tree.root:
            return descendant_taxa
        if not parent_is_root:
            split_lengths[canonical_bipartition(descendant_taxa, taxa)] = (
                node.branch_length
            )
        return descendant_taxa

    child_taxa = [
        visit(child, parent_is_root=suppress_binary_root)
        for child in tree.root.children
    ]
    if tree.root.branch_length is not None:
        split_lengths[frozenset(taxa)] = tree.root.branch_length
    if not suppress_binary_root:
        return split_lengths

    if len(tree.root.children) != 2:
        return split_lengths

    left_child, right_child = tree.root.children
    left_taxa, right_taxa = child_taxa
    if not left_taxa or not right_taxa:
        return split_lengths

    merged_split = canonical_bipartition(left_taxa, taxa)
    split_lengths[merged_split] = _sum_branch_lengths(
        left_child.branch_length,
        right_child.branch_length,
    )
    return split_lengths


def _build_branch_score_report(
    left_path: Path,
    right_path: Path,
    left: PhyloTree,
    right: PhyloTree,
    *,
    taxon_overlap_policy: TaxonOverlapPolicy,
) -> BranchScoreComparisonReport:
    left_taxa = set(left.tip_names)
    right_taxa = set(right.tip_names)
    shared_taxa, left_only_taxa, right_only_taxa = _resolve_shared_taxa(
        left_taxa,
        right_taxa,
        taxon_overlap_policy=taxon_overlap_policy,
    )
    same_taxon_set = not left_only_taxa and not right_only_taxa
    if same_taxon_set:
        left_pruned = left
        right_pruned = right
    else:
        left_pruned, _, _ = _prune_tree_against_taxa(
            left,
            shared_taxa,
            clear_root_branch_length=False,
        )
        right_pruned, _, _ = _prune_tree_against_taxa(
            right,
            shared_taxa,
            clear_root_branch_length=False,
        )

    left_lengths = _unrooted_branch_score_lengths(left_pruned)
    right_lengths = _unrooted_branch_score_lengths(right_pruned)
    split_ids = sorted(
        set(left_lengths) | set(right_lengths),
        key=split_sort_key,
    )
    rows: list[BranchScoreSplit] = []
    sum_of_squares = 0.0
    missing_length_split_count = 0
    shared_split_count = 0
    left_only_split_count = 0
    right_only_split_count = 0
    for signature in split_ids:
        left_present = signature in left_lengths
        right_present = signature in right_lengths
        if left_present and right_present:
            status = "shared"
            shared_split_count += 1
        elif left_present:
            status = "left_only"
            left_only_split_count += 1
        else:
            status = "right_only"
            right_only_split_count += 1
        left_length = left_lengths.get(signature)
        right_length = right_lengths.get(signature)
        branch_score_difference: float | None
        squared_difference: float | None
        if (left_present and left_length is None) or (
            right_present and right_length is None
        ):
            branch_score_difference = None
            squared_difference = None
            missing_length_split_count += 1
        else:
            numeric_left = 0.0 if not left_present else float(left_length)
            numeric_right = 0.0 if not right_present else float(right_length)
            branch_score_difference = numeric_right - numeric_left
            squared_difference = branch_score_difference**2
            sum_of_squares += branch_score_difference**2
        rows.append(
            BranchScoreSplit(
                split_id=_split_id(signature),
                comparison_status=status,
                left_length=left_length,
                right_length=right_length,
                branch_score_difference=branch_score_difference,
                squared_difference=squared_difference,
            )
        )
    branch_score_distance = None if missing_length_split_count else sum_of_squares**0.5
    return BranchScoreComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=sorted(shared_taxa),
        left_only_taxa=left_only_taxa,
        right_only_taxa=right_only_taxa,
        taxon_overlap_policy=taxon_overlap_policy,
        same_taxon_set=same_taxon_set,
        branch_score_distance=branch_score_distance,
        split_count=len(rows),
        shared_split_count=shared_split_count,
        left_only_split_count=left_only_split_count,
        right_only_split_count=right_only_split_count,
        missing_length_split_count=missing_length_split_count,
        splits=rows,
    )


def _compare_branch_lengths_for_trees(
    left_path: Path,
    right_path: Path,
    left: PhyloTree,
    right: PhyloTree,
    *,
    taxon_overlap_policy: TaxonOverlapPolicy,
) -> InMemoryBranchLengthComparison:
    left_taxa = set(left.tip_names)
    right_taxa = set(right.tip_names)
    shared_taxa, left_only_taxa, right_only_taxa = _resolve_shared_taxa(
        left_taxa,
        right_taxa,
        taxon_overlap_policy=taxon_overlap_policy,
    )

    left_clades = informative_rooted_clade_nodes(left, shared_taxa)
    right_clades = informative_rooted_clade_nodes(right, shared_taxa)
    shared_clade_ids = sorted(
        left_clades.keys() & right_clades.keys(), key=lambda item: sorted(item)
    )
    pairs: list[BranchLengthPair] = []
    for clade_id in shared_clade_ids:
        left_length = left_clades[clade_id].branch_length
        right_length = right_clades[clade_id].branch_length
        delta = (
            None
            if left_length is None or right_length is None
            else right_length - left_length
        )
        ratio = None
        if left_length not in {None, 0} and right_length is not None:
            ratio = right_length / left_length
        pairs.append(
            BranchLengthPair(
                split_id="|".join(sorted(clade_id)),
                left_length=left_length,
                right_length=right_length,
                delta=delta,
                ratio=ratio,
            )
        )
    branch_score = _build_branch_score_report(
        left_path,
        right_path,
        left,
        right,
        taxon_overlap_policy=taxon_overlap_policy,
    )
    return InMemoryBranchLengthComparison(
        shared_taxa=sorted(shared_taxa),
        left_only_taxa=left_only_taxa,
        right_only_taxa=right_only_taxa,
        taxon_overlap_policy=taxon_overlap_policy,
        shared_splits=pairs,
        branch_score_distance=branch_score.branch_score_distance,
        branch_score_splits=branch_score.splits,
        missing_length_split_count=branch_score.missing_length_split_count,
    )

def _resolve_shared_taxa_for_many_trees(
    tree_paths: list[Path],
) -> tuple[list[PhyloTree], set[str], list[list[str]]]:
    if len(tree_paths) < 2:
        raise ValueError("clade-overlap comparison requires at least two trees")
    trees = [_load_tree(path) for path in tree_paths]
    taxon_sets = [set(tree.tip_names) for tree in trees]
    shared_taxa = set.intersection(*taxon_sets)
    if len(shared_taxa) < 2:
        raise ValueError("clade-overlap comparison requires at least two shared taxa")
    excluded_taxa = [sorted(taxa - shared_taxa) for taxa in taxon_sets]
    return trees, shared_taxa, excluded_taxa


def prune_trees_to_shared_taxa(
    left_path: Path,
    right_path: Path,
) -> tuple[PhyloTree, PhyloTree, SharedTaxaPruningReport]:
    """Prune two trees to the exact shared taxon set."""
    left = _load_tree(left_path)
    right = _load_tree(right_path)
    left_taxa = set(left.tip_names)
    right_taxa = set(right.tip_names)
    shared_taxa = sorted(left_taxa & right_taxa)
    if len(shared_taxa) < 2:
        raise ValueError("shared-taxon pruning requires at least two shared taxa")

    pruned_left, left_pruning = prune_tree_to_requested_taxa(left_path, shared_taxa)
    pruned_right, right_pruning = prune_tree_to_requested_taxa(right_path, shared_taxa)
    post_pruning_comparison = _build_tree_comparison_report(
        left_path,
        right_path,
        pruned_left,
        pruned_right,
        rf_mode="rooted",
        taxon_overlap_policy="require-identical",
    )
    return (
        pruned_left,
        pruned_right,
        SharedTaxaPruningReport(
            left_path=left_path,
            right_path=right_path,
            shared_taxa=shared_taxa,
            left_only_taxa=sorted(left_taxa - right_taxa),
            right_only_taxa=sorted(right_taxa - left_taxa),
            left_pruning=left_pruning,
            right_pruning=right_pruning,
            post_pruning_comparison=post_pruning_comparison,
        ),
    )


def compare_clade_sets(left_path: Path, right_path: Path) -> CladeSetComparisonReport:
    """Compare rooted informative clade sets across two trees."""
    overlap = compare_clade_overlap([left_path, right_path])
    left_summary, right_summary = overlap.tree_summaries
    return CladeSetComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=overlap.shared_taxa,
        shared_clades=overlap.shared_clades,
        left_only_clades=left_summary.unique_clades,
        right_only_clades=right_summary.unique_clades,
    )


def compare_clade_overlap(tree_paths: list[Path]) -> CladeOverlapComparisonReport:
    """Compare rooted clade overlap across two or more trees."""
    trees, shared_taxa, excluded_taxa = _resolve_shared_taxa_for_many_trees(tree_paths)
    clade_maps = [informative_rooted_clades(tree, shared_taxa) for tree in trees]
    clade_node_maps = [
        informative_rooted_clade_nodes(tree, shared_taxa) for tree in trees
    ]
    all_clades = sorted(
        set().union(*clade_maps),
        key=lambda signature: (len(signature), tuple(sorted(signature))),
    )
    shared_clades = [
        _split_id(clade)
        for clade in all_clades
        if all(clade in clade_map for clade_map in clade_maps)
    ]
    conflicting_clades = [
        _split_id(clade)
        for clade in all_clades
        if not all(clade in clade_map for clade_map in clade_maps)
    ]
    tree_summaries: list[TreeCladeOverlapSummary] = []
    for path, clade_map, clade_node_map, tree_excluded_taxa in zip(
        tree_paths,
        clade_maps,
        clade_node_maps,
        excluded_taxa,
        strict=True,
    ):
        unique_clades = [
            _split_id(clade)
            for clade in sorted(
                (
                    clade
                    for clade in clade_map
                    if sum(clade in other_map for other_map in clade_maps) == 1
                ),
                key=lambda signature: (len(signature), tuple(sorted(signature))),
            )
        ]
        support_clade_count = sum(
            1
            for clade in clade_map
            if node_support_value(clade_node_map[clade]) is not None
        )
        tree_summaries.append(
            TreeCladeOverlapSummary(
                tree_path=path,
                clade_count=len(clade_map),
                support_clade_count=support_clade_count,
                unique_clades=unique_clades,
                excluded_taxa=tree_excluded_taxa,
            )
        )
    clade_rows: list[CladeOverlapRow] = []
    for clade in all_clades:
        observations: list[CladeOverlapObservation] = []
        present_tree_count = 0
        for path, clade_map, clade_node_map in zip(
            tree_paths, clade_maps, clade_node_maps, strict=True
        ):
            present = clade in clade_map
            if present:
                present_tree_count += 1
            support = None
            if present:
                support = node_support_value(clade_node_map[clade])
            observations.append(
                CladeOverlapObservation(
                    tree_path=path,
                    present=present,
                    support=support,
                )
            )
        clade_rows.append(
            CladeOverlapRow(
                clade_id=_split_id(clade),
                present_in_all_trees=present_tree_count == len(tree_paths),
                present_tree_count=present_tree_count,
                absent_tree_count=len(tree_paths) - present_tree_count,
                observations=observations,
            )
        )
    return CladeOverlapComparisonReport(
        tree_paths=tree_paths,
        shared_taxa=sorted(shared_taxa),
        shared_clades=shared_clades,
        conflicting_clades=conflicting_clades,
        tree_summaries=tree_summaries,
        clade_rows=clade_rows,
    )


def detect_clade_changes(left_path: Path, right_path: Path) -> CladeChangeReport:
    """Report clades lost from the left tree and gained in the right tree."""
    report = compare_clade_sets(left_path, right_path)
    return CladeChangeReport(
        left_path=left_path,
        right_path=right_path,
        lost_clades=report.left_only_clades,
        gained_clades=report.right_only_clades,
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


def _pipe_join(values: list[str]) -> str:
    return "|".join(values)


def write_shared_taxa_pruning_table(
    path: Path,
    left_path: Path,
    right_path: Path,
) -> Path:
    """Write one row per tree summarizing shared-taxon pruning evidence."""
    _, _, report = prune_trees_to_shared_taxa(left_path, right_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "tree_side",
                "tree_path",
                "original_tip_count",
                "retained_tip_count",
                "removed_tip_count",
                "requested_taxa",
                "kept_taxa",
                "removed_taxa",
                "absent_requested_taxa",
                "removed_taxa_with_reasons",
                "transformation",
                "root_to_tip_complete",
                "min_root_to_tip",
                "max_root_to_tip",
                "unary_internal_nodes",
                "original_total_branch_length",
                "pruned_total_branch_length",
                "branch_length_delta",
                "lost_taxa_count",
                "lost_taxa_fraction",
                "lost_clade_count",
                "lost_clade_fraction",
                "lost_branch_length",
                "lost_branch_length_fraction",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for tree_side, pruning in (
            ("left", report.left_pruning),
            ("right", report.right_pruning),
        ):
            writer.writerow(
                {
                    "tree_side": tree_side,
                    "tree_path": str(pruning.tree_path),
                    "original_tip_count": pruning.original_tip_count,
                    "retained_tip_count": len(pruning.kept_taxa),
                    "removed_tip_count": len(pruning.removed_taxa),
                    "requested_taxa": _pipe_join(pruning.requested_taxa),
                    "kept_taxa": _pipe_join(pruning.kept_taxa),
                    "removed_taxa": _pipe_join(pruning.removed_taxa),
                    "absent_requested_taxa": _pipe_join(pruning.absent_requested_taxa),
                    "removed_taxa_with_reasons": _pipe_join(
                        [
                            f"{row.taxon}:{row.reason}"
                            for row in pruning.removed_taxa_with_reasons
                        ]
                    ),
                    "transformation": pruning.summary.transformation,
                    "root_to_tip_complete": str(
                        pruning.pruning_audit.root_to_tip_complete
                    ).lower(),
                    "min_root_to_tip": pruning.pruning_audit.min_root_to_tip,
                    "max_root_to_tip": pruning.pruning_audit.max_root_to_tip,
                    "unary_internal_nodes": _pipe_join(
                        pruning.pruning_audit.unary_internal_nodes
                    ),
                    "original_total_branch_length": (
                        pruning.pruning_audit.original_total_branch_length
                    ),
                    "pruned_total_branch_length": (
                        pruning.pruning_audit.pruned_total_branch_length
                    ),
                    "branch_length_delta": pruning.pruning_audit.branch_length_delta,
                    "lost_taxa_count": pruning.information_loss.lost_taxa_count,
                    "lost_taxa_fraction": pruning.information_loss.lost_taxa_fraction,
                    "lost_clade_count": pruning.information_loss.lost_clade_count,
                    "lost_clade_fraction": pruning.information_loss.lost_clade_fraction,
                    "lost_branch_length": pruning.information_loss.lost_branch_length,
                    "lost_branch_length_fraction": (
                        pruning.information_loss.lost_branch_length_fraction
                    ),
                }
            )
    return path


def write_shared_taxa_removed_taxa_table(
    path: Path,
    left_path: Path,
    right_path: Path,
) -> Path:
    """Write one row per removed taxon from shared-taxon pruning."""
    _, _, report = prune_trees_to_shared_taxa(left_path, right_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["tree_side", "tree_path", "taxon", "reason"],
            delimiter="\t",
        )
        writer.writeheader()
        for tree_side, pruning in (
            ("left", report.left_pruning),
            ("right", report.right_pruning),
        ):
            for removed in pruning.removed_taxa_with_reasons:
                writer.writerow(
                    {
                        "tree_side": tree_side,
                        "tree_path": str(pruning.tree_path),
                        "taxon": removed.taxon,
                        "reason": removed.reason,
                    }
                )
    return path


def compare_support_values(
    left_path: Path,
    right_path: Path,
    *,
    strong_support_threshold: float = _STRONG_SUPPORT_THRESHOLD,
    weak_support_threshold: float = _WEAK_SUPPORT_THRESHOLD,
    support_disagreement_threshold: float = _SUPPORT_DISAGREEMENT_THRESHOLD,
) -> SupportComparisonReport:
    """Compare clade support values and support-aware conflicts across two trees."""
    left = _load_tree(left_path)
    right = _load_tree(right_path)
    return _build_support_comparison_report(
        left_path,
        right_path,
        left,
        right,
        strong_support_threshold=strong_support_threshold,
        weak_support_threshold=weak_support_threshold,
        support_disagreement_threshold=support_disagreement_threshold,
    )


def _build_support_comparison_report(
    left_path: Path,
    right_path: Path,
    left: PhyloTree,
    right: PhyloTree,
    *,
    strong_support_threshold: float,
    weak_support_threshold: float,
    support_disagreement_threshold: float,
) -> SupportComparisonReport:
    left_taxa = set(left.tip_names)
    right_taxa = set(right.tip_names)
    shared_taxa = left_taxa & right_taxa
    if len(shared_taxa) < 2:
        raise ValueError("support comparison requires at least two shared taxa")

    left_clades = informative_rooted_clade_nodes(left, shared_taxa)
    right_clades = informative_rooted_clade_nodes(right, shared_taxa)
    shared_clade_ids = sorted(left_clades.keys() & right_clades.keys(), key=split_sort_key)
    all_clade_ids = sorted(
        left_clades.keys() | right_clades.keys(),
        key=split_sort_key,
    )
    shared_clades: list[CladeSupportPair] = []
    conflicting_clades: list[SupportConflictRow] = []
    for clade_id in shared_clade_ids:
        left_support = node_support_value(left_clades[clade_id])
        right_support = node_support_value(right_clades[clade_id])
        left_fraction = support_fraction(left_support)
        right_fraction = support_fraction(right_support)
        shared_clades.append(
            CladeSupportPair(
                split_id=_split_id(clade_id),
                left_support=left_support,
                right_support=right_support,
                left_support_fraction=left_fraction,
                right_support_fraction=right_fraction,
                support_fraction_delta=(
                    None
                    if left_fraction is None or right_fraction is None
                    else abs(left_fraction - right_fraction)
                ),
                support_disagreement=(
                    left_fraction is not None
                    and right_fraction is not None
                    and abs(left_fraction - right_fraction)
                    >= support_disagreement_threshold
                ),
            )
        )
    for clade_id in all_clade_ids:
        left_present = clade_id in left_clades
        right_present = clade_id in right_clades
        if left_present and right_present:
            continue
        left_support = (
            node_support_value(left_clades[clade_id])
            if left_present
            else None
        )
        right_support = (
            node_support_value(right_clades[clade_id])
            if right_present
            else None
        )
        left_fraction = support_fraction(left_support)
        right_fraction = support_fraction(right_support)
        strongest_support_fraction = (
            max(value for value in (left_fraction, right_fraction) if value is not None)
            if left_fraction is not None or right_fraction is not None
            else None
        )
        support_strength = _support_strength(
            left_support if left_present else right_support,
            strong_support_threshold=strong_support_threshold,
            weak_support_threshold=weak_support_threshold,
        )
        if strongest_support_fraction is None:
            conflict_classification = "support_unavailable"
            detail = "clade conflict could not be ranked because no branch support was available"
        elif strongest_support_fraction >= strong_support_threshold:
            conflict_classification = "high_support_conflict"
            detail = "conflicting clade carried strong branch support in the tree where it was present"
        elif strongest_support_fraction >= weak_support_threshold:
            conflict_classification = "moderate_support_disagreement"
            detail = "conflicting clade carried moderate branch support in the tree where it was present"
        else:
            conflict_classification = "low_support_disagreement"
            detail = "conflicting clade was only weakly supported in the tree where it was present"
        conflicting_clades.append(
            SupportConflictRow(
                split_id=_split_id(clade_id),
                comparison_status="left_only" if left_present else "right_only",
                left_present=left_present,
                right_present=right_present,
                left_support=left_support,
                right_support=right_support,
                left_support_fraction=left_fraction,
                right_support_fraction=right_fraction,
                strongest_support_fraction=strongest_support_fraction,
                support_strength=support_strength,
                conflict_classification=conflict_classification,
                detail=detail,
            )
        )
    return SupportComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=sorted(shared_taxa),
        strong_support_threshold=strong_support_threshold,
        weak_support_threshold=weak_support_threshold,
        support_disagreement_threshold=support_disagreement_threshold,
        shared_clades=shared_clades,
        conflicting_clades=conflicting_clades,
    )


def compare_branch_lengths(
    left_path: Path,
    right_path: Path,
    *,
    taxon_overlap_policy: TaxonOverlapPolicy = "prune-to-shared",
) -> BranchLengthComparisonReport:
    """Compare branch lengths across matching informative splits."""
    left = _load_tree(left_path)
    right = _load_tree(right_path)
    comparison = _compare_branch_lengths_for_trees(
        left_path,
        right_path,
        left,
        right,
        taxon_overlap_policy=taxon_overlap_policy,
    )
    branch_score = _build_branch_score_report(
        left_path,
        right_path,
        left,
        right,
        taxon_overlap_policy=taxon_overlap_policy,
    )
    return BranchLengthComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=comparison.shared_taxa,
        left_only_taxa=comparison.left_only_taxa,
        right_only_taxa=comparison.right_only_taxa,
        taxon_overlap_policy=comparison.taxon_overlap_policy,
        same_taxon_set=not comparison.left_only_taxa and not comparison.right_only_taxa,
        shared_splits=comparison.shared_splits,
        branch_score=branch_score,
    )


def compare_branch_score_distance(
    left_path: Path,
    right_path: Path,
    *,
    taxon_overlap_policy: TaxonOverlapPolicy = "prune-to-shared",
) -> BranchScoreComparisonReport:
    """Compare two trees using Felsenstein branch-score distance."""
    left = _load_tree(left_path)
    right = _load_tree(right_path)
    return _build_branch_score_report(
        left_path,
        right_path,
        left,
        right,
        taxon_overlap_policy=taxon_overlap_policy,
    )


def write_clade_overlap_table(path: Path, tree_paths: list[Path]) -> Path:
    """Write one row per clade-per-tree overlap observation."""
    report = compare_clade_overlap(tree_paths)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "clade_id",
                "tree_path",
                "present",
                "support",
                "present_in_all_trees",
                "present_tree_count",
                "absent_tree_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.clade_rows:
            for observation in row.observations:
                writer.writerow(
                    {
                        "clade_id": row.clade_id,
                        "tree_path": str(observation.tree_path),
                        "present": str(observation.present).lower(),
                        "support": ""
                        if observation.support is None
                        else observation.support,
                        "present_in_all_trees": str(row.present_in_all_trees).lower(),
                        "present_tree_count": row.present_tree_count,
                        "absent_tree_count": row.absent_tree_count,
                    }
                )
    return path


def write_support_comparison_table(
    path: Path, left_path: Path, right_path: Path
) -> Path:
    """Write a flat TSV ledger for support-aware clade comparison."""
    report = compare_support_values(left_path, right_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "split_id",
                "row_kind",
                "comparison_status",
                "left_present",
                "right_present",
                "left_support",
                "right_support",
                "left_support_fraction",
                "right_support_fraction",
                "support_fraction_delta",
                "support_disagreement",
                "strongest_support_fraction",
                "support_strength",
                "conflict_classification",
                "detail",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.shared_clades:
            writer.writerow(
                {
                    "split_id": row.split_id,
                    "row_kind": "shared_clade",
                    "comparison_status": "shared",
                    "left_present": "true",
                    "right_present": "true",
                    "left_support": ""
                    if row.left_support is None
                    else row.left_support,
                    "right_support": ""
                    if row.right_support is None
                    else row.right_support,
                    "left_support_fraction": (
                        ""
                        if row.left_support_fraction is None
                        else row.left_support_fraction
                    ),
                    "right_support_fraction": (
                        ""
                        if row.right_support_fraction is None
                        else row.right_support_fraction
                    ),
                    "support_fraction_delta": (
                        ""
                        if row.support_fraction_delta is None
                        else row.support_fraction_delta
                    ),
                    "support_disagreement": str(row.support_disagreement).lower(),
                    "strongest_support_fraction": "",
                    "support_strength": "",
                    "conflict_classification": "",
                    "detail": (
                        "shared clade support differs across trees"
                        if row.support_disagreement
                        else "shared clade support is aligned across trees"
                    ),
                }
            )
        for row in report.conflicting_clades:
            writer.writerow(
                {
                    "split_id": row.split_id,
                    "row_kind": "conflicting_clade",
                    "comparison_status": row.comparison_status,
                    "left_present": str(row.left_present).lower(),
                    "right_present": str(row.right_present).lower(),
                    "left_support": ""
                    if row.left_support is None
                    else row.left_support,
                    "right_support": ""
                    if row.right_support is None
                    else row.right_support,
                    "left_support_fraction": (
                        ""
                        if row.left_support_fraction is None
                        else row.left_support_fraction
                    ),
                    "right_support_fraction": (
                        ""
                        if row.right_support_fraction is None
                        else row.right_support_fraction
                    ),
                    "support_fraction_delta": "",
                    "support_disagreement": "false",
                    "strongest_support_fraction": (
                        ""
                        if row.strongest_support_fraction is None
                        else row.strongest_support_fraction
                    ),
                    "support_strength": row.support_strength,
                    "conflict_classification": row.conflict_classification,
                    "detail": row.detail,
                }
            )
    return path


def write_tree_comparison_table(path: Path, left_path: Path, right_path: Path) -> Path:
    """Write a flat TSV table covering the compared clade and split surfaces."""
    clades = compare_clade_sets(left_path, right_path)
    support = compare_support_values(left_path, right_path)
    branch_lengths = compare_branch_lengths(left_path, right_path)
    support_by_id = {row.split_id: row for row in support.shared_clades}
    support_conflict_by_id = {row.split_id: row for row in support.conflicting_clades}
    branch_by_id = {row.split_id: row for row in branch_lengths.shared_splits}
    branch_score_by_id = {
        row.split_id: row for row in branch_lengths.branch_score.splits
    }
    all_split_ids = sorted(
        set(clades.shared_clades)
        | set(clades.left_only_clades)
        | set(clades.right_only_clades)
        | set(support_by_id)
        | set(support_conflict_by_id)
        | set(branch_by_id)
        | set(branch_score_by_id)
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "split_id",
                "comparison_status",
                "shared_clade",
                "left_support",
                "right_support",
                "left_support_fraction",
                "right_support_fraction",
                "support_fraction_delta",
                "support_disagreement",
                "support_conflict_classification",
                "support_conflict_strength",
                "left_length",
                "right_length",
                "length_delta",
                "length_ratio",
                "branch_score_status",
                "branch_score_difference",
                "branch_score_squared_difference",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for split_id in all_split_ids:
            support_row = support_by_id.get(split_id)
            support_conflict_row = support_conflict_by_id.get(split_id)
            branch_row = branch_by_id.get(split_id)
            branch_score_row = branch_score_by_id.get(split_id)
            if split_id in clades.shared_clades:
                status = "shared"
            elif split_id in clades.left_only_clades:
                status = "left_only"
            else:
                status = "right_only"
            writer.writerow(
                {
                    "split_id": split_id,
                    "comparison_status": status,
                    "shared_clade": str(split_id in clades.shared_clades).lower(),
                    "left_support": ""
                    if (
                        (support_row is None or support_row.left_support is None)
                        and (
                            support_conflict_row is None
                            or support_conflict_row.left_support is None
                        )
                    )
                    else (
                        support_row.left_support
                        if support_row is not None
                        else support_conflict_row.left_support
                    ),
                    "right_support": ""
                    if (
                        (support_row is None or support_row.right_support is None)
                        and (
                            support_conflict_row is None
                            or support_conflict_row.right_support is None
                        )
                    )
                    else (
                        support_row.right_support
                        if support_row is not None
                        else support_conflict_row.right_support
                    ),
                    "left_support_fraction": ""
                    if (
                        (
                            support_row is None
                            or support_row.left_support_fraction is None
                        )
                        and (
                            support_conflict_row is None
                            or support_conflict_row.left_support_fraction is None
                        )
                    )
                    else (
                        support_row.left_support_fraction
                        if support_row is not None
                        else support_conflict_row.left_support_fraction
                    ),
                    "right_support_fraction": ""
                    if (
                        (
                            support_row is None
                            or support_row.right_support_fraction is None
                        )
                        and (
                            support_conflict_row is None
                            or support_conflict_row.right_support_fraction is None
                        )
                    )
                    else (
                        support_row.right_support_fraction
                        if support_row is not None
                        else support_conflict_row.right_support_fraction
                    ),
                    "support_fraction_delta": ""
                    if support_row is None or support_row.support_fraction_delta is None
                    else support_row.support_fraction_delta,
                    "support_disagreement": (
                        "false"
                        if support_row is None
                        else str(support_row.support_disagreement).lower()
                    ),
                    "support_conflict_classification": (
                        ""
                        if support_conflict_row is None
                        else support_conflict_row.conflict_classification
                    ),
                    "support_conflict_strength": (
                        ""
                        if support_conflict_row is None
                        else support_conflict_row.support_strength
                    ),
                    "left_length": ""
                    if branch_row is None or branch_row.left_length is None
                    else branch_row.left_length,
                    "right_length": ""
                    if branch_row is None or branch_row.right_length is None
                    else branch_row.right_length,
                    "length_delta": ""
                    if branch_row is None or branch_row.delta is None
                    else branch_row.delta,
                    "length_ratio": ""
                    if branch_row is None or branch_row.ratio is None
                    else branch_row.ratio,
                    "branch_score_status": ""
                    if branch_score_row is None
                    else branch_score_row.comparison_status,
                    "branch_score_difference": ""
                    if branch_score_row is None
                    or branch_score_row.branch_score_difference is None
                    else branch_score_row.branch_score_difference,
                    "branch_score_squared_difference": ""
                    if branch_score_row is None
                    or branch_score_row.squared_difference is None
                    else branch_score_row.squared_difference,
                }
            )
    return path
