from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.diagnostics.validation.structure import _load_tree
from bijux_phylogenetics.phylo.pruning import _prune_tree_against_taxa
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_bipartition,
    informative_rooted_clade_nodes,
    split_sort_key,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .comparison import _resolve_shared_taxa
from .models import (
    BranchLengthComparisonReport,
    BranchLengthPair,
    BranchScoreComparisonReport,
    BranchScoreSplit,
    InMemoryBranchLengthComparison,
    TaxonOverlapPolicy,
)


def _split_id(signature: frozenset[str]) -> str:
    return "|".join(sorted(signature))


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
