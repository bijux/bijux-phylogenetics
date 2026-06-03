from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment import AlignmentRecord
from bijux_phylogenetics.phylo.topology import (
    TreeTransformationSummary,
    _summarize_transformation,
)
from bijux_phylogenetics.phylo.topology.clades import informative_rooted_clades
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import MetadataJoinError


@dataclass(slots=True)
class TreePruningReport:
    """Explicit record of tree pruning against a taxon table."""

    tree_path: Path
    keep_from_path: Path
    taxon_column: str
    original_tip_count: int
    kept_taxa: list[str]
    removed_taxa: list[str]
    removed_taxa_with_reasons: list[RemovedTaxonReason]
    pruning_audit: PruningArtifactAudit
    information_loss: PruningInformationLoss
    summary: TreeTransformationSummary


@dataclass(slots=True)
class RequestedTaxaPruningReport:
    """Explicit record of tree pruning against a caller-provided taxon list."""

    tree_path: Path
    original_tip_count: int
    requested_taxa: list[str]
    kept_taxa: list[str]
    removed_taxa: list[str]
    absent_requested_taxa: list[str]
    removed_taxa_with_reasons: list[RemovedTaxonReason]
    pruning_audit: PruningArtifactAudit
    information_loss: PruningInformationLoss
    summary: TreeTransformationSummary


@dataclass(slots=True)
class AlignmentPruningReport:
    """Explicit record of alignment pruning against a tree tip set."""

    alignment_path: Path
    tree_path: Path
    original_sequence_count: int
    kept_ids: list[str]
    removed_ids: list[str]
    removed_ids_with_reasons: list[RemovedTaxonReason]


@dataclass(slots=True)
class RemovedTaxonReason:
    """One removed taxon and the explicit reason it was excluded."""

    taxon: str
    reason: str


@dataclass(slots=True)
class PruningArtifactAudit:
    """Post-pruning structural audit for retained branch lengths and unary cleanup."""

    rooted: bool | None
    root_to_tip_complete: bool
    min_root_to_tip: float | None
    max_root_to_tip: float | None
    unary_internal_nodes: list[str]
    original_total_branch_length: float
    pruned_total_branch_length: float
    branch_length_delta: float


@dataclass(slots=True)
class PruningInformationLoss:
    """Quantified information loss after a pruning transform."""

    lost_taxa_count: int
    lost_taxa_fraction: float
    lost_clade_count: int
    lost_clade_fraction: float
    lost_branch_length: float
    lost_branch_length_fraction: float
    lost_metadata_count: int | None
    lost_metadata_fraction: float | None


def _merge_branch_lengths(left: float | None, right: float | None) -> float | None:
    if left is None:
        return right
    if right is None:
        return left
    return left + right


def _ape_style_rooted_state(root: TreeNode) -> bool:
    """Infer the rooted state ape exposes after valid tip pruning."""
    return root.is_leaf() or len(root.children) <= 2


def _node_signature(node: TreeNode) -> str:
    if node.is_leaf():
        return node.name or "<unnamed>"
    taxa: list[str] = []
    for child in node.children:
        signature = _node_signature(child)
        if "|" in signature:
            taxa.extend(signature.split("|"))
        elif signature != "<unnamed>":
            taxa.append(signature)
    return "|".join(sorted(taxa)) if taxa else (node.name or "<unnamed>")


def _singleton_internal_nodes(tree: PhyloTree) -> list[str]:
    return sorted(
        _node_signature(node)
        for node in tree.iter_nodes()
        if not node.is_leaf() and len(node.children) == 1
    )


def _pruning_artifact_audit(
    original: PhyloTree, pruned: PhyloTree
) -> PruningArtifactAudit:
    root_to_tip = pruned.root_to_tip_lengths()
    numeric = [float(length) for length in root_to_tip if length is not None]
    root_state = (
        pruned.rooted if pruned.rooted is not None else (len(pruned.root.children) == 2)
    )
    original_total = round(original.total_branch_length(), 15)
    pruned_total = round(pruned.total_branch_length(), 15)
    return PruningArtifactAudit(
        rooted=root_state,
        root_to_tip_complete=all(length is not None for length in root_to_tip),
        min_root_to_tip=min(numeric) if numeric else None,
        max_root_to_tip=max(numeric) if numeric else None,
        unary_internal_nodes=_singleton_internal_nodes(pruned),
        original_total_branch_length=original_total,
        pruned_total_branch_length=pruned_total,
        branch_length_delta=round(pruned_total - original_total, 15),
    )


def _pruning_information_loss(
    original: PhyloTree,
    pruned: PhyloTree,
    *,
    metadata_taxa: set[str] | None = None,
) -> PruningInformationLoss:
    original_clades = informative_rooted_clades(original)
    pruned_clades = informative_rooted_clades(pruned)
    lost_taxa_count = max(0, original.tip_count - pruned.tip_count)
    lost_taxa_fraction = (
        0.0 if original.tip_count == 0 else lost_taxa_count / original.tip_count
    )
    lost_clade_count = max(0, len(original_clades) - len(pruned_clades))
    lost_clade_fraction = (
        0.0 if not original_clades else lost_clade_count / len(original_clades)
    )
    original_total = round(original.total_branch_length(), 15)
    pruned_total = round(pruned.total_branch_length(), 15)
    lost_branch_length = max(0.0, round(original_total - pruned_total, 15))
    lost_branch_length_fraction = (
        0.0 if original_total == 0.0 else lost_branch_length / original_total
    )
    lost_metadata_count = None
    lost_metadata_fraction = None
    if metadata_taxa is not None:
        original_tree_taxa = set(original.tip_names)
        lost_metadata_taxa = original_tree_taxa - metadata_taxa
        lost_metadata_count = len(lost_metadata_taxa)
        lost_metadata_fraction = (
            0.0
            if not original_tree_taxa
            else lost_metadata_count / len(original_tree_taxa)
        )
    return PruningInformationLoss(
        lost_taxa_count=lost_taxa_count,
        lost_taxa_fraction=lost_taxa_fraction,
        lost_clade_count=lost_clade_count,
        lost_clade_fraction=lost_clade_fraction,
        lost_branch_length=lost_branch_length,
        lost_branch_length_fraction=lost_branch_length_fraction,
        lost_metadata_count=lost_metadata_count,
        lost_metadata_fraction=lost_metadata_fraction,
    )


def _collapse_unary(node: TreeNode) -> TreeNode:
    if len(node.children) != 1:
        return node
    only_child = _collapse_unary(node.children[0])
    return TreeNode(
        name=only_child.name,
        branch_length=_merge_branch_lengths(
            node.branch_length, only_child.branch_length
        ),
        children=only_child.children,
        metadata=deepcopy(only_child.metadata),
        edge_metadata=deepcopy(only_child.edge_metadata),
    )


def _prune_node(node: TreeNode, keep_taxa: set[str]) -> TreeNode | None:
    if node.is_leaf():
        if node.name is None or node.name not in keep_taxa:
            return None
        return TreeNode(
            name=node.name,
            branch_length=node.branch_length,
            children=[],
            metadata=deepcopy(node.metadata),
            edge_metadata=deepcopy(node.edge_metadata),
        )

    children = [
        child
        for child in (_prune_node(child, keep_taxa) for child in node.children)
        if child is not None
    ]
    if not children:
        return None
    pruned_node = TreeNode(
        name=node.name,
        branch_length=node.branch_length,
        children=children,
        metadata=deepcopy(node.metadata),
        edge_metadata=deepcopy(node.edge_metadata),
    )
    if len(children) == 1:
        return _collapse_unary(pruned_node)
    return pruned_node


def _prune_tree_against_taxa(
    tree: PhyloTree,
    keep_taxa: set[str],
    *,
    clear_root_branch_length: bool = True,
) -> tuple[PhyloTree, list[str], list[str]]:
    retained_tips = sorted(name for name in tree.tip_names if name in keep_taxa)
    removed_tips = sorted(name for name in tree.tip_names if name not in keep_taxa)

    pruned_root = _prune_node(tree.root, keep_taxa)
    if pruned_root is None or not retained_tips:
        raise MetadataJoinError("no overlapping taxa remain after pruning request")
    if len(retained_tips) < 2:
        raise ValueError("tree pruning requires at least two retained taxa")

    if clear_root_branch_length and not pruned_root.is_leaf():
        pruned_root.branch_length = None
    pruned_tree = PhyloTree(
        root=pruned_root,
        source_format=tree.source_format,
        rooted=_ape_style_rooted_state(pruned_root),
    )
    return pruned_tree, retained_tips, removed_tips


def prune_tree_to_taxa(
    tree_path: Path,
    keep_from_path: Path,
    *,
    taxon_column: str | None = None,
) -> tuple[PhyloTree, TreePruningReport]:
    """Prune a tree to the taxa present in a metadata or traits table."""
    tree = load_tree(tree_path)
    keep_table = load_taxon_table(keep_from_path, taxon_column=taxon_column)
    pruned_tree, retained_tips, removed_tips = _prune_tree_against_taxa(
        tree, set(keep_table.taxa)
    )
    summary = _summarize_transformation(
        tree, pruned_tree, transformation="prune-tree-to-table"
    )
    return pruned_tree, TreePruningReport(
        tree_path=tree_path,
        keep_from_path=keep_from_path,
        taxon_column=keep_table.taxon_column,
        original_tip_count=tree.tip_count,
        kept_taxa=retained_tips,
        removed_taxa=removed_tips,
        removed_taxa_with_reasons=[
            RemovedTaxonReason(taxon=taxon, reason="absent_from_keep_table")
            for taxon in removed_tips
        ],
        pruning_audit=_pruning_artifact_audit(tree, pruned_tree),
        information_loss=_pruning_information_loss(
            tree, pruned_tree, metadata_taxa=set(keep_table.taxa)
        ),
        summary=summary,
    )


def prune_tree_to_requested_taxa(
    tree_path: Path,
    requested_taxa: list[str],
) -> tuple[PhyloTree, RequestedTaxaPruningReport]:
    """Prune a tree to an explicit requested taxon list."""
    tree = load_tree(tree_path)
    requested_set = set(requested_taxa)
    pruned_tree, retained_tips, removed_tips = _prune_tree_against_taxa(
        tree, requested_set
    )
    summary = _summarize_transformation(
        tree, pruned_tree, transformation="prune-tree-to-requested-taxa"
    )
    return pruned_tree, RequestedTaxaPruningReport(
        tree_path=tree_path,
        original_tip_count=tree.tip_count,
        requested_taxa=sorted(requested_set),
        kept_taxa=retained_tips,
        removed_taxa=removed_tips,
        absent_requested_taxa=sorted(requested_set - set(tree.tip_names)),
        removed_taxa_with_reasons=[
            RemovedTaxonReason(taxon=taxon, reason="not_requested")
            for taxon in removed_tips
        ],
        pruning_audit=_pruning_artifact_audit(tree, pruned_tree),
        information_loss=_pruning_information_loss(tree, pruned_tree),
        summary=summary,
    )


def prune_tree_object_to_requested_taxa(
    tree: PhyloTree,
    requested_taxa: list[str],
) -> PhyloTree:
    """Prune one in-memory tree to an explicit requested taxon list."""
    requested_set = set(requested_taxa)
    pruned_tree, _retained_tips, _removed_tips = _prune_tree_against_taxa(
        tree,
        requested_set,
    )
    return pruned_tree


def drop_tree_taxa(
    tree_path: Path,
    excluded_taxa: list[str],
) -> tuple[PhyloTree, RequestedTaxaPruningReport]:
    """Drop a caller-provided exclusion list from a tree."""
    tree = load_tree(tree_path)
    excluded_set = set(excluded_taxa)
    retained_taxa = [name for name in tree.tip_names if name not in excluded_set]
    pruned_tree, retained_tips, removed_tips = _prune_tree_against_taxa(
        tree, set(retained_taxa)
    )
    summary = _summarize_transformation(
        tree, pruned_tree, transformation="drop-tree-taxa"
    )
    return pruned_tree, RequestedTaxaPruningReport(
        tree_path=tree_path,
        original_tip_count=tree.tip_count,
        requested_taxa=sorted(excluded_set),
        kept_taxa=retained_tips,
        removed_taxa=removed_tips,
        absent_requested_taxa=sorted(excluded_set - set(tree.tip_names)),
        removed_taxa_with_reasons=[
            RemovedTaxonReason(
                taxon=taxon,
                reason="excluded_explicitly"
                if taxon in excluded_set
                else "pruned_after_exclusion",
            )
            for taxon in removed_tips
        ],
        pruning_audit=_pruning_artifact_audit(tree, pruned_tree),
        information_loss=_pruning_information_loss(tree, pruned_tree),
        summary=summary,
    )


def write_pruned_taxa(path: Path, removed_taxa: list[str]) -> Path:
    """Write the removed taxa for an explicit pruning run."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["taxon"]
    lines.extend(removed_taxa)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def prune_alignment_to_tree(
    alignment_path: Path,
    tree_path: Path,
) -> tuple[list[AlignmentRecord], AlignmentPruningReport]:
    """Prune an alignment to the taxa present in a tree."""
    records = load_fasta_alignment(alignment_path)
    tree = load_tree(tree_path)
    tree_taxa = set(tree.tip_names)
    kept_records = [record for record in records if record.identifier in tree_taxa]
    if not kept_records:
        raise MetadataJoinError(
            "no overlapping taxa remain after alignment pruning request"
        )

    kept_ids = [record.identifier for record in kept_records]
    removed_ids = [
        record.identifier for record in records if record.identifier not in tree_taxa
    ]
    return kept_records, AlignmentPruningReport(
        alignment_path=alignment_path,
        tree_path=tree_path,
        original_sequence_count=len(records),
        kept_ids=kept_ids,
        removed_ids=removed_ids,
        removed_ids_with_reasons=[
            RemovedTaxonReason(taxon=identifier, reason="absent_from_tree")
            for identifier in removed_ids
        ],
    )


def prune_tree_to_alignment(
    tree_path: Path,
    alignment_path: Path,
) -> tuple[PhyloTree, TreePruningReport]:
    """Prune a tree to the identifiers present in an alignment."""
    tree = load_tree(tree_path)
    records = load_fasta_alignment(alignment_path)
    pruned_tree, retained_tips, removed_tips = _prune_tree_against_taxa(
        tree,
        {record.identifier for record in records},
    )
    summary = _summarize_transformation(
        tree, pruned_tree, transformation="prune-tree-to-alignment"
    )
    return pruned_tree, TreePruningReport(
        tree_path=tree_path,
        keep_from_path=alignment_path,
        taxon_column="identifier",
        original_tip_count=tree.tip_count,
        kept_taxa=retained_tips,
        removed_taxa=removed_tips,
        removed_taxa_with_reasons=[
            RemovedTaxonReason(taxon=taxon, reason="absent_from_alignment")
            for taxon in removed_tips
        ],
        pruning_audit=_pruning_artifact_audit(tree, pruned_tree),
        information_loss=_pruning_information_loss(tree, pruned_tree),
        summary=summary,
    )
