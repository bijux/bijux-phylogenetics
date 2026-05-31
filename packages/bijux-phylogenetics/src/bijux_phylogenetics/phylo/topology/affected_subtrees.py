from __future__ import annotations

from bijux_phylogenetics.phylo.topology.clades import canonical_clade_id, split_sort_key
from bijux_phylogenetics.phylo.topology.models import AffectedSubtreeReport

from .tree import PhyloTree, descendant_taxa


def summarize_affected_subtrees(
    original_tree: PhyloTree,
    moved_tree: PhyloTree,
) -> AffectedSubtreeReport:
    """Summarize changed and unchanged non-root branch clades across one move."""
    original_signatures = _nonroot_branch_signatures(original_tree)
    moved_signatures = _nonroot_branch_signatures(moved_tree)
    retired_signatures = original_signatures - moved_signatures
    introduced_signatures = moved_signatures - original_signatures
    unaffected_signatures = original_signatures & moved_signatures
    return AffectedSubtreeReport(
        original_branch_clade_ids=_ordered_clade_ids(original_signatures),
        moved_branch_clade_ids=_ordered_clade_ids(moved_signatures),
        retired_branch_clade_ids=_ordered_clade_ids(retired_signatures),
        introduced_branch_clade_ids=_ordered_clade_ids(introduced_signatures),
        affected_branch_clade_ids=_ordered_clade_ids(
            retired_signatures | introduced_signatures
        ),
        unaffected_branch_clade_ids=_ordered_clade_ids(unaffected_signatures),
    )


def _nonroot_branch_signatures(tree: PhyloTree) -> set[frozenset[str]]:
    return {
        frozenset(descendant_taxa(node))
        for node in tree.iter_nodes(order="preorder")
        if node is not tree.root
    }


def _ordered_clade_ids(signatures: set[frozenset[str]]) -> list[str]:
    return [
        canonical_clade_id(signature)
        for signature in sorted(signatures, key=split_sort_key)
    ]
