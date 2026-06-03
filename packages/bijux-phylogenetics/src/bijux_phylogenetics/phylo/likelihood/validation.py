from __future__ import annotations

from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import (
    AlignmentTaxonMismatchError,
    InvalidBranchLengthError,
)


def validate_explicit_branch_lengths(
    tree: PhyloTree,
    *,
    model_name: str,
) -> None:
    """Require explicit nonnegative branch lengths on every edge."""
    for _parent, child in tree.iter_edges():
        if child.branch_length is None:
            raise InvalidBranchLengthError(
                f"{model_name} fixed-topology likelihood requires explicit branch lengths on every edge"
            )
        if child.branch_length < 0.0:
            raise InvalidBranchLengthError(
                f"{model_name} likelihood does not accept negative branch lengths"
            )


def validate_tree_taxa_against_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    model_name: str,
) -> None:
    """Require identical named taxon sets between one tree and one alignment."""
    tree_taxa = [leaf.name for leaf in tree.iter_leaves()]
    if any(name is None for name in tree_taxa):
        raise AlignmentTaxonMismatchError(
            f"{model_name} likelihood requires every tree tip to have a matching alignment identifier"
        )
    observed_tree_taxa = [name for name in tree_taxa if name is not None]
    if len(set(observed_tree_taxa)) != len(observed_tree_taxa):
        raise AlignmentTaxonMismatchError(
            f"{model_name} likelihood requires uniquely named tree tips"
        )
    expected_taxa = compressed_patterns.taxon_order
    if set(observed_tree_taxa) != set(expected_taxa):
        missing_from_alignment = sorted(set(observed_tree_taxa) - set(expected_taxa))
        missing_from_tree = sorted(set(expected_taxa) - set(observed_tree_taxa))
        details: list[str] = []
        if missing_from_alignment:
            details.append(f"tree-only taxa: {', '.join(missing_from_alignment)}")
        if missing_from_tree:
            details.append(f"alignment-only taxa: {', '.join(missing_from_tree)}")
        detail_suffix = f" ({'; '.join(details)})" if details else ""
        raise AlignmentTaxonMismatchError(
            f"{model_name} likelihood requires identical tree and alignment taxon sets"
            f"{detail_suffix}"
        )
