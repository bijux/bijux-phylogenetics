from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .topology_search import (
    normalize_nucleotide_topology_search_records,
    resolve_nucleotide_topology_search_records,
    resolve_nucleotide_topology_search_tree,
    validate_nucleotide_topology_search_tree,
)
from .validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)


def validate_nucleotide_likelihood_starting_tree(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    model_name: str,
    workflow_name: str = "nucleotide likelihood starting-tree validation",
) -> None:
    """Require one structurally valid scored-start tree before likelihood search."""
    model_label = f"{model_name.strip().upper()} starting-tree"
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name=model_label,
    )
    validate_explicit_branch_lengths(tree, model_name=model_label)
    validate_nucleotide_topology_search_tree(tree, workflow_name=workflow_name)


def validate_nucleotide_likelihood_starting_tree_from_alignment(
    tree: PhyloTree | Path,
    records: list[AlignmentRecord] | Path,
    *,
    model_name: str,
    workflow_name: str = "nucleotide likelihood starting-tree validation",
) -> None:
    """Validate one likelihood start tree against one DNA alignment input."""
    resolved_tree, _resolved_tree_path = resolve_nucleotide_topology_search_tree(tree)
    resolved_records, _resolved_alignment_path = (
        resolve_nucleotide_topology_search_records(records)
    )
    _normalized_records, compressed_patterns = (
        normalize_nucleotide_topology_search_records(
            resolved_records,
            owner_name=workflow_name,
        )
    )
    validate_nucleotide_likelihood_starting_tree(
        resolved_tree,
        compressed_patterns,
        model_name=model_name,
        workflow_name=workflow_name,
    )
