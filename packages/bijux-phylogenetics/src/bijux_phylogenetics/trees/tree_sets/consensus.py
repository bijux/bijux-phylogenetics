from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick, write_newick
from bijux_phylogenetics.phylo.topology.clades import informative_rooted_clades
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .contracts import ConsensusTreeReport
from .inventory import (
    _analyze_tree_set,
    _clade_branch_lengths,
    _require_exact_taxa,
    _terminal_branch_lengths,
    _TreeSetAnalysis,
    _validate_same_taxa,
)


def _maximal_nested_clades(
    parent: frozenset[str], clades: set[frozenset[str]]
) -> list[frozenset[str]]:
    nested = [clade for clade in clades if clade < parent]
    return sorted(
        [
            clade
            for clade in nested
            if not any(clade < other < parent for other in nested)
        ],
        key=lambda clade: (len(clade), sorted(clade)),
    )


def _mean(values: list[float]) -> float:
    return round(sum(values) / len(values), 15)


def _build_consensus_node(
    taxa: frozenset[str],
    *,
    included_clades: set[frozenset[str]],
    clade_support: dict[frozenset[str], float],
    clade_lengths: dict[frozenset[str], float],
    terminal_lengths: dict[str, float],
    is_root: bool = False,
) -> TreeNode:
    child_clades = _maximal_nested_clades(taxa, included_clades)
    covered: set[str] = set()
    children: list[TreeNode] = []
    for child_clade in child_clades:
        covered.update(child_clade)
        children.append(
            _build_consensus_node(
                child_clade,
                included_clades=included_clades,
                clade_support=clade_support,
                clade_lengths=clade_lengths,
                terminal_lengths=terminal_lengths,
            )
        )
    for taxon in sorted(taxa - covered):
        children.append(TreeNode(name=taxon, branch_length=terminal_lengths.get(taxon)))
    if len(children) == 1:
        return children[0]
    label = None if is_root else format(clade_support[taxa], ".15g")
    return TreeNode(
        name=label,
        branch_length=None if is_root else clade_lengths.get(taxa),
        children=children,
    )


def _build_consensus_tree_from_selected_clades(
    *,
    shared_taxa: list[str],
    source_format: str,
    included_clades: set[frozenset[str]],
    clade_support: dict[frozenset[str], float],
    clade_lengths: dict[frozenset[str], float],
    terminal_lengths: dict[str, float],
) -> PhyloTree:
    return PhyloTree(
        root=_build_consensus_node(
            frozenset(shared_taxa),
            included_clades=included_clades,
            clade_support=clade_support,
            clade_lengths=clade_lengths,
            terminal_lengths=terminal_lengths,
            is_root=True,
        ),
        source_format=source_format,
        rooted=True,
    )


def _build_consensus_tree_with_threshold(
    analysis: _TreeSetAnalysis,
    *,
    threshold: float,
) -> tuple[PhyloTree, ConsensusTreeReport]:
    shared_taxa = _require_exact_taxa(analysis)
    counts = analysis.clade_counts or {}
    majority_clades = {
        clade
        for clade, count in counts.items()
        if count / len(analysis.trees) >= threshold
    }
    clade_support = {
        clade: round((counts[clade] / len(analysis.trees)) * 100.0, 15)
        for clade in majority_clades
    }
    clade_lengths = {
        clade: _mean(lengths)
        for clade, lengths in analysis.clade_branch_lengths.items()
        if clade in majority_clades and lengths
    }
    terminal_length_means = {
        taxon: _mean(lengths)
        for taxon, lengths in analysis.terminal_lengths.items()
        if lengths
    }
    tree = _build_consensus_tree_from_selected_clades(
        shared_taxa=shared_taxa,
        source_format=analysis.source_format,
        included_clades=majority_clades,
        clade_support=clade_support,
        clade_lengths=clade_lengths,
        terminal_lengths=terminal_length_means,
    )
    if math.isclose(threshold, 1.0):
        consensus_method = "strict"
    elif math.isclose(threshold, 0.5):
        consensus_method = "majority-rule"
    else:
        consensus_method = "thresholded"
    return tree, ConsensusTreeReport(
        path=analysis.path,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=shared_taxa,
        consensus_method=consensus_method,
        consensus_threshold=threshold,
        included_clade_count=len(majority_clades),
        consensus_newick=dumps_newick(tree),
    )


def _build_consensus_tree_with_threshold_from_trees(
    trees: list[PhyloTree],
    *,
    threshold: float,
) -> tuple[PhyloTree, int]:
    shared_taxa = _validate_same_taxa(trees)
    shared_taxa_set = set(shared_taxa)
    counts: dict[frozenset[str], int] = {}
    clade_branch_lengths: dict[frozenset[str], list[float]] = {}
    terminal_lengths: dict[str, list[float]] = {}
    for tree in trees:
        for clade in informative_rooted_clades(tree, shared_taxa_set):
            counts[clade] = counts.get(clade, 0) + 1
        for clade, length in _clade_branch_lengths(tree, shared_taxa_set).items():
            if length is not None:
                clade_branch_lengths.setdefault(clade, []).append(float(length))
        for taxon, length in _terminal_branch_lengths(tree).items():
            if length is not None:
                terminal_lengths.setdefault(taxon, []).append(float(length))
    majority_clades = {
        clade for clade, count in counts.items() if count / len(trees) >= threshold
    }
    clade_support = {
        clade: round((counts[clade] / len(trees)) * 100.0, 15)
        for clade in majority_clades
    }
    clade_lengths = {
        clade: _mean(lengths)
        for clade, lengths in clade_branch_lengths.items()
        if clade in majority_clades and lengths
    }
    terminal_length_means = {
        taxon: _mean(lengths) for taxon, lengths in terminal_lengths.items() if lengths
    }
    tree = _build_consensus_tree_from_selected_clades(
        shared_taxa=shared_taxa,
        source_format=trees[0].source_format,
        included_clades=majority_clades,
        clade_support=clade_support,
        clade_lengths=clade_lengths,
        terminal_lengths=terminal_length_means,
    )
    return tree, len(majority_clades)


def compute_consensus_tree(path: Path) -> tuple[PhyloTree, ConsensusTreeReport]:
    """Compute a majority-rule consensus tree from a tree set."""
    return compute_consensus_tree_with_threshold(path, threshold=0.5)


def compute_strict_consensus_tree(path: Path) -> tuple[PhyloTree, ConsensusTreeReport]:
    """Compute a strict consensus tree from a tree set."""
    return compute_consensus_tree_with_threshold(path, threshold=1.0)


def compute_consensus_tree_with_threshold(
    path: Path,
    *,
    threshold: float,
) -> tuple[PhyloTree, ConsensusTreeReport]:
    """Compute a deterministic consensus tree at a caller-supplied clade frequency threshold."""
    if not 0.0 < threshold <= 1.0:
        raise ValueError(
            f"consensus threshold must be greater than 0 and at most 1, got {threshold}"
        )
    return _build_consensus_tree_with_threshold(
        _analyze_tree_set(path), threshold=threshold
    )


def write_consensus_tree(path: Path, tree: PhyloTree) -> Path:
    """Write a consensus tree as canonical Newick."""
    return write_newick(path, tree)
