from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
import tracemalloc

from bijux_phylogenetics.io.newick import (
    dumps_newick,
    iter_newick_tree_records_from_path,
    load_newick_tree_set,
    loads_newick,
)
from bijux_phylogenetics.io.trees import detect_tree_format
from bijux_phylogenetics.phylo.topology.clades import informative_rooted_clades
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidAlignmentError,
    UnsupportedTreeFormatError,
)

from .contracts import TreeSetProcessingSummary, TreeSetRecord, TreeSetReport
from .topology import _rooted_topology_id, _unrooted_topology_id


@dataclass(slots=True)
class _TreeSetAnalysis:
    path: Path
    source_format: str
    processing: TreeSetProcessingSummary
    trees: list[PhyloTree]
    shared_taxa: list[str]
    taxa_union: list[str]
    exact_taxa: list[str] | None
    records: list[TreeSetRecord]
    rooted_topology_counts: dict[str, int]
    unrooted_topology_counts: dict[str, int]
    rooted_representatives: dict[str, tuple[int, str, PhyloTree]]
    clade_counts: dict[frozenset[str], int] | None
    clade_branch_lengths: dict[frozenset[str], list[float]]
    terminal_lengths: dict[str, list[float]]


def _shared_taxa(trees: list[PhyloTree]) -> set[str]:
    shared = set(trees[0].tip_names)
    for tree in trees[1:]:
        shared &= set(tree.tip_names)
    return shared


def _taxa_union(trees: list[PhyloTree]) -> set[str]:
    taxa: set[str] = set()
    for tree in trees:
        taxa.update(tree.tip_names)
    return taxa


def _iter_tree_set(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"tree-set file not found: {path}")
    source_format = detect_tree_format(path)
    if source_format != "newick":
        raise UnsupportedTreeFormatError(
            f"tree-set workflows require Newick tree-set records, got {source_format} for {path}"
        )
    for source_index, statement in iter_newick_tree_records_from_path(path):
        try:
            tree = loads_newick(statement)
        except Exception as error:  # pragma: no cover - parser exceptions vary
            yield source_format, source_index, None, str(error)
            continue
        yield source_format, source_index, tree, None


def _exact_taxa_or_none(trees: list[PhyloTree]) -> list[str] | None:
    first = sorted(trees[0].tip_names)
    for tree in trees[1:]:
        if sorted(tree.tip_names) != first:
            return None
    return first


def _processing_summary(
    *, started: float, started_tracing: bool
) -> TreeSetProcessingSummary:
    _current, peak = tracemalloc.get_traced_memory()
    if not started_tracing:
        tracemalloc.stop()
    return TreeSetProcessingSummary(
        runtime_seconds=round(perf_counter() - started, 6),
        peak_memory_bytes=peak,
        skipped_malformed_tree_count=0,
    )


def _clade_branch_lengths(
    tree: PhyloTree, shared_taxa: set[str]
) -> dict[frozenset[str], float | None]:
    lengths: dict[frozenset[str], float | None] = {}

    def visit(node: TreeNode) -> set[str]:
        if node.is_leaf():
            return (
                {node.name}
                if node.name is not None and node.name in shared_taxa
                else set()
            )
        taxa: set[str] = set()
        for child in node.children:
            taxa.update(visit(child))
        if 1 < len(taxa) < len(shared_taxa):
            lengths[frozenset(taxa)] = node.branch_length
        return taxa

    visit(tree.root)
    return lengths


def _terminal_branch_lengths(tree: PhyloTree) -> dict[str, float | None]:
    return {
        name: length
        for name, length in tree.terminal_branch_lengths()
        if name is not None
    }


def _analyze_tree_set(path: Path) -> _TreeSetAnalysis:
    started = perf_counter()
    started_tracing = tracemalloc.is_tracing()
    if not started_tracing:
        tracemalloc.start()

    source_format: str | None = None
    skipped_malformed_tree_count = 0
    trees: list[PhyloTree] = []
    try:
        for parsed_format, _source_index, tree, error_message in _iter_tree_set(path):
            source_format = parsed_format
            if error_message is not None:
                skipped_malformed_tree_count += 1
                continue
            if tree is None:
                continue
            trees.append(tree)
        if not trees:
            raise InvalidAlignmentError(f"tree set contains no trees: {path}")
        shared_taxa = sorted(_shared_taxa(trees))
        shared_taxa_set = set(shared_taxa)
        taxa_union = sorted(_taxa_union(trees))
        records: list[TreeSetRecord] = []
        rooted_topology_counts: dict[str, int] = {}
        unrooted_topology_counts: dict[str, int] = {}
        rooted_representatives: dict[str, tuple[int, str, PhyloTree]] = {}
        for index, tree in enumerate(trees, start=1):
            rooted_topology_id = _rooted_topology_id(tree, shared_taxa_set)
            unrooted_topology_id = _unrooted_topology_id(tree, shared_taxa_set)
            records.append(
                TreeSetRecord(
                    index=index,
                    tip_count=tree.tip_count,
                    taxa=sorted(tree.tip_names),
                    rooted_topology_id=rooted_topology_id,
                    unrooted_topology_id=unrooted_topology_id,
                )
            )
            rooted_topology_counts[rooted_topology_id] = (
                rooted_topology_counts.get(rooted_topology_id, 0) + 1
            )
            unrooted_topology_counts[unrooted_topology_id] = (
                unrooted_topology_counts.get(unrooted_topology_id, 0) + 1
            )
            rooted_representatives.setdefault(
                rooted_topology_id, (index, dumps_newick(tree), tree)
            )

        exact_taxa = _exact_taxa_or_none(trees)
        clade_counts: dict[frozenset[str], int] | None = None
        clade_branch_lengths: dict[frozenset[str], list[float]] = {}
        terminal_lengths: dict[str, list[float]] = {}
        if exact_taxa is not None:
            exact_taxa_set = set(exact_taxa)
            clade_counts = {}
            for tree in trees:
                for clade in informative_rooted_clades(tree, exact_taxa_set):
                    clade_counts[clade] = clade_counts.get(clade, 0) + 1
                for clade, length in _clade_branch_lengths(
                    tree, exact_taxa_set
                ).items():
                    if length is not None:
                        clade_branch_lengths.setdefault(clade, []).append(float(length))
                for taxon, length in _terminal_branch_lengths(tree).items():
                    if length is not None:
                        terminal_lengths.setdefault(taxon, []).append(float(length))
    finally:
        processing = _processing_summary(
            started=started, started_tracing=started_tracing
        )
    processing = TreeSetProcessingSummary(
        runtime_seconds=processing.runtime_seconds,
        peak_memory_bytes=processing.peak_memory_bytes,
        skipped_malformed_tree_count=skipped_malformed_tree_count,
    )
    return _TreeSetAnalysis(
        path=path,
        source_format=source_format or detect_tree_format(path),
        processing=processing,
        trees=trees,
        shared_taxa=shared_taxa,
        taxa_union=taxa_union,
        exact_taxa=exact_taxa,
        records=records,
        rooted_topology_counts=rooted_topology_counts,
        unrooted_topology_counts=unrooted_topology_counts,
        rooted_representatives=rooted_representatives,
        clade_counts=clade_counts,
        clade_branch_lengths=clade_branch_lengths,
        terminal_lengths=terminal_lengths,
    )


def _require_exact_taxa(analysis: _TreeSetAnalysis) -> list[str]:
    if analysis.exact_taxa is None:
        raise InvalidAlignmentError(
            "tree-set analysis requires all trees to share the exact same taxon set"
        )
    return analysis.exact_taxa


def _require_tree_set(path: Path) -> tuple[str, list[PhyloTree]]:
    if not path.exists():
        raise FileNotFoundError(f"tree-set file not found: {path}")
    source_format = detect_tree_format(path)
    if source_format != "newick":
        raise UnsupportedTreeFormatError(
            f"tree-set workflows require Newick tree-set records, got {source_format} for {path}"
        )
    trees = load_newick_tree_set(path)
    if not trees:
        raise InvalidAlignmentError(f"tree set contains no trees: {path}")
    return source_format, trees


def _validate_same_taxa(trees: list[PhyloTree]) -> list[str]:
    first = sorted(trees[0].tip_names)
    for tree in trees[1:]:
        if sorted(tree.tip_names) != first:
            raise InvalidAlignmentError(
                "tree-set analysis requires all trees to share the exact same taxon set"
            )
    return first


def load_tree_set(path: Path) -> TreeSetReport:
    """Read a set of trees and summarize their topology diversity over shared taxa."""
    analysis = _analyze_tree_set(path)
    return TreeSetReport(
        path=path,
        source_format=analysis.source_format,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=analysis.shared_taxa,
        taxa_union=analysis.taxa_union,
        rooted_topology_count=len(analysis.rooted_topology_counts),
        unrooted_topology_count=len(analysis.unrooted_topology_counts),
        records=analysis.records,
    )
