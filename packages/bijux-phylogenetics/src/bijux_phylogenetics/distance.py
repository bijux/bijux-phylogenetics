from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from bijux_phylogenetics.core.alignment import AlignmentRecord
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.errors import InvalidAlignmentError
from bijux_phylogenetics.io.fasta import infer_alignment_alphabet, load_fasta_alignment
from bijux_phylogenetics.io.biopython import tree_from_biophylo

DistanceModel = str
GapHandlingMode = str

_COMPARISON_NUCLEOTIDES = {"A", "C", "G", "T"}
_MISSING_OR_GAP = {"-", "?", "N", "X"}


@dataclass(frozen=True, slots=True)
class PairwiseGeneticDistance:
    """One pairwise genetic distance entry for an aligned DNA dataset."""

    left_identifier: str
    right_identifier: str
    distance: float | None
    comparable_sites: int


@dataclass(slots=True)
class GeneticDistanceMatrix:
    """Deterministic pairwise genetic distance matrix for one alignment."""

    path: Path
    model: DistanceModel
    gap_handling: GapHandlingMode
    identifiers: list[str]
    pairs: list[PairwiseGeneticDistance]


@dataclass(slots=True)
class DistanceTreeBuildReport:
    """Explicit report for a distance-based tree build."""

    alignment_path: Path
    model: DistanceModel
    gap_handling: GapHandlingMode
    method: str
    taxon_count: int
    pair_count: int


@dataclass(slots=True)
class DistanceTreeTopologyComparison:
    """Topology comparison between NJ and UPGMA trees built from one alignment."""

    alignment_path: Path
    model: DistanceModel
    gap_handling: GapHandlingMode
    shared_taxa: list[str]
    nj_informative_clades: int
    upgma_informative_clades: int
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    topology_equal: bool
    same_unrooted_topology: bool
    same_taxa_different_rooting: bool


def _normalize_residue(residue: str) -> str:
    upper = residue.upper()
    if upper == "U":
        return "T"
    return upper


def _load_dna_alignment(path: Path) -> list[AlignmentRecord]:
    records = load_fasta_alignment(path)
    alphabet = infer_alignment_alphabet(records)
    if alphabet not in {"dna", "rna"}:
        raise InvalidAlignmentError(
            f"genetic distance analysis requires a nucleotide alignment, got inferred alphabet '{alphabet}'"
        )
    return records


def _pairwise_distance(left: str, right: str) -> tuple[float | None, int]:
    comparable_pairs = [
        (_normalize_residue(left_residue), _normalize_residue(right_residue))
        for left_residue, right_residue in zip(left, right, strict=True)
        if _normalize_residue(left_residue) in _COMPARISON_NUCLEOTIDES
        and _normalize_residue(right_residue) in _COMPARISON_NUCLEOTIDES
        and _normalize_residue(left_residue) not in _MISSING_OR_GAP
        and _normalize_residue(right_residue) not in _MISSING_OR_GAP
    ]
    if not comparable_pairs:
        return None, 0
    mismatches = sum(1 for left_residue, right_residue in comparable_pairs if left_residue != right_residue)
    comparable_sites = len(comparable_pairs)
    return round(mismatches / comparable_sites, 15), comparable_sites


def _complete_deletion_positions(records: list[AlignmentRecord]) -> list[int]:
    return [
        position
        for position in range(len(records[0].sequence))
        if all(_normalize_residue(record.sequence[position]) in _COMPARISON_NUCLEOTIDES for record in records)
    ]


def _distance_over_positions(left: str, right: str, positions: list[int]) -> tuple[float | None, int]:
    if not positions:
        return None, 0
    mismatches = sum(1 for position in positions if _normalize_residue(left[position]) != _normalize_residue(right[position]))
    comparable_sites = len(positions)
    return round(mismatches / comparable_sites, 15), comparable_sites


def _jukes_cantor_distance(p_distance: float | None) -> float | None:
    if p_distance is None:
        return None
    if p_distance == 0.0:
        return 0.0
    if p_distance >= 0.75:
        return None
    return round((-3.0 / 4.0) * math.log(1.0 - (4.0 * p_distance / 3.0)), 15)


def _bio_distance_matrix(report: GeneticDistanceMatrix) -> DistanceMatrix:
    undefined_pairs = [
        f"{pair.left_identifier}/{pair.right_identifier}"
        for pair in report.pairs
        if pair.distance is None
    ]
    if undefined_pairs:
        raise InvalidAlignmentError(
            "distance matrix contains undefined entries for: " + ", ".join(undefined_pairs)
        )
    rows: list[list[float]] = []
    for row_index, left_identifier in enumerate(report.identifiers):
        row: list[float] = []
        for right_identifier in report.identifiers[: row_index + 1]:
            if left_identifier == right_identifier:
                row.append(0.0)
                continue
            pair = next(
                pair
                for pair in report.pairs
                if {
                    pair.left_identifier,
                    pair.right_identifier,
                } == {left_identifier, right_identifier}
            )
            row.append(float(pair.distance))
        rows.append(row)
    return DistanceMatrix(report.identifiers, rows)


def _informative_clades(tree: PhyloTree, shared_taxa: set[str]) -> set[frozenset[str]]:
    clades: set[frozenset[str]] = set()

    def visit(node: TreeNode) -> set[str]:
        if node.is_leaf():
            return {node.name} if node.name in shared_taxa else set()

        taxa: set[str] = set()
        for child in node.children:
            taxa.update(visit(child))

        if 1 < len(taxa) < len(shared_taxa):
            clades.add(frozenset(taxa))
        return taxa

    visit(tree.root)
    return clades


def _canonical_bipartition(taxa: set[str], universe: set[str]) -> frozenset[str]:
    complement = universe - taxa
    left = sorted(taxa)
    right = sorted(complement)
    if (len(left), left) <= (len(right), right):
        return frozenset(taxa)
    return frozenset(complement)


def _unrooted_splits(tree: PhyloTree, shared_taxa: set[str]) -> set[frozenset[str]]:
    splits: set[frozenset[str]] = set()

    def visit(node: TreeNode) -> set[str]:
        if node.is_leaf():
            return {node.name} if node.name in shared_taxa else set()

        taxa: set[str] = set()
        for child in node.children:
            taxa.update(visit(child))

        if node is not tree.root and 1 < len(taxa) < len(shared_taxa) - 1:
            splits.add(_canonical_bipartition(taxa, shared_taxa))
        return taxa

    visit(tree.root)
    return splits


def compute_pairwise_genetic_distance_matrix(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
) -> GeneticDistanceMatrix:
    """Compute a deterministic pairwise genetic distance matrix for a DNA alignment."""
    if model not in {"p-distance", "jukes-cantor"}:
        raise ValueError(f"unsupported distance model: {model}")
    if gap_handling not in {"pairwise-deletion", "complete-deletion"}:
        raise ValueError(f"unsupported gap handling mode: {gap_handling}")

    records = _load_dna_alignment(path)
    retained_positions = _complete_deletion_positions(records) if gap_handling == "complete-deletion" else None
    pairs: list[PairwiseGeneticDistance] = []
    for left_index, left in enumerate(records):
        for right_index, right in enumerate(records):
            if right_index < left_index:
                continue
            if left_index == right_index:
                if gap_handling == "complete-deletion":
                    comparable_sites = len(retained_positions or [])
                else:
                    comparable_sites = sum(
                        1
                        for residue in left.sequence
                        if _normalize_residue(residue) in _COMPARISON_NUCLEOTIDES
                    )
                pairs.append(
                    PairwiseGeneticDistance(
                        left_identifier=left.identifier,
                        right_identifier=right.identifier,
                        distance=0.0,
                        comparable_sites=comparable_sites,
                    )
                )
                continue
            if gap_handling == "complete-deletion":
                distance, comparable_sites = _distance_over_positions(
                    left.sequence,
                    right.sequence,
                    retained_positions or [],
                )
            else:
                distance, comparable_sites = _pairwise_distance(left.sequence, right.sequence)
            transformed_distance = _jukes_cantor_distance(distance) if model == "jukes-cantor" else distance
            pairs.append(
                PairwiseGeneticDistance(
                    left_identifier=left.identifier,
                    right_identifier=right.identifier,
                    distance=transformed_distance,
                    comparable_sites=comparable_sites,
                )
            )
    return GeneticDistanceMatrix(
        path=path,
        model=model,
        gap_handling=gap_handling,
        identifiers=[record.identifier for record in records],
        pairs=pairs,
    )


def write_genetic_distance_matrix(path: Path, report: GeneticDistanceMatrix) -> Path:
    """Write a pairwise genetic distance matrix as a deterministic TSV."""
    rows = {(pair.left_identifier, pair.right_identifier): pair for pair in report.pairs}
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["left_identifier\tright_identifier\tdistance\tcomparable_sites"]
    for left in report.identifiers:
        for right in report.identifiers:
            pair = rows.get((left, right)) or rows.get((right, left))
            if pair is None:
                continue
            distance = "" if pair.distance is None else format(pair.distance, ".15g")
            lines.append(f"{left}\t{right}\t{distance}\t{pair.comparable_sites}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def build_distance_tree(
    path: Path,
    *,
    method: str,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
) -> tuple[PhyloTree, DistanceTreeBuildReport]:
    """Build a distance-based tree from an aligned DNA dataset."""
    report = compute_pairwise_genetic_distance_matrix(
        path,
        model=model,
        gap_handling=gap_handling,
    )
    if len(report.identifiers) < 2:
        raise InvalidAlignmentError("distance tree building requires at least two taxa")

    constructor = DistanceTreeConstructor()
    distance_matrix = _bio_distance_matrix(report)
    if method == "neighbor-joining":
        tree = constructor.nj(distance_matrix)
    elif method == "upgma":
        tree = constructor.upgma(distance_matrix)
    else:
        raise ValueError(f"unsupported tree-building method: {method}")

    return tree_from_biophylo(tree, source_format="newick"), DistanceTreeBuildReport(
        alignment_path=path,
        model=model,
        gap_handling=gap_handling,
        method=method,
        taxon_count=len(report.identifiers),
        pair_count=len(report.pairs),
    )


def compare_distance_tree_topologies(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
) -> DistanceTreeTopologyComparison:
    """Compare NJ and UPGMA topologies built from the same alignment."""
    nj_tree, _ = build_distance_tree(
        path,
        method="neighbor-joining",
        model=model,
        gap_handling=gap_handling,
    )
    upgma_tree, _ = build_distance_tree(
        path,
        method="upgma",
        model=model,
        gap_handling=gap_handling,
    )
    shared_taxa = set(nj_tree.tip_names) & set(upgma_tree.tip_names)
    nj_clades = _informative_clades(nj_tree, shared_taxa)
    upgma_clades = _informative_clades(upgma_tree, shared_taxa)
    symmetric_difference = nj_clades.symmetric_difference(upgma_clades)
    denominator = len(nj_clades) + len(upgma_clades)
    topology_equal = len(symmetric_difference) == 0
    same_unrooted_topology = _unrooted_splits(nj_tree, shared_taxa) == _unrooted_splits(upgma_tree, shared_taxa)
    return DistanceTreeTopologyComparison(
        alignment_path=path,
        model=model,
        gap_handling=gap_handling,
        shared_taxa=sorted(shared_taxa),
        nj_informative_clades=len(nj_clades),
        upgma_informative_clades=len(upgma_clades),
        robinson_foulds_distance=len(symmetric_difference),
        normalized_robinson_foulds=0.0 if denominator == 0 else len(symmetric_difference) / denominator,
        topology_equal=topology_equal,
        same_unrooted_topology=same_unrooted_topology,
        same_taxa_different_rooting=topology_equal is False and same_unrooted_topology,
    )
