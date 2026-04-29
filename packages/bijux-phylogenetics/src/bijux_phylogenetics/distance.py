from __future__ import annotations

import csv
from dataclasses import dataclass
import math
from pathlib import Path

from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from bijux_phylogenetics.core.alignment import AlignmentRecord
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.errors import InvalidAlignmentError, InvalidDistanceMatrixError
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


@dataclass(frozen=True, slots=True)
class ImportedDistanceEntry:
    """One directional entry from an imported long-form distance matrix table."""

    left_identifier: str
    right_identifier: str
    distance: float
    comparable_sites: int | None


@dataclass(frozen=True, slots=True)
class DistanceMatrixAsymmetry:
    """Two directional entries that disagree numerically."""

    left_identifier: str
    right_identifier: str
    left_to_right_distance: float
    right_to_left_distance: float


@dataclass(frozen=True, slots=True)
class NonMetricDistanceObservation:
    """One triangle-inequality violation within an imported distance matrix."""

    left_identifier: str
    middle_identifier: str
    right_identifier: str
    direct_distance: float
    indirect_distance: float


@dataclass(slots=True)
class ImportedDistanceMatrixReport:
    """Validation report for an imported long-form distance matrix table."""

    path: Path
    identifiers: list[str]
    pair_count: int
    complete: bool
    zero_diagonal: bool
    symmetric: bool
    nonnegative: bool
    missing_pairs: list[str]
    diagonal_problems: list[str]
    negative_distance_pairs: list[str]
    asymmetric_pairs: list[DistanceMatrixAsymmetry]
    nonmetric_observations: list[NonMetricDistanceObservation]
    warnings: list[str]


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


def _pair_key(left_identifier: str, right_identifier: str) -> tuple[str, str]:
    return tuple(sorted((left_identifier, right_identifier)))


def load_imported_distance_matrix(path: Path) -> list[ImportedDistanceEntry]:
    """Load a long-form imported distance matrix table."""
    if not path.exists():
        raise FileNotFoundError(f"distance matrix file not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = {"left_identifier", "right_identifier", "distance"}
        if reader.fieldnames is None or not required_columns <= set(reader.fieldnames):
            raise InvalidDistanceMatrixError(
                "distance matrix must contain left_identifier, right_identifier, and distance columns"
            )
        comparable_sites_column = "comparable_sites" if "comparable_sites" in reader.fieldnames else None
        entries: list[ImportedDistanceEntry] = []
        seen_directional_pairs: set[tuple[str, str]] = set()
        for row_index, row in enumerate(reader, start=2):
            left_identifier = str(row.get("left_identifier", "")).strip()
            right_identifier = str(row.get("right_identifier", "")).strip()
            raw_distance = str(row.get("distance", "")).strip()
            if not left_identifier or not right_identifier:
                raise InvalidDistanceMatrixError(
                    f"row {row_index} in {path} is missing a left_identifier or right_identifier value"
                )
            if not raw_distance:
                raise InvalidDistanceMatrixError(f"row {row_index} in {path} is missing a distance value")
            try:
                distance = float(raw_distance)
            except ValueError as error:
                raise InvalidDistanceMatrixError(
                    f"row {row_index} in {path} contains a non-numeric distance value '{raw_distance}'"
                ) from error
            comparable_sites: int | None = None
            if comparable_sites_column is not None:
                raw_comparable_sites = str(row.get(comparable_sites_column, "")).strip()
                if raw_comparable_sites:
                    try:
                        comparable_sites = int(raw_comparable_sites)
                    except ValueError as error:
                        raise InvalidDistanceMatrixError(
                            f"row {row_index} in {path} contains a non-integer comparable_sites value '{raw_comparable_sites}'"
                        ) from error
            directional_pair = (left_identifier, right_identifier)
            if directional_pair in seen_directional_pairs:
                raise InvalidDistanceMatrixError(
                    f"distance matrix contains duplicate directional entry {left_identifier}/{right_identifier}"
                )
            seen_directional_pairs.add(directional_pair)
            entries.append(
                ImportedDistanceEntry(
                    left_identifier=left_identifier,
                    right_identifier=right_identifier,
                    distance=round(distance, 15),
                    comparable_sites=comparable_sites,
                )
            )

    if not entries:
        raise InvalidDistanceMatrixError(f"distance matrix contains no rows: {path}")
    return entries


def _symmetric_distances(entries: list[ImportedDistanceEntry]) -> dict[tuple[str, str], float]:
    distances: dict[tuple[str, str], float] = {}
    by_direction = {(entry.left_identifier, entry.right_identifier): entry.distance for entry in entries}
    identifiers = sorted({entry.left_identifier for entry in entries} | {entry.right_identifier for entry in entries})
    for left_identifier in identifiers:
        for right_identifier in identifiers:
            pair_key = _pair_key(left_identifier, right_identifier)
            if pair_key in distances:
                continue
            if left_identifier == right_identifier:
                if (left_identifier, right_identifier) in by_direction:
                    distances[pair_key] = by_direction[(left_identifier, right_identifier)]
                continue
            if (left_identifier, right_identifier) in by_direction:
                distances[pair_key] = by_direction[(left_identifier, right_identifier)]
            elif (right_identifier, left_identifier) in by_direction:
                distances[pair_key] = by_direction[(right_identifier, left_identifier)]
    return distances


def validate_imported_distance_matrix(path: Path) -> ImportedDistanceMatrixReport:
    """Validate a long-form imported distance matrix table."""
    entries = load_imported_distance_matrix(path)
    identifiers = sorted({entry.left_identifier for entry in entries} | {entry.right_identifier for entry in entries})
    by_direction = {(entry.left_identifier, entry.right_identifier): entry for entry in entries}

    missing_pairs: list[str] = []
    diagonal_problems: list[str] = []
    negative_distance_pairs: list[str] = []
    asymmetric_pairs: list[DistanceMatrixAsymmetry] = []
    symmetric_distances = _symmetric_distances(entries)

    for left_identifier in identifiers:
        diagonal = by_direction.get((left_identifier, left_identifier))
        if diagonal is None:
            missing_pairs.append(f"{left_identifier}/{left_identifier}")
        elif diagonal.distance != 0.0:
            diagonal_problems.append(f"{left_identifier}/{left_identifier} has diagonal distance {diagonal.distance:g}")
        if diagonal is not None and diagonal.distance < 0:
            negative_distance_pairs.append(f"{left_identifier}/{left_identifier}")

        for right_identifier in identifiers:
            if left_identifier >= right_identifier:
                continue
            left_to_right = by_direction.get((left_identifier, right_identifier))
            right_to_left = by_direction.get((right_identifier, left_identifier))
            if left_to_right is None and right_to_left is None:
                missing_pairs.append(f"{left_identifier}/{right_identifier}")
                continue
            if left_to_right is not None and left_to_right.distance < 0:
                negative_distance_pairs.append(f"{left_identifier}/{right_identifier}")
            if right_to_left is not None and right_to_left.distance < 0:
                negative_distance_pairs.append(f"{right_identifier}/{left_identifier}")
            if left_to_right is not None and right_to_left is not None and not math.isclose(
                left_to_right.distance,
                right_to_left.distance,
                rel_tol=1e-12,
                abs_tol=1e-12,
            ):
                asymmetric_pairs.append(
                    DistanceMatrixAsymmetry(
                        left_identifier=left_identifier,
                        right_identifier=right_identifier,
                        left_to_right_distance=left_to_right.distance,
                        right_to_left_distance=right_to_left.distance,
                    )
                )

    nonmetric_observations: list[NonMetricDistanceObservation] = []
    if not missing_pairs and not diagonal_problems and not negative_distance_pairs:
        for left_index, left_identifier in enumerate(identifiers):
            for middle_index, middle_identifier in enumerate(identifiers):
                if middle_index == left_index:
                    continue
                for right_index, right_identifier in enumerate(identifiers):
                    if len({left_index, middle_index, right_index}) < 3:
                        continue
                    if left_identifier > right_identifier:
                        continue
                    direct_distance = symmetric_distances.get(_pair_key(left_identifier, right_identifier))
                    left_middle = symmetric_distances.get(_pair_key(left_identifier, middle_identifier))
                    middle_right = symmetric_distances.get(_pair_key(middle_identifier, right_identifier))
                    if direct_distance is None or left_middle is None or middle_right is None:
                        continue
                    indirect_distance = round(left_middle + middle_right, 15)
                    if direct_distance > indirect_distance + 1e-12:
                        nonmetric_observations.append(
                            NonMetricDistanceObservation(
                                left_identifier=left_identifier,
                                middle_identifier=middle_identifier,
                                right_identifier=right_identifier,
                                direct_distance=direct_distance,
                                indirect_distance=indirect_distance,
                            )
                        )

    warnings: list[str] = []
    if missing_pairs:
        warnings.append("distance matrix is missing one or more required pairs")
    if diagonal_problems:
        warnings.append("distance matrix contains nonzero diagonal entries")
    if negative_distance_pairs:
        warnings.append("distance matrix contains negative distances")
    if asymmetric_pairs:
        warnings.append("distance matrix contains asymmetric directional entries")
    if nonmetric_observations:
        warnings.append("distance matrix violates triangle inequality for one or more taxon triples")

    return ImportedDistanceMatrixReport(
        path=path,
        identifiers=identifiers,
        pair_count=len(entries),
        complete=not missing_pairs,
        zero_diagonal=not diagonal_problems,
        symmetric=not asymmetric_pairs,
        nonnegative=not negative_distance_pairs,
        missing_pairs=missing_pairs,
        diagonal_problems=diagonal_problems,
        negative_distance_pairs=sorted(set(negative_distance_pairs)),
        asymmetric_pairs=sorted(asymmetric_pairs, key=lambda row: (row.left_identifier, row.right_identifier)),
        nonmetric_observations=sorted(
            nonmetric_observations,
            key=lambda row: (row.left_identifier, row.middle_identifier, row.right_identifier),
        ),
        warnings=warnings,
    )


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
