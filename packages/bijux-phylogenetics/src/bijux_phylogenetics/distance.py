from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.alignment import AlignmentRecord
from bijux_phylogenetics.errors import InvalidAlignmentError
from bijux_phylogenetics.io.fasta import infer_alignment_alphabet, load_fasta_alignment

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


def compute_pairwise_genetic_distance_matrix(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
) -> GeneticDistanceMatrix:
    """Compute a deterministic pairwise genetic distance matrix for a DNA alignment."""
    if model != "p-distance":
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
            pairs.append(
                PairwiseGeneticDistance(
                    left_identifier=left.identifier,
                    right_identifier=right.identifier,
                    distance=distance,
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
