from __future__ import annotations

from pathlib import Path

from Bio.Phylo.TreeConstruction import DistanceMatrix

from bijux_phylogenetics.io.fasta import (
    infer_alignment_alphabet,
    load_fasta_alignment,
)
from bijux_phylogenetics.io.fasta.matrix import load_dna_bin_alignment
from bijux_phylogenetics.phylo.alignment import AlignmentRecord, DnaBinAlignment
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

from ..models import (
    AmbiguityPolicy,
    DistanceModel,
    GapHandlingMode,
    GeneticDistanceMatrix,
    GeneticDistanceModelParameters,
    PairwiseGeneticDistance,
)
from ..shared import (
    _allowed_models_for_alphabet,
    _normalize_distance_model,
    _pair_key,
)
from .correction_models import _estimate_nucleotide_model_parameters
from .pairwise_analysis import _pair_summary
from .site_policy import _complete_deletion_positions


def _load_alignment_for_model(
    path: Path, *, model: DistanceModel
) -> tuple[list[AlignmentRecord], str]:
    normalized_model = _normalize_distance_model(model)
    if normalized_model != "amino-acid-p-distance":
        matrix = load_dna_bin_alignment(path, normalize_uracil=True)
        if normalized_model not in _allowed_models_for_alphabet(matrix.source_alphabet):
            raise InvalidAlignmentError(
                f"distance model '{normalized_model}' is not supported for inferred alphabet '{matrix.source_alphabet}'"
            )
        records = [
            AlignmentRecord(
                identifier=record.identifier,
                sequence=record.sequence.upper(),
            )
            for record in matrix.records
        ]
        return records, matrix.source_alphabet

    records = load_fasta_alignment(path)
    alphabet = infer_alignment_alphabet(records)
    if normalized_model not in _allowed_models_for_alphabet(alphabet):
        raise InvalidAlignmentError(
            f"distance model '{normalized_model}' is not supported for inferred alphabet '{alphabet}'"
        )
    return records, alphabet


def _validate_gap_and_ambiguity_policy(
    *,
    gap_handling: GapHandlingMode,
    ambiguity_policy: AmbiguityPolicy,
) -> None:
    if gap_handling not in {"pairwise-deletion", "complete-deletion"}:
        raise ValueError(f"unsupported gap handling mode: {gap_handling}")
    if ambiguity_policy not in {
        "ignore",
        "partial-match",
        "strict-mismatch",
        "report-only",
    }:
        raise ValueError(f"unsupported ambiguity policy: {ambiguity_policy}")


def _build_pairwise_distance_rows(
    records: list[AlignmentRecord],
    *,
    alphabet: str,
    model: DistanceModel,
    ambiguity_policy: AmbiguityPolicy,
    retained_positions: list[int] | None,
    model_parameters: GeneticDistanceModelParameters | None,
) -> list[PairwiseGeneticDistance]:
    pairs: list[PairwiseGeneticDistance] = []
    for left_index, left in enumerate(records):
        for right_index, right in enumerate(records):
            if right_index < left_index:
                continue
            summary = _pair_summary(
                left.sequence,
                right.sequence,
                alphabet=alphabet,
                ambiguity_policy=ambiguity_policy,
                retained_positions=retained_positions,
                model=model,
                model_parameters=model_parameters,
            )
            pairs.append(
                PairwiseGeneticDistance(
                    left_identifier=left.identifier,
                    right_identifier=right.identifier,
                    distance=summary.distance if left_index != right_index else 0.0,
                    comparable_sites=summary.comparable_sites,
                    mismatch_sites=summary.mismatch_sites
                    if left_index != right_index
                    else 0.0,
                    transition_sites=summary.transition_sites
                    if left_index != right_index
                    else 0.0,
                    ag_transition_sites=summary.ag_transition_sites
                    if left_index != right_index
                    else 0.0,
                    ct_transition_sites=summary.ct_transition_sites
                    if left_index != right_index
                    else 0.0,
                    transversion_sites=summary.transversion_sites
                    if left_index != right_index
                    else 0.0,
                    ambiguity_sites=summary.ambiguity_sites,
                    skipped_sites=summary.skipped_sites,
                    saturated=False if left_index == right_index else summary.saturated,
                    saturation_reason=None
                    if left_index == right_index
                    else summary.saturation_reason,
                )
            )
    return pairs


def compute_pairwise_genetic_distance_matrix_from_dna_bin_alignment(
    alignment: DnaBinAlignment,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> GeneticDistanceMatrix:
    """Compute a deterministic nucleotide distance matrix from one DNAbin-compatible alignment."""
    model = _normalize_distance_model(model)
    if model not in _allowed_models_for_alphabet(alignment.source_alphabet):
        raise InvalidAlignmentError(
            f"distance model '{model}' is not supported for inferred alphabet '{alignment.source_alphabet}'"
        )
    if model == "amino-acid-p-distance":
        raise InvalidAlignmentError(
            "dnabin-compatible nucleotide distance loading does not support amino-acid distances"
        )
    _validate_gap_and_ambiguity_policy(
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )

    records = [
        AlignmentRecord(
            identifier=record.identifier,
            sequence=record.sequence.upper(),
        )
        for record in alignment.records
    ]
    alphabet = alignment.source_alphabet
    model_parameters: GeneticDistanceModelParameters | None = None
    warnings: list[str] = []
    if model in {"felsenstein-81", "tamura-nei-93"}:
        model_parameters, warnings = _estimate_nucleotide_model_parameters(records)
    retained_positions = (
        _complete_deletion_positions(
            records, alphabet=alphabet, ambiguity_policy=ambiguity_policy
        )
        if gap_handling == "complete-deletion"
        else None
    )
    return GeneticDistanceMatrix(
        path=alignment.path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        inferred_alphabet=alphabet,
        alignment_length=alignment.alignment_length,
        identifiers=[record.identifier for record in records],
        model_parameters=model_parameters,
        warnings=warnings,
        pairs=_build_pairwise_distance_rows(
            records,
            alphabet=alphabet,
            model=model,
            ambiguity_policy=ambiguity_policy,
            retained_positions=retained_positions,
            model_parameters=model_parameters,
        ),
    )


def _build_alignment_distance_lookup(
    report: GeneticDistanceMatrix,
) -> dict[tuple[str, str], float]:
    distances: dict[tuple[str, str], float] = {}
    for pair in report.pairs:
        if pair.left_identifier == pair.right_identifier or pair.distance is None:
            continue
        distances[_pair_key(pair.left_identifier, pair.right_identifier)] = float(
            pair.distance
        )
    return distances


def _bio_distance_matrix(report: GeneticDistanceMatrix) -> DistanceMatrix:
    undefined_pairs = [
        f"{pair.left_identifier}/{pair.right_identifier}"
        for pair in report.pairs
        if pair.distance is None
    ]
    if undefined_pairs:
        raise InvalidAlignmentError(
            "distance matrix contains undefined entries for: "
            + ", ".join(undefined_pairs)
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
                if {pair.left_identifier, pair.right_identifier}
                == {left_identifier, right_identifier}
            )
            row.append(float(pair.distance))
        rows.append(row)
    return DistanceMatrix(report.identifiers, rows)


def _distance_lookup(report: GeneticDistanceMatrix) -> dict[tuple[str, str], float]:
    lookup: dict[tuple[str, str], float] = {}
    for identifier in report.identifiers:
        lookup[(identifier, identifier)] = 0.0
    for pair in report.pairs:
        if pair.distance is None:
            raise InvalidAlignmentError(
                "distance matrix contains undefined entries for: "
                f"{pair.left_identifier}/{pair.right_identifier}"
            )
        lookup[(pair.left_identifier, pair.right_identifier)] = float(pair.distance)
        lookup[(pair.right_identifier, pair.left_identifier)] = float(pair.distance)
    return lookup


def compute_pairwise_genetic_distance_matrix(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> GeneticDistanceMatrix:
    """Compute a deterministic pairwise genetic distance matrix for an aligned dataset."""
    model = _normalize_distance_model(model)
    if model != "amino-acid-p-distance":
        matrix = load_dna_bin_alignment(path, normalize_uracil=True)
        return compute_pairwise_genetic_distance_matrix_from_dna_bin_alignment(
            matrix,
            model=model,
            gap_handling=gap_handling,
            ambiguity_policy=ambiguity_policy,
        )

    _validate_gap_and_ambiguity_policy(
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    records, alphabet = _load_alignment_for_model(path, model=model)
    retained_positions = (
        _complete_deletion_positions(
            records, alphabet=alphabet, ambiguity_policy=ambiguity_policy
        )
        if gap_handling == "complete-deletion"
        else None
    )
    return GeneticDistanceMatrix(
        path=path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        inferred_alphabet=alphabet,
        alignment_length=len(records[0].sequence),
        identifiers=[record.identifier for record in records],
        model_parameters=None,
        warnings=[],
        pairs=_build_pairwise_distance_rows(
            records,
            alphabet=alphabet,
            model=model,
            ambiguity_policy=ambiguity_policy,
            retained_positions=retained_positions,
            model_parameters=None,
        ),
    )
