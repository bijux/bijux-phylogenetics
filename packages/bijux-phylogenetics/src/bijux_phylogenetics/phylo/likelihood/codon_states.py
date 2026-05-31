from __future__ import annotations

from dataclasses import dataclass

from Bio.Data import CodonTable
import numpy

from bijux_phylogenetics.phylo.likelihood.ctmc import (
    normalize_ctmc_rate_matrix_by_expected_substitution_rate,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError, PhylogeneticsError

_DNA_BASE_ORDER = ("A", "C", "G", "T")


@dataclass(frozen=True, slots=True)
class CodonStateSpace:
    """Resolved non-stop codon state space for one genetic code."""

    genetic_code_id: int
    genetic_code_name: str
    state_order: tuple[str, ...]
    stop_codons: frozenset[str]
    state_index: dict[str, int]


def resolve_codon_state_space(
    genetic_code: int | str | None = None,
) -> CodonStateSpace:
    """Resolve one stable lexicographic sense-codon state order."""
    code_id, code_name, stop_codons = _resolve_genetic_code(genetic_code)
    state_order = tuple(
        codon
        for codon in _iter_lexicographic_codons()
        if codon not in stop_codons
    )
    return CodonStateSpace(
        genetic_code_id=code_id,
        genetic_code_name=code_name,
        state_order=state_order,
        stop_codons=frozenset(stop_codons),
        state_index={codon: index for index, codon in enumerate(state_order)},
    )


def validate_codon_frequency_vector(
    codon_frequencies: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None,
    *,
    state_space: CodonStateSpace,
    owner_name: str,
) -> tuple[numpy.ndarray, str]:
    """Normalize one codon-frequency vector in the resolved state order."""
    if codon_frequencies is None:
        state_count = len(state_space.state_order)
        return (
            numpy.full(state_count, 1.0 / state_count, dtype=float),
            "uniform",
        )

    if isinstance(codon_frequencies, dict):
        unexpected_codons = sorted(
            set(codon_frequencies) - set(state_space.state_order)
        )
        missing_codons = sorted(set(state_space.state_order) - set(codon_frequencies))
        if unexpected_codons or missing_codons:
            details: list[str] = []
            if unexpected_codons:
                details.append(
                    "unexpected codons: " + ", ".join(unexpected_codons[:8])
                )
            if missing_codons:
                details.append("missing codons: " + ", ".join(missing_codons[:8]))
            raise InvalidAlignmentError(
                f"{owner_name} requires frequencies for exactly the resolved sense-codon state order"
                + (f" ({'; '.join(details)})" if details else "")
            )
        vector = numpy.asarray(
            [codon_frequencies[codon] for codon in state_space.state_order],
            dtype=float,
        )
        return _normalize_codon_frequency_vector(
            vector,
            state_space=state_space,
            owner_name=owner_name,
            source="provided",
        )

    vector = numpy.asarray(codon_frequencies, dtype=float)
    if vector.shape != (len(state_space.state_order),):
        raise InvalidAlignmentError(
            f"{owner_name} requires exactly {len(state_space.state_order)} codon frequencies in the resolved state order"
        )
    return _normalize_codon_frequency_vector(
        vector,
        state_space=state_space,
        owner_name=owner_name,
        source="provided",
    )


def build_equal_rate_codon_ctmc_rate_matrix(
    codon_frequencies: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    *,
    genetic_code: int | str | None = None,
) -> tuple[CodonStateSpace, numpy.ndarray, numpy.ndarray, str]:
    """Build one normalized 61-state single-nucleotide codon CTMC generator."""
    state_space = resolve_codon_state_space(genetic_code)
    frequencies, frequency_source = validate_codon_frequency_vector(
        codon_frequencies,
        state_space=state_space,
        owner_name="codon CTMC likelihood",
    )
    state_count = len(state_space.state_order)
    rate_matrix = numpy.zeros((state_count, state_count), dtype=float)
    for row_index, source_codon in enumerate(state_space.state_order):
        for column_index, target_codon in enumerate(state_space.state_order):
            if row_index == column_index:
                continue
            if _single_nucleotide_difference(source_codon, target_codon):
                rate_matrix[row_index, column_index] = frequencies[column_index]
        rate_matrix[row_index, row_index] = -float(rate_matrix[row_index, :].sum())
    try:
        normalized_rate_matrix = normalize_ctmc_rate_matrix_by_expected_substitution_rate(
            rate_matrix,
            frequencies,
            state_labels=state_space.state_order,
        )
    except PhylogeneticsError as error:
        raise InvalidAlignmentError(
            "codon CTMC likelihood requires a frequency vector with positive expected substitution rate"
        ) from error
    return (
        state_space,
        normalized_rate_matrix,
        frequencies,
        frequency_source,
    )


def _normalize_codon_frequency_vector(
    vector: numpy.ndarray,
    *,
    state_space: CodonStateSpace,
    owner_name: str,
    source: str,
) -> tuple[numpy.ndarray, str]:
    if not numpy.all(numpy.isfinite(vector)):
        raise InvalidAlignmentError(
            f"{owner_name} codon frequencies must contain only finite values"
        )
    if numpy.any(vector < 0.0):
        raise InvalidAlignmentError(
            f"{owner_name} codon frequencies must be nonnegative"
        )
    total = float(vector.sum())
    if total <= 0.0:
        raise InvalidAlignmentError(
            f"{owner_name} codon frequencies must sum to a positive value"
        )
    normalized = vector / total
    positive_frequency_count = int(numpy.count_nonzero(normalized > 0.0))
    if positive_frequency_count < 2:
        raise InvalidAlignmentError(
            f"{owner_name} codon frequencies must assign positive weight to at least two sense codons"
        )
    if normalized.shape != (len(state_space.state_order),):
        raise InvalidAlignmentError(
            f"{owner_name} codon frequencies must match the resolved sense-codon state order"
        )
    return normalized, source


def _resolve_genetic_code(
    genetic_code: int | str | None,
) -> tuple[int, str, set[str]]:
    if genetic_code is None:
        table = CodonTable.unambiguous_dna_by_id[1]
        return 1, table.names[0], set(table.stop_codons)
    if isinstance(genetic_code, int):
        try:
            table = CodonTable.unambiguous_dna_by_id[genetic_code]
        except KeyError as error:
            raise InvalidAlignmentError(
                f"unsupported genetic code id for codon likelihood: {genetic_code}"
            ) from error
        return genetic_code, table.names[0], set(table.stop_codons)
    text = genetic_code.strip()
    if not text:
        raise InvalidAlignmentError(
            "codon likelihood genetic code name must not be empty"
        )
    if text.isdigit():
        return _resolve_genetic_code(int(text))
    lowered = text.lower()
    for table in CodonTable.unambiguous_dna_by_id.values():
        if lowered in {name.lower() for name in table.names}:
            return int(table.id), table.names[0], set(table.stop_codons)
    raise InvalidAlignmentError(
        f"unsupported genetic code for codon likelihood: {genetic_code}"
    )


def _iter_lexicographic_codons() -> tuple[str, ...]:
    return tuple(
        first + second + third
        for first in _DNA_BASE_ORDER
        for second in _DNA_BASE_ORDER
        for third in _DNA_BASE_ORDER
    )


def _single_nucleotide_difference(left: str, right: str) -> bool:
    return sum(1 for left_base, right_base in zip(left, right, strict=True) if left_base != right_base) == 1
