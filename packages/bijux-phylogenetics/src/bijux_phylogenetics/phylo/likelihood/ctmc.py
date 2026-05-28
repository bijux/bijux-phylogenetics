from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import math

import numpy

from bijux_phylogenetics.runtime.errors import PhylogeneticsError


@dataclass(frozen=True, slots=True)
class ValidatedCtmcRateMatrix:
    """One validated finite-state CTMC generator matrix."""

    rate_matrix: numpy.ndarray
    state_labels: tuple[str, ...]
    row_sum_tolerance: float

    @property
    def state_count(self) -> int:
        return len(self.state_labels)

    def index_for_state(self, state_label: str) -> int:
        try:
            return self.state_labels.index(state_label)
        except ValueError as error:
            raise KeyError(state_label) from error


def validate_ctmc_rate_matrix(
    rate_matrix: numpy.ndarray | Sequence[Sequence[float]],
    *,
    state_labels: Sequence[str] | None = None,
    row_sum_tolerance: float = 1e-9,
) -> ValidatedCtmcRateMatrix:
    validated_tolerance = _validate_row_sum_tolerance(row_sum_tolerance)
    candidate = numpy.asarray(rate_matrix, dtype=float)
    if candidate.ndim != 2 or candidate.shape[0] != candidate.shape[1]:
        raise PhylogeneticsError(
            "ctmc rate matrix must be square",
            code="ctmc_rate_matrix_not_square",
            details={"shape": list(candidate.shape)},
        )
    if candidate.shape[0] < 1:
        raise PhylogeneticsError(
            "ctmc rate matrix must contain at least one state",
            code="ctmc_rate_matrix_empty",
        )
    resolved_labels = _resolve_state_labels(
        state_count=candidate.shape[0],
        state_labels=state_labels,
    )
    if not numpy.all(numpy.isfinite(candidate)):
        row_index, column_index = _first_invalid_entry(
            candidate,
            predicate=lambda value: not math.isfinite(value),
        )
        raise PhylogeneticsError(
            "ctmc rate matrix must contain only finite values",
            code="ctmc_rate_matrix_value_not_finite",
            details={
                "row_index": row_index,
                "column_index": column_index,
                "state_label": resolved_labels[row_index],
                "value": float(candidate[row_index, column_index]),
            },
        )

    off_diagonal = candidate.copy()
    numpy.fill_diagonal(off_diagonal, 0.0)
    if numpy.any(off_diagonal < 0.0):
        row_index, column_index = _first_invalid_entry(
            off_diagonal,
            predicate=lambda value: value < 0.0,
        )
        raise PhylogeneticsError(
            "ctmc rate matrix must not contain negative off-diagonal rates",
            code="ctmc_rate_matrix_off_diagonal_negative",
            details={
                "row_index": row_index,
                "column_index": column_index,
                "source_state": resolved_labels[row_index],
                "target_state": resolved_labels[column_index],
                "value": float(candidate[row_index, column_index]),
            },
        )

    row_sums = candidate.sum(axis=1)
    if not numpy.allclose(
        row_sums,
        numpy.zeros(candidate.shape[0], dtype=float),
        rtol=0.0,
        atol=validated_tolerance,
    ):
        offending_rows = [
            {
                "row_index": row_index,
                "state_label": resolved_labels[row_index],
                "row_sum": float(row_sum),
            }
            for row_index, row_sum in enumerate(row_sums)
            if not math.isclose(
                float(row_sum),
                0.0,
                rel_tol=0.0,
                abs_tol=validated_tolerance,
            )
        ]
        raise PhylogeneticsError(
            "ctmc rate matrix rows must sum to zero",
            code="ctmc_rate_matrix_row_sums_nonzero",
            details={
                "row_sum_tolerance": validated_tolerance,
                "offending_rows": offending_rows,
            },
        )

    diagonal = numpy.diag(candidate)
    if numpy.any(diagonal > validated_tolerance):
        row_index = _first_positive_diagonal_index(diagonal, validated_tolerance)
        raise PhylogeneticsError(
            "ctmc rate matrix diagonal entries must be zero or negative",
            code="ctmc_rate_matrix_diagonal_positive",
            details={
                "row_index": row_index,
                "state_label": resolved_labels[row_index],
                "value": float(diagonal[row_index]),
                "row_sum_tolerance": validated_tolerance,
            },
        )

    return ValidatedCtmcRateMatrix(
        rate_matrix=candidate.copy(),
        state_labels=resolved_labels,
        row_sum_tolerance=validated_tolerance,
    )


def _resolve_state_labels(
    *,
    state_count: int,
    state_labels: Sequence[str] | None,
) -> tuple[str, ...]:
    if state_labels is None:
        return tuple(f"state_{index}" for index in range(state_count))
    resolved_labels = tuple(state_labels)
    if len(resolved_labels) != state_count:
        raise PhylogeneticsError(
            "ctmc rate matrix state-label count must match matrix dimensions",
            code="ctmc_rate_matrix_state_label_count_mismatch",
            details={
                "state_count": state_count,
                "state_label_count": len(resolved_labels),
            },
        )
    if len(set(resolved_labels)) != len(resolved_labels):
        duplicate_labels = sorted(
            {label for label in resolved_labels if resolved_labels.count(label) > 1}
        )
        raise PhylogeneticsError(
            "ctmc rate matrix state labels must be unique",
            code="ctmc_rate_matrix_state_labels_duplicate",
            details={"duplicate_state_labels": duplicate_labels},
        )
    return resolved_labels


def _validate_row_sum_tolerance(row_sum_tolerance: float) -> float:
    if (
        math.isnan(row_sum_tolerance)
        or math.isinf(row_sum_tolerance)
        or row_sum_tolerance < 0.0
    ):
        raise PhylogeneticsError(
            "ctmc rate matrix row-sum tolerance must be finite and non-negative",
            code="ctmc_rate_matrix_row_sum_tolerance_invalid",
            details={"row_sum_tolerance": row_sum_tolerance},
        )
    return float(row_sum_tolerance)


def _first_invalid_entry(
    matrix: numpy.ndarray,
    *,
    predicate,
) -> tuple[int, int]:
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            if predicate(float(matrix[row_index, column_index])):
                return row_index, column_index
    raise ValueError("matrix does not contain an invalid entry")


def _first_positive_diagonal_index(
    diagonal: numpy.ndarray,
    tolerance: float,
) -> int:
    for row_index, value in enumerate(diagonal):
        if float(value) > tolerance:
            return row_index
    raise ValueError("diagonal does not contain a positive entry")


__all__ = [
    "ValidatedCtmcRateMatrix",
    "validate_ctmc_rate_matrix",
]
