from __future__ import annotations

from collections.abc import Mapping, Sequence
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


@dataclass(frozen=True, slots=True)
class SolvedCtmcStationaryDistribution:
    """One stationary distribution solved or verified for a CTMC generator."""

    state_labels: tuple[str, ...]
    probabilities: tuple[float, ...]
    residual_by_state: tuple[float, ...]
    normalization_error: float
    probability_tolerance: float

    @property
    def state_count(self) -> int:
        return len(self.state_labels)

    def probability_for(self, state_label: str) -> float:
        try:
            index = self.state_labels.index(state_label)
        except ValueError as error:
            raise KeyError(state_label) from error
        return self.probabilities[index]

    def as_mapping(self) -> dict[str, float]:
        return dict(zip(self.state_labels, self.probabilities, strict=True))


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


def solve_ctmc_stationary_distribution(
    rate_matrix: numpy.ndarray | Sequence[Sequence[float]] | ValidatedCtmcRateMatrix,
    *,
    state_labels: Sequence[str] | None = None,
    row_sum_tolerance: float = 1e-9,
    probability_tolerance: float = 1e-9,
) -> SolvedCtmcStationaryDistribution:
    validated_rate_matrix = _resolve_validated_rate_matrix(
        rate_matrix,
        state_labels=state_labels,
        row_sum_tolerance=row_sum_tolerance,
    )
    validated_probability_tolerance = _validate_probability_tolerance(
        probability_tolerance
    )
    _validate_stationary_distribution_uniqueness(
        validated_rate_matrix=validated_rate_matrix,
        probability_tolerance=validated_probability_tolerance,
    )
    stationary_vector = _solve_stationary_vector(validated_rate_matrix.rate_matrix)
    _validate_stationary_vector(
        stationary_vector,
        validated_rate_matrix=validated_rate_matrix,
        probability_tolerance=validated_probability_tolerance,
    )
    residual_vector = stationary_vector @ validated_rate_matrix.rate_matrix
    return SolvedCtmcStationaryDistribution(
        state_labels=validated_rate_matrix.state_labels,
        probabilities=tuple(float(value) for value in stationary_vector),
        residual_by_state=tuple(float(value) for value in residual_vector),
        normalization_error=abs(float(numpy.sum(stationary_vector)) - 1.0),
        probability_tolerance=validated_probability_tolerance,
    )


def compute_ctmc_expected_substitution_rate(
    rate_matrix: numpy.ndarray | Sequence[Sequence[float]] | ValidatedCtmcRateMatrix,
    stationary_distribution: Mapping[str, float] | Sequence[float] | numpy.ndarray,
    *,
    state_labels: Sequence[str] | None = None,
    row_sum_tolerance: float = 1e-9,
    probability_tolerance: float = 1e-9,
) -> float:
    validated_rate_matrix = _resolve_validated_rate_matrix(
        rate_matrix,
        state_labels=state_labels,
        row_sum_tolerance=row_sum_tolerance,
    )
    validated_probability_tolerance = _validate_probability_tolerance(
        probability_tolerance
    )
    probability_vector = _resolve_probability_vector_for_expected_rate(
        stationary_distribution,
        validated_rate_matrix=validated_rate_matrix,
        probability_tolerance=validated_probability_tolerance,
    )
    expected_rate = -float(
        numpy.sum(probability_vector * numpy.diag(validated_rate_matrix.rate_matrix))
    )
    if expected_rate <= 0.0 or not math.isfinite(expected_rate):
        raise PhylogeneticsError(
            "ctmc rate matrix requires a positive finite expected substitution rate",
            code="ctmc_expected_substitution_rate_nonpositive",
            details={
                "expected_substitution_rate": expected_rate,
                "state_labels": list(validated_rate_matrix.state_labels),
            },
        )
    return expected_rate


def normalize_ctmc_rate_matrix_by_expected_substitution_rate(
    rate_matrix: numpy.ndarray | Sequence[Sequence[float]] | ValidatedCtmcRateMatrix,
    stationary_distribution: Mapping[str, float] | Sequence[float] | numpy.ndarray,
    *,
    state_labels: Sequence[str] | None = None,
    row_sum_tolerance: float = 1e-9,
    probability_tolerance: float = 1e-9,
) -> numpy.ndarray:
    validated_rate_matrix = _resolve_validated_rate_matrix(
        rate_matrix,
        state_labels=state_labels,
        row_sum_tolerance=row_sum_tolerance,
    )
    expected_rate = compute_ctmc_expected_substitution_rate(
        validated_rate_matrix,
        stationary_distribution,
        probability_tolerance=probability_tolerance,
    )
    return validated_rate_matrix.rate_matrix / expected_rate


def verify_ctmc_stationary_distribution(
    rate_matrix: numpy.ndarray | Sequence[Sequence[float]] | ValidatedCtmcRateMatrix,
    stationary_distribution: Mapping[str, float] | Sequence[float] | numpy.ndarray,
    *,
    state_labels: Sequence[str] | None = None,
    row_sum_tolerance: float = 1e-9,
    probability_tolerance: float = 1e-9,
) -> SolvedCtmcStationaryDistribution:
    validated_rate_matrix = _resolve_validated_rate_matrix(
        rate_matrix,
        state_labels=state_labels,
        row_sum_tolerance=row_sum_tolerance,
    )
    validated_probability_tolerance = _validate_probability_tolerance(
        probability_tolerance
    )
    _validate_stationary_distribution_uniqueness(
        validated_rate_matrix=validated_rate_matrix,
        probability_tolerance=validated_probability_tolerance,
    )
    candidate_vector = _resolve_stationary_distribution_vector(
        stationary_distribution,
        validated_rate_matrix=validated_rate_matrix,
    )
    _validate_stationary_vector(
        candidate_vector,
        validated_rate_matrix=validated_rate_matrix,
        probability_tolerance=validated_probability_tolerance,
    )
    residual_vector = candidate_vector @ validated_rate_matrix.rate_matrix
    return SolvedCtmcStationaryDistribution(
        state_labels=validated_rate_matrix.state_labels,
        probabilities=tuple(float(value) for value in candidate_vector),
        residual_by_state=tuple(float(value) for value in residual_vector),
        normalization_error=abs(float(numpy.sum(candidate_vector)) - 1.0),
        probability_tolerance=validated_probability_tolerance,
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


def _resolve_validated_rate_matrix(
    rate_matrix: numpy.ndarray | Sequence[Sequence[float]] | ValidatedCtmcRateMatrix,
    *,
    state_labels: Sequence[str] | None,
    row_sum_tolerance: float,
) -> ValidatedCtmcRateMatrix:
    if isinstance(rate_matrix, ValidatedCtmcRateMatrix):
        if state_labels is not None and tuple(state_labels) != rate_matrix.state_labels:
            raise PhylogeneticsError(
                "ctmc stationary-distribution state labels must match the validated rate matrix",
                code="ctmc_stationary_distribution_state_labels_mismatch",
                details={
                    "validated_state_labels": list(rate_matrix.state_labels),
                    "requested_state_labels": list(state_labels),
                },
            )
        return rate_matrix
    return validate_ctmc_rate_matrix(
        rate_matrix,
        state_labels=state_labels,
        row_sum_tolerance=row_sum_tolerance,
    )


def _solve_stationary_vector(rate_matrix: numpy.ndarray) -> numpy.ndarray:
    state_count = rate_matrix.shape[0]
    linear_system = rate_matrix.T.copy()
    linear_system[-1, :] = 1.0
    targets = numpy.zeros(state_count, dtype=float)
    targets[-1] = 1.0
    try:
        stationary_vector = numpy.linalg.solve(linear_system, targets)
    except numpy.linalg.LinAlgError:
        stationary_vector, *_ = numpy.linalg.lstsq(linear_system, targets, rcond=None)
    return numpy.real_if_close(stationary_vector, tol=1000).astype(float)


def _validate_stationary_distribution_uniqueness(
    *,
    validated_rate_matrix: ValidatedCtmcRateMatrix,
    probability_tolerance: float,
) -> None:
    singular_values = numpy.linalg.svd(
        validated_rate_matrix.rate_matrix.T,
        compute_uv=False,
    )
    rank_tolerance = _stationary_rank_tolerance(
        rate_matrix=validated_rate_matrix.rate_matrix,
        row_sum_tolerance=validated_rate_matrix.row_sum_tolerance,
        probability_tolerance=probability_tolerance,
    )
    nullity = sum(
        1
        for singular_value in singular_values
        if float(singular_value) <= rank_tolerance
    )
    if nullity == 1:
        return
    raise PhylogeneticsError(
        "ctmc stationary distribution is not unique for this rate matrix",
        code="ctmc_stationary_distribution_not_unique",
        details={
            "state_labels": list(validated_rate_matrix.state_labels),
            "nullity": nullity,
            "rank_tolerance": rank_tolerance,
            "singular_values": [float(value) for value in singular_values],
        },
    )


def _resolve_stationary_distribution_vector(
    stationary_distribution: Mapping[str, float] | Sequence[float] | numpy.ndarray,
    *,
    validated_rate_matrix: ValidatedCtmcRateMatrix,
) -> numpy.ndarray:
    if isinstance(stationary_distribution, Mapping):
        unexpected_states = sorted(
            state
            for state in stationary_distribution
            if state not in set(validated_rate_matrix.state_labels)
        )
        if unexpected_states:
            raise PhylogeneticsError(
                "ctmc stationary distribution contains unexpected states",
                code="ctmc_stationary_distribution_unexpected_states",
                details={"unexpected_states": unexpected_states},
            )
        missing_states = [
            state
            for state in validated_rate_matrix.state_labels
            if state not in stationary_distribution
        ]
        if missing_states:
            raise PhylogeneticsError(
                "ctmc stationary distribution is missing required states",
                code="ctmc_stationary_distribution_missing_states",
                details={"missing_states": missing_states},
            )
        return numpy.array(
            [
                float(stationary_distribution[state])
                for state in validated_rate_matrix.state_labels
            ],
            dtype=float,
        )
    vector = numpy.asarray(stationary_distribution, dtype=float)
    if vector.ndim != 1 or vector.shape[0] != validated_rate_matrix.state_count:
        raise PhylogeneticsError(
            "ctmc stationary distribution must match the CTMC state count",
            code="ctmc_stationary_distribution_length_mismatch",
            details={
                "state_count": validated_rate_matrix.state_count,
                "distribution_shape": list(vector.shape),
            },
        )
    return vector


def _resolve_probability_vector_for_expected_rate(
    stationary_distribution: Mapping[str, float] | Sequence[float] | numpy.ndarray,
    *,
    validated_rate_matrix: ValidatedCtmcRateMatrix,
    probability_tolerance: float,
) -> numpy.ndarray:
    candidate = _resolve_stationary_distribution_vector(
        stationary_distribution,
        validated_rate_matrix=validated_rate_matrix,
    ).astype(float, copy=True)
    if not numpy.all(numpy.isfinite(candidate)):
        raise PhylogeneticsError(
            "ctmc expected-rate probabilities must contain only finite values",
            code="ctmc_expected_substitution_rate_probability_not_finite",
        )
    if numpy.any(candidate < -probability_tolerance):
        offending_indices = [
            index
            for index, value in enumerate(candidate)
            if float(value) < -probability_tolerance
        ]
        raise PhylogeneticsError(
            "ctmc expected-rate probabilities must not contain negative values",
            code="ctmc_expected_substitution_rate_probability_negative",
            details={
                "offending_states": [
                    validated_rate_matrix.state_labels[index]
                    for index in offending_indices
                ],
                "probabilities": [
                    float(candidate[index]) for index in offending_indices
                ],
                "probability_tolerance": probability_tolerance,
            },
        )
    candidate[candidate < 0.0] = 0.0
    total_probability = float(numpy.sum(candidate))
    normalization_error = abs(total_probability - 1.0)
    if normalization_error > probability_tolerance:
        raise PhylogeneticsError(
            "ctmc expected-rate probabilities must sum to one within tolerance",
            code="ctmc_expected_substitution_rate_probability_not_normalized",
            details={
                "total_probability": total_probability,
                "expected_total": 1.0,
                "absolute_error": normalization_error,
                "probability_tolerance": probability_tolerance,
            },
        )
    return candidate / total_probability


def _validate_stationary_vector(
    stationary_vector: numpy.ndarray,
    *,
    validated_rate_matrix: ValidatedCtmcRateMatrix,
    probability_tolerance: float,
) -> None:
    if not numpy.all(numpy.isfinite(stationary_vector)):
        raise PhylogeneticsError(
            "ctmc stationary distribution must contain only finite values",
            code="ctmc_stationary_distribution_value_not_finite",
        )
    if numpy.any(stationary_vector < -probability_tolerance):
        offending_indices = [
            index
            for index, value in enumerate(stationary_vector)
            if float(value) < -probability_tolerance
        ]
        raise PhylogeneticsError(
            "ctmc stationary distribution must not contain negative probabilities",
            code="ctmc_stationary_distribution_negative",
            details={
                "offending_states": [
                    validated_rate_matrix.state_labels[index]
                    for index in offending_indices
                ],
                "probabilities": [
                    float(stationary_vector[index]) for index in offending_indices
                ],
                "probability_tolerance": probability_tolerance,
            },
        )
    clipped_vector = stationary_vector.copy()
    clipped_vector[clipped_vector < 0.0] = 0.0
    normalization_error = abs(float(numpy.sum(clipped_vector)) - 1.0)
    if normalization_error > probability_tolerance:
        raise PhylogeneticsError(
            "ctmc stationary distribution must sum to one within tolerance",
            code="ctmc_stationary_distribution_not_normalized",
            details={
                "total_probability": float(numpy.sum(clipped_vector)),
                "expected_total": 1.0,
                "absolute_error": normalization_error,
                "probability_tolerance": probability_tolerance,
            },
        )
    residual_vector = clipped_vector @ validated_rate_matrix.rate_matrix
    max_abs_residual = max(abs(float(value)) for value in residual_vector)
    if max_abs_residual > max(
        validated_rate_matrix.row_sum_tolerance,
        probability_tolerance,
    ):
        raise PhylogeneticsError(
            "ctmc stationary distribution does not satisfy pi Q = 0 within tolerance",
            code="ctmc_stationary_distribution_residual_nonzero",
            details={
                "max_absolute_residual": max_abs_residual,
                "residual_by_state": {
                    state_label: float(residual)
                    for state_label, residual in zip(
                        validated_rate_matrix.state_labels,
                        residual_vector,
                        strict=True,
                    )
                },
                "probability_tolerance": probability_tolerance,
                "row_sum_tolerance": validated_rate_matrix.row_sum_tolerance,
            },
        )
    stationary_vector[:] = clipped_vector / float(numpy.sum(clipped_vector))


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


def _validate_probability_tolerance(probability_tolerance: float) -> float:
    if (
        math.isnan(probability_tolerance)
        or math.isinf(probability_tolerance)
        or probability_tolerance < 0.0
    ):
        raise PhylogeneticsError(
            "ctmc stationary-distribution tolerance must be finite and non-negative",
            code="ctmc_stationary_distribution_tolerance_invalid",
            details={"probability_tolerance": probability_tolerance},
        )
    return float(probability_tolerance)


def _stationary_rank_tolerance(
    *,
    rate_matrix: numpy.ndarray,
    row_sum_tolerance: float,
    probability_tolerance: float,
) -> float:
    matrix_scale = max(
        1.0,
        float(numpy.linalg.norm(rate_matrix, ord=2)),
        float(numpy.max(numpy.abs(rate_matrix))),
    )
    return (
        max(row_sum_tolerance, probability_tolerance)
        * matrix_scale
        * rate_matrix.shape[0]
    )


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
    "compute_ctmc_expected_substitution_rate",
    "normalize_ctmc_rate_matrix_by_expected_substitution_rate",
    "SolvedCtmcStationaryDistribution",
    "ValidatedCtmcRateMatrix",
    "solve_ctmc_stationary_distribution",
    "validate_ctmc_rate_matrix",
    "verify_ctmc_stationary_distribution",
]
