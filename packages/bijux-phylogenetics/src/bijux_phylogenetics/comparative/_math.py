from __future__ import annotations

import math


def identity_matrix(size: int) -> list[list[float]]:
    """Return an identity matrix of the requested size."""
    return [
        [1.0 if row_index == column_index else 0.0 for column_index in range(size)]
        for row_index in range(size)
    ]


def matrix_copy(matrix: list[list[float]]) -> list[list[float]]:
    """Copy a matrix without sharing nested row lists."""
    return [list(row) for row in matrix]


def transpose(matrix: list[list[float]]) -> list[list[float]]:
    """Transpose a rectangular matrix."""
    if not matrix:
        return []
    return [list(column) for column in zip(*matrix, strict=False)]


def dot(left: list[float], right: list[float]) -> float:
    """Return the dot product of two vectors."""
    return sum(
        left_value * right_value
        for left_value, right_value in zip(left, right, strict=True)
    )


def matrix_vector_multiply(
    matrix: list[list[float]], vector: list[float]
) -> list[float]:
    """Multiply a matrix by a vector."""
    return [dot(row, vector) for row in matrix]


def matrix_multiply(
    left: list[list[float]], right: list[list[float]]
) -> list[list[float]]:
    """Multiply two matrices."""
    right_transposed = transpose(right)
    return [
        [dot(left_row, right_column) for right_column in right_transposed]
        for left_row in left
    ]


def quadratic_form(vector: list[float], matrix: list[list[float]]) -> float:
    """Evaluate x' A x for a vector x and matrix A."""
    return dot(vector, matrix_vector_multiply(matrix, vector))


def invert_matrix(matrix: list[list[float]]) -> list[list[float]]:
    """Invert a square matrix with Gauss-Jordan elimination and partial pivoting."""
    size = len(matrix)
    if size == 0:
        return []
    augmented = [
        [float(value) for value in row] + identity_row
        for row, identity_row in zip(
            matrix_copy(matrix), identity_matrix(size), strict=True
        )
    ]
    for pivot_index in range(size):
        pivot_row = max(
            range(pivot_index, size),
            key=lambda row_index: abs(augmented[row_index][pivot_index]),
        )
        pivot_value = augmented[pivot_row][pivot_index]
        if math.isclose(pivot_value, 0.0, abs_tol=1e-12):
            raise ValueError("matrix is singular and cannot be inverted")
        if pivot_row != pivot_index:
            augmented[pivot_index], augmented[pivot_row] = (
                augmented[pivot_row],
                augmented[pivot_index],
            )
        pivot_value = augmented[pivot_index][pivot_index]
        augmented[pivot_index] = [
            value / pivot_value for value in augmented[pivot_index]
        ]
        for row_index in range(size):
            if row_index == pivot_index:
                continue
            factor = augmented[row_index][pivot_index]
            if math.isclose(factor, 0.0, abs_tol=1e-15):
                continue
            augmented[row_index] = [
                row_value - factor * pivot_value
                for row_value, pivot_value in zip(
                    augmented[row_index], augmented[pivot_index], strict=True
                )
            ]
    return [row[size:] for row in augmented]


def log_determinant(matrix: list[list[float]]) -> float:
    """Return the natural log determinant of a square matrix."""
    size = len(matrix)
    if size == 0:
        return 0.0
    working = matrix_copy(matrix)
    sign = 1.0
    log_abs_det = 0.0
    for pivot_index in range(size):
        pivot_row = max(
            range(pivot_index, size),
            key=lambda row_index: abs(working[row_index][pivot_index]),
        )
        pivot_value = working[pivot_row][pivot_index]
        if math.isclose(pivot_value, 0.0, abs_tol=1e-12):
            raise ValueError("matrix determinant is zero")
        if pivot_row != pivot_index:
            working[pivot_index], working[pivot_row] = (
                working[pivot_row],
                working[pivot_index],
            )
            sign *= -1.0
        pivot_value = working[pivot_index][pivot_index]
        if pivot_value < 0:
            sign *= -1.0
        log_abs_det += math.log(abs(pivot_value))
        for row_index in range(pivot_index + 1, size):
            factor = working[row_index][pivot_index] / pivot_value
            if math.isclose(factor, 0.0, abs_tol=1e-15):
                continue
            for column_index in range(pivot_index, size):
                working[row_index][column_index] -= (
                    factor * working[pivot_index][column_index]
                )
    if sign <= 0:
        raise ValueError("matrix determinant is not positive")
    return log_abs_det


def stable_covariance(
    matrix: list[list[float]], *, epsilon: float = 1e-8
) -> list[list[float]]:
    """Return a lightly regularized covariance matrix for numerical stability."""
    stabilized = matrix_copy(matrix)
    for index in range(len(stabilized)):
        stabilized[index][index] += epsilon
    return stabilized
