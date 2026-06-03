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
    diagonal_scale = max(
        (abs(matrix[index][index]) for index in range(len(matrix))),
        default=0.0,
    )
    stabilizer = max(1e-12, min(epsilon, diagonal_scale * 5e-8))
    stabilized = matrix_copy(matrix)
    for index in range(len(stabilized)):
        stabilized[index][index] += stabilizer
    return stabilized


def matrix_infinity_norm(matrix: list[list[float]]) -> float:
    """Return the infinity norm of a rectangular matrix."""
    if not matrix:
        return 0.0
    return max(sum(abs(value) for value in row) for row in matrix)


def matrix_condition_number(matrix: list[list[float]]) -> float:
    """Return the infinity-norm condition number of an invertible square matrix."""
    inverse = invert_matrix(matrix)
    return matrix_infinity_norm(matrix) * matrix_infinity_norm(inverse)


def symmetric_matrix_eigenvalues(
    matrix: list[list[float]],
    *,
    tolerance: float = 1e-15,
    max_iterations: int = 10_000,
) -> list[float]:
    """Return the eigenvalues of a real symmetric matrix via Jacobi rotations."""
    size = len(matrix)
    if size == 0:
        return []
    if size == 1:
        return [float(matrix[0][0])]
    working = matrix_copy(matrix)
    for _ in range(max_iterations):
        pivot_row = 0
        pivot_column = 1
        pivot_value = 0.0
        for row_index in range(size):
            for column_index in range(row_index + 1, size):
                candidate = abs(working[row_index][column_index])
                if candidate > pivot_value:
                    pivot_row = row_index
                    pivot_column = column_index
                    pivot_value = candidate
        if pivot_value <= tolerance:
            return [working[index][index] for index in range(size)]
        app = working[pivot_row][pivot_row]
        aqq = working[pivot_column][pivot_column]
        apq = working[pivot_row][pivot_column]
        tau = (aqq - app) / (2.0 * apq)
        tangent = (
            math.copysign(1.0, tau) / (abs(tau) + math.sqrt(1.0 + tau * tau))
            if not math.isclose(tau, 0.0, abs_tol=tolerance)
            else 1.0
        )
        cosine = 1.0 / math.sqrt(1.0 + tangent * tangent)
        sine = tangent * cosine
        for index in range(size):
            if index in (pivot_row, pivot_column):
                continue
            left = working[index][pivot_row]
            right = working[index][pivot_column]
            working[index][pivot_row] = working[pivot_row][index] = (
                cosine * left - sine * right
            )
            working[index][pivot_column] = working[pivot_column][index] = (
                sine * left + cosine * right
            )
        working[pivot_row][pivot_row] = (
            cosine * cosine * app - 2.0 * sine * cosine * apq + sine * sine * aqq
        )
        working[pivot_column][pivot_column] = (
            sine * sine * app + 2.0 * sine * cosine * apq + cosine * cosine * aqq
        )
        working[pivot_row][pivot_column] = 0.0
        working[pivot_column][pivot_row] = 0.0
    raise ValueError("symmetric eigenvalue iteration did not converge")


def symmetric_matrix_condition_number(
    matrix: list[list[float]], *, tolerance: float = 1e-12
) -> float:
    """Return the exact singular-value condition number of a symmetric matrix."""
    singular_values = sorted(
        abs(value)
        for value in symmetric_matrix_eigenvalues(matrix, tolerance=tolerance)
    )
    if not singular_values:
        return 0.0
    if math.isclose(singular_values[0], 0.0, abs_tol=tolerance):
        return math.inf
    return singular_values[-1] / singular_values[0]


def matrix_rank(matrix: list[list[float]], *, tolerance: float) -> int:
    """Return the numeric rank of a matrix under one absolute pivot tolerance."""
    if not matrix:
        return 0
    working = [list(map(float, row)) for row in matrix]
    row_count = len(working)
    column_count = len(working[0])
    rank = 0
    pivot_row = 0
    for pivot_column in range(column_count):
        candidate_row = max(
            range(pivot_row, row_count),
            key=lambda index: abs(working[index][pivot_column]),
        )
        pivot_value = working[candidate_row][pivot_column]
        if math.isclose(pivot_value, 0.0, abs_tol=tolerance):
            continue
        working[pivot_row], working[candidate_row] = (
            working[candidate_row],
            working[pivot_row],
        )
        pivot = working[pivot_row][pivot_column]
        working[pivot_row] = [value / pivot for value in working[pivot_row]]
        for row_index in range(row_count):
            if row_index == pivot_row:
                continue
            factor = working[row_index][pivot_column]
            if math.isclose(factor, 0.0, abs_tol=tolerance):
                continue
            working[row_index] = [
                row_value - factor * pivot_value
                for row_value, pivot_value in zip(
                    working[row_index], working[pivot_row], strict=True
                )
            ]
        rank += 1
        pivot_row += 1
        if pivot_row == row_count:
            break
    return rank


def _beta_continued_fraction(
    a: float,
    b: float,
    x: float,
    *,
    max_iterations: int = 200,
    tolerance: float = 3e-14,
) -> float:
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - (qab * x / qap)
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    fraction = d
    for iteration in range(1, max_iterations + 1):
        even_step = 2 * iteration
        numerator = (
            iteration * (b - iteration) * x / ((qam + even_step) * (a + even_step))
        )
        d = 1.0 + numerator * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + numerator / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        fraction *= d * c

        numerator = (
            -(a + iteration)
            * (qab + iteration)
            * x
            / ((a + even_step) * (qap + even_step))
        )
        d = 1.0 + numerator * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + numerator / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = d * c
        fraction *= delta
        if abs(delta - 1.0) <= tolerance:
            return fraction
    raise ValueError("beta continued fraction did not converge")


def regularized_incomplete_beta(a: float, b: float, x: float) -> float:
    """Return the regularized incomplete beta function I_x(a, b)."""
    if a <= 0.0 or b <= 0.0:
        raise ValueError("incomplete beta parameters must be positive")
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_beta_term = (
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log1p(-x)
    )
    front = math.exp(log_beta_term)
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _beta_continued_fraction(a, b, x) / a
    return 1.0 - front * _beta_continued_fraction(b, a, 1.0 - x) / b


def student_t_cdf(value: float, degrees_of_freedom: float) -> float:
    """Return the Student-t cumulative distribution function."""
    if degrees_of_freedom <= 0.0:
        raise ValueError("degrees_of_freedom must be positive")
    if math.isclose(value, 0.0, abs_tol=1e-15):
        return 0.5
    x = degrees_of_freedom / (degrees_of_freedom + value * value)
    tail_mass = 0.5 * regularized_incomplete_beta(degrees_of_freedom / 2.0, 0.5, x)
    if value > 0.0:
        return 1.0 - tail_mass
    return tail_mass


def student_t_two_sided_p_value(
    statistic: float,
    degrees_of_freedom: float,
) -> float:
    """Return the two-sided Student-t p-value for a test statistic."""
    tail_probability = 1.0 - student_t_cdf(abs(statistic), degrees_of_freedom)
    return min(max(2.0 * tail_probability, 0.0), 1.0)


def student_t_quantile(probability: float, degrees_of_freedom: float) -> float:
    """Return the Student-t quantile for one cumulative probability."""
    if degrees_of_freedom <= 0.0:
        raise ValueError("degrees_of_freedom must be positive")
    if not 0.0 < probability < 1.0:
        raise ValueError("probability must fall strictly between 0 and 1")
    if math.isclose(probability, 0.5, abs_tol=1e-15):
        return 0.0
    if probability < 0.5:
        return -student_t_quantile(1.0 - probability, degrees_of_freedom)
    lower = 0.0
    upper = 1.0
    while student_t_cdf(upper, degrees_of_freedom) < probability:
        upper *= 2.0
    for _ in range(100):
        midpoint = (lower + upper) / 2.0
        if student_t_cdf(midpoint, degrees_of_freedom) < probability:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2.0
