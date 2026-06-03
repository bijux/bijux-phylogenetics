from __future__ import annotations

import math

from bijux_phylogenetics.ancestral.common import (
    AncestralContinuousDataset,
    stable_value,
)
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    summarize_ultrametric_tip_depths,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .models import ContinuousAncestralBrownianFitDiagnostics

_BROWNIAN_CONDITION_THRESHOLD = 1e12
_SOLVER_REGULARIZATION_EPSILON = 1e-8


def _apply_anc_ml_profile_to_brownian_diagnostics(
    diagnostics: ContinuousAncestralBrownianFitDiagnostics | None,
    *,
    log_likelihood: float,
    sigma_squared: float,
) -> ContinuousAncestralBrownianFitDiagnostics | None:
    if diagnostics is None:
        return None
    return ContinuousAncestralBrownianFitDiagnostics(
        covariance_model=diagnostics.covariance_model,
        tree_is_ultrametric=diagnostics.tree_is_ultrametric,
        minimum_root_to_tip_depth=diagnostics.minimum_root_to_tip_depth,
        maximum_root_to_tip_depth=diagnostics.maximum_root_to_tip_depth,
        minimum_branch_length=diagnostics.minimum_branch_length,
        maximum_branch_length=diagnostics.maximum_branch_length,
        covariance_matrix_dimension=diagnostics.covariance_matrix_dimension,
        covariance_matrix_rank=diagnostics.covariance_matrix_rank,
        covariance_singular=diagnostics.covariance_singular,
        covariance_near_singular=diagnostics.covariance_near_singular,
        covariance_positive_definite=diagnostics.covariance_positive_definite,
        covariance_condition_number=diagnostics.covariance_condition_number,
        covariance_log_determinant=diagnostics.covariance_log_determinant,
        solver_name=diagnostics.solver_name,
        solver_regularized=diagnostics.solver_regularized,
        solver_regularization_epsilon=diagnostics.solver_regularization_epsilon,
        log_likelihood=log_likelihood,
        residual_sigma_squared=sigma_squared,
    )


def _summarize_brownian_fit_diagnostics(
    dataset: AncestralContinuousDataset,
) -> ContinuousAncestralBrownianFitDiagnostics:
    trait_values = [dataset.values_by_taxon[taxon] for taxon in dataset.taxa]
    covariance_matrix = _build_brownian_covariance_matrix(dataset.tree, dataset.taxa)
    covariance_matrix_rank = _matrix_rank(covariance_matrix, tolerance=1e-12)
    covariance_singular = covariance_matrix_rank < len(dataset.taxa)
    covariance_condition_number = math.inf
    if not covariance_singular:
        covariance_condition_number = _symmetric_matrix_condition_number(
            covariance_matrix
        )
    covariance_near_singular = (
        covariance_singular
        or covariance_condition_number >= _BROWNIAN_CONDITION_THRESHOLD
    )
    covariance_positive_definite, covariance_log_determinant = (
        _matrix_positive_definite_diagnostics(covariance_matrix)
    )
    solver_regularized = False
    solver_regularization_epsilon: float | None = None
    solver_covariance = covariance_matrix
    try:
        inverse_covariance = _invert_matrix(solver_covariance)
        solver_log_determinant = _log_determinant(solver_covariance)
    except ValueError:
        solver_regularized = True
        solver_regularization_epsilon = _SOLVER_REGULARIZATION_EPSILON
        solver_covariance = _stable_covariance(
            covariance_matrix,
            epsilon=_SOLVER_REGULARIZATION_EPSILON,
        )
        inverse_covariance = _invert_matrix(solver_covariance)
        solver_log_determinant = _log_determinant(solver_covariance)
    ones = [1.0] * len(trait_values)
    denominator = _quadratic_form(ones, inverse_covariance)
    root_estimate = (
        sum(
            ones[row_index]
            * sum(
                inverse_covariance[row_index][column_index] * trait_values[column_index]
                for column_index in range(len(trait_values))
            )
            for row_index in range(len(trait_values))
        )
        / denominator
    )
    residuals = [value - root_estimate for value in trait_values]
    residual_sigma_squared = max(
        _quadratic_form(residuals, inverse_covariance) / len(trait_values),
        1e-12,
    )
    log_likelihood = -0.5 * (
        len(trait_values) * math.log(2.0 * math.pi * residual_sigma_squared)
        + solver_log_determinant
        + len(trait_values)
    )
    root_depths = _tip_root_depths(dataset.tree, dataset.taxa)
    ultrametric_summary = summarize_ultrametric_tip_depths(root_depths, tolerance=1e-12)
    minimum_branch_length, maximum_branch_length = _branch_length_range(dataset.tree)
    return ContinuousAncestralBrownianFitDiagnostics(
        covariance_model="brownian-shared-path",
        tree_is_ultrametric=ultrametric_summary.ultrametric,
        minimum_root_to_tip_depth=stable_value(ultrametric_summary.minimum_tip_depth),
        maximum_root_to_tip_depth=stable_value(ultrametric_summary.maximum_tip_depth),
        minimum_branch_length=stable_value(minimum_branch_length),
        maximum_branch_length=stable_value(maximum_branch_length),
        covariance_matrix_dimension=len(dataset.taxa),
        covariance_matrix_rank=covariance_matrix_rank,
        covariance_singular=covariance_singular,
        covariance_near_singular=covariance_near_singular,
        covariance_positive_definite=covariance_positive_definite,
        covariance_condition_number=stable_value(covariance_condition_number),
        covariance_log_determinant=(
            None
            if covariance_log_determinant is None
            else stable_value(covariance_log_determinant)
        ),
        solver_name="gauss-jordan-inverse",
        solver_regularized=solver_regularized,
        solver_regularization_epsilon=solver_regularization_epsilon,
        log_likelihood=stable_value(log_likelihood),
        residual_sigma_squared=stable_value(residual_sigma_squared),
    )


def _matrix_positive_definite_diagnostics(
    covariance_matrix: list[list[float]],
) -> tuple[bool, float | None]:
    try:
        _invert_matrix(covariance_matrix)
        return True, _log_determinant(covariance_matrix)
    except ValueError:
        return False, None


def _branch_length_range(tree: PhyloTree) -> tuple[float, float]:
    branch_lengths = [
        float(node.branch_length)
        for node in tree.iter_nodes()
        if node is not tree.root and node.branch_length is not None
    ]
    if not branch_lengths:
        return 0.0, 0.0
    return min(branch_lengths), max(branch_lengths)


def _build_brownian_covariance_matrix(
    tree: PhyloTree, taxa: list[str]
) -> list[list[float]]:
    leaf_paths = _leaf_ancestor_depths(tree)
    matrix: list[list[float]] = []
    for left_taxon in taxa:
        left_path = leaf_paths[left_taxon]
        row: list[float] = []
        for right_taxon in taxa:
            right_path = leaf_paths[right_taxon]
            shared_ancestor_ids = set(left_path) & set(right_path)
            shared_depth = max(left_path[node_id] for node_id in shared_ancestor_ids)
            row.append(shared_depth)
        matrix.append(row)
    return matrix


def _build_anc_ml_covariance_matrix(
    tree: PhyloTree,
    taxa: list[str],
) -> tuple[list[TreeNode], list[list[float]]]:
    ancestor_depths = _node_ancestor_depths(tree)
    leaves = {node.name: node for node in tree.iter_leaves()}
    ordered_tip_nodes = [leaves[taxon] for taxon in taxa]
    ordered_internal_nodes = [
        node
        for node in tree.iter_internal_nodes(order="preorder")
        if node is not tree.root
    ]
    ordered_nodes = ordered_tip_nodes + ordered_internal_nodes
    covariance_matrix: list[list[float]] = []
    for left_node in ordered_nodes:
        left_path = ancestor_depths[left_node.node_id or ""]
        row: list[float] = []
        for right_node in ordered_nodes:
            right_path = ancestor_depths[right_node.node_id or ""]
            shared_ancestor_ids = set(left_path) & set(right_path)
            shared_depth = max(left_path[node_id] for node_id in shared_ancestor_ids)
            row.append(shared_depth)
        covariance_matrix.append(row)
    return ordered_nodes, covariance_matrix


def _tip_root_depths(tree: PhyloTree, taxa: list[str]) -> dict[str, float]:
    leaf_paths = _leaf_ancestor_depths(tree)
    return {taxon: max(leaf_paths[taxon].values()) for taxon in taxa}


def _leaf_ancestor_depths(tree: PhyloTree) -> dict[str, dict[str, float]]:
    depths_by_leaf: dict[str, dict[str, float]] = {}

    def visit(node: TreeNode, depth: float, path: dict[str, float]) -> None:
        branch_length = 0.0 if node is tree.root else float(node.branch_length or 0.0)
        current_depth = depth + branch_length
        current_path = dict(path)
        current_path[node.node_id or ""] = current_depth
        if node.is_leaf():
            if node.name is None:
                raise ValueError(
                    "leaf taxon name is required for ancestral reconstruction"
                )
            depths_by_leaf[node.name] = current_path
            return
        for child in node.children:
            visit(child, current_depth, current_path)

    visit(tree.root, 0.0, {})
    return depths_by_leaf


def _node_ancestor_depths(tree: PhyloTree) -> dict[str, dict[str, float]]:
    depths_by_node: dict[str, dict[str, float]] = {}

    def visit(node: TreeNode, depth: float, path: dict[str, float]) -> None:
        branch_length = 0.0 if node is tree.root else float(node.branch_length or 0.0)
        current_depth = depth + branch_length
        current_path = dict(path)
        current_path[node.node_id or ""] = current_depth
        depths_by_node[node.node_id or ""] = current_path
        for child in node.children:
            visit(child, current_depth, current_path)

    visit(tree.root, 0.0, {})
    return depths_by_node


def _build_anc_ml_parameter_design(
    tip_count: int,
    internal_count: int,
) -> list[list[float]]:
    design_matrix: list[list[float]] = []
    for _ in range(tip_count):
        design_matrix.append([-1.0] + [0.0] * internal_count)
    for internal_index in range(internal_count):
        row = [-1.0] + [0.0] * internal_count
        row[internal_index + 1] = 1.0
        design_matrix.append(row)
    return design_matrix


def _design_quadratic_form(
    design_matrix: list[list[float]],
    precision_matrix: list[list[float]],
) -> list[list[float]]:
    row_count = len(design_matrix)
    parameter_count = len(design_matrix[0]) if design_matrix else 0
    return [
        [
            sum(
                design_matrix[row_index][left_parameter]
                * sum(
                    precision_matrix[row_index][column_index]
                    * design_matrix[column_index][right_parameter]
                    for column_index in range(row_count)
                )
                for row_index in range(row_count)
            )
            for right_parameter in range(parameter_count)
        ]
        for left_parameter in range(parameter_count)
    ]


def _quadratic_form(vector: list[float], matrix: list[list[float]]) -> float:
    return sum(
        vector[row_index]
        * sum(
            matrix[row_index][column_index] * vector[column_index]
            for column_index in range(len(vector))
        )
        for row_index in range(len(vector))
    )


def _invert_matrix(matrix: list[list[float]]) -> list[list[float]]:
    size = len(matrix)
    if size == 0:
        return []
    augmented = [
        [float(value) for value in row]
        + [1.0 if row_index == column_index else 0.0 for column_index in range(size)]
        for row_index, row in enumerate(matrix)
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


def _log_determinant(matrix: list[list[float]]) -> float:
    size = len(matrix)
    if size == 0:
        return 0.0
    working = [list(row) for row in matrix]
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
        if pivot_value < 0.0:
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
    if sign <= 0.0:
        raise ValueError("matrix determinant is not positive")
    return log_abs_det


def _stable_covariance(
    matrix: list[list[float]], *, epsilon: float
) -> list[list[float]]:
    stabilized = [list(row) for row in matrix]
    for index in range(len(stabilized)):
        stabilized[index][index] += epsilon
    return stabilized


def _symmetric_matrix_eigenvalues(
    matrix: list[list[float]],
    *,
    tolerance: float = 1e-15,
    max_iterations: int = 10_000,
) -> list[float]:
    size = len(matrix)
    if size == 0:
        return []
    if size == 1:
        return [float(matrix[0][0])]
    working = [list(row) for row in matrix]
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


def _symmetric_matrix_condition_number(
    matrix: list[list[float]], *, tolerance: float = 1e-12
) -> float:
    singular_values = sorted(
        abs(value)
        for value in _symmetric_matrix_eigenvalues(matrix, tolerance=tolerance)
    )
    if not singular_values:
        return 0.0
    if math.isclose(singular_values[0], 0.0, abs_tol=tolerance):
        return math.inf
    return singular_values[-1] / singular_values[0]


def _matrix_rank(matrix: list[list[float]], *, tolerance: float) -> int:
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
        for column_index in range(pivot_column, column_count):
            working[pivot_row][column_index] /= pivot
        for row_index in range(row_count):
            if row_index == pivot_row:
                continue
            factor = working[row_index][pivot_column]
            if math.isclose(factor, 0.0, abs_tol=tolerance):
                continue
            for column_index in range(pivot_column, column_count):
                working[row_index][column_index] -= (
                    factor * working[pivot_row][column_index]
                )
        rank += 1
        pivot_row += 1
        if pivot_row == row_count:
            break
    return rank
