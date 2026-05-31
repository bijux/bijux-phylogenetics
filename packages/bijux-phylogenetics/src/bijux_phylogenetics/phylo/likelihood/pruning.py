from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import math

import numpy

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .ctmc import ValidatedCtmcRateMatrix, validate_ctmc_rate_matrix


@dataclass(slots=True)
class FiniteStatePruningPass:
    """Scaled postorder conditional likelihood vectors for one rooted tree."""

    conditional_by_node_id: dict[str, numpy.ndarray]
    subtree_log_scaling_by_node_id: dict[str, float]
    state_count: int
    normalized_internal_conditionals: bool

    def conditional_for_node(self, node: TreeNode) -> numpy.ndarray:
        if node.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        return self.conditional_by_node_id[node.node_id]

    def subtree_log_scaling_for_node(self, node: TreeNode) -> float:
        if node.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        return self.subtree_log_scaling_by_node_id[node.node_id]


@dataclass(slots=True)
class FiniteStateTransitionMatrixEvaluator:
    """Branch-length transition evaluator for one validated finite-state model."""

    validated_rate_matrix: ValidatedCtmcRateMatrix
    eigenvalues: numpy.ndarray
    eigenvectors: numpy.ndarray
    inverse_vectors: numpy.ndarray
    cache_matrices: bool
    transition_matrix_by_branch_length: dict[float, numpy.ndarray]
    matrix_exponential_evaluation_count: int = 0

    @property
    def cached_branch_length_count(self) -> int:
        return len(self.transition_matrix_by_branch_length)

    @property
    def eigendecomposition(self) -> "FiniteStateRateMatrixEigendecomposition":
        return FiniteStateRateMatrixEigendecomposition(
            eigenvalues=self.eigenvalues,
            eigenvectors=self.eigenvectors,
            inverse_vectors=self.inverse_vectors,
        )

    def transition_probability_matrix(self, branch_length: float) -> numpy.ndarray:
        normalized_branch_length = _normalize_branch_length(branch_length)
        if not self.cache_matrices:
            self.matrix_exponential_evaluation_count += 1
            return _compute_transition_probability_matrix(
                self.validated_rate_matrix,
                normalized_branch_length,
                eigendecomposition=self.eigendecomposition,
            )
        cached_matrix = self.transition_matrix_by_branch_length.get(
            normalized_branch_length
        )
        if cached_matrix is None:
            cached_matrix = _compute_transition_probability_matrix(
                self.validated_rate_matrix,
                normalized_branch_length,
                eigendecomposition=self.eigendecomposition,
            )
            cached_matrix.setflags(write=False)
            self.transition_matrix_by_branch_length[normalized_branch_length] = (
                cached_matrix
            )
            self.matrix_exponential_evaluation_count += 1
        return cached_matrix.copy()


def build_transition_matrix_evaluator(
    rate_matrix: numpy.ndarray | ValidatedCtmcRateMatrix,
    *,
    cache_matrices: bool = True,
) -> FiniteStateTransitionMatrixEvaluator:
    """Build one reusable transition evaluator for a single finite-state model."""
    validated_rate_matrix = _validated_rate_matrix(rate_matrix)
    eigendecomposition = _compute_rate_matrix_eigendecomposition(validated_rate_matrix)
    return FiniteStateTransitionMatrixEvaluator(
        validated_rate_matrix=validated_rate_matrix,
        eigenvalues=eigendecomposition.eigenvalues,
        eigenvectors=eigendecomposition.eigenvectors,
        inverse_vectors=eigendecomposition.inverse_vectors,
        cache_matrices=cache_matrices,
        transition_matrix_by_branch_length={},
    )


def transition_probability_matrix(
    rate_matrix: numpy.ndarray,
    branch_length: float,
) -> numpy.ndarray:
    """Exponentiate one finite-state rate matrix into a branch transition matrix."""
    validated_rate_matrix = _validated_rate_matrix(rate_matrix)
    return _compute_transition_probability_matrix(
        validated_rate_matrix,
        _normalize_branch_length(branch_length),
    )


def _validated_rate_matrix(
    rate_matrix: numpy.ndarray | ValidatedCtmcRateMatrix,
) -> ValidatedCtmcRateMatrix:
    if isinstance(rate_matrix, ValidatedCtmcRateMatrix):
        return rate_matrix
    return validate_ctmc_rate_matrix(rate_matrix)


def _normalize_branch_length(branch_length: float) -> float:
    if branch_length <= 0.0:
        return 0.0
    return float(branch_length)


@dataclass(slots=True)
class FiniteStateRateMatrixEigendecomposition:
    """Cached eigendecomposition for one validated CTMC rate matrix."""

    eigenvalues: numpy.ndarray
    eigenvectors: numpy.ndarray
    inverse_vectors: numpy.ndarray


def _compute_rate_matrix_eigendecomposition(
    validated_rate_matrix: ValidatedCtmcRateMatrix,
) -> FiniteStateRateMatrixEigendecomposition:
    eigenvalues, eigenvectors = numpy.linalg.eig(validated_rate_matrix.rate_matrix)
    inverse_vectors = numpy.linalg.inv(eigenvectors)
    return FiniteStateRateMatrixEigendecomposition(
        eigenvalues=eigenvalues,
        eigenvectors=eigenvectors,
        inverse_vectors=inverse_vectors,
    )


def _compute_transition_probability_matrix(
    validated_rate_matrix: ValidatedCtmcRateMatrix,
    branch_length: float,
    *,
    eigendecomposition: FiniteStateRateMatrixEigendecomposition | None = None,
) -> numpy.ndarray:
    if branch_length == 0.0:
        return numpy.eye(validated_rate_matrix.state_count, dtype=float)
    if eigendecomposition is None:
        eigendecomposition = _compute_rate_matrix_eigendecomposition(
            validated_rate_matrix
        )
    eigenvalues = eigendecomposition.eigenvalues
    eigenvectors = eigendecomposition.eigenvectors
    inverse_vectors = eigendecomposition.inverse_vectors
    diagonal = numpy.diag(numpy.exp(eigenvalues * branch_length))
    transition = eigenvectors @ diagonal @ inverse_vectors
    transition = numpy.real_if_close(transition, tol=1000).astype(float)
    transition[transition < 0.0] = 0.0
    row_sums = transition.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0.0] = 1.0
    return transition / row_sums


def postorder_conditional_likelihoods(
    tree: PhyloTree,
    *,
    state_count: int,
    leaf_likelihood: Callable[[TreeNode], numpy.ndarray],
    transition_matrix_for_child: Callable[[TreeNode], numpy.ndarray],
    normalize_internal_conditionals: bool = True,
) -> FiniteStatePruningPass:
    """Run Felsenstein's postorder pruning recursion for one finite-state model."""
    if state_count < 1:
        raise ValueError("state_count must be positive")

    conditional_by_node_id: dict[str, numpy.ndarray] = {}
    subtree_log_scaling_by_node_id: dict[str, float] = {}

    for node in tree.iter_nodes(order="postorder"):
        if node.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        if node.is_leaf():
            conditional = _validated_probability_vector(
                leaf_likelihood(node),
                state_count=state_count,
                context=f"leaf '{node.name or node.node_id}'",
            )
            conditional_by_node_id[node.node_id] = conditional
            subtree_log_scaling_by_node_id[node.node_id] = 0.0
            continue

        conditional = numpy.ones(state_count, dtype=float)
        child_log_scaling = 0.0
        child_impossible = False
        for child in node.children:
            if child.node_id is None:
                raise ValueError("tree child node is missing a stable node_id")
            child_transition = _validated_transition_matrix(
                transition_matrix_for_child(child),
                state_count=state_count,
                context=f"branch '{child.node_id}'",
            )
            child_conditional = conditional_by_node_id[child.node_id]
            conditional *= child_transition @ child_conditional
            child_scale = subtree_log_scaling_by_node_id[child.node_id]
            if math.isinf(child_scale) and child_scale < 0.0:
                child_impossible = True
            else:
                child_log_scaling += child_scale

        scale = float(conditional.sum())
        if child_impossible or scale <= 0.0 or not math.isfinite(scale):
            conditional_by_node_id[node.node_id] = conditional
            subtree_log_scaling_by_node_id[node.node_id] = float("-inf")
            continue
        if normalize_internal_conditionals:
            conditional /= scale
            child_log_scaling += math.log(scale)
        conditional_by_node_id[node.node_id] = conditional
        subtree_log_scaling_by_node_id[node.node_id] = child_log_scaling

    return FiniteStatePruningPass(
        conditional_by_node_id=conditional_by_node_id,
        subtree_log_scaling_by_node_id=subtree_log_scaling_by_node_id,
        state_count=state_count,
        normalized_internal_conditionals=normalize_internal_conditionals,
    )


def log_likelihood_from_root_prior(
    tree: PhyloTree,
    pruning_pass: FiniteStatePruningPass,
    *,
    root_prior: numpy.ndarray,
) -> float:
    """Combine one scaled pruning pass with one explicit root prior."""
    weighted_root = _validated_probability_vector(
        root_prior,
        state_count=pruning_pass.state_count,
        context="root prior",
    ) * pruning_pass.conditional_for_node(tree.root)
    root_scale = float(weighted_root.sum())
    if root_scale <= 0.0 or not math.isfinite(root_scale):
        return float("-inf")
    subtree_log_scaling = pruning_pass.subtree_log_scaling_for_node(tree.root)
    if math.isinf(subtree_log_scaling) and subtree_log_scaling < 0.0:
        return float("-inf")
    return subtree_log_scaling + math.log(root_scale)


def _validated_probability_vector(
    vector: numpy.ndarray,
    *,
    state_count: int,
    context: str,
) -> numpy.ndarray:
    candidate = numpy.asarray(vector, dtype=float)
    if candidate.shape != (state_count,):
        raise ValueError(
            f"{context} must be a one-dimensional vector with {state_count} states"
        )
    if not numpy.all(numpy.isfinite(candidate)):
        raise ValueError(f"{context} must contain only finite values")
    if numpy.any(candidate < 0.0):
        raise ValueError(f"{context} must not contain negative probabilities")
    return candidate.copy()


def _validated_transition_matrix(
    matrix: numpy.ndarray,
    *,
    state_count: int,
    context: str,
) -> numpy.ndarray:
    candidate = numpy.asarray(matrix, dtype=float)
    if candidate.shape != (state_count, state_count):
        raise ValueError(
            f"{context} must be a {state_count}x{state_count} transition matrix"
        )
    if not numpy.all(numpy.isfinite(candidate)):
        raise ValueError(f"{context} must contain only finite values")
    if numpy.any(candidate < 0.0):
        raise ValueError(f"{context} must not contain negative probabilities")
    return candidate.copy()
