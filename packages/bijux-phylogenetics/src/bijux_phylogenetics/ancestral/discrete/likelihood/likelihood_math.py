from __future__ import annotations

import math

import numpy


def tree_log_likelihood(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray | None,
    root_prior_mode: str = "given",
) -> float:
    state_index = {state: index for index, state in enumerate(state_order)}
    transition_cache: dict[float, numpy.ndarray] = {}

    def transition(cached_branch_length: float) -> numpy.ndarray:
        cached = transition_cache.get(cached_branch_length)
        if cached is None:
            cached = transition_probability_matrix(rate_matrix, cached_branch_length)
            transition_cache[cached_branch_length] = cached
        return cached

    def visit(node) -> tuple[numpy.ndarray, float]:
        if node.is_leaf():
            likelihood = numpy.zeros(len(state_order), dtype=float)
            likelihood[state_index[states_by_taxon[node.name]]] = 1.0
            return likelihood, 0.0
        partial = numpy.ones(len(state_order), dtype=float)
        log_scale = 0.0
        for child in node.children:
            child_partial, child_scale = visit(child)
            branch_transition = transition(branch_length(child))
            partial *= branch_transition @ child_partial
            log_scale += child_scale
        scale = float(partial.sum())
        if scale <= 0.0:
            return partial, float("-inf")
        partial /= scale
        return partial, log_scale + math.log(scale)

    root_partial, subtree_log_scale = visit(tree.root)
    if root_prior_mode == "observed":
        observed_root_total = float(root_partial.sum())
        if observed_root_total <= 0.0:
            return float("-inf")
        root_scale = float((root_partial @ root_partial) / observed_root_total)
        if root_scale <= 0.0:
            return float("-inf")
        return subtree_log_scale + math.log(root_scale)
    if root_prior is None:
        raise ValueError("root_prior is required unless root_prior_mode is 'observed'")
    root_weight = root_prior * root_partial
    root_scale = float(root_weight.sum())
    if root_scale <= 0.0:
        return float("-inf")
    return subtree_log_scale + math.log(root_scale)


def transition_probability_matrix(
    rate_matrix: numpy.ndarray,
    branch_length: float,
) -> numpy.ndarray:
    if branch_length <= 0.0:
        return numpy.eye(rate_matrix.shape[0], dtype=float)
    eigenvalues, eigenvectors = numpy.linalg.eig(rate_matrix)
    inverse_vectors = numpy.linalg.inv(eigenvectors)
    diagonal = numpy.diag(numpy.exp(eigenvalues * branch_length))
    transition = eigenvectors @ diagonal @ inverse_vectors
    transition = numpy.real_if_close(transition, tol=1000).astype(float)
    transition[transition < 0.0] = 0.0
    row_sums = transition.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0.0] = 1.0
    return transition / row_sums


def branch_length(node) -> float:
    if node.branch_length is None:
        return 1.0
    return max(float(node.branch_length), 0.0)
