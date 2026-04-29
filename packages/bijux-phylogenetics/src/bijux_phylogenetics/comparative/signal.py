from __future__ import annotations

from dataclasses import dataclass
import math
import random
from pathlib import Path

from bijux_phylogenetics.comparative._math import invert_matrix, log_determinant, quadratic_form
from bijux_phylogenetics.comparative.common import (
    ComparativeDataset,
    lambda_transform_covariance,
    load_comparative_dataset,
)


@dataclass(slots=True)
class IndependentContrast:
    """One phylogenetic independent contrast at an internal node."""

    node: str
    left_taxa: list[str]
    right_taxa: list[str]
    contrast: float
    expected_variance: float
    ancestral_value: float


@dataclass(slots=True)
class IndependentContrastReport:
    """Independent contrasts across a rooted binary tree."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    root_estimate: float
    contrasts: list[IndependentContrast]


@dataclass(slots=True)
class BlombergKReport:
    """Blomberg's K estimate for one numeric trait."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    k: float
    generalized_mean: float
    observed_mean_square: float
    phylogenetic_mean_square: float
    expected_mean_square_ratio: float


@dataclass(slots=True)
class PagelLambdaReport:
    """Pagel's lambda estimate for one numeric trait."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    lambda_value: float
    log_likelihood: float
    null_log_likelihood: float
    brownian_log_likelihood: float


@dataclass(slots=True)
class PhylogeneticSignalTestReport:
    """Permutation-based phylogenetic signal test."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    observed_k: float
    estimated_lambda: float
    p_value: float
    permutations: int
    permuted_k_at_or_above_observed: int


def compute_phylogenetic_independent_contrasts(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> IndependentContrastReport:
    """Compute Felsenstein independent contrasts for a rooted binary tree."""
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=2,
        require_rooted=True,
        require_binary=True,
    )
    lookup = {taxon: value for taxon, value in zip(dataset.taxa, dataset.trait_values, strict=True)}
    contrasts, root_estimate, _ = _compute_node_contrasts(dataset.tree.root, lookup)
    return IndependentContrastReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=len(dataset.taxa),
        root_estimate=root_estimate,
        contrasts=contrasts,
    )


def compute_blombergs_k(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> BlombergKReport:
    """Compute Blomberg's K under a Brownian covariance model."""
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=3,
        require_rooted=True,
    )
    inverse_covariance = invert_matrix(dataset.covariance_matrix)
    generalized_mean = _generalized_mean(dataset.trait_values, inverse_covariance)
    residuals = [value - generalized_mean for value in dataset.trait_values]
    observed_mean_square = sum(value * value for value in residuals) / (len(residuals) - 1)
    phylogenetic_mean_square = quadratic_form(residuals, inverse_covariance) / (len(residuals) - 1)
    sum_inverse = sum(sum(row) for row in inverse_covariance)
    expected_mean_square_ratio = (
        (
            sum(dataset.covariance_matrix[index][index] for index in range(len(dataset.covariance_matrix)))
            - (len(dataset.trait_values) / sum_inverse)
        )
        / (len(dataset.trait_values) - 1)
    )
    k = (observed_mean_square / phylogenetic_mean_square) / expected_mean_square_ratio
    return BlombergKReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=len(dataset.taxa),
        k=k,
        generalized_mean=generalized_mean,
        observed_mean_square=observed_mean_square,
        phylogenetic_mean_square=phylogenetic_mean_square,
        expected_mean_square_ratio=expected_mean_square_ratio,
    )


def estimate_pagels_lambda(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    coarse_step: float = 0.05,
    fine_step: float = 0.005,
) -> PagelLambdaReport:
    """Estimate Pagel's lambda by likelihood search over [0, 1]."""
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=3,
        require_rooted=True,
    )
    coarse_values = _grid_values(0.0, 1.0, coarse_step)
    coarse_best_lambda, coarse_best_log_likelihood = max(
        ((_lambda, _lambda_log_likelihood(dataset, _lambda)) for _lambda in coarse_values),
        key=lambda item: item[1],
    )
    fine_start = max(0.0, coarse_best_lambda - coarse_step)
    fine_stop = min(1.0, coarse_best_lambda + coarse_step)
    fine_values = _grid_values(fine_start, fine_stop, fine_step)
    lambda_value, log_likelihood = max(
        ((_lambda, _lambda_log_likelihood(dataset, _lambda)) for _lambda in fine_values),
        key=lambda item: item[1],
    )
    return PagelLambdaReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=len(dataset.taxa),
        lambda_value=lambda_value,
        log_likelihood=log_likelihood,
        null_log_likelihood=_lambda_log_likelihood(dataset, 0.0),
        brownian_log_likelihood=_lambda_log_likelihood(dataset, 1.0),
    )


def compute_phylogenetic_signal_test(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    permutations: int = 199,
    seed: int = 1,
) -> PhylogeneticSignalTestReport:
    """Test phylogenetic signal with a permutation distribution of Blomberg's K."""
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=3,
        require_rooted=True,
    )
    observed_k = _compute_blombergs_k_from_dataset(dataset)
    estimated_lambda = estimate_pagels_lambda(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    ).lambda_value
    randomizer = random.Random(seed)
    exceed_count = 0
    permuted_values = list(dataset.trait_values)
    for _ in range(permutations):
        randomizer.shuffle(permuted_values)
        permuted_dataset = ComparativeDataset(
            tree_path=dataset.tree_path,
            traits_path=dataset.traits_path,
            tree=dataset.tree,
            taxon_column=dataset.taxon_column,
            trait=dataset.trait,
            taxa=dataset.taxa,
            trait_values=list(permuted_values),
            covariance_matrix=dataset.covariance_matrix,
            readiness=dataset.readiness,
        )
        if _compute_blombergs_k_from_dataset(permuted_dataset) >= observed_k:
            exceed_count += 1
    p_value = (exceed_count + 1) / (permutations + 1)
    return PhylogeneticSignalTestReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=len(dataset.taxa),
        observed_k=observed_k,
        estimated_lambda=estimated_lambda,
        p_value=p_value,
        permutations=permutations,
        permuted_k_at_or_above_observed=exceed_count,
    )


def _compute_node_contrasts(node, values_by_taxon: dict[str, float]) -> tuple[list[IndependentContrast], float, float]:
    if node.is_leaf():
        if node.name is None:
            raise ValueError("leaf taxon name is required for comparative contrasts")
        if node.branch_length is None:
            raise ValueError("branch lengths are required for comparative contrasts")
        return [], values_by_taxon[node.name], node.branch_length
    if len(node.children) != 2:
        raise ValueError("independent contrasts require a strictly binary tree")

    left_contrasts, left_value, left_variance = _compute_node_contrasts(node.children[0], values_by_taxon)
    right_contrasts, right_value, right_variance = _compute_node_contrasts(node.children[1], values_by_taxon)
    expected_variance = left_variance + right_variance
    contrast = (left_value - right_value) / math.sqrt(expected_variance)
    ancestral_value = (
        (left_value / left_variance) + (right_value / right_variance)
    ) / ((1.0 / left_variance) + (1.0 / right_variance))
    propagated_variance = (left_variance * right_variance) / (left_variance + right_variance)
    if node.branch_length is not None:
        propagated_variance += node.branch_length

    report = IndependentContrast(
        node="|".join(sorted(set(_leaf_names(node.children[0])) | set(_leaf_names(node.children[1])))),
        left_taxa=_leaf_names(node.children[0]),
        right_taxa=_leaf_names(node.children[1]),
        contrast=contrast,
        expected_variance=expected_variance,
        ancestral_value=ancestral_value,
    )
    return left_contrasts + right_contrasts + [report], ancestral_value, propagated_variance


def _leaf_names(node) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_leaf_names(child))
    return sorted(taxa)


def _generalized_mean(values: list[float], inverse_covariance: list[list[float]]) -> float:
    ones = [1.0] * len(values)
    numerator = sum(
        left * value
        for left, value in zip(
            [sum(row) for row in inverse_covariance],
            values,
            strict=True,
        )
    )
    denominator = sum(sum(row) for row in inverse_covariance)
    return numerator / denominator if denominator else sum(values) / len(values)


def _compute_blombergs_k_from_dataset(dataset: ComparativeDataset) -> float:
    inverse_covariance = invert_matrix(dataset.covariance_matrix)
    generalized_mean = _generalized_mean(dataset.trait_values, inverse_covariance)
    residuals = [value - generalized_mean for value in dataset.trait_values]
    observed_mean_square = sum(value * value for value in residuals) / (len(residuals) - 1)
    phylogenetic_mean_square = quadratic_form(residuals, inverse_covariance) / (len(residuals) - 1)
    sum_inverse = sum(sum(row) for row in inverse_covariance)
    expected_mean_square_ratio = (
        (
            sum(dataset.covariance_matrix[index][index] for index in range(len(dataset.covariance_matrix)))
            - (len(dataset.trait_values) / sum_inverse)
        )
        / (len(dataset.trait_values) - 1)
    )
    return (observed_mean_square / phylogenetic_mean_square) / expected_mean_square_ratio


def _lambda_log_likelihood(dataset: ComparativeDataset, lambda_value: float) -> float:
    covariance = lambda_transform_covariance(dataset.covariance_matrix, lambda_value)
    inverse_covariance = invert_matrix(covariance)
    generalized_mean = _generalized_mean(dataset.trait_values, inverse_covariance)
    residuals = [value - generalized_mean for value in dataset.trait_values]
    sigma_squared = quadratic_form(residuals, inverse_covariance) / len(dataset.trait_values)
    return -0.5 * (
        len(dataset.trait_values) * math.log(2.0 * math.pi * sigma_squared)
        + log_determinant(covariance)
        + len(dataset.trait_values)
    )


def _grid_values(start: float, stop: float, step: float) -> list[float]:
    values: list[float] = []
    current = start
    while current <= stop + (step / 2):
        values.append(round(current, 6))
        current += step
    return values
