from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import random

from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    log_determinant,
    quadratic_form,
)
from bijux_phylogenetics.comparative.common import (
    ComparativeDataset,
    lambda_transform_covariance,
    load_comparative_dataset,
    tip_root_depths,
)
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    summarize_ultrametric_tip_depths,
)
from bijux_phylogenetics.phylo.topology.node_identity import ape_node_id_for_node
from bijux_phylogenetics.runtime.errors import ComparativeMethodError


@dataclass(slots=True)
class IndependentContrast:
    """One phylogenetic independent contrast at an internal node."""

    node_id: int
    node: str
    left_taxa: list[str]
    right_taxa: list[str]
    contrast: float
    expected_variance: float
    ancestral_value: float


@dataclass(slots=True)
class IndependentContrastInputAudit:
    """Owned input-policy audit for one phylogenetic independent-contrast analysis."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    taxa: list[str]
    tree_is_ultrametric: bool
    minimum_root_to_tip_depth: float
    maximum_root_to_tip_depth: float
    ultrametric_policy: str
    missing_value_policy: str
    pruned_missing_value_taxa: list[str]
    warnings: list[str]


@dataclass(slots=True)
class IndependentContrastReport:
    """Independent contrasts across a rooted binary tree."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    input_audit: IndependentContrastInputAudit
    root_estimate: float
    contrasts: list[IndependentContrast]


@dataclass(slots=True)
class PhylogeneticSignalInputAudit:
    """Owned input-policy audit for one phylogenetic signal analysis."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    taxa: list[str]
    tree_is_ultrametric: bool
    minimum_root_to_tip_depth: float
    maximum_root_to_tip_depth: float
    ultrametric_policy: str
    missing_value_policy: str
    pruned_missing_value_taxa: list[str]
    warnings: list[str]


@dataclass(slots=True)
class BlombergKReport:
    """Blomberg's K estimate for one numeric trait."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    input_audit: PhylogeneticSignalInputAudit
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
    input_audit: PhylogeneticSignalInputAudit
    lambda_value: float
    log_likelihood: float
    null_log_likelihood: float
    brownian_log_likelihood: float
    likelihood_ratio_statistic: float
    likelihood_ratio_p_value: float
    p_value_method: str
    optimizer_diagnostics: PagelLambdaOptimizerDiagnostics
    profile_rows: list[PagelLambdaProfileRow]


@dataclass(slots=True)
class PagelLambdaLikelihoodReport:
    """Fixed-lambda likelihood evaluation for one numeric trait."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    input_audit: PhylogeneticSignalInputAudit
    lambda_value: float
    log_likelihood: float


@dataclass(slots=True)
class PagelLambdaOptimizerDiagnostics:
    """Optimizer diagnostics for one Pagel's-lambda likelihood search."""

    optimizer_name: str
    coarse_step: float
    fine_step: float
    coarse_grid_point_count: int
    fine_grid_point_count: int
    function_evaluation_count: int
    coarse_best_lambda: float
    coarse_best_log_likelihood: float
    fine_search_start: float
    fine_search_stop: float
    converged: bool
    hit_lower_boundary: bool
    hit_upper_boundary: bool


@dataclass(slots=True)
class PagelLambdaProfileRow:
    """One evaluated lambda profile point from the governed grid search."""

    lambda_value: float
    log_likelihood: float
    delta_log_likelihood: float
    within_95_confidence_interval: bool


@dataclass(slots=True)
class PhylogeneticSignalTestReport:
    """Permutation-based phylogenetic signal test."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    input_audit: PhylogeneticSignalInputAudit
    observed_k: float
    estimated_lambda: float
    p_value: float
    permutations: int
    seed: int
    permuted_k_at_or_above_observed: int
    null_distribution_minimum: float
    null_distribution_mean: float
    null_distribution_maximum: float
    permutation_rows: list[PhylogeneticSignalPermutation]


@dataclass(slots=True)
class PhylogeneticSignalPermutation:
    """One permutation row for a Blomberg-K signal test."""

    permutation_index: int
    permuted_k: float
    at_or_above_observed: bool


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
    return compute_phylogenetic_independent_contrasts_from_dataset(dataset)


def compute_phylogenetic_independent_contrasts_from_dataset(
    dataset: ComparativeDataset,
) -> IndependentContrastReport:
    """Compute Felsenstein independent contrasts from one native comparative dataset."""
    lookup = dict(zip(dataset.taxa, dataset.trait_values, strict=True))
    input_audit = _build_independent_contrast_input_audit(dataset)
    contrasts, root_estimate, _ = _compute_node_contrasts(
        dataset.tree,
        dataset.tree.root,
        lookup,
    )
    return IndependentContrastReport(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        trait=dataset.trait,
        taxon_count=len(dataset.taxa),
        input_audit=input_audit,
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
    dataset = _load_signal_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    input_audit = _build_signal_input_audit(dataset)
    inverse_covariance = invert_matrix(dataset.covariance_matrix)
    generalized_mean = _generalized_mean(dataset.trait_values, inverse_covariance)
    residuals = [value - generalized_mean for value in dataset.trait_values]
    observed_mean_square = sum(value * value for value in residuals) / (
        len(residuals) - 1
    )
    phylogenetic_mean_square = quadratic_form(residuals, inverse_covariance) / (
        len(residuals) - 1
    )
    sum_inverse = sum(sum(row) for row in inverse_covariance)
    expected_mean_square_ratio = (
        sum(
            dataset.covariance_matrix[index][index]
            for index in range(len(dataset.covariance_matrix))
        )
        - (len(dataset.trait_values) / sum_inverse)
    ) / (len(dataset.trait_values) - 1)
    k = (observed_mean_square / phylogenetic_mean_square) / expected_mean_square_ratio
    return BlombergKReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=len(dataset.taxa),
        input_audit=input_audit,
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
    dataset = _load_signal_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    input_audit = _build_signal_input_audit(dataset)
    fit_result = _estimate_pagels_lambda_from_dataset(
        dataset,
        coarse_step=coarse_step,
        fine_step=fine_step,
    )
    null_log_likelihood = _lambda_log_likelihood(dataset, 0.0)
    brownian_log_likelihood = _lambda_log_likelihood(dataset, 1.0)
    likelihood_ratio_statistic = max(
        0.0,
        2.0 * (fit_result.log_likelihood - null_log_likelihood),
    )
    return PagelLambdaReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=len(dataset.taxa),
        input_audit=input_audit,
        lambda_value=fit_result.lambda_value,
        log_likelihood=fit_result.log_likelihood,
        null_log_likelihood=null_log_likelihood,
        brownian_log_likelihood=brownian_log_likelihood,
        likelihood_ratio_statistic=likelihood_ratio_statistic,
        likelihood_ratio_p_value=_likelihood_ratio_p_value(likelihood_ratio_statistic),
        p_value_method="chi-square-approximation",
        optimizer_diagnostics=fit_result.optimizer_diagnostics,
        profile_rows=fit_result.profile_rows,
    )


def evaluate_pagels_lambda_likelihood(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    lambda_value: float,
    taxon_column: str | None = None,
) -> PagelLambdaLikelihoodReport:
    """Evaluate the Gaussian log-likelihood at one fixed Pagel's lambda value."""
    dataset = _load_signal_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    input_audit = _build_signal_input_audit(dataset)
    return evaluate_pagels_lambda_likelihood_from_dataset(
        dataset,
        lambda_value=lambda_value,
        input_audit=input_audit,
    )


def evaluate_pagels_lambda_likelihood_from_dataset(
    dataset: ComparativeDataset,
    *,
    lambda_value: float,
    input_audit: PhylogeneticSignalInputAudit | None = None,
) -> PagelLambdaLikelihoodReport:
    """Evaluate the Gaussian log-likelihood at one fixed Pagel's lambda value."""
    if not 0.0 <= lambda_value <= 1.0:
        raise ValueError(f"lambda_value must be within [0, 1], got {lambda_value}")
    audit = input_audit or _build_signal_input_audit(dataset)
    return PagelLambdaLikelihoodReport(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        trait=dataset.trait,
        taxon_count=len(dataset.taxa),
        input_audit=audit,
        lambda_value=lambda_value,
        log_likelihood=_lambda_log_likelihood(dataset, lambda_value),
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
    dataset = _load_signal_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    input_audit = _build_signal_input_audit(dataset)
    observed_k = _compute_blombergs_k_from_dataset(dataset)
    estimated_lambda = _estimate_pagels_lambda_from_dataset(dataset).lambda_value
    randomizer = random.Random(seed)  # nosec B311
    exceed_count = 0
    permuted_values = list(dataset.trait_values)
    permutation_rows: list[PhylogeneticSignalPermutation] = []
    for permutation_index in range(1, permutations + 1):
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
        permuted_k = _compute_blombergs_k_from_dataset(permuted_dataset)
        at_or_above_observed = permuted_k >= observed_k
        if at_or_above_observed:
            exceed_count += 1
        permutation_rows.append(
            PhylogeneticSignalPermutation(
                permutation_index=permutation_index,
                permuted_k=permuted_k,
                at_or_above_observed=at_or_above_observed,
            )
        )
    permuted_k_values = [row.permuted_k for row in permutation_rows]
    p_value = (exceed_count + 1) / (permutations + 1)
    return PhylogeneticSignalTestReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=len(dataset.taxa),
        input_audit=input_audit,
        observed_k=observed_k,
        estimated_lambda=estimated_lambda,
        p_value=p_value,
        permutations=permutations,
        seed=seed,
        permuted_k_at_or_above_observed=exceed_count,
        null_distribution_minimum=min(permuted_k_values),
        null_distribution_mean=sum(permuted_k_values) / len(permuted_k_values),
        null_distribution_maximum=max(permuted_k_values),
        permutation_rows=permutation_rows,
    )


def _load_signal_dataset(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None,
) -> ComparativeDataset:
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=3,
        require_rooted=True,
    )
    _require_signal_variation(dataset)
    return dataset


def _build_signal_input_audit(
    dataset: ComparativeDataset,
) -> PhylogeneticSignalInputAudit:
    ultrametric_summary = summarize_ultrametric_tip_depths(
        tip_root_depths(dataset.tree, dataset.taxa),
        tolerance=1e-12,
    )
    return PhylogeneticSignalInputAudit(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        trait=dataset.trait,
        taxon_count=len(dataset.taxa),
        taxa=list(dataset.taxa),
        tree_is_ultrametric=ultrametric_summary.ultrametric,
        minimum_root_to_tip_depth=ultrametric_summary.minimum_tip_depth,
        maximum_root_to_tip_depth=ultrametric_summary.maximum_tip_depth,
        ultrametric_policy="accept-rooted-trees-and-report-ultrametricity",
        missing_value_policy="prune-overlapping-missing-values",
        pruned_missing_value_taxa=list(dataset.readiness.pruned_missing_value_taxa),
        warnings=list(dataset.readiness.warnings),
    )


def _build_independent_contrast_input_audit(
    dataset: ComparativeDataset,
) -> IndependentContrastInputAudit:
    ultrametric_summary = summarize_ultrametric_tip_depths(
        tip_root_depths(dataset.tree, dataset.taxa),
        tolerance=1e-12,
    )
    return IndependentContrastInputAudit(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        trait=dataset.trait,
        taxon_count=len(dataset.taxa),
        taxa=list(dataset.taxa),
        tree_is_ultrametric=ultrametric_summary.ultrametric,
        minimum_root_to_tip_depth=ultrametric_summary.minimum_tip_depth,
        maximum_root_to_tip_depth=ultrametric_summary.maximum_tip_depth,
        ultrametric_policy="accept-rooted-trees-and-report-ultrametricity",
        missing_value_policy="prune-overlapping-missing-values",
        pruned_missing_value_taxa=list(dataset.readiness.pruned_missing_value_taxa),
        warnings=list(dataset.readiness.warnings),
    )


def _require_signal_variation(dataset: ComparativeDataset) -> None:
    minimum_value = min(dataset.trait_values)
    maximum_value = max(dataset.trait_values)
    if math.isclose(minimum_value, maximum_value, abs_tol=1e-12):
        raise ComparativeMethodError(
            "phylogenetic signal requires at least two distinct numeric trait values after pruning"
        )


def _compute_node_contrasts(
    tree,
    node,
    values_by_taxon: dict[str, float],
) -> tuple[list[IndependentContrast], float, float]:
    if node.is_leaf():
        if node.name is None:
            raise ValueError("leaf taxon name is required for comparative contrasts")
        if node.branch_length is None:
            raise ValueError("branch lengths are required for comparative contrasts")
        return [], values_by_taxon[node.name], node.branch_length
    if len(node.children) != 2:
        raise ValueError("independent contrasts require a strictly binary tree")

    left_contrasts, left_value, left_variance = _compute_node_contrasts(
        tree,
        node.children[0],
        values_by_taxon,
    )
    right_contrasts, right_value, right_variance = _compute_node_contrasts(
        tree,
        node.children[1],
        values_by_taxon,
    )
    expected_variance = left_variance + right_variance
    contrast = (left_value - right_value) / math.sqrt(expected_variance)
    ancestral_value = (
        (left_value / left_variance) + (right_value / right_variance)
    ) / ((1.0 / left_variance) + (1.0 / right_variance))
    propagated_variance = (left_variance * right_variance) / (
        left_variance + right_variance
    )
    if node.branch_length is not None:
        propagated_variance += node.branch_length

    report = IndependentContrast(
        node_id=ape_node_id_for_node(tree, node),
        node="|".join(
            sorted(
                set(_leaf_names(node.children[0])) | set(_leaf_names(node.children[1]))
            )
        ),
        left_taxa=_leaf_names(node.children[0]),
        right_taxa=_leaf_names(node.children[1]),
        contrast=contrast,
        expected_variance=expected_variance,
        ancestral_value=ancestral_value,
    )
    return (
        left_contrasts + right_contrasts + [report],
        ancestral_value,
        propagated_variance,
    )


def _leaf_names(node) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_leaf_names(child))
    return sorted(taxa)


def _generalized_mean(
    values: list[float], inverse_covariance: list[list[float]]
) -> float:
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
    observed_mean_square = sum(value * value for value in residuals) / (
        len(residuals) - 1
    )
    phylogenetic_mean_square = quadratic_form(residuals, inverse_covariance) / (
        len(residuals) - 1
    )
    sum_inverse = sum(sum(row) for row in inverse_covariance)
    expected_mean_square_ratio = (
        sum(
            dataset.covariance_matrix[index][index]
            for index in range(len(dataset.covariance_matrix))
        )
        - (len(dataset.trait_values) / sum_inverse)
    ) / (len(dataset.trait_values) - 1)
    return (
        observed_mean_square / phylogenetic_mean_square
    ) / expected_mean_square_ratio


def _estimate_pagels_lambda_from_dataset(
    dataset: ComparativeDataset,
    *,
    coarse_step: float = 0.05,
    fine_step: float = 0.005,
) -> _PagelLambdaFitResult:
    coarse_values = _grid_values(0.0, 1.0, coarse_step)
    coarse_rows = [
        (_lambda, _lambda_log_likelihood(dataset, _lambda)) for _lambda in coarse_values
    ]
    coarse_best_lambda, coarse_best_log_likelihood = max(
        coarse_rows,
        key=lambda item: item[1],
    )
    fine_start = max(0.0, coarse_best_lambda - coarse_step)
    fine_stop = min(1.0, coarse_best_lambda + coarse_step)
    fine_values = _grid_values(fine_start, fine_stop, fine_step)
    fine_rows = [
        (_lambda, _lambda_log_likelihood(dataset, _lambda)) for _lambda in fine_values
    ]
    lambda_value, log_likelihood = max(
        fine_rows,
        key=lambda item: item[1],
    )
    profile_threshold = 0.5 * 3.841458820694124
    profile_rows = [
        PagelLambdaProfileRow(
            lambda_value=row_lambda,
            log_likelihood=row_log_likelihood,
            delta_log_likelihood=max(0.0, log_likelihood - row_log_likelihood),
            within_95_confidence_interval=(
                (log_likelihood - row_log_likelihood) <= profile_threshold
            ),
        )
        for row_lambda, row_log_likelihood in fine_rows
    ]
    diagnostics = PagelLambdaOptimizerDiagnostics(
        optimizer_name="two-stage-grid-search",
        coarse_step=coarse_step,
        fine_step=fine_step,
        coarse_grid_point_count=len(coarse_rows),
        fine_grid_point_count=len(fine_rows),
        function_evaluation_count=len(coarse_rows) + len(fine_rows),
        coarse_best_lambda=coarse_best_lambda,
        coarse_best_log_likelihood=coarse_best_log_likelihood,
        fine_search_start=fine_start,
        fine_search_stop=fine_stop,
        converged=True,
        hit_lower_boundary=math.isclose(lambda_value, 0.0, abs_tol=fine_step / 2.0),
        hit_upper_boundary=math.isclose(lambda_value, 1.0, abs_tol=fine_step / 2.0),
    )
    return _PagelLambdaFitResult(
        lambda_value=lambda_value,
        log_likelihood=log_likelihood,
        optimizer_diagnostics=diagnostics,
        profile_rows=profile_rows,
    )


def _lambda_log_likelihood(dataset: ComparativeDataset, lambda_value: float) -> float:
    covariance = lambda_transform_covariance(dataset.covariance_matrix, lambda_value)
    inverse_covariance = invert_matrix(covariance)
    generalized_mean = _generalized_mean(dataset.trait_values, inverse_covariance)
    residuals = [value - generalized_mean for value in dataset.trait_values]
    sigma_squared = quadratic_form(residuals, inverse_covariance) / len(
        dataset.trait_values
    )
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


@dataclass(slots=True)
class _PagelLambdaFitResult:
    lambda_value: float
    log_likelihood: float
    optimizer_diagnostics: PagelLambdaOptimizerDiagnostics
    profile_rows: list[PagelLambdaProfileRow]


def _likelihood_ratio_p_value(likelihood_ratio_statistic: float) -> float:
    return math.erfc(math.sqrt(likelihood_ratio_statistic / 2.0))
