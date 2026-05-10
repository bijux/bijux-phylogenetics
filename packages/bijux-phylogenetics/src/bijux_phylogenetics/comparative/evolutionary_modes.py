from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.comparative.common import (
    ComparativeDataset,
    build_brownian_covariance_matrix,
    load_comparative_dataset,
    node_signature,
    stable_covariance,
)
from bijux_phylogenetics.comparative.models import (
    ComparativeModelComparisonRow,
    ComparativeResidualSummary,
    _brownian_parameter_intervals,
    _build_residual_diagnostics,
    _comparison_row,
    _fit_intercept_only_model,
)
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.errors import ComparativeMethodError
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree


ALLOWED_EVOLUTIONARY_MODES = {
    "brownian",
    "ornstein-uhlenbeck",
    "early-burst",
}


@dataclass(slots=True)
class EvolutionaryModeBranchLengthRow:
    """One deterministic branch-length change under a governed tree rescaling."""

    node: str
    descendant_taxa: list[str]
    original_branch_length: float
    transformed_branch_length: float
    parent_depth: float
    child_depth: float


@dataclass(slots=True)
class ComparativeTreeRescalingReport:
    """Canonical summary of one OU or early-burst tree-rescaling surface."""

    tree_path: Path
    mode: str
    parameter_name: str
    parameter_value: float
    tip_count: int
    original_total_branch_length: float
    transformed_total_branch_length: float
    transformed_tree_newick: str
    branch_rows: list[EvolutionaryModeBranchLengthRow]


@dataclass(slots=True)
class ContinuousEvolutionaryModeFitReport:
    """Intercept-only continuous-trait fit under one governed evolutionary mode."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    taxon_count: int
    taxa: list[str]
    mode: str
    parameter_name: str | None
    parameter_value: float | None
    root_state: float
    rate: float
    log_likelihood: float
    aic: float
    aicc: float
    fitted_values: list[float]
    residuals: list[float]
    transformed_tree_newick: str
    confidence_intervals: list[object]
    residual_diagnostics: ComparativeResidualSummary
    assumptions: list[str]


@dataclass(slots=True)
class LikelihoodRatioTestResult:
    """Likelihood-ratio comparison between two governed evolutionary-mode fits."""

    comparison_id: str
    left_mode: str
    right_mode: str
    statistic: float
    degrees_of_freedom: int
    p_value: float


@dataclass(slots=True)
class ContinuousEvolutionaryModeComparisonReport:
    """Model-comparison summary over Brownian, OU, and early-burst fits."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    rows: list[ComparativeModelComparisonRow]
    better_model: str
    likelihood_ratio_tests: list[LikelihoodRatioTestResult]


def rescale_tree_ornstein_uhlenbeck(
    tree_path: Path,
    *,
    alpha: float,
    sigsq: float = 1.0,
) -> ComparativeTreeRescalingReport:
    """Apply the geiger-style OU branch rescaling to a rooted tree."""
    return _build_tree_rescaling_report(
        load_tree(tree_path),
        tree_path,
        mode="ornstein-uhlenbeck",
        parameter_name="alpha",
        parameter_value=alpha,
        sigsq=sigsq,
    )


def rescale_tree_early_burst(
    tree_path: Path,
    *,
    rate_change: float,
    sigsq: float = 1.0,
) -> ComparativeTreeRescalingReport:
    """Apply the geiger-style early-burst branch rescaling to a rooted tree."""
    return _build_tree_rescaling_report(
        load_tree(tree_path),
        tree_path,
        mode="early-burst",
        parameter_name="rate_change",
        parameter_value=rate_change,
        sigsq=sigsq,
    )


def transform_tree_for_evolutionary_mode(
    tree: PhyloTree,
    *,
    mode: str,
    parameter_value: float,
    sigsq: float = 1.0,
) -> PhyloTree:
    """Transform an in-memory tree under the governed OU or early-burst branch rule."""
    return _transform_tree(
        tree,
        mode=mode,
        parameter_value=parameter_value,
        sigsq=sigsq,
    )


def fit_continuous_evolutionary_mode(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    mode: str,
    taxon_column: str | None = None,
    ou_bounds: tuple[float, float] = (0.0, 10.0),
    early_burst_bounds: tuple[float, float] = (0.0, 50.0),
) -> ContinuousEvolutionaryModeFitReport:
    """Fit a Brownian, OU-rescaled, or early-burst-rescaled intercept-only trait model."""
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=3,
        require_rooted=True,
        require_binary=False,
    )
    return _fit_evolutionary_mode_from_dataset(
        dataset,
        mode=mode,
        ou_bounds=ou_bounds,
        early_burst_bounds=early_burst_bounds,
    )


def compare_continuous_evolutionary_modes(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    ou_bounds: tuple[float, float] = (0.0, 10.0),
    early_burst_bounds: tuple[float, float] = (0.0, 50.0),
) -> ContinuousEvolutionaryModeComparisonReport:
    """Compare Brownian, OU-rescaled, and early-burst-rescaled intercept-only fits."""
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=3,
        require_rooted=True,
        require_binary=False,
    )
    fits = [
        _fit_evolutionary_mode_from_dataset(
            dataset,
            mode=mode,
            ou_bounds=ou_bounds,
            early_burst_bounds=early_burst_bounds,
        )
        for mode in ("brownian", "ornstein-uhlenbeck", "early-burst")
    ]
    rows = [
        _comparison_row(
            fit.mode,
            2 if fit.parameter_value is None else 3,
            fit.log_likelihood,
            fit.taxon_count,
        )
        for fit in fits
    ]
    best_aicc = min(row.aicc for row in rows)
    for row in rows:
        row.selected = math.isclose(row.aicc, best_aicc, rel_tol=0.0, abs_tol=1e-12)
    likelihood_ratio_tests = [
        _likelihood_ratio_test(
            comparison_id="brownian-vs-ornstein-uhlenbeck",
            left_fit=fits[0],
            right_fit=fits[1],
        ),
        _likelihood_ratio_test(
            comparison_id="brownian-vs-early-burst",
            left_fit=fits[0],
            right_fit=fits[2],
        ),
        _likelihood_ratio_test(
            comparison_id="ornstein-uhlenbeck-vs-early-burst",
            left_fit=fits[1],
            right_fit=fits[2],
        ),
    ]
    better_model = next(row.model for row in rows if row.selected)
    return ContinuousEvolutionaryModeComparisonReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=dataset.readiness.tree_taxa,
        rows=rows,
        better_model=better_model,
        likelihood_ratio_tests=likelihood_ratio_tests,
    )


def _fit_evolutionary_mode_from_dataset(
    dataset: ComparativeDataset,
    *,
    mode: str,
    ou_bounds: tuple[float, float],
    early_burst_bounds: tuple[float, float],
) -> ContinuousEvolutionaryModeFitReport:
    if mode not in ALLOWED_EVOLUTIONARY_MODES:
        raise ComparativeMethodError(
            "unsupported evolutionary mode; expected one of: "
            + ", ".join(sorted(ALLOWED_EVOLUTIONARY_MODES))
        )

    if mode == "brownian":
        transformed_tree = _clone_tree(dataset.tree)
        covariance = stable_covariance(
            build_brownian_covariance_matrix(transformed_tree, dataset.taxa)
        )
        fit = _fit_intercept_only_model(dataset, covariance)
        parameter_name = None
        parameter_value = None
        assumptions = [
            "Brownian mode retains the original rooted branch lengths.",
            "Trait variance accumulates proportionally with shared branch length.",
        ]
    elif mode == "ornstein-uhlenbeck":
        parameter_name = "alpha"
        parameter_value, transformed_tree, covariance = _best_transformed_mode_fit(
            dataset,
            mode=mode,
            bounds=ou_bounds,
        )
        fit = _fit_intercept_only_model(dataset, covariance)
        assumptions = [
            "OU mode follows the lecture-style tree rescaling before Brownian intercept fitting.",
            "Alpha is selected by maximizing log likelihood over a governed bounded search grid.",
        ]
    else:
        parameter_name = "rate_change"
        parameter_value, transformed_tree, covariance = _best_transformed_mode_fit(
            dataset,
            mode=mode,
            bounds=early_burst_bounds,
        )
        fit = _fit_intercept_only_model(dataset, covariance)
        assumptions = [
            "Early-burst mode follows the lecture-style tree rescaling before Brownian intercept fitting.",
            "The rate-change parameter is selected by maximizing log likelihood over a governed bounded search grid.",
        ]

    row = _comparison_row(
        mode,
        2 if parameter_value is None else 3,
        fit.log_likelihood,
        len(dataset.taxa),
    )
    intervals = _brownian_parameter_intervals(
        fit.theta,
        fit.sigma_squared,
        covariance,
    )
    if parameter_name is not None and parameter_value is not None:
        intervals = [
            *[
                interval
                for interval in intervals
                if getattr(interval, "name", None) != "root_state"
            ],
        ]
    residual_diagnostics = _build_residual_diagnostics(
        dataset,
        covariance,
        fit.residuals,
        fit.sigma_squared,
    )
    return ContinuousEvolutionaryModeFitReport(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        taxon_column=dataset.taxon_column,
        trait=dataset.trait,
        taxon_count=len(dataset.taxa),
        taxa=list(dataset.taxa),
        mode=mode,
        parameter_name=parameter_name,
        parameter_value=_stable_value(parameter_value)
        if parameter_value is not None
        else None,
        root_state=_stable_value(fit.theta),
        rate=_stable_value(fit.sigma_squared),
        log_likelihood=_stable_value(fit.log_likelihood),
        aic=_stable_value(row.aic),
        aicc=_stable_value(row.aicc),
        fitted_values=[_stable_value(value) for value in fit.fitted_values],
        residuals=[_stable_value(value) for value in fit.residuals],
        transformed_tree_newick=dumps_newick(transformed_tree),
        confidence_intervals=intervals,
        residual_diagnostics=residual_diagnostics,
        assumptions=assumptions,
    )


def _best_transformed_mode_fit(
    dataset: ComparativeDataset,
    *,
    mode: str,
    bounds: tuple[float, float],
) -> tuple[float, PhyloTree, list[list[float]]]:
    lower, upper = bounds
    if upper <= lower:
        raise ComparativeMethodError("parameter bounds must be strictly increasing")
    if mode == "ornstein-uhlenbeck":
        lower = max(lower, 1e-6)

    coarse = _linspace(lower, upper, 81)
    best_parameter = coarse[0]
    best_tree = _transform_tree(dataset.tree, mode=mode, parameter_value=best_parameter)
    best_covariance = stable_covariance(
        build_brownian_covariance_matrix(best_tree, dataset.taxa)
    )
    best_fit = _fit_intercept_only_model(dataset, best_covariance)
    best_index = 0
    for index, candidate in enumerate(coarse[1:], start=1):
        transformed_tree = _transform_tree(
            dataset.tree,
            mode=mode,
            parameter_value=candidate,
        )
        covariance = stable_covariance(
            build_brownian_covariance_matrix(transformed_tree, dataset.taxa)
        )
        fit = _fit_intercept_only_model(dataset, covariance)
        if fit.log_likelihood > best_fit.log_likelihood:
            best_parameter = candidate
            best_tree = transformed_tree
            best_covariance = covariance
            best_fit = fit
            best_index = index

    left = coarse[max(0, best_index - 1)]
    right = coarse[min(len(coarse) - 1, best_index + 1)]
    for candidate in _linspace(left, right, 81):
        if math.isclose(candidate, best_parameter, rel_tol=0.0, abs_tol=1e-12):
            continue
        transformed_tree = _transform_tree(
            dataset.tree,
            mode=mode,
            parameter_value=candidate,
        )
        covariance = stable_covariance(
            build_brownian_covariance_matrix(transformed_tree, dataset.taxa)
        )
        fit = _fit_intercept_only_model(dataset, covariance)
        if fit.log_likelihood > best_fit.log_likelihood:
            best_parameter = candidate
            best_tree = transformed_tree
            best_covariance = covariance
            best_fit = fit
    return best_parameter, best_tree, best_covariance


def _linspace(start: float, stop: float, count: int) -> list[float]:
    if count < 2:
        return [start]
    step = (stop - start) / float(count - 1)
    return [start + (step * index) for index in range(count)]


def _build_tree_rescaling_report(
    tree: PhyloTree,
    tree_path: Path,
    *,
    mode: str,
    parameter_name: str,
    parameter_value: float,
    sigsq: float,
) -> ComparativeTreeRescalingReport:
    transformed_tree = _transform_tree(
        tree,
        mode=mode,
        parameter_value=parameter_value,
        sigsq=sigsq,
    )
    branch_rows = _branch_length_rows(
        original_tree=tree,
        transformed_tree=transformed_tree,
    )
    return ComparativeTreeRescalingReport(
        tree_path=tree_path,
        mode=mode,
        parameter_name=parameter_name,
        parameter_value=_stable_value(parameter_value),
        tip_count=tree.tip_count,
        original_total_branch_length=_stable_value(tree.total_branch_length()),
        transformed_total_branch_length=_stable_value(
            transformed_tree.total_branch_length()
        ),
        transformed_tree_newick=dumps_newick(transformed_tree),
        branch_rows=branch_rows,
    )


def _transform_tree(
    tree: PhyloTree,
    *,
    mode: str,
    parameter_value: float,
    sigsq: float = 1.0,
) -> PhyloTree:
    if mode not in {"ornstein-uhlenbeck", "early-burst"}:
        raise ComparativeMethodError(
            "tree transformation mode must be 'ornstein-uhlenbeck' or 'early-burst'"
        )
    if parameter_value < 0.0:
        raise ComparativeMethodError("evolutionary mode parameter must be non-negative")
    cloned_root = _clone_node(tree.root)
    total_depth = _max_tip_depth(tree.root, depth=0.0)

    def visit(node: TreeNode, depth: float) -> None:
        for child in node.children:
            original_length = float(child.branch_length or 0.0)
            child_depth = depth + original_length
            if mode == "ornstein-uhlenbeck":
                child.branch_length = _ou_branch_length(
                    parent_depth=depth,
                    child_depth=child_depth,
                    total_depth=total_depth,
                    alpha=parameter_value,
                    sigsq=sigsq,
                )
            else:
                child.branch_length = _early_burst_branch_length(
                    parent_depth=depth,
                    child_depth=child_depth,
                    rate_change=parameter_value,
                    sigsq=sigsq,
                )
            visit(child, child_depth)

    visit(cloned_root, 0.0)
    return PhyloTree(
        root=cloned_root,
        source_format=tree.source_format,
        rooted=tree.rooted,
    )


def _branch_length_rows(
    *,
    original_tree: PhyloTree,
    transformed_tree: PhyloTree,
) -> list[EvolutionaryModeBranchLengthRow]:
    original_rows = _tree_branch_lookup(original_tree)
    transformed_rows = _tree_branch_lookup(transformed_tree)
    rows: list[EvolutionaryModeBranchLengthRow] = []
    for node_id in sorted(original_rows):
        original = original_rows[node_id]
        transformed = transformed_rows[node_id]
        rows.append(
            EvolutionaryModeBranchLengthRow(
                node=node_id,
                descendant_taxa=list(original["descendant_taxa"]),
                original_branch_length=_stable_value(original["branch_length"]),
                transformed_branch_length=_stable_value(transformed["branch_length"]),
                parent_depth=_stable_value(original["parent_depth"]),
                child_depth=_stable_value(original["child_depth"]),
            )
        )
    return rows


def _tree_branch_lookup(tree: PhyloTree) -> dict[str, dict[str, object]]:
    rows: dict[str, dict[str, object]] = {}

    def visit(node: TreeNode, depth: float) -> None:
        for child in node.children:
            branch_length = float(child.branch_length or 0.0)
            child_depth = depth + branch_length
            branch_id = node_signature(child)
            rows[branch_id] = {
                "branch_length": branch_length,
                "parent_depth": depth,
                "child_depth": child_depth,
                "descendant_taxa": _descendant_taxa(child),
            }
            visit(child, child_depth)

    visit(tree.root, 0.0)
    return rows


def _ou_branch_length(
    *,
    parent_depth: float,
    child_depth: float,
    total_depth: float,
    alpha: float,
    sigsq: float,
) -> float:
    if alpha <= 0.0:
        raise ComparativeMethodError("OU alpha must be positive")

    def _term(depth: float) -> float:
        return (
            (1.0 / (2.0 * alpha))
            * math.exp(-2.0 * alpha * (total_depth - depth))
            * (1.0 - math.exp(-2.0 * alpha * depth))
        )

    return max(0.0, (_term(child_depth) - _term(parent_depth)) * sigsq)


def _early_burst_branch_length(
    *,
    parent_depth: float,
    child_depth: float,
    rate_change: float,
    sigsq: float,
) -> float:
    if math.isclose(rate_change, 0.0, rel_tol=0.0, abs_tol=1e-12):
        return max(0.0, (child_depth - parent_depth) * sigsq)
    transformed = (
        (math.exp(rate_change * child_depth) - math.exp(rate_change * parent_depth))
        / rate_change
    )
    return max(0.0, transformed * sigsq)


def _likelihood_ratio_test(
    *,
    comparison_id: str,
    left_fit: ContinuousEvolutionaryModeFitReport,
    right_fit: ContinuousEvolutionaryModeFitReport,
) -> LikelihoodRatioTestResult:
    statistic = max(0.0, -2.0 * (left_fit.log_likelihood - right_fit.log_likelihood))
    return LikelihoodRatioTestResult(
        comparison_id=comparison_id,
        left_mode=left_fit.mode,
        right_mode=right_fit.mode,
        statistic=_stable_value(statistic),
        degrees_of_freedom=1,
        p_value=_stable_value(math.erfc(math.sqrt(statistic / 2.0))),
    )


def _clone_node(node: TreeNode) -> TreeNode:
    return TreeNode(
        name=node.name,
        branch_length=node.branch_length,
        children=[_clone_node(child) for child in node.children],
    )


def _clone_tree(tree: PhyloTree) -> PhyloTree:
    return PhyloTree(
        root=_clone_node(tree.root),
        source_format=tree.source_format,
        rooted=tree.rooted,
    )


def _max_tip_depth(node: TreeNode, *, depth: float) -> float:
    if node.is_leaf():
        return depth
    return max(
        _max_tip_depth(
            child,
            depth=depth + float(child.branch_length or 0.0),
        )
        for child in node.children
    )


def _descendant_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_descendant_taxa(child))
    return sorted(taxa)


def _stable_value(value: float | None) -> float:
    if value is None:
        raise ValueError("expected a float value")
    return float(format(round(float(value), 15), ".15g"))
