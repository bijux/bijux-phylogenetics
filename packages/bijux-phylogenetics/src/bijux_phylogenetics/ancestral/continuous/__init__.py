from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    AncestralContinuousDataset,
    dump_pruned_tree,
    load_continuous_dataset,
    node_descendant_taxa,
    node_signature,
    stable_value,
)
from bijux_phylogenetics.core.topology import _root_tree_by_outgroup_node
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.core.ultrametric import summarize_ultrametric_tip_depths
from .models import (
    ContinuousAncestralBrownianFitDiagnostics as ContinuousAncestralBrownianFitDiagnostics,
    ContinuousAncestralEstimate as ContinuousAncestralEstimate,
    ContinuousAncestralExclusion as ContinuousAncestralExclusion,
    ContinuousAncestralOptimizerDiagnostics as ContinuousAncestralOptimizerDiagnostics,
    ContinuousAncestralReport as ContinuousAncestralReport,
    ContinuousAncestralSummary as ContinuousAncestralSummary,
)
from .reporting import (
    continuous_ancestral_exclusions as continuous_ancestral_exclusions,
    summarize_continuous_ancestral_report as summarize_continuous_ancestral_report,
    write_continuous_ancestral_exclusion_table as write_continuous_ancestral_exclusion_table,
    write_continuous_ancestral_summary_table as write_continuous_ancestral_summary_table,
    write_continuous_ancestral_uncertainty_table as write_continuous_ancestral_uncertainty_table,
)

_NORMAL_95_CRITICAL = 1.959963984540054
_BROWNIAN_CONDITION_THRESHOLD = 1e12
_SOLVER_REGULARIZATION_EPSILON = 1e-8


def reconstruct_continuous_ancestral_states(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "brownian",
    estimator: str | None = None,
    alpha: float = 1.0,
) -> ContinuousAncestralReport:
    """Reconstruct continuous ancestral states under a Brownian or OU-style model."""
    if model not in {"brownian", "ou"}:
        raise ValueError(f"unsupported continuous ancestral model: {model}")
    if alpha <= 0:
        raise ValueError(
            f"alpha must be positive for continuous ancestral reconstruction, got {alpha}"
        )
    dataset = load_continuous_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    return reconstruct_continuous_ancestral_states_from_dataset(
        dataset,
        model=model,
        estimator=estimator,
        alpha=alpha,
    )


def reconstruct_continuous_ancestral_states_from_dataset(
    dataset: AncestralContinuousDataset,
    *,
    model: str = "brownian",
    estimator: str | None = None,
    alpha: float = 1.0,
) -> ContinuousAncestralReport:
    """Reconstruct continuous ancestral states from one native ancestral dataset."""
    if model not in {"brownian", "ou"}:
        raise ValueError(f"unsupported continuous ancestral model: {model}")
    if alpha <= 0:
        raise ValueError(
            f"alpha must be positive for continuous ancestral reconstruction, got {alpha}"
        )
    resolved_estimator = _resolve_continuous_estimator(model, estimator)
    brownian_fit_diagnostics = (
        _summarize_brownian_fit_diagnostics(dataset) if model == "brownian" else None
    )
    optimizer_diagnostics: ContinuousAncestralOptimizerDiagnostics | None = None
    anc_ml_profile_fit: _ContinuousAncMlProfileFit | None = None
    if model == "brownian" and resolved_estimator == "anc-ml":
        anc_ml_profile_fit = _fit_continuous_anc_ml_profile(dataset)
        brownian_fit_diagnostics = _apply_anc_ml_profile_to_brownian_diagnostics(
            brownian_fit_diagnostics,
            anc_ml_profile_fit,
        )
        optimizer_diagnostics = anc_ml_profile_fit.optimizer_diagnostics
    return _reconstruct_continuous_from_dataset(
        dataset,
        working_tree=dataset.tree,
        model=model,
        estimator=resolved_estimator,
        alpha=alpha,
        brownian_fit_diagnostics=brownian_fit_diagnostics,
        optimizer_diagnostics=optimizer_diagnostics,
        anc_ml_profile_fit=anc_ml_profile_fit,
    )


def _reconstruct_continuous_from_dataset(
    dataset: AncestralContinuousDataset,
    *,
    working_tree: PhyloTree,
    model: str,
    estimator: str,
    alpha: float,
    brownian_fit_diagnostics: ContinuousAncestralBrownianFitDiagnostics | None,
    optimizer_diagnostics: ContinuousAncestralOptimizerDiagnostics | None,
    anc_ml_profile_fit: _ContinuousAncMlProfileFit | None,
) -> ContinuousAncestralReport:
    global_mean = sum(dataset.values_by_taxon[taxon] for taxon in dataset.taxa) / len(
        dataset.taxa
    )
    sigma = _sample_standard_deviation(
        [dataset.values_by_taxon[taxon] for taxon in dataset.taxa]
    )
    trait_range = (
        max(dataset.values_by_taxon.values()) - min(dataset.values_by_taxon.values())
        if dataset.values_by_taxon
        else 0.0
    )
    if estimator == "fast-anc":
        estimates = _build_fast_anc_estimates(
            dataset,
            trait_range=trait_range,
        )
    elif estimator == "anc-ml":
        if anc_ml_profile_fit is None:
            raise ValueError(
                "anc-ml continuous ancestral reconstruction requires one anc-ml profile fit"
            )
        estimates = _build_anc_ml_estimates(
            dataset,
            trait_range=trait_range,
            profile_fit=anc_ml_profile_fit,
        )
    else:
        estimates = _build_local_continuous_estimates(
            dataset,
            working_tree=working_tree,
            model=model,
            alpha=alpha,
            global_mean=global_mean,
            sigma=sigma,
            trait_range=trait_range,
        )
    ordered_estimates = _ordered_estimates(dataset, estimates)
    unstable_nodes = [
        estimate.node
        for estimate in ordered_estimates
        if not estimate.is_tip and estimate.unstable
    ]
    weak_support_nodes = [
        estimate.node
        for estimate in ordered_estimates
        if not estimate.is_tip and estimate.confidence < 0.75
    ]
    warnings = list(dataset.warnings)
    if unstable_nodes:
        warnings.append(
            "one or more continuous ancestral estimates have broad uncertainty intervals"
        )
    if weak_support_nodes:
        warnings.append(
            "low-confidence ancestral estimates should not be overinterpreted for evolutionary timing or trait polarity"
        )
    if (
        brownian_fit_diagnostics is not None
        and brownian_fit_diagnostics.covariance_near_singular
    ):
        warnings.append(
            "Brownian covariance diagnostics indicate a singular or ill-conditioned fit, so ancestral uncertainty should be interpreted cautiously"
        )
    if (
        brownian_fit_diagnostics is not None
        and brownian_fit_diagnostics.solver_regularized
    ):
        warnings.append(
            "Brownian covariance inversion required light diagonal regularization for numerical stability"
        )
    return ContinuousAncestralReport(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        taxon_column=dataset.taxon_column,
        trait=dataset.trait,
        model=model,
        estimator=estimator,
        alpha=stable_value(alpha),
        taxon_count=len(dataset.taxa),
        analysis_tree_newick=dump_pruned_tree(working_tree),
        dropped_missing_taxa=dataset.dropped_missing_taxa,
        dropped_non_numeric_taxa=dataset.dropped_non_numeric_taxa,
        warnings=warnings,
        unstable_nodes=unstable_nodes,
        weak_support_nodes=weak_support_nodes,
        brownian_fit_diagnostics=brownian_fit_diagnostics,
        optimizer_diagnostics=optimizer_diagnostics,
        estimates=ordered_estimates,
    )


def _sample_standard_deviation(values: list[float]) -> float:
    if len(values) < 2:
        return 1.0
    mean_value = sum(values) / len(values)
    variance = sum((value - mean_value) ** 2 for value in values) / (len(values) - 1)
    return math.sqrt(max(variance, 1e-12))


def _ordered_estimates(
    dataset: AncestralContinuousDataset,
    estimates: list[ContinuousAncestralEstimate],
) -> list[ContinuousAncestralEstimate]:
    node_order = {
        node_signature(node): index
        for index, node in enumerate(dataset.tree.iter_nodes())
    }
    return sorted(estimates, key=lambda estimate: node_order[estimate.node])


def _resolve_continuous_estimator(model: str, estimator: str | None) -> str:
    default_estimators = {
        "brownian": "ace-pic",
        "ou": "generalized-least-squares",
    }
    allowed_estimators = {
        "brownian": {"ace-pic", "anc-ml", "fast-anc"},
        "ou": {"generalized-least-squares"},
    }
    resolved = default_estimators[model] if estimator is None else estimator
    if resolved not in allowed_estimators[model]:
        supported = ", ".join(sorted(allowed_estimators[model]))
        raise ValueError(
            f"unsupported continuous ancestral estimator '{resolved}' for model '{model}'; expected one of: {supported}"
        )
    return resolved


def _build_local_continuous_estimates(
    dataset: AncestralContinuousDataset,
    *,
    working_tree: PhyloTree,
    model: str,
    alpha: float,
    global_mean: float,
    sigma: float,
    trait_range: float,
) -> list[ContinuousAncestralEstimate]:
    estimates: list[ContinuousAncestralEstimate] = []

    def visit(node: TreeNode) -> tuple[float, float]:
        if node.is_leaf():
            estimate = dataset.values_by_taxon[node.name]
            estimates.append(
                ContinuousAncestralEstimate(
                    node=node_signature(node),
                    node_name=node.name,
                    is_tip=True,
                    descendant_taxa=node_descendant_taxa(node),
                    estimate=stable_value(estimate),
                    standard_error=0.0,
                    lower_95_interval=stable_value(estimate),
                    upper_95_interval=stable_value(estimate),
                    uncertainty_width=0.0,
                    confidence=1.0,
                    interpretation="observed tip value",
                    unstable=False,
                    downstream_risks=[],
                )
            )
            return estimate, float(node.branch_length or 0.0)

        if len(node.children) != 2:
            raise ValueError(
                "continuous ancestral reconstruction requires a fully dichotomous rooted tree"
            )

        left_child, right_child = node.children
        left_estimate, left_working_length = visit(left_child)
        right_estimate, right_working_length = visit(right_child)

        if model == "brownian":
            sum_working_lengths = max(
                left_working_length + right_working_length,
                1e-12,
            )
            estimate = (
                left_estimate * right_working_length
                + right_estimate * left_working_length
            ) / sum_working_lengths
            variance = sum_working_lengths
            returned_length = (
                float(node.branch_length or 0.0)
                + (left_working_length * right_working_length) / sum_working_lengths
            )
        else:
            child_payloads: list[tuple[float, float]] = []
            for child, child_estimate, child_variance in (
                (left_child, left_estimate, left_working_length),
                (right_child, right_estimate, right_working_length),
            ):
                branch_length = float(child.branch_length or 0.0)
                shrink = math.exp(-alpha * branch_length)
                transformed_estimate = (
                    shrink * child_estimate + (1.0 - shrink) * global_mean
                )
                stationary_variance = ((sigma**2) / (2.0 * alpha)) * (
                    1.0 - math.exp(-2.0 * alpha * branch_length)
                )
                propagated_variance = (
                    child_variance * math.exp(-2.0 * alpha * branch_length)
                ) + stationary_variance
                child_payloads.append(
                    (transformed_estimate, max(propagated_variance, 1e-12))
                )
            weight_sum = sum(
                1.0 / child_variance for _, child_variance in child_payloads
            )
            estimate = (
                sum(
                    (value / child_variance) for value, child_variance in child_payloads
                )
                / weight_sum
            )
            variance = 1.0 / weight_sum
            returned_length = variance

        standard_error = math.sqrt(max(variance, 0.0))
        lower = estimate - _NORMAL_95_CRITICAL * standard_error
        upper = estimate + _NORMAL_95_CRITICAL * standard_error
        uncertainty_width = max(0.0, upper - lower)
        confidence, unstable = _continuous_confidence(uncertainty_width, trait_range)
        interpretation = _continuous_interpretation(
            uncertainty_width,
            trait_range,
            unstable=unstable,
        )
        estimates.append(
            ContinuousAncestralEstimate(
                node=node_signature(node),
                node_name=node.name,
                is_tip=False,
                descendant_taxa=node_descendant_taxa(node),
                estimate=stable_value(estimate),
                standard_error=stable_value(standard_error),
                lower_95_interval=stable_value(lower),
                upper_95_interval=stable_value(upper),
                uncertainty_width=stable_value(uncertainty_width),
                confidence=stable_value(confidence),
                interpretation=interpretation,
                unstable=unstable,
                downstream_risks=_continuous_downstream_risks(unstable),
            )
        )
        return estimate, returned_length

    visit(working_tree.root)
    return estimates


def _build_fast_anc_estimates(
    dataset: AncestralContinuousDataset,
    *,
    trait_range: float,
) -> list[ContinuousAncestralEstimate]:
    estimates: list[ContinuousAncestralEstimate] = []
    for node in dataset.tree.iter_leaves():
        estimate = dataset.values_by_taxon[node.name]
        estimates.append(
            ContinuousAncestralEstimate(
                node=node_signature(node),
                node_name=node.name,
                is_tip=True,
                descendant_taxa=node_descendant_taxa(node),
                estimate=stable_value(estimate),
                standard_error=0.0,
                lower_95_interval=stable_value(estimate),
                upper_95_interval=stable_value(estimate),
                uncertainty_width=0.0,
                confidence=1.0,
                interpretation="observed tip value",
                unstable=False,
                downstream_risks=[],
            )
        )
    for node in dataset.tree.iter_internal_nodes(order="preorder"):
        rerooted_tree = _root_tree_by_outgroup_node(dataset.tree, outgroup_node=node)
        if len(rerooted_tree.root.children) != 2:
            raise ValueError(
                "fast ancestral reconstruction requires a fully dichotomous rerooted tree"
            )
        left_child, right_child = rerooted_tree.root.children
        left_contrasts, left_value, left_variance = _compute_fast_anc_pic_payload(
            left_child,
            dataset.values_by_taxon,
        )
        right_contrasts, right_value, right_variance = _compute_fast_anc_pic_payload(
            right_child,
            dataset.values_by_taxon,
        )
        root_expected_variance = max(left_variance + right_variance, 1e-12)
        estimate = ((left_value / left_variance) + (right_value / right_variance)) / (
            (1.0 / left_variance) + (1.0 / right_variance)
        )
        root_contrast = (left_value - right_value) / math.sqrt(root_expected_variance)
        contrast_sum = sum(
            contrast * contrast for contrast in left_contrasts + right_contrasts
        ) + (root_contrast * root_contrast)
        contrast_count = len(left_contrasts) + len(right_contrasts) + 1
        contrast_mean_square = contrast_sum / contrast_count
        variance = (
            ((1.0 / left_variance) + (1.0 / right_variance)) ** -1
        ) * contrast_mean_square
        standard_error = math.sqrt(max(variance, 0.0))
        lower = estimate - _NORMAL_95_CRITICAL * standard_error
        upper = estimate + _NORMAL_95_CRITICAL * standard_error
        uncertainty_width = max(0.0, upper - lower)
        confidence, unstable = _continuous_confidence(uncertainty_width, trait_range)
        interpretation = _continuous_interpretation(
            uncertainty_width,
            trait_range,
            unstable=unstable,
        )
        estimates.append(
            ContinuousAncestralEstimate(
                node=node_signature(node),
                node_name=node.name,
                is_tip=False,
                descendant_taxa=node_descendant_taxa(node),
                estimate=stable_value(estimate),
                standard_error=stable_value(standard_error),
                lower_95_interval=stable_value(lower),
                upper_95_interval=stable_value(upper),
                uncertainty_width=stable_value(uncertainty_width),
                confidence=stable_value(confidence),
                interpretation=interpretation,
                unstable=unstable,
                downstream_risks=_continuous_downstream_risks(unstable),
            )
        )
    return estimates


@dataclass(slots=True)
class _ContinuousAncMlProfileFit:
    root_state: float
    internal_states_by_node: dict[str, float]
    state_variances_by_node: dict[str, float]
    log_likelihood: float
    sigma_squared: float
    optimizer_diagnostics: ContinuousAncestralOptimizerDiagnostics


def _build_anc_ml_estimates(
    dataset: AncestralContinuousDataset,
    *,
    trait_range: float,
    profile_fit: _ContinuousAncMlProfileFit,
) -> list[ContinuousAncestralEstimate]:
    anc_ml_interval_critical = 1.96
    estimates: list[ContinuousAncestralEstimate] = []
    for node in dataset.tree.iter_leaves():
        estimate = dataset.values_by_taxon[node.name]
        estimates.append(
            ContinuousAncestralEstimate(
                node=node_signature(node),
                node_name=node.name,
                is_tip=True,
                descendant_taxa=node_descendant_taxa(node),
                estimate=stable_value(estimate),
                standard_error=0.0,
                lower_95_interval=stable_value(estimate),
                upper_95_interval=stable_value(estimate),
                uncertainty_width=0.0,
                confidence=1.0,
                interpretation="observed tip value",
                unstable=False,
                downstream_risks=[],
            )
        )
    for node in dataset.tree.iter_internal_nodes(order="preorder"):
        node_key = node_signature(node)
        estimate = (
            profile_fit.root_state
            if node is dataset.tree.root
            else profile_fit.internal_states_by_node[node_key]
        )
        variance = profile_fit.state_variances_by_node[node_key]
        standard_error = math.sqrt(max(variance, 0.0))
        lower = estimate - anc_ml_interval_critical * standard_error
        upper = estimate + anc_ml_interval_critical * standard_error
        uncertainty_width = max(0.0, upper - lower)
        confidence, unstable = _continuous_confidence(uncertainty_width, trait_range)
        interpretation = _continuous_interpretation(
            uncertainty_width,
            trait_range,
            unstable=unstable,
        )
        estimates.append(
            ContinuousAncestralEstimate(
                node=node_key,
                node_name=node.name,
                is_tip=False,
                descendant_taxa=node_descendant_taxa(node),
                estimate=stable_value(estimate),
                standard_error=stable_value(standard_error),
                lower_95_interval=stable_value(lower),
                upper_95_interval=stable_value(upper),
                uncertainty_width=stable_value(uncertainty_width),
                confidence=stable_value(confidence),
                interpretation=interpretation,
                unstable=unstable,
                downstream_risks=_continuous_downstream_risks(unstable),
            )
        )
    return estimates


def _fit_continuous_anc_ml_profile(
    dataset: AncestralContinuousDataset,
) -> _ContinuousAncMlProfileFit:
    ordered_nodes, covariance_matrix = _build_anc_ml_covariance_matrix(
        dataset.tree,
        dataset.taxa,
    )
    try:
        inverse_covariance = _invert_matrix(covariance_matrix)
    except ValueError as error:
        raise ValueError(
            "continuous ancestral anc-ml reconstruction requires a nonsingular full Brownian covariance matrix"
        ) from error
    tip_count = len(dataset.taxa)
    internal_nodes = ordered_nodes[tip_count:]
    internal_count = len(internal_nodes)
    precision_tips = [row[:tip_count] for row in inverse_covariance[:tip_count]]
    precision_tip_internal = [row[tip_count:] for row in inverse_covariance[:tip_count]]
    precision_internal_tip = [row[:tip_count] for row in inverse_covariance[tip_count:]]
    precision_internal = [row[tip_count:] for row in inverse_covariance[tip_count:]]
    inverse_internal_precision = (
        _invert_matrix(precision_internal) if internal_count > 0 else []
    )
    schur_projection = [
        [
            sum(
                precision_tip_internal[row_index][shared_index]
                * inverse_internal_precision[shared_index][column_index]
                for shared_index in range(internal_count)
            )
            for column_index in range(internal_count)
        ]
        for row_index in range(tip_count)
    ]
    schur_precision = [
        [
            precision_tips[row_index][column_index]
            - sum(
                schur_projection[row_index][shared_index]
                * precision_internal_tip[shared_index][column_index]
                for shared_index in range(internal_count)
            )
            for column_index in range(tip_count)
        ]
        for row_index in range(tip_count)
    ]
    tip_values = [dataset.values_by_taxon[taxon] for taxon in dataset.taxa]
    unit_vector = [1.0] * tip_count
    denominator = _quadratic_form(unit_vector, schur_precision)
    root_state = (
        sum(
            unit_vector[row_index]
            * sum(
                schur_precision[row_index][column_index] * tip_values[column_index]
                for column_index in range(tip_count)
            )
            for row_index in range(tip_count)
        )
        / denominator
    )
    tip_residuals = [value - root_state for value in tip_values]
    internal_states_by_node: dict[str, float] = {}
    centered_internal_states: list[float] = []
    if internal_count > 0:
        internal_rhs = [
            sum(
                precision_internal_tip[row_index][column_index]
                * tip_residuals[column_index]
                for column_index in range(tip_count)
            )
            for row_index in range(internal_count)
        ]
        centered_internal_states = [
            -sum(
                inverse_internal_precision[row_index][shared_index]
                * internal_rhs[shared_index]
                for shared_index in range(internal_count)
            )
            for row_index in range(internal_count)
        ]
        internal_states_by_node = {
            node_signature(node): stable_value(
                root_state + centered_internal_states[index]
            )
            for index, node in enumerate(internal_nodes)
        }
    completed_residuals = tip_residuals + centered_internal_states
    sigma_squared = max(
        _quadratic_form(completed_residuals, inverse_covariance)
        / len(covariance_matrix),
        1e-12,
    )
    log_likelihood = -0.5 * (
        len(covariance_matrix) * math.log(2.0 * math.pi * sigma_squared)
        + _log_determinant(covariance_matrix)
        + len(covariance_matrix)
    )
    parameter_design = _build_anc_ml_parameter_design(tip_count, internal_count)
    design_precision = _design_quadratic_form(parameter_design, inverse_covariance)
    try:
        state_covariance = [
            [sigma_squared * value for value in row]
            for row in _invert_matrix(design_precision)
        ]
    except ValueError as error:
        raise ValueError(
            "continuous ancestral anc-ml state uncertainty requires an invertible parameter-information matrix"
        ) from error
    state_variances_by_node = {
        node_signature(dataset.tree.root): stable_value(state_covariance[0][0])
    }
    for index, node in enumerate(internal_nodes):
        state_variances_by_node[node_signature(node)] = stable_value(
            state_covariance[index + 1][index + 1]
        )
    return _ContinuousAncMlProfileFit(
        root_state=stable_value(root_state),
        internal_states_by_node=internal_states_by_node,
        state_variances_by_node=state_variances_by_node,
        log_likelihood=stable_value(log_likelihood),
        sigma_squared=stable_value(sigma_squared),
        optimizer_diagnostics=ContinuousAncestralOptimizerDiagnostics(
            optimizer_name="closed-form-profile-solution",
            converged=True,
            iteration_count=0,
            function_evaluation_count=1,
            convergence_status="profile-solved",
            message="closed-form Brownian maximum-likelihood profile solved without iterative optimization",
        ),
    )


def _apply_anc_ml_profile_to_brownian_diagnostics(
    diagnostics: ContinuousAncestralBrownianFitDiagnostics | None,
    profile_fit: _ContinuousAncMlProfileFit,
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
        log_likelihood=profile_fit.log_likelihood,
        residual_sigma_squared=profile_fit.sigma_squared,
    )


def _compute_fast_anc_pic_payload(
    node: TreeNode,
    values_by_taxon: dict[str, float],
) -> tuple[list[float], float, float]:
    if node.is_leaf():
        if node.name is None:
            raise ValueError(
                "leaf taxon name is required for fast ancestral reconstruction"
            )
        if node.branch_length is None:
            raise ValueError(
                "branch lengths are required for fast ancestral reconstruction"
            )
        return [], values_by_taxon[node.name], node.branch_length
    if len(node.children) != 2:
        raise ValueError(
            "fast ancestral reconstruction requires a strictly binary tree"
        )
    left_contrasts, left_value, left_variance = _compute_fast_anc_pic_payload(
        node.children[0], values_by_taxon
    )
    right_contrasts, right_value, right_variance = _compute_fast_anc_pic_payload(
        node.children[1], values_by_taxon
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
    return (
        left_contrasts + right_contrasts + [contrast],
        ancestral_value,
        propagated_variance,
    )


def _continuous_confidence(
    uncertainty_width: float, trait_range: float
) -> tuple[float, bool]:
    if uncertainty_width == 0.0:
        return 1.0, False
    scale = max(trait_range, 1e-12)
    relative_width = uncertainty_width / scale
    confidence = max(0.0, min(1.0, 1.0 - min(relative_width, 1.0)))
    return stable_value(confidence), relative_width > 0.6 or confidence < 0.55


def _continuous_interpretation(
    uncertainty_width: float, trait_range: float, *, unstable: bool
) -> str:
    if uncertainty_width == 0.0:
        return "observed tip value"
    scale = max(trait_range, 1e-12)
    relative_width = uncertainty_width / scale
    if unstable:
        return "unstable node estimate"
    if relative_width <= 0.25:
        return "narrow uncertainty"
    if relative_width <= 0.6:
        return "moderate uncertainty"
    return "broad uncertainty"


def _continuous_downstream_risks(unstable: bool) -> list[str]:
    if not unstable:
        return []
    return [
        "node ordering and trait-polarity interpretations may change across alternative trees or models",
        "publication claims about deep ancestral values should be treated as provisional",
    ]


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
