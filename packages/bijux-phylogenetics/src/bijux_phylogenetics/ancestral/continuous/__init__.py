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
from .diagnostics import (
    _apply_anc_ml_profile_to_brownian_diagnostics,
    _build_anc_ml_covariance_matrix,
    _build_anc_ml_parameter_design,
    _design_quadratic_form,
    _invert_matrix,
    _log_determinant,
    _quadratic_form,
    _summarize_brownian_fit_diagnostics,
)
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
            log_likelihood=anc_ml_profile_fit.log_likelihood,
            sigma_squared=anc_ml_profile_fit.sigma_squared,
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
