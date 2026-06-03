from __future__ import annotations

import math

from bijux_phylogenetics.comparative.common import summarize_numeric_trait_readiness
from bijux_phylogenetics.comparative.evolutionary_modes import (
    CONTINUOUS_GAUSSIAN_LIKELIHOOD_CONSTANT_POLICY,
    FITCONTINUOUS_MODEL_RANKING_POLICY,
    ContinuousModeSearchControls,
    compare_fitcontinuous_model_ranking,
    fit_continuous_evolutionary_mode,
)
from bijux_phylogenetics.parity.geiger.registry import GeigerParityCase

from .comparison import comparison_rows, parameter_rows
from .payload_policy import (
    bijux_optimizer_result,
    comparison_modes,
    missing_value_policy,
    parameter_bound_policy,
    standard_error_policy,
)


def build_bijux_continuous_case_payload(
    case: GeigerParityCase,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    tree_path, traits_path = case.input_fixtures
    readiness = summarize_numeric_trait_readiness(
        tree_path,
        traits_path,
        trait=case.trait_name,
        taxon_column=case.taxon_column,
    )
    search_controls = None
    if (
        case.coarse_grid_point_count is not None
        or case.fine_grid_point_count is not None
        or case.initial_parameter_value is not None
    ):
        search_controls = ContinuousModeSearchControls(
            coarse_grid_point_count=(
                81
                if case.coarse_grid_point_count is None
                else case.coarse_grid_point_count
            ),
            fine_grid_point_count=(
                81 if case.fine_grid_point_count is None else case.fine_grid_point_count
            ),
            initial_parameter_value=case.initial_parameter_value,
        )
    report = fit_continuous_evolutionary_mode(
        tree_path,
        traits_path,
        trait=case.trait_name,
        mode=case.python_mode,
        taxon_column=case.taxon_column,
        search_controls=search_controls,
        lambda_bounds=(0.0, 1.0) if case.lambda_bounds is None else case.lambda_bounds,
        kappa_bounds=(0.0, 3.0) if case.kappa_bounds is None else case.kappa_bounds,
        delta_bounds=(0.0, 3.0) if case.delta_bounds is None else case.delta_bounds,
        ou_bounds=(0.0, 10.0) if case.ou_bounds is None else case.ou_bounds,
        early_burst_bounds=(0.0, 50.0)
        if case.early_burst_bounds is None
        else case.early_burst_bounds,
    )
    excluded_taxa = sorted(
        {
            *readiness.missing_from_traits,
            *readiness.pruned_missing_value_taxa,
            *readiness.pruned_non_numeric_taxa,
        }
    )
    summary = {
        "taxon_count": report.taxon_count,
        "trait_name": report.trait,
        "model_name": case.model_name,
        "excluded_taxon_count": len(excluded_taxa),
        "excluded_taxa": excluded_taxa,
        "missing_value_taxa": list(readiness.pruned_missing_value_taxa),
        "non_numeric_taxa": list(readiness.pruned_non_numeric_taxa),
        "missing_from_traits": list(readiness.missing_from_traits),
        "extra_trait_taxa": list(readiness.extra_trait_taxa),
        "missing_value_policy": missing_value_policy(),
        "standard_error_policy": standard_error_policy(),
        "parameter_bound_policy": parameter_bound_policy(case),
        "hit_lower_parameter_boundary": (
            False
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.hit_lower_boundary
        ),
        "hit_upper_parameter_boundary": (
            False
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.hit_upper_boundary
        ),
        "identifiability_warning_kinds": [
            warning.kind for warning in report.identifiability_warnings
        ],
        "identifiability_warning_count": len(report.identifiability_warnings),
        "root_state": report.root_state,
        "rate": report.rate,
        "log_likelihood": report.log_likelihood,
        "aic": report.aic,
        "aicc": report.aicc,
        "likelihood_constant_policy": report.likelihood_constant_policy,
        "likelihood_comparison_policy": report.likelihood_comparison_policy,
        "parameter_name": report.parameter_name,
        "parameter_value": report.parameter_value,
        "optimizer_settings": case.optimizer_settings,
        "optimizer_result": bijux_optimizer_result(case, report),
    }
    return summary, parameter_rows(summary)


def build_bijux_model_comparison_payload(
    case: GeigerParityCase,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    tree_path, traits_path = case.input_fixtures
    report = compare_fitcontinuous_model_ranking(
        tree_path,
        traits_path,
        trait=case.trait_name,
        taxon_column=case.taxon_column,
        modes=comparison_modes(case.candidate_model_names),
        lambda_bounds=(0.0, 1.0) if case.lambda_bounds is None else case.lambda_bounds,
        kappa_bounds=(0.0, 3.0) if case.kappa_bounds is None else case.kappa_bounds,
        delta_bounds=(0.0, 3.0) if case.delta_bounds is None else case.delta_bounds,
        ou_bounds=(0.0, 10.0) if case.ou_bounds is None else case.ou_bounds,
        early_burst_bounds=(0.0, 50.0)
        if case.early_burst_bounds is None
        else case.early_burst_bounds,
    )
    runner_up_rows = [
        row
        for row in report.rows
        if row.model != report.better_model and row.comparable and row.rank is not None
    ]
    runner_up_row = runner_up_rows[0] if runner_up_rows else None
    summary = {
        "taxon_count": report.taxon_count,
        "trait_name": report.trait,
        "model_name": case.model_name,
        "selected_model": report.better_model,
        "model_ranking": [row.model for row in report.rows],
        "comparable_model_count": sum(1 for row in report.rows if row.comparable),
        "noncomparable_model_count": sum(
            1 for row in report.rows if not row.comparable
        ),
        "runner_up_model": None if runner_up_row is None else runner_up_row.model,
        "runner_up_aicc_delta": (
            math.nan if runner_up_row is None else runner_up_row.delta_aicc
        ),
        "warning_count": len(report.warnings),
        "likelihood_constant_policy": (
            report.likelihood_constant_policy
            if report.likelihood_constant_policy is not None
            else CONTINUOUS_GAUSSIAN_LIKELIHOOD_CONSTANT_POLICY
        ),
        "likelihood_comparison_policy": FITCONTINUOUS_MODEL_RANKING_POLICY,
        "noncomparable_likelihood_models": list(report.noncomparable_likelihood_models),
        "optimizer_settings": case.optimizer_settings,
    }
    return summary, comparison_rows(report)
