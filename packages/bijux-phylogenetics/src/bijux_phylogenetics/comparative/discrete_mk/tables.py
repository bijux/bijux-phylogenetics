from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import write_ancestral_rows

from .models import DiscreteMkFitReport


def write_discrete_mk_summary_table(path: Path, report: DiscreteMkFitReport) -> Path:
    """Write one flat summary ledger for a discrete Mk fit."""
    baseline = report.baseline_comparison
    diagnostics = report.optimizer_diagnostics
    transform_fit = report.transform_fit
    transform_baseline = report.transform_baseline_comparison
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "ascertainment_policy",
            "ascertainment_conditioning_log_probability",
            "invariant_pattern_log_probability",
            "transform",
            "transform_parameter_name",
            "transform_parameter_value",
            "transform_lower_bound",
            "transform_upper_bound",
            "transform_starting_parameter_policy",
            "transform_starting_parameter_value",
            "transform_starting_parameter_log_likelihood",
            "transform_coarse_grid_point_count",
            "transform_fine_grid_point_count",
            "transform_refinement_start_count",
            "transform_function_evaluation_count",
            "transform_hit_lower_parameter_boundary",
            "transform_hit_upper_parameter_boundary",
            "transform_warning_count",
            "transform_tree_is_ultrametric",
            "transform_tree_minimum_tip_depth",
            "transform_tree_maximum_tip_depth",
            "transform_baseline",
            "transform_baseline_parameter_count",
            "transform_baseline_log_likelihood",
            "transform_baseline_aic",
            "transform_delta_log_likelihood",
            "transform_delta_aic",
            "preferred_transform_by_aic",
            "state_ordering",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "observed_state_count",
            "sparse_state_count",
            "log_likelihood",
            "parameter_count",
            "aic",
            "aicc",
            "likelihood_constant_policy",
            "likelihood_comparison_policy",
            "optimizer_name",
            "optimizer_converged",
            "optimizer_iteration_count",
            "optimizer_function_evaluation_count",
            "optimizer_simplex_shrink_count",
            "optimizer_initial_candidate_count",
            "optimizer_best_initial_scale",
            "optimizer_hit_lower_parameter_bound",
            "optimizer_hit_upper_parameter_bound",
            "overparameterized",
            "baseline_model",
            "baseline_parameter_count",
            "baseline_log_likelihood",
            "baseline_aic",
            "delta_log_likelihood",
            "delta_aic",
            "preferred_model_by_aic",
        ],
        rows=[
            {
                "trait": report.trait,
                "taxon_column": report.taxon_column,
                "model": report.model,
                "ascertainment_policy": report.ascertainment_policy,
                "ascertainment_conditioning_log_probability": (
                    ""
                    if report.ascertainment_conditioning_log_probability is None
                    else format(
                        report.ascertainment_conditioning_log_probability,
                        ".15g",
                    )
                ),
                "invariant_pattern_log_probability": (
                    ""
                    if report.invariant_pattern_log_probability is None
                    else format(report.invariant_pattern_log_probability, ".15g")
                ),
                "transform": ""
                if transform_fit is None
                else transform_fit.transform_name,
                "transform_parameter_name": (
                    "" if transform_fit is None else transform_fit.parameter_name
                ),
                "transform_parameter_value": (
                    ""
                    if transform_fit is None
                    else format(transform_fit.parameter_value, ".15g")
                ),
                "transform_lower_bound": (
                    ""
                    if transform_fit is None
                    else format(transform_fit.lower_bound, ".15g")
                ),
                "transform_upper_bound": (
                    ""
                    if transform_fit is None
                    else format(transform_fit.upper_bound, ".15g")
                ),
                "transform_starting_parameter_policy": (
                    ""
                    if transform_fit is None
                    else transform_fit.starting_parameter_policy
                ),
                "transform_starting_parameter_value": (
                    ""
                    if transform_fit is None
                    else format(transform_fit.starting_parameter_value, ".15g")
                ),
                "transform_starting_parameter_log_likelihood": (
                    ""
                    if transform_fit is None
                    else format(transform_fit.starting_parameter_log_likelihood, ".15g")
                ),
                "transform_coarse_grid_point_count": (
                    ""
                    if transform_fit is None
                    else str(transform_fit.coarse_grid_point_count)
                ),
                "transform_fine_grid_point_count": (
                    ""
                    if transform_fit is None
                    else str(transform_fit.fine_grid_point_count)
                ),
                "transform_refinement_start_count": (
                    ""
                    if transform_fit is None
                    else str(transform_fit.refinement_start_count)
                ),
                "transform_function_evaluation_count": (
                    ""
                    if transform_fit is None
                    else str(transform_fit.function_evaluation_count)
                ),
                "transform_hit_lower_parameter_boundary": (
                    ""
                    if transform_fit is None
                    else str(transform_fit.hit_lower_parameter_boundary).lower()
                ),
                "transform_hit_upper_parameter_boundary": (
                    ""
                    if transform_fit is None
                    else str(transform_fit.hit_upper_parameter_boundary).lower()
                ),
                "transform_warning_count": (
                    "" if transform_fit is None else str(len(transform_fit.warnings))
                ),
                "transform_tree_is_ultrametric": (
                    ""
                    if transform_fit is None
                    else str(transform_fit.transformed_tree_is_ultrametric).lower()
                ),
                "transform_tree_minimum_tip_depth": (
                    ""
                    if transform_fit is None
                    else format(
                        transform_fit.transformed_tree_minimum_tip_depth, ".15g"
                    )
                ),
                "transform_tree_maximum_tip_depth": (
                    ""
                    if transform_fit is None
                    else format(
                        transform_fit.transformed_tree_maximum_tip_depth, ".15g"
                    )
                ),
                "transform_baseline": (
                    ""
                    if transform_baseline is None
                    else transform_baseline.baseline_transform
                ),
                "transform_baseline_parameter_count": ""
                if transform_baseline is None
                else str(transform_baseline.baseline_parameter_count),
                "transform_baseline_log_likelihood": ""
                if transform_baseline is None
                else format(transform_baseline.baseline_log_likelihood, ".15g"),
                "transform_baseline_aic": ""
                if transform_baseline is None
                else format(transform_baseline.baseline_aic, ".15g"),
                "transform_delta_log_likelihood": ""
                if transform_baseline is None
                else format(transform_baseline.delta_log_likelihood, ".15g"),
                "transform_delta_aic": ""
                if transform_baseline is None
                else format(transform_baseline.delta_aic, ".15g"),
                "preferred_transform_by_aic": ""
                if transform_baseline is None
                else transform_baseline.preferred_transform_by_aic,
                "state_ordering": report.state_ordering,
                "analyzed_taxon_count": str(report.taxon_count),
                "excluded_taxon_count": str(
                    len(report.input_audit.pruned_missing_value_taxa)
                ),
                "observed_state_count": str(len(report.input_audit.observed_states)),
                "sparse_state_count": str(len(report.input_audit.sparse_states)),
                "log_likelihood": format(report.log_likelihood, ".15g"),
                "parameter_count": str(report.parameter_count),
                "aic": format(report.aic, ".15g"),
                "aicc": "inf"
                if math.isinf(report.aicc)
                else format(report.aicc, ".15g"),
                "likelihood_constant_policy": report.likelihood_constant_policy,
                "likelihood_comparison_policy": report.likelihood_comparison_policy,
                "optimizer_name": diagnostics.optimizer_name,
                "optimizer_converged": str(diagnostics.converged).lower(),
                "optimizer_iteration_count": str(diagnostics.iteration_count),
                "optimizer_function_evaluation_count": str(
                    diagnostics.function_evaluation_count
                ),
                "optimizer_simplex_shrink_count": str(diagnostics.simplex_shrink_count),
                "optimizer_initial_candidate_count": str(
                    diagnostics.initial_candidate_count
                ),
                "optimizer_best_initial_scale": format(
                    diagnostics.best_initial_scale,
                    ".15g",
                ),
                "optimizer_hit_lower_parameter_bound": str(
                    diagnostics.hit_lower_parameter_bound
                ).lower(),
                "optimizer_hit_upper_parameter_bound": str(
                    diagnostics.hit_upper_parameter_bound
                ).lower(),
                "overparameterized": str(report.overparameterized).lower(),
                "baseline_model": "" if baseline is None else baseline.baseline_model,
                "baseline_parameter_count": ""
                if baseline is None
                else str(baseline.baseline_parameter_count),
                "baseline_log_likelihood": ""
                if baseline is None
                else format(baseline.baseline_log_likelihood, ".15g"),
                "baseline_aic": ""
                if baseline is None
                else format(baseline.baseline_aic, ".15g"),
                "delta_log_likelihood": ""
                if baseline is None
                else format(baseline.delta_log_likelihood, ".15g"),
                "delta_aic": ""
                if baseline is None
                else format(baseline.delta_aic, ".15g"),
                "preferred_model_by_aic": (
                    "" if baseline is None else baseline.preferred_model_by_aic
                ),
            }
        ],
    )


def write_discrete_mk_rate_table(path: Path, report: DiscreteMkFitReport) -> Path:
    """Write one directed rate-matrix ledger for a discrete Mk fit."""
    return write_ancestral_rows(
        path,
        columns=[
            "source_state",
            "target_state",
            "transition_allowed",
            "step_distance",
            "rate",
        ],
        rows=[
            {
                "source_state": row.source_state,
                "target_state": row.target_state,
                "transition_allowed": str(row.transition_allowed).lower(),
                "step_distance": str(row.step_distance),
                "rate": format(row.rate, ".15g"),
            }
            for row in report.transition_rate_rows
        ],
    )


def write_discrete_mk_pattern_likelihood_table(
    path: Path,
    report: DiscreteMkFitReport,
) -> Path:
    """Write the observed Mk trait-pattern likelihood rows used to reconstruct totals."""
    return write_ancestral_rows(
        path,
        columns=[
            "pattern_id",
            "pattern_weight",
            "tip_states",
            "raw_log_likelihood",
            "ascertainment_conditioning_log_probability",
            "log_likelihood",
        ],
        rows=[
            {
                "pattern_id": row.pattern_id,
                "pattern_weight": str(row.pattern_weight),
                "tip_states": "|".join(row.tip_states),
                "raw_log_likelihood": format(row.raw_log_likelihood, ".15g"),
                "ascertainment_conditioning_log_probability": (
                    ""
                    if row.ascertainment_conditioning_log_probability is None
                    else format(
                        row.ascertainment_conditioning_log_probability,
                        ".15g",
                    )
                ),
                "log_likelihood": format(row.log_likelihood, ".15g"),
            }
            for row in report.pattern_likelihood_rows
        ],
    )
