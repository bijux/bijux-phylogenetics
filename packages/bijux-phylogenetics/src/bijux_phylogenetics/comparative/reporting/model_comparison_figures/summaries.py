from __future__ import annotations

import math

from ...continuous import BrownianMotionFitReport, OUTraitModelReport
from ...model_selection import ComparativeModelComparisonReport
from .contracts import (
    ComparativeModelCriteriaRow,
    ComparativeModelFigureAudit,
    ComparativeModelFigureCaptionDraft,
    ComparativeModelFigureLegendEntry,
    ComparativeModelFitRow,
    ComparativeModelLikelihoodRow,
    ComparativeModelParameterRow,
)


def format_model_figure_number(value: float | None) -> str:
    if value is None:
        return "n/a"
    if math.isnan(value):
        return "nan"
    if math.isinf(value):
        return "inf" if value > 0.0 else "-inf"
    return format(value, ".15g")


def delta_from_best(value: float, best: float) -> float:
    if math.isfinite(value) and math.isfinite(best):
        return value - best
    if math.isinf(value) and math.isinf(best):
        return 0.0
    if math.isinf(value):
        return math.inf
    if math.isinf(best):
        return -math.inf
    return math.nan


def model_figure_color(model: str, *, selected_model: str) -> str:
    if model == selected_model:
        return "#0f766e"
    if model == "brownian":
        return "#0369a1"
    return "#9333ea"


def parameter_label(parameter: str) -> str:
    return {
        "root_state": "root state",
        "rate": "brownian rate",
        "alpha": "ou alpha",
        "theta": "ou theta",
        "sigma_squared": "ou sigma^2",
    }.get(parameter, parameter.replace("_", " "))


def build_model_criteria_rows(
    report: ComparativeModelComparisonReport,
) -> list[ComparativeModelCriteriaRow]:
    best_aic = min(row.aic for row in report.rows)
    best_aicc = min(row.aicc for row in report.rows)
    return [
        ComparativeModelCriteriaRow(
            model=row.model,
            parameter_count=row.parameter_count,
            log_likelihood=row.log_likelihood,
            aic=row.aic,
            aicc=row.aicc,
            delta_aic=delta_from_best(row.aic, best_aic),
            delta_aicc=delta_from_best(row.aicc, best_aicc),
            selected=row.selected,
        )
        for row in report.rows
    ]


def build_model_likelihood_rows(
    report: ComparativeModelComparisonReport,
) -> list[ComparativeModelLikelihoodRow]:
    best_log_likelihood = max(row.log_likelihood for row in report.rows)
    return [
        ComparativeModelLikelihoodRow(
            model=row.model,
            log_likelihood=row.log_likelihood,
            delta_log_likelihood=best_log_likelihood - row.log_likelihood,
            selected=row.selected,
        )
        for row in report.rows
    ]


def build_model_parameter_rows(
    brownian: BrownianMotionFitReport,
    ou: OUTraitModelReport,
) -> list[ComparativeModelParameterRow]:
    rows: list[ComparativeModelParameterRow] = []
    for interval in brownian.confidence_intervals:
        rows.append(
            ComparativeModelParameterRow(
                model="brownian",
                parameter=interval.name,
                estimate=interval.estimate,
                lower_95=interval.lower_95,
                upper_95=interval.upper_95,
                interval_method=interval.method,
            )
        )
    for interval in ou.confidence_intervals:
        rows.append(
            ComparativeModelParameterRow(
                model="ou",
                parameter=interval.name,
                estimate=interval.estimate,
                lower_95=interval.lower_95,
                upper_95=interval.upper_95,
                interval_method=interval.method,
            )
        )
    return rows


def build_model_fit_rows(
    brownian: BrownianMotionFitReport,
    ou: OUTraitModelReport,
    *,
    selected_model: str,
) -> list[ComparativeModelFitRow]:
    return [
        ComparativeModelFitRow(
            model="brownian",
            taxon_count=brownian.taxon_count,
            residual_variance=brownian.residual_diagnostics.residual_variance,
            max_abs_standardized_residual=(
                brownian.residual_diagnostics.max_abs_standardized_residual
            ),
            phylogenetic_residual_lambda=(
                brownian.residual_diagnostics.phylogenetic_residual_lambda
            ),
            outlier_taxon_count=len(brownian.residual_diagnostics.outlier_taxa),
            warning_count=len(brownian.residual_diagnostics.warnings),
            convergence_status="analytic",
            selected=selected_model == "brownian",
        ),
        ComparativeModelFitRow(
            model="ou",
            taxon_count=ou.taxon_count,
            residual_variance=ou.residual_diagnostics.residual_variance,
            max_abs_standardized_residual=(
                ou.residual_diagnostics.max_abs_standardized_residual
            ),
            phylogenetic_residual_lambda=(
                ou.residual_diagnostics.phylogenetic_residual_lambda
            ),
            outlier_taxon_count=len(ou.residual_diagnostics.outlier_taxa),
            warning_count=len(ou.residual_diagnostics.warnings)
            + len(ou.identifiability_warnings),
            convergence_status=ou.convergence_status,
            selected=selected_model == "ou",
        ),
    ]


def build_model_figure_legend() -> list[ComparativeModelFigureLegendEntry]:
    return [
        ComparativeModelFigureLegendEntry(
            surface="information-criteria",
            label="selected model",
            swatch="#0f766e",
            detail="the highlighted row marks the model selected by AICc when the candidate models retain finite information-criterion values",
        ),
        ComparativeModelFigureLegendEntry(
            surface="likelihood",
            label="log-likelihood bar",
            swatch="#0369a1",
            detail="bar length tracks relative log-likelihood so the gap between candidate models remains visible instead of being hidden inside a summary sentence",
        ),
        ComparativeModelFigureLegendEntry(
            surface="parameters",
            label="95% parameter interval",
            swatch="#9333ea",
            detail="horizontal whiskers retain the estimated parameter value with its fitted 95% interval for each continuous-trait model parameter",
        ),
        ComparativeModelFigureLegendEntry(
            surface="fit-summary",
            label="diagnostic row",
            swatch="#475569",
            detail="the fit-summary surface keeps residual scale, standardized outlier burden, phylogenetic residual lambda, and warning counts visible for both models",
        ),
    ]


def build_model_figure_audit(
    *,
    comparison_report: ComparativeModelComparisonReport,
    brownian: BrownianMotionFitReport,
    ou: OUTraitModelReport,
    criteria_rows: list[ComparativeModelCriteriaRow],
    likelihood_rows: list[ComparativeModelLikelihoodRow],
    parameter_rows: list[ComparativeModelParameterRow],
    fit_rows: list[ComparativeModelFitRow],
    legend_entries: list[ComparativeModelFigureLegendEntry],
) -> ComparativeModelFigureAudit:
    criteria_surface_visible = len(criteria_rows) == len(comparison_report.rows)
    likelihood_surface_visible = len(likelihood_rows) == len(comparison_report.rows)
    parameter_surface_visible = len(parameter_rows) >= 5
    fit_surface_visible = len(fit_rows) == len(comparison_report.rows)
    legend_complete = {entry.surface for entry in legend_entries} == {
        "information-criteria",
        "likelihood",
        "parameters",
        "fit-summary",
    }
    caption_ready = (
        criteria_surface_visible
        and likelihood_surface_visible
        and parameter_surface_visible
        and fit_surface_visible
    )
    finite_rows = [row for row in criteria_rows if math.isfinite(row.aicc)]
    selected_row = next(
        row for row in criteria_rows if row.model == comparison_report.better_model
    )
    runner_up_candidates = [
        row.aicc for row in criteria_rows if row.model != comparison_report.better_model
    ]
    runner_up_aicc = min(runner_up_candidates) if runner_up_candidates else math.nan
    aicc_delta = (
        runner_up_aicc - selected_row.aicc
        if math.isfinite(runner_up_aicc) and math.isfinite(selected_row.aicc)
        else None
    )
    support_distinct = aicc_delta is not None and aicc_delta >= 2.0
    warning_count = (
        len(brownian.residual_diagnostics.warnings)
        + len(ou.residual_diagnostics.warnings)
        + len(ou.identifiability_warnings)
    )
    publication_ready = (
        criteria_surface_visible
        and likelihood_surface_visible
        and parameter_surface_visible
        and fit_surface_visible
        and legend_complete
        and caption_ready
        and len(finite_rows) == len(criteria_rows)
        and support_distinct
    )
    reviewer_summary = [
        f"criteria rows rendered: {len(criteria_rows)}/{len(comparison_report.rows)}",
        f"likelihood rows rendered: {len(likelihood_rows)}/{len(comparison_report.rows)}",
        f"parameter rows rendered: {len(parameter_rows)}",
        f"fit rows rendered: {len(fit_rows)}/{len(comparison_report.rows)}",
        f"selected model: {comparison_report.better_model}",
        f"aicc delta to runner-up: {format_model_figure_number(aicc_delta)}",
    ]
    limitations: list[str] = []
    if len(finite_rows) != len(criteria_rows):
        limitations.append(
            "one or more candidate models do not retain finite AICc values, so the model-comparison package stays reviewable but not publication-ready"
        )
    if not support_distinct:
        limitations.append(
            "the AICc separation between the selected model and the runner-up stays below the publication threshold of 2.0, so process support remains ambiguous"
        )
    if len(ou.identifiability_warnings) > 0:
        limitations.append(
            "the OU fit retains one or more identifiability warnings, so alpha- and theta-level interpretation remains fragile"
        )
    return ComparativeModelFigureAudit(
        publication_ready=publication_ready,
        criteria_surface_visible=criteria_surface_visible,
        likelihood_surface_visible=likelihood_surface_visible,
        parameter_surface_visible=parameter_surface_visible,
        fit_surface_visible=fit_surface_visible,
        legend_complete=legend_complete,
        caption_ready=caption_ready,
        finite_aicc_model_count=len(finite_rows),
        support_distinct=support_distinct,
        selected_model=comparison_report.better_model,
        aicc_delta=aicc_delta,
        plotted_model_count=len(comparison_report.rows),
        rendered_parameter_count=len(parameter_rows),
        rendered_fit_row_count=len(fit_rows),
        warning_count=warning_count,
        reviewer_summary=reviewer_summary,
        limitations=limitations,
    )


def build_model_figure_caption_draft(
    *,
    comparison_report: ComparativeModelComparisonReport,
    audit: ComparativeModelFigureAudit,
    fit_rows: list[ComparativeModelFitRow],
) -> ComparativeModelFigureCaptionDraft:
    selected_fit = next(row for row in fit_rows if row.model == audit.selected_model)
    return ComparativeModelFigureCaptionDraft(
        title="Comparative continuous-trait model review across information criteria, likelihood, parameters, and fit diagnostics",
        lead_sentence=(
            f"This package compares Brownian-motion and Ornstein-Uhlenbeck fits for `{comparison_report.trait}` across {comparison_report.taxon_count} overlapping taxa with reviewer-facing figure surfaces instead of leaving model choice buried in one JSON payload."
        ),
        criteria_sentence=(
            f"The information-criterion surface retains AIC and AICc for both candidate models and currently selects `{audit.selected_model}` with an AICc delta of {format_model_figure_number(audit.aicc_delta)} relative to the runner-up."
        ),
        likelihood_sentence=(
            "The likelihood surface keeps the relative log-likelihood gap explicit so likelihood fit can be reviewed separately from the information-criterion penalty."
        ),
        parameter_sentence=(
            f"The parameter surface renders {audit.rendered_parameter_count} fitted BM and OU parameter intervals, including Brownian root-state and rate estimates plus OU alpha, theta, and sigma-squared."
        ),
        fit_sentence=(
            f"The fit-summary surface keeps residual variance, standardized outlier burden, phylogenetic residual lambda, and warning counts visible for the selected `{audit.selected_model}` model, which currently carries {selected_fit.warning_count} warning rows."
        ),
        limitation_sentence=audit.limitations[0],
        caption_ready=audit.caption_ready,
    )
