from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative.common import (
    ComparativeReadinessReport,
    summarize_numeric_trait_readiness,
)
from bijux_phylogenetics.comparative.continuous.model_fitting import (
    ComparativeParameterInterval,
    ComparativeResidualSummary,
    OUIdentifiabilityWarning,
    fit_ornstein_uhlenbeck_model,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows


@dataclass(slots=True)
class OUTraitEvolutionExclusion:
    """One taxon excluded before OU trait-evolution fitting."""

    taxon: str
    reason: str


@dataclass(slots=True)
class OUTraitEvolutionSummaryReport:
    """Reviewer-facing OU trait-evolution fit with explicit pruning context."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    tree_taxon_count: int
    analyzed_taxa: list[str]
    analyzed_taxon_count: int
    excluded_taxa: list[OUTraitEvolutionExclusion]
    alpha: float
    theta: float
    sigma_squared: float
    log_likelihood: float
    aic: float
    aicc: float
    convergence_status: str
    confidence_intervals: list[ComparativeParameterInterval]
    identifiability_warnings: list[OUIdentifiabilityWarning]
    residual_diagnostics: ComparativeResidualSummary
    assumptions: list[str]
    warnings: list[str]
    readiness: ComparativeReadinessReport


def summarize_ou_trait_evolution(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> OUTraitEvolutionSummaryReport:
    """Summarize one OU continuous-trait fit with explicit exclusion tracking."""
    readiness = summarize_numeric_trait_readiness(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    fit = fit_ornstein_uhlenbeck_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    excluded_taxa = _build_excluded_taxa(readiness)
    parameter_count = 3
    sample_size = fit.taxon_count
    aic = (2.0 * parameter_count) - (2.0 * fit.log_likelihood)
    if sample_size <= parameter_count + 1:
        aicc = float("inf")
    else:
        aicc = aic + (
            (2.0 * parameter_count * (parameter_count + 1))
            / (sample_size - parameter_count - 1)
        )
    warnings = list(
        dict.fromkeys(
            [
                *readiness.warnings,
                *fit.residual_diagnostics.warnings,
                *[warning.message for warning in fit.identifiability_warnings],
            ]
        )
    )
    return OUTraitEvolutionSummaryReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=fit.taxon_column,
        trait=trait,
        tree_taxon_count=readiness.tree_taxa,
        analyzed_taxa=list(fit.taxa),
        analyzed_taxon_count=fit.taxon_count,
        excluded_taxa=excluded_taxa,
        alpha=fit.alpha,
        theta=fit.theta,
        sigma_squared=fit.sigma_squared,
        log_likelihood=fit.log_likelihood,
        aic=aic,
        aicc=aicc,
        convergence_status=fit.convergence_status,
        confidence_intervals=list(fit.confidence_intervals),
        identifiability_warnings=list(fit.identifiability_warnings),
        residual_diagnostics=fit.residual_diagnostics,
        assumptions=list(fit.assumptions),
        warnings=warnings,
        readiness=readiness,
    )


def write_ou_trait_evolution_summary_table(
    path: Path, report: OUTraitEvolutionSummaryReport
) -> Path:
    """Write one summary ledger for an OU trait-evolution fit."""
    interval_by_name = {
        interval.name: interval for interval in report.confidence_intervals
    }
    alpha_interval = interval_by_name.get("alpha")
    theta_interval = interval_by_name.get("theta")
    sigma_interval = interval_by_name.get("sigma_squared")
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "tree_taxon_count",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "alpha",
            "alpha_lower_95",
            "alpha_upper_95",
            "theta",
            "theta_lower_95",
            "theta_upper_95",
            "sigma_squared",
            "sigma_squared_lower_95",
            "sigma_squared_upper_95",
            "log_likelihood",
            "aic",
            "aicc",
            "convergence_status",
            "identifiability_warning_count",
            "residual_variance",
            "max_abs_standardized_residual",
            "phylogenetic_residual_lambda",
        ],
        rows=[
            {
                "trait": report.trait,
                "taxon_column": report.taxon_column,
                "tree_taxon_count": report.tree_taxon_count,
                "analyzed_taxon_count": report.analyzed_taxon_count,
                "excluded_taxon_count": len(report.excluded_taxa),
                "alpha": format(report.alpha, ".15g"),
                "alpha_lower_95": _format_interval_bound(alpha_interval, "lower_95"),
                "alpha_upper_95": _format_interval_bound(alpha_interval, "upper_95"),
                "theta": format(report.theta, ".15g"),
                "theta_lower_95": _format_interval_bound(theta_interval, "lower_95"),
                "theta_upper_95": _format_interval_bound(theta_interval, "upper_95"),
                "sigma_squared": format(report.sigma_squared, ".15g"),
                "sigma_squared_lower_95": _format_interval_bound(
                    sigma_interval, "lower_95"
                ),
                "sigma_squared_upper_95": _format_interval_bound(
                    sigma_interval, "upper_95"
                ),
                "log_likelihood": format(report.log_likelihood, ".15g"),
                "aic": format(report.aic, ".15g"),
                "aicc": format(report.aicc, ".15g"),
                "convergence_status": report.convergence_status,
                "identifiability_warning_count": len(report.identifiability_warnings),
                "residual_variance": format(
                    report.residual_diagnostics.residual_variance, ".15g"
                ),
                "max_abs_standardized_residual": format(
                    report.residual_diagnostics.max_abs_standardized_residual, ".15g"
                ),
                "phylogenetic_residual_lambda": format(
                    report.residual_diagnostics.phylogenetic_residual_lambda, ".15g"
                ),
            }
        ],
    )


def write_ou_trait_evolution_exclusion_table(
    path: Path, report: OUTraitEvolutionSummaryReport
) -> Path:
    """Write one excluded-taxon ledger for an OU trait-evolution fit."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {"taxon": row.taxon, "reason": row.reason} for row in report.excluded_taxa
        ],
    )


def _build_excluded_taxa(
    readiness: ComparativeReadinessReport,
) -> list[OUTraitEvolutionExclusion]:
    rows: list[OUTraitEvolutionExclusion] = []
    rows.extend(
        OUTraitEvolutionExclusion(taxon=taxon, reason="missing_from_trait_table")
        for taxon in readiness.missing_from_traits
    )
    rows.extend(
        OUTraitEvolutionExclusion(taxon=taxon, reason="missing_trait_value")
        for taxon in readiness.pruned_missing_value_taxa
    )
    rows.extend(
        OUTraitEvolutionExclusion(taxon=taxon, reason="non_numeric_trait_value")
        for taxon in readiness.pruned_non_numeric_taxa
    )
    rows.extend(
        OUTraitEvolutionExclusion(taxon=taxon, reason="absent_from_tree")
        for taxon in readiness.extra_trait_taxa
    )
    return rows


def _format_interval_bound(
    interval: ComparativeParameterInterval | None,
    field_name: str,
) -> str:
    if interval is None:
        return ""
    return format(getattr(interval, field_name), ".15g")
