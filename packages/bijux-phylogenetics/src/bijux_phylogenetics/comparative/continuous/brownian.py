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
)
from bijux_phylogenetics.comparative.evolutionary_modes import (
    fit_continuous_evolutionary_mode,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows


@dataclass(slots=True)
class BrownianTraitEvolutionExclusion:
    """One taxon excluded before Brownian trait-evolution fitting."""

    taxon: str
    reason: str


@dataclass(slots=True)
class BrownianTraitEvolutionSummaryReport:
    """Reviewer-facing Brownian trait-evolution fit with explicit pruning context."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    tree_taxon_count: int
    analyzed_taxa: list[str]
    analyzed_taxon_count: int
    excluded_taxa: list[BrownianTraitEvolutionExclusion]
    root_state: float
    sigma_squared: float
    log_likelihood: float
    aic: float
    aicc: float
    confidence_intervals: list[ComparativeParameterInterval]
    residual_diagnostics: ComparativeResidualSummary
    assumptions: list[str]
    warnings: list[str]
    readiness: ComparativeReadinessReport


def summarize_brownian_trait_evolution(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> BrownianTraitEvolutionSummaryReport:
    """Summarize one Brownian continuous-trait fit with explicit exclusion tracking."""
    readiness = summarize_numeric_trait_readiness(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    fit = fit_continuous_evolutionary_mode(
        tree_path,
        traits_path,
        trait=trait,
        mode="brownian",
        taxon_column=taxon_column,
    )
    excluded_taxa = _build_excluded_taxa(readiness)
    warnings = list(
        dict.fromkeys([*readiness.warnings, *fit.residual_diagnostics.warnings])
    )
    return BrownianTraitEvolutionSummaryReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=fit.taxon_column,
        trait=trait,
        tree_taxon_count=readiness.tree_taxa,
        analyzed_taxa=list(fit.taxa),
        analyzed_taxon_count=fit.taxon_count,
        excluded_taxa=excluded_taxa,
        root_state=fit.root_state,
        sigma_squared=fit.rate,
        log_likelihood=fit.log_likelihood,
        aic=fit.aic,
        aicc=fit.aicc,
        confidence_intervals=list(fit.confidence_intervals),
        residual_diagnostics=fit.residual_diagnostics,
        assumptions=list(fit.assumptions),
        warnings=warnings,
        readiness=readiness,
    )


def write_brownian_trait_evolution_summary_table(
    path: Path, report: BrownianTraitEvolutionSummaryReport
) -> Path:
    """Write one summary ledger for a Brownian trait-evolution fit."""
    interval_by_name = {
        interval.name: interval for interval in report.confidence_intervals
    }
    root_state_interval = interval_by_name.get("root_state")
    sigma_interval = interval_by_name.get("rate")
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "tree_taxon_count",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "root_state",
            "root_state_lower_95",
            "root_state_upper_95",
            "sigma_squared",
            "sigma_squared_lower_95",
            "sigma_squared_upper_95",
            "log_likelihood",
            "aic",
            "aicc",
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
                "root_state": format(report.root_state, ".15g"),
                "root_state_lower_95": _format_interval_bound(
                    root_state_interval, "lower_95"
                ),
                "root_state_upper_95": _format_interval_bound(
                    root_state_interval, "upper_95"
                ),
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


def write_brownian_trait_evolution_exclusion_table(
    path: Path, report: BrownianTraitEvolutionSummaryReport
) -> Path:
    """Write one excluded-taxon ledger for a Brownian trait-evolution fit."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
            }
            for row in report.excluded_taxa
        ],
    )


def _build_excluded_taxa(
    readiness: ComparativeReadinessReport,
) -> list[BrownianTraitEvolutionExclusion]:
    rows: list[BrownianTraitEvolutionExclusion] = []
    rows.extend(
        BrownianTraitEvolutionExclusion(
            taxon=taxon,
            reason="missing_from_trait_table",
        )
        for taxon in readiness.missing_from_traits
    )
    rows.extend(
        BrownianTraitEvolutionExclusion(
            taxon=taxon,
            reason="missing_trait_value",
        )
        for taxon in readiness.pruned_missing_value_taxa
    )
    rows.extend(
        BrownianTraitEvolutionExclusion(
            taxon=taxon,
            reason="non_numeric_trait_value",
        )
        for taxon in readiness.pruned_non_numeric_taxa
    )
    rows.extend(
        BrownianTraitEvolutionExclusion(
            taxon=taxon,
            reason="absent_from_tree",
        )
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
