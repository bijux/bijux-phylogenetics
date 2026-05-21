from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.common import write_ancestral_rows

from .models import (
    ContinuousAncestralExclusion,
    ContinuousAncestralReport,
    ContinuousAncestralSummary,
)


def summarize_continuous_ancestral_report(
    report: ContinuousAncestralReport,
) -> ContinuousAncestralSummary:
    """Summarize the main review facts for one continuous ancestral report."""
    internal_estimates = [
        estimate for estimate in report.estimates if not estimate.is_tip
    ]
    if not internal_estimates:
        raise ValueError(
            "continuous ancestral summary requires at least one internal-node estimate"
        )
    root_estimate = max(
        internal_estimates,
        key=lambda estimate: (
            len(estimate.descendant_taxa),
            estimate.node,
        ),
    )
    return ContinuousAncestralSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        estimator=report.estimator,
        alpha=report.alpha,
        analyzed_taxon_count=report.taxon_count,
        excluded_taxon_count=len(report.missing_from_traits_taxa)
        + len(report.dropped_missing_taxa)
        + len(report.dropped_non_numeric_taxa),
        missing_tip_taxon_count=len(report.missing_from_traits_taxa)
        + len(report.dropped_missing_taxa),
        non_numeric_tip_taxon_count=len(report.dropped_non_numeric_taxa),
        internal_node_count=len(internal_estimates),
        unstable_node_count=len(report.unstable_nodes),
        weak_support_node_count=len(report.weak_support_nodes),
        root_node=root_estimate.node,
        root_estimate=root_estimate.estimate,
        root_standard_error=root_estimate.standard_error,
        root_lower_95_interval=root_estimate.lower_95_interval,
        root_upper_95_interval=root_estimate.upper_95_interval,
        tree_is_ultrametric=(
            None
            if report.brownian_fit_diagnostics is None
            else report.brownian_fit_diagnostics.tree_is_ultrametric
        ),
        covariance_near_singular=(
            None
            if report.brownian_fit_diagnostics is None
            else report.brownian_fit_diagnostics.covariance_near_singular
        ),
        covariance_condition_number=(
            None
            if report.brownian_fit_diagnostics is None
            else report.brownian_fit_diagnostics.covariance_condition_number
        ),
        log_likelihood=(
            None
            if report.brownian_fit_diagnostics is None
            else report.brownian_fit_diagnostics.log_likelihood
        ),
        residual_sigma_squared=(
            None
            if report.brownian_fit_diagnostics is None
            else report.brownian_fit_diagnostics.residual_sigma_squared
        ),
        optimizer_name=(
            None
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.optimizer_name
        ),
        optimizer_converged=(
            None
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.converged
        ),
        optimizer_iteration_count=(
            None
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.iteration_count
        ),
        optimizer_function_evaluation_count=(
            None
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.function_evaluation_count
        ),
        optimizer_convergence_status=(
            None
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.convergence_status
        ),
        warning_count=len(report.warnings),
    )


def continuous_ancestral_exclusions(
    report: ContinuousAncestralReport,
) -> list[ContinuousAncestralExclusion]:
    """Return one explicit exclusion row per dropped tip taxon."""
    rows = [
        ContinuousAncestralExclusion(
            taxon=taxon,
            reason="missing_trait_value",
        )
        for taxon in report.missing_from_traits_taxa + report.dropped_missing_taxa
    ]
    rows.extend(
        ContinuousAncestralExclusion(
            taxon=taxon,
            reason="non_numeric_trait_value",
        )
        for taxon in report.dropped_non_numeric_taxa
    )
    return rows


def write_continuous_ancestral_summary_table(
    path: Path, report: ContinuousAncestralReport
) -> Path:
    """Write one summary ledger for a continuous ancestral reconstruction."""
    summary = summarize_continuous_ancestral_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "estimator",
            "alpha",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "missing_tip_taxon_count",
            "non_numeric_tip_taxon_count",
            "internal_node_count",
            "unstable_node_count",
            "weak_support_node_count",
            "root_node",
            "root_estimate",
            "root_standard_error",
            "root_lower_95_interval",
            "root_upper_95_interval",
            "tree_is_ultrametric",
            "covariance_near_singular",
            "covariance_condition_number",
            "log_likelihood",
            "residual_sigma_squared",
            "optimizer_name",
            "optimizer_converged",
            "optimizer_iteration_count",
            "optimizer_function_evaluation_count",
            "optimizer_convergence_status",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "estimator": summary.estimator,
                "alpha": str(summary.alpha),
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "missing_tip_taxon_count": str(summary.missing_tip_taxon_count),
                "non_numeric_tip_taxon_count": str(summary.non_numeric_tip_taxon_count),
                "internal_node_count": str(summary.internal_node_count),
                "unstable_node_count": str(summary.unstable_node_count),
                "weak_support_node_count": str(summary.weak_support_node_count),
                "root_node": summary.root_node,
                "root_estimate": str(summary.root_estimate),
                "root_standard_error": str(summary.root_standard_error),
                "root_lower_95_interval": str(summary.root_lower_95_interval),
                "root_upper_95_interval": str(summary.root_upper_95_interval),
                "tree_is_ultrametric": (
                    ""
                    if summary.tree_is_ultrametric is None
                    else str(summary.tree_is_ultrametric).lower()
                ),
                "covariance_near_singular": (
                    ""
                    if summary.covariance_near_singular is None
                    else str(summary.covariance_near_singular).lower()
                ),
                "covariance_condition_number": (
                    ""
                    if summary.covariance_condition_number is None
                    else str(summary.covariance_condition_number)
                ),
                "log_likelihood": (
                    ""
                    if summary.log_likelihood is None
                    else str(summary.log_likelihood)
                ),
                "residual_sigma_squared": (
                    ""
                    if summary.residual_sigma_squared is None
                    else str(summary.residual_sigma_squared)
                ),
                "optimizer_name": ""
                if summary.optimizer_name is None
                else summary.optimizer_name,
                "optimizer_converged": (
                    ""
                    if summary.optimizer_converged is None
                    else str(summary.optimizer_converged).lower()
                ),
                "optimizer_iteration_count": (
                    ""
                    if summary.optimizer_iteration_count is None
                    else str(summary.optimizer_iteration_count)
                ),
                "optimizer_function_evaluation_count": (
                    ""
                    if summary.optimizer_function_evaluation_count is None
                    else str(summary.optimizer_function_evaluation_count)
                ),
                "optimizer_convergence_status": (
                    ""
                    if summary.optimizer_convergence_status is None
                    else summary.optimizer_convergence_status
                ),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_continuous_ancestral_uncertainty_table(
    path: Path, report: ContinuousAncestralReport
) -> Path:
    """Write one internal-node uncertainty ledger for a continuous reconstruction."""
    return write_ancestral_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "estimate",
            "standard_error",
            "lower_95_interval",
            "upper_95_interval",
            "uncertainty_width",
            "confidence",
            "interpretation",
            "unstable",
        ],
        rows=[
            {
                "node": estimate.node,
                "node_name": estimate.node_name or "",
                "descendant_taxa": ",".join(estimate.descendant_taxa),
                "estimate": str(estimate.estimate),
                "standard_error": str(estimate.standard_error),
                "lower_95_interval": str(estimate.lower_95_interval),
                "upper_95_interval": str(estimate.upper_95_interval),
                "uncertainty_width": str(estimate.uncertainty_width),
                "confidence": str(estimate.confidence),
                "interpretation": estimate.interpretation,
                "unstable": str(estimate.unstable).lower(),
            }
            for estimate in report.estimates
            if not estimate.is_tip
        ],
    )


def write_continuous_ancestral_exclusion_table(
    path: Path, report: ContinuousAncestralReport
) -> Path:
    """Write one explicit excluded-tip ledger for a continuous reconstruction."""
    exclusions = continuous_ancestral_exclusions(report)
    return write_ancestral_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
            }
            for row in exclusions
        ],
    )
