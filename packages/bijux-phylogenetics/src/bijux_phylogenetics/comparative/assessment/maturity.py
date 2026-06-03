from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative.assessment.sensitivity import (
    run_comparative_sensitivity_analysis,
)
from bijux_phylogenetics.comparative.common import (
    ComparativeDataset,
    load_comparative_dataset,
)
from bijux_phylogenetics.comparative.continuous import (
    compare_brownian_and_ou_models,
    fit_brownian_motion_model,
    fit_ornstein_uhlenbeck_model,
)
from bijux_phylogenetics.comparative.continuous.model_fitting import (
    _estimate_lambda_for_values,
)
from bijux_phylogenetics.comparative.pgls import inspect_pgls_inputs, run_pgls
from bijux_phylogenetics.comparative.validation.reference_examples import (
    validate_comparative_reference_examples,
)


@dataclass(slots=True)
class ComparativeResidualDiagnosticSurface:
    """Reviewer-facing residual diagnostics for one comparative analysis surface."""

    analysis: str
    residual_variance: float
    max_abs_standardized_residual: float
    phylogenetic_residual_lambda: float | None
    max_leverage: float | None
    outlier_taxa: list[str]
    high_leverage_taxa: list[str]
    warnings: list[str]


@dataclass(slots=True)
class ComparativeSensitivitySummary:
    """Compact summary of leave-one-taxon-out comparative sensitivity."""

    model: str
    influential_taxa: list[str]
    max_abs_delta_log_likelihood: float
    max_abs_delta_primary_parameter: float


@dataclass(slots=True)
class ComparativeMethodMaturityReport:
    """Integrated comparative audit over one user-supplied response trait workflow."""

    tree_path: Path
    traits_path: Path
    trait: str
    selected_model: str
    reference_validation_passed: bool
    residual_diagnostics: list[ComparativeResidualDiagnosticSurface]
    sensitivity: ComparativeSensitivitySummary
    warnings: list[str]


def assess_comparative_method_maturity(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> ComparativeMethodMaturityReport:
    """Summarize residual and sensitivity trust signals for one comparative workflow."""
    pgls_inputs = inspect_pgls_inputs(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
    )
    brownian = fit_brownian_motion_model(
        tree_path, traits_path, trait=pgls_inputs.response, taxon_column=taxon_column
    )
    ou = fit_ornstein_uhlenbeck_model(
        tree_path, traits_path, trait=pgls_inputs.response, taxon_column=taxon_column
    )
    model_comparison = compare_brownian_and_ou_models(
        tree_path,
        traits_path,
        trait=pgls_inputs.response,
        taxon_column=taxon_column,
    )
    pgls = run_pgls(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    pgls_dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=pgls_inputs.response,
        taxon_column=taxon_column,
        minimum_taxa=max(3, len(pgls.taxa)),
        require_rooted=True,
        require_binary=False,
    )
    residual_lambda = _estimate_lambda_for_values(
        ComparativeDataset(
            tree_path=pgls_dataset.tree_path,
            traits_path=pgls_dataset.traits_path,
            tree=pgls_dataset.tree,
            taxon_column=pgls_dataset.taxon_column,
            trait=pgls_dataset.trait,
            taxa=pgls.taxa,
            trait_values=[0.0] * len(pgls.taxa),
            covariance_matrix=_subset_covariance(
                pgls_dataset.covariance_matrix, pgls_dataset.taxa, pgls.taxa
            ),
            readiness=pgls_dataset.readiness,
        ),
        pgls.residuals,
    )
    pgls_leverage_cutoff = (2.0 * len(pgls.encoded_columns)) / len(pgls.taxa)
    pgls_high_leverage = [
        row.taxon
        for row in pgls.diagnostics.leverage_rows
        if row.leverage >= pgls_leverage_cutoff
    ]
    residual_surfaces = [
        ComparativeResidualDiagnosticSurface(
            analysis="brownian",
            residual_variance=brownian.residual_diagnostics.residual_variance,
            max_abs_standardized_residual=brownian.residual_diagnostics.max_abs_standardized_residual,
            phylogenetic_residual_lambda=brownian.residual_diagnostics.phylogenetic_residual_lambda,
            max_leverage=None,
            outlier_taxa=[
                row.taxon for row in brownian.residual_diagnostics.outlier_taxa
            ],
            high_leverage_taxa=[],
            warnings=list(brownian.residual_diagnostics.warnings),
        ),
        ComparativeResidualDiagnosticSurface(
            analysis="ou",
            residual_variance=ou.residual_diagnostics.residual_variance,
            max_abs_standardized_residual=ou.residual_diagnostics.max_abs_standardized_residual,
            phylogenetic_residual_lambda=ou.residual_diagnostics.phylogenetic_residual_lambda,
            max_leverage=None,
            outlier_taxa=[row.taxon for row in ou.residual_diagnostics.outlier_taxa],
            high_leverage_taxa=[],
            warnings=[
                *ou.residual_diagnostics.warnings,
                *[warning.message for warning in ou.identifiability_warnings],
            ],
        ),
        ComparativeResidualDiagnosticSurface(
            analysis="pgls",
            residual_variance=pgls.residual_variance,
            max_abs_standardized_residual=max(
                abs(row.standardized_residual) for row in pgls.diagnostics.leverage_rows
            ),
            phylogenetic_residual_lambda=residual_lambda,
            max_leverage=max(row.leverage for row in pgls.diagnostics.leverage_rows),
            outlier_taxa=[row.taxon for row in pgls.diagnostics.outlier_taxa],
            high_leverage_taxa=pgls_high_leverage,
            warnings=_build_pgls_residual_warnings(
                residual_lambda,
                outlier_taxa=[row.taxon for row in pgls.diagnostics.outlier_taxa],
                high_leverage_taxa=pgls_high_leverage,
            ),
        ),
    ]
    sensitivity = run_comparative_sensitivity_analysis(
        tree_path,
        traits_path,
        trait=pgls_inputs.response,
        model=model_comparison.better_model,
        taxon_column=taxon_column,
    )
    sensitivity_summary = ComparativeSensitivitySummary(
        model=sensitivity.model,
        influential_taxa=sensitivity.most_influential_taxa,
        max_abs_delta_log_likelihood=max(
            abs(row.delta_log_likelihood) for row in sensitivity.rows
        ),
        max_abs_delta_primary_parameter=max(
            abs(row.delta_primary_parameter) for row in sensitivity.rows
        ),
    )
    warnings = sorted(
        {
            *brownian.residual_diagnostics.warnings,
            *ou.residual_diagnostics.warnings,
            *[warning.message for warning in ou.identifiability_warnings],
            *pgls_inputs.warnings,
            *[warning for surface in residual_surfaces for warning in surface.warnings],
        }
    )
    return ComparativeMethodMaturityReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=pgls_inputs.response,
        selected_model=model_comparison.better_model,
        reference_validation_passed=validate_comparative_reference_examples().all_passed,
        residual_diagnostics=residual_surfaces,
        sensitivity=sensitivity_summary,
        warnings=warnings,
    )


def _subset_covariance(
    covariance_matrix: list[list[float]],
    source_taxa: list[str],
    target_taxa: list[str],
) -> list[list[float]]:
    positions = {taxon: index for index, taxon in enumerate(source_taxa)}
    return [
        [
            covariance_matrix[positions[left_taxon]][positions[right_taxon]]
            for right_taxon in target_taxa
        ]
        for left_taxon in target_taxa
    ]


def _build_pgls_residual_warnings(
    residual_lambda: float,
    *,
    outlier_taxa: list[str],
    high_leverage_taxa: list[str],
) -> list[str]:
    warnings: list[str] = []
    if residual_lambda > 0.5:
        warnings.append("PGLS residuals retain moderate phylogenetic structure")
    if outlier_taxa:
        warnings.append("PGLS residual diagnostics identify one or more outlier taxa")
    if high_leverage_taxa:
        warnings.append(
            "PGLS leverage diagnostics identify one or more high-leverage taxa"
        )
    return warnings
