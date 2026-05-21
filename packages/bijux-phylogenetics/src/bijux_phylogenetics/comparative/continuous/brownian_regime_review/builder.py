from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.common import (
    load_comparative_dataset,
    summarize_numeric_trait_readiness,
)
from bijux_phylogenetics.comparative.continuous.model_fitting import (
    fit_brownian_motion_model,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .contracts import BrownianRegimeFitSummaryReport
from .covariance_fitting import fit_multirate_brownian_model
from .regime_mapping import (
    build_excluded_taxa,
    build_regime_covariance_components,
    load_branch_regime_rows,
)
from .review_metrics import (
    build_comparison_rows,
    build_identifiability_warnings,
    build_multirate_residual_diagnostics,
    build_profile_rows,
    build_regime_rows,
    chi_square_survival,
)


def summarize_brownian_regime_rates(
    tree_path: Path,
    traits_path: Path,
    regime_map_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    branch_id_column: str | None = None,
    regime_column: str = "regime",
) -> BrownianRegimeFitSummaryReport:
    readiness = summarize_numeric_trait_readiness(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=3,
        require_rooted=True,
        require_binary=False,
    )
    source_tree = load_tree(tree_path)
    branch_rows, resolved_branch_id_column = load_branch_regime_rows(
        regime_map_path,
        tree=source_tree,
        analyzed_taxa=dataset.taxa,
        branch_id_column=branch_id_column,
        regime_column=regime_column,
    )
    regime_matrices = build_regime_covariance_components(dataset.taxa, branch_rows)
    if len(regime_matrices) < 2:
        raise ComparativeMethodError(
            "multi-rate Brownian fitting requires at least two contributing regimes"
        )
    baseline = fit_brownian_motion_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    fit = fit_multirate_brownian_model(
        dataset,
        regime_matrices=regime_matrices,
        baseline_rate=baseline.rate,
    )
    comparison_rows = build_comparison_rows(
        taxon_count=dataset.taxon_count
        if hasattr(dataset, "taxon_count")
        else len(dataset.taxa),
        baseline_log_likelihood=baseline.log_likelihood,
        multirate_log_likelihood=fit.log_likelihood,
        regime_count=len(regime_matrices),
    )
    better_model = next(row.model for row in comparison_rows if row.selected)
    likelihood_ratio_degrees_of_freedom = len(regime_matrices) - 1
    likelihood_ratio_statistic = max(
        0.0,
        -2.0 * (baseline.log_likelihood - fit.log_likelihood),
    )
    likelihood_ratio_p_value = chi_square_survival(
        likelihood_ratio_statistic,
        likelihood_ratio_degrees_of_freedom,
    )
    profile_rows = build_profile_rows(
        dataset,
        regime_matrices=regime_matrices,
        fitted_rates=fit.regime_rates,
        baseline_rate=baseline.rate,
    )
    regime_rows = build_regime_rows(branch_rows, fit.regime_rates, profile_rows)
    identifiability_warnings = build_identifiability_warnings(
        regime_rows,
        profile_rows,
        better_model=better_model,
    )
    residual_diagnostics = build_multirate_residual_diagnostics(
        dataset,
        fit.covariance,
        fit.residuals,
    )
    warnings = list(
        dict.fromkeys(
            [
                *readiness.warnings,
                *residual_diagnostics.warnings,
                *[warning.message for warning in identifiability_warnings],
            ]
        )
    )
    return BrownianRegimeFitSummaryReport(
        tree_path=tree_path,
        traits_path=traits_path,
        regime_map_path=regime_map_path,
        taxon_column=dataset.taxon_column,
        branch_id_column=resolved_branch_id_column,
        regime_column=regime_column,
        trait=trait,
        tree_taxon_count=readiness.tree_taxa,
        analyzed_taxa=list(dataset.taxa),
        analyzed_taxon_count=len(dataset.taxa),
        excluded_taxa=build_excluded_taxa(readiness),
        branch_rows=branch_rows,
        regime_rows=regime_rows,
        profile_rows=profile_rows,
        root_state=fit.root_state,
        root_state_interval=fit.root_state_interval,
        log_likelihood=fit.log_likelihood,
        aic=next(row.aic for row in comparison_rows if row.model == "brownian-regimes"),
        aicc=next(
            row.aicc for row in comparison_rows if row.model == "brownian-regimes"
        ),
        comparison_rows=comparison_rows,
        better_model=better_model,
        likelihood_ratio_statistic=likelihood_ratio_statistic,
        likelihood_ratio_degrees_of_freedom=likelihood_ratio_degrees_of_freedom,
        likelihood_ratio_p_value=likelihood_ratio_p_value,
        likelihood_ratio_p_value_method="wilson-hilferty-approximation",
        identifiability_warnings=identifiability_warnings,
        residual_diagnostics=residual_diagnostics,
        assumptions=[
            "Every non-root branch must be assigned to one user-provided regime.",
            "Each regime contributes its own Brownian sigma-squared rate to shared-path covariance.",
            "Regime-specific uncertainty is reported from conditional likelihood profiles with other regimes fixed at their best-supported values.",
        ],
        warnings=warnings,
        readiness=readiness,
    )
