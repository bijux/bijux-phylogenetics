from __future__ import annotations

from importlib import import_module
from typing import Any

SUPPORTED_EVIDENCE_API_MODULES = ("bijux_phylogenetics.comparative",)
SUPPORTED_EVIDENCE_API_LOCATORS = (
    "bijux_phylogenetics.comparative:inspect_pgls_inputs",
    "bijux_phylogenetics.comparative:build_pgls_model_matrix",
    "bijux_phylogenetics.comparative:summarize_brownian_covariance_pgls",
    "bijux_phylogenetics.comparative:summarize_brownian_trait_evolution",
    "bijux_phylogenetics.comparative:summarize_ou_covariance_pgls",
    "bijux_phylogenetics.comparative:summarize_ou_trait_evolution",
    "bijux_phylogenetics.comparative:summarize_pgls_lambda_fit",
    "bijux_phylogenetics.comparative:summarize_pgls_categorical_contrasts",
    "bijux_phylogenetics.comparative:summarize_pgls_interaction_coefficients",
    "bijux_phylogenetics.comparative:summarize_independent_contrast_regression",
    "bijux_phylogenetics.comparative:summarize_phylogenetic_signal",
    "bijux_phylogenetics.comparative:summarize_phylogenetic_logistic",
    "bijux_phylogenetics.comparative:analyze_comparative_clade_stability",
    "bijux_phylogenetics.comparative:analyze_comparative_residual_clades",
    "bijux_phylogenetics.comparative:run_posterior_tree_pgls",
    "bijux_phylogenetics.comparative:compare_comparative_regression_models",
    "bijux_phylogenetics.comparative:run_multivariate_comparative_regression",
    "bijux_phylogenetics.comparative:run_pgls",
    "bijux_phylogenetics.comparative:write_brownian_covariance_table",
    "bijux_phylogenetics.comparative:write_brownian_trait_evolution_exclusion_table",
    "bijux_phylogenetics.comparative:write_brownian_trait_evolution_summary_table",
    "bijux_phylogenetics.comparative:write_ou_trait_evolution_exclusion_table",
    "bijux_phylogenetics.comparative:write_ou_trait_evolution_summary_table",
    "bijux_phylogenetics.comparative:write_comparative_clade_coefficient_change_table",
    "bijux_phylogenetics.comparative:write_comparative_clade_stability_table",
    "bijux_phylogenetics.comparative:write_comparative_residual_clade_table",
    "bijux_phylogenetics.comparative:write_comparative_residual_taxon_table",
    "bijux_phylogenetics.comparative:write_posterior_tree_pgls_coefficient_table",
    "bijux_phylogenetics.comparative:write_posterior_tree_pgls_summary_table",
    "bijux_phylogenetics.comparative:write_posterior_tree_pgls_tree_table",
    "bijux_phylogenetics.comparative:write_comparative_regression_excluded_taxa_table",
    "bijux_phylogenetics.comparative:write_comparative_regression_model_ranking_table",
    "bijux_phylogenetics.comparative:write_comparative_regression_pairwise_table",
    "bijux_phylogenetics.comparative:write_phylogenetic_logistic_coefficient_table",
    "bijux_phylogenetics.comparative:write_phylogenetic_logistic_excluded_taxa_table",
    "bijux_phylogenetics.comparative:write_phylogenetic_logistic_fitted_table",
    "bijux_phylogenetics.comparative:write_multivariate_excluded_taxa_table",
    "bijux_phylogenetics.comparative:write_multivariate_residual_association_table",
    "bijux_phylogenetics.comparative:write_multivariate_residual_covariance_table",
    "bijux_phylogenetics.comparative:write_independent_contrast_regression_table",
    "bijux_phylogenetics.comparative:write_independent_contrast_table",
    "bijux_phylogenetics.comparative:write_ou_alpha_profile_table",
    "bijux_phylogenetics.comparative:write_ou_covariance_table",
    "bijux_phylogenetics.comparative:write_phylogenetic_signal_permutation_table",
    "bijux_phylogenetics.comparative:write_phylogenetic_signal_summary_table",
    "bijux_phylogenetics.comparative:write_pgls_categorical_contrast_table",
    "bijux_phylogenetics.comparative:write_pgls_interaction_coefficient_table",
    "bijux_phylogenetics.comparative:write_pgls_lambda_profile_table",
    "bijux_phylogenetics.comparative:write_pgls_model_matrix_table",
    "bijux_phylogenetics.comparative:run_pgls_multiple_testing",
    "bijux_phylogenetics.comparative:compute_phylogenetic_signal_test",
    "bijux_phylogenetics.comparative:estimate_pagels_lambda",
    "bijux_phylogenetics.comparative:audit_ou_identifiability_reference_examples",
    "bijux_phylogenetics.comparative:compare_brownian_and_ou_models",
    "bijux_phylogenetics.comparative:fit_brownian_motion_model",
    "bijux_phylogenetics.comparative:fit_ornstein_uhlenbeck_model",
    "bijux_phylogenetics.comparative:transform_tree_for_evolutionary_mode",
)


def resolve_supported_evidence_api(locator: str) -> Any:
    """Resolve one governed evidence-consumer locator against the public runtime API."""
    if locator not in SUPPORTED_EVIDENCE_API_LOCATORS:
        raise KeyError(locator)
    module_name, export_name = locator.split(":", maxsplit=1)
    module = import_module(module_name)
    return getattr(module, export_name)
