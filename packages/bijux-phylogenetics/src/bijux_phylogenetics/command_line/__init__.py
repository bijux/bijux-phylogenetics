from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Any

__all__ = ["build_parser", "main", "run_command"]

from bijux_phylogenetics import __version__
from bijux_phylogenetics.ancestral.confidence import (
    build_continuous_ancestral_confidence_rows,
    build_continuous_ancestral_tree_set_confidence_rows,
    build_discrete_ancestral_confidence_rows,
    build_discrete_ancestral_tree_set_confidence_rows,
    summarize_continuous_ancestral_confidence,
    summarize_continuous_ancestral_tree_set_confidence,
    summarize_discrete_ancestral_confidence,
    summarize_discrete_ancestral_tree_set_confidence,
    write_ancestral_confidence_summary_table,
    write_continuous_ancestral_confidence_table,
    write_continuous_ancestral_tree_set_confidence_table,
    write_discrete_ancestral_confidence_table,
    write_discrete_ancestral_tree_set_confidence_table,
)
from bijux_phylogenetics.ancestral.continuous import (
    continuous_ancestral_exclusions,
    reconstruct_continuous_ancestral_states,
    summarize_continuous_ancestral_report,
    write_continuous_ancestral_exclusion_table,
    write_continuous_ancestral_summary_table,
    write_continuous_ancestral_uncertainty_table,
)
from bijux_phylogenetics.ancestral.discrete import (
    discrete_ancestral_exclusions,
    reconstruct_discrete_ancestral_states,
    summarize_discrete_ancestral_report,
    write_discrete_ancestral_exclusion_table,
    write_discrete_ancestral_fit_table,
    write_discrete_ancestral_probability_table,
    write_discrete_ancestral_summary_table,
    write_discrete_ancestral_transition_table,
)
from bijux_phylogenetics.ancestral.discrete_reference import (
    validate_discrete_ancestral_reference_examples,
)
from bijux_phylogenetics.ancestral.irreversible_discrete import (
    summarize_irreversible_discrete_reconstruction,
    summarize_irreversible_discrete_report,
    write_irreversible_discrete_fit_table,
    write_irreversible_discrete_node_table,
    write_irreversible_discrete_summary_table,
    write_irreversible_discrete_transition_table,
)
from bijux_phylogenetics.ancestral.ordered_discrete import (
    summarize_ordered_discrete_reconstruction,
    summarize_ordered_discrete_report,
    write_ordered_discrete_fit_table,
    write_ordered_discrete_node_table,
    write_ordered_discrete_summary_table,
    write_ordered_discrete_transition_table,
)
from bijux_phylogenetics.ancestral.package import build_ancestral_figure_package
from bijux_phylogenetics.ancestral.report_package import (
    build_ancestral_report_package,
)
from bijux_phylogenetics.ancestral.root_sensitivity import (
    summarize_ancestral_root_sensitivity,
    summarize_ancestral_root_sensitivity_report,
    write_ancestral_root_assumption_table,
    write_ancestral_root_sensitivity_node_table,
    write_ancestral_root_sensitivity_summary_table,
)
from bijux_phylogenetics.ancestral.sensitivity import build_ancestral_sensitivity_report
from bijux_phylogenetics.ancestral.service import (
    compare_continuous_ancestral_models,
    compare_discrete_ancestral_reconstructions,
    render_ancestral_state_report,
    write_ancestral_state_table,
    write_discrete_ancestral_comparison_table,
)
from bijux_phylogenetics.ancestral.transitions import (
    summarize_ancestral_transition_report,
    summarize_ancestral_transition_tree_set,
    summarize_ancestral_transition_tree_set_report,
    summarize_ancestral_transitions,
    write_ancestral_transition_branch_table,
    write_ancestral_transition_count_table,
    write_ancestral_transition_exclusion_table,
    write_ancestral_transition_summary_table,
    write_ancestral_transition_tree_set_branch_table,
    write_ancestral_transition_tree_set_count_table,
    write_ancestral_transition_tree_set_summary_table,
    write_ancestral_transition_tree_set_tree_table,
)
from bijux_phylogenetics.ancestral.tree_set import (
    summarize_continuous_ancestral_tree_set,
    summarize_continuous_ancestral_tree_set_report,
    summarize_discrete_ancestral_tree_set,
    summarize_discrete_ancestral_tree_set_report,
    write_ancestral_tree_set_exclusion_table,
    write_ancestral_tree_set_tree_table,
    write_continuous_ancestral_tree_set_clade_table,
    write_continuous_ancestral_tree_set_node_table,
    write_continuous_ancestral_tree_set_summary_table,
    write_discrete_ancestral_tree_set_clade_table,
    write_discrete_ancestral_tree_set_node_table,
    write_discrete_ancestral_tree_set_summary_table,
)
from bijux_phylogenetics.ancestral.visualization import (
    render_ancestral_state_visualization,
)
from bijux_phylogenetics.benchmark import (
    benchmark_alignment_diagnostics,
    benchmark_large_alignment_scaling,
    benchmark_large_dataset_stress_suite,
    benchmark_large_tree_model_fitting,
    benchmark_large_tree_scaling,
    benchmark_large_tree_set_scaling,
    benchmark_real_dataset_macroevolution,
    benchmark_tree_comparison,
    benchmark_tree_validation,
    benchmark_workflow_practical_limits,
)
from bijux_phylogenetics.parity import build_generated_geiger_parity_report
from bijux_phylogenetics.biogeography import (
    summarize_constrained_geographic_model,
    summarize_constrained_geographic_report,
    summarize_geographic_migration_event_tree_set,
    summarize_geographic_migration_events,
    summarize_geographic_sampling_bias,
    summarize_geographic_state_model,
    summarize_time_stratified_geographic_transitions,
    write_constrained_geographic_exclusion_table,
    write_constrained_geographic_fit_table,
    write_constrained_geographic_summary_table,
    write_constrained_geographic_transition_table,
    write_geographic_exclusion_table,
    write_geographic_migration_event_summary_table,
    write_geographic_migration_event_table,
    write_geographic_migration_exclusion_table,
    write_geographic_migration_tree_set_event_summary_table,
    write_geographic_migration_tree_set_event_table,
    write_geographic_migration_tree_set_exclusion_table,
    write_geographic_migration_tree_set_summary_table,
    write_geographic_migration_tree_set_tree_table,
    write_geographic_region_probability_table,
    write_geographic_sampling_bias_exclusion_table,
    write_geographic_sampling_bias_node_table,
    write_geographic_sampling_bias_summary_table,
    write_geographic_sampling_bias_transition_table,
    write_geographic_sampling_count_table,
    write_geographic_state_summary_table,
    write_geographic_transition_event_table,
    write_geographic_transition_rate_table,
    write_time_stratified_branch_table,
    write_time_stratified_exclusion_table,
    write_time_stratified_transition_matrix_table,
    write_time_stratified_transition_summary_table,
    write_unsupported_geographic_transition_claim_table,
)
from bijux_phylogenetics.biogeography.report_package import (
    build_biogeography_report_package,
)
from bijux_phylogenetics.biogeography.transition_chronology import (
    summarize_biogeographic_transition_chronology,
    write_dated_biogeography_event_table,
    write_dated_biogeography_exclusion_table,
    write_dated_biogeography_node_table,
    write_dated_biogeography_summary_table,
    write_dated_biogeography_time_bin_table,
)
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.output import (
    _print_commands,
    _print_result,
)
from bijux_phylogenetics.command_line.routing import (
    _command_inputs,
    _finalize_outputs,
)
from bijux_phylogenetics.command_line.engines import run_phylo_command
from bijux_phylogenetics.command_line.metadata import (
    add_metadata_commands,
    run_metadata_command,
)
from bijux_phylogenetics.command_line.prune import (
    add_prune_command,
    run_prune_command,
)
from bijux_phylogenetics.command_line.alignment import (
    add_alignment_commands,
    run_alignment_command,
)
from bijux_phylogenetics.command_line.annotate import (
    add_annotate_command,
    run_annotate_command,
)
from bijux_phylogenetics.command_line.comparative_continuous import (
    add_comparative_continuous_commands,
    run_comparative_continuous_command,
)
from bijux_phylogenetics.command_line.comparative_evolution import (
    add_comparative_evolution_commands,
    run_comparative_evolution_command,
)
from bijux_phylogenetics.command_line.parity import (
    add_parity_command,
    run_parity_command,
)
from bijux_phylogenetics.command_line.benchmark import (
    add_benchmark_commands,
    run_benchmark_command,
)
from bijux_phylogenetics.command_line.simulate import (
    add_simulate_command,
    run_simulate_command,
)
from bijux_phylogenetics.command_line.biogeography import (
    add_biogeography_commands,
    run_biogeography_command,
)
from bijux_phylogenetics.command_line.phylogeography import (
    add_phylogeography_commands,
    run_phylogeography_command,
)
from bijux_phylogenetics.command_line.ancestral import (
    add_ancestral_commands,
    run_ancestral_command,
)
from bijux_phylogenetics.command_line.host_association import (
    add_host_association_commands,
    run_host_association_command,
)
from bijux_phylogenetics.command_line.ecological_niche import (
    add_ecological_niche_commands,
    run_ecological_niche_command,
)
from bijux_phylogenetics.command_line.discrete_evolution import (
    add_discrete_evolution_commands,
    run_discrete_evolution_command,
)
from bijux_phylogenetics.command_line.diversification import (
    add_diversification_commands,
    run_diversification_command,
)
from bijux_phylogenetics.command_line.compare import (
    add_compare_command,
    run_compare_command,
)
from bijux_phylogenetics.command_line.adapter_beast import (
    add_beast_adapter_commands,
    run_beast_adapter_command,
)
from bijux_phylogenetics.command_line.adapter_bayesian import (
    add_bayesian_adapter_commands,
    run_bayesian_adapter_command,
)
from bijux_phylogenetics.command_line.adapter_inference import (
    add_inference_adapter_commands,
    run_inference_adapter_command,
)
from bijux_phylogenetics.command_line.adapter_mrbayes import (
    add_mrbayes_adapter_commands,
    run_mrbayes_adapter_command,
)
from bijux_phylogenetics.command_line.distance import (
    add_distance_commands,
    run_distance_command,
)
from bijux_phylogenetics.command_line.demo import (
    add_demo_command,
    run_demo_command,
)
from bijux_phylogenetics.command_line.evidence import (
    add_evidence_command,
    run_evidence_command,
)
from bijux_phylogenetics.command_line.report import (
    add_report_command,
    run_report_command,
)
from bijux_phylogenetics.command_line.render import (
    add_render_command,
    run_render_command,
)
from bijux_phylogenetics.command_line.diagnose import (
    add_diagnose_command,
    run_diagnose_command,
)
from bijux_phylogenetics.command_line.tree_set import (
    add_tree_set_commands,
    run_tree_set_command,
)
from bijux_phylogenetics.command_line.tree_inspection import (
    add_tree_inspection_commands,
    run_inspect_command,
    run_validate_command,
)
from bijux_phylogenetics.command_line.topology import (
    add_topology_commands,
    run_topology_command,
)
from bijux_phylogenetics.command_line.taxonomy import (
    add_taxonomy_commands,
    run_taxonomy_command,
)
from bijux_phylogenetics.command_line.traits import (
    add_traits_commands,
    run_traits_command,
)
from bijux_phylogenetics.command_line.arguments import (
    _adapter_version_args,
    _add_external_adapter_execution_arguments,
    _add_manifest_argument,
    _add_preflight_executable_arguments,
    _json_requested,
    _parse_assignment_map,
    _parse_time_bin_definition,
    _parse_transition_pairs,
    _validate_ancestral_discrete_model_arguments,
)
from bijux_phylogenetics.comparative.clade_residuals import (
    analyze_comparative_residual_clades,
    write_comparative_residual_clade_table,
    write_comparative_residual_taxon_table,
)
from bijux_phylogenetics.comparative.clade_stability import (
    analyze_comparative_clade_stability,
    write_comparative_clade_coefficient_change_table,
    write_comparative_clade_stability_table,
)
from bijux_phylogenetics.comparative.clade_traits import (
    summarize_clade_traits,
    write_clade_trait_clade_table,
    write_clade_trait_exclusion_table,
    write_clade_trait_summary_table,
)
from bijux_phylogenetics.comparative.covariance_audit import (
    summarize_comparative_covariance_audit,
    write_comparative_covariance_audit_candidate_table,
    write_comparative_covariance_audit_excluded_taxa_table,
    write_comparative_covariance_audit_summary_table,
)
from bijux_phylogenetics.comparative.models import (
    assess_comparative_method_maturity,
)
from bijux_phylogenetics.comparative.multivariate_regression import (
    run_multivariate_comparative_regression,
    write_multivariate_excluded_taxa_table,
    write_multivariate_residual_association_table,
    write_multivariate_residual_correlation_table,
    write_multivariate_residual_covariance_table,
    write_multivariate_response_coefficient_table,
    write_multivariate_response_model_table,
)
from bijux_phylogenetics.comparative.pgls import (
    inspect_pgls_inputs,
    run_pgls,
    run_pgls_multiple_testing,
    write_pgls_model_matrix_table,
)
from bijux_phylogenetics.comparative.pgls_brownian_covariance import (
    summarize_brownian_covariance_pgls,
    write_brownian_covariance_table,
)
from bijux_phylogenetics.comparative.pgls_categorical_contrasts import (
    summarize_pgls_categorical_contrasts,
    write_pgls_categorical_contrast_table,
)
from bijux_phylogenetics.comparative.pgls_interaction_coefficients import (
    summarize_pgls_interaction_coefficients,
    write_pgls_interaction_coefficient_table,
)
from bijux_phylogenetics.comparative.pgls_lambda_fit import (
    write_pgls_lambda_profile_table,
)
from bijux_phylogenetics.comparative.pgls_ou_covariance import (
    summarize_ou_covariance_pgls,
    write_ou_alpha_profile_table,
    write_ou_covariance_table,
)
from bijux_phylogenetics.comparative.phylogenetic_logistic import (
    summarize_phylogenetic_logistic,
    write_phylogenetic_logistic_coefficient_table,
    write_phylogenetic_logistic_excluded_taxa_table,
    write_phylogenetic_logistic_fitted_table,
)
from bijux_phylogenetics.comparative.phylogenetic_anova import (
    summarize_phylogenetic_anova,
    write_phylogenetic_anova_exclusion_table,
    write_phylogenetic_anova_group_table,
    write_phylogenetic_anova_pairwise_table,
    write_phylogenetic_anova_simulation_table,
    write_phylogenetic_anova_summary_table,
)
from bijux_phylogenetics.comparative.phylogenetic_residuals import (
    summarize_phylogenetic_residuals,
    write_phylogenetic_residual_coefficient_table,
    write_phylogenetic_residual_exclusion_table,
    write_phylogenetic_residual_summary_table,
    write_phylogenetic_residual_taxon_table,
)
from bijux_phylogenetics.comparative.posterior_tree_pgls import (
    run_posterior_tree_pgls,
    write_posterior_tree_pgls_coefficient_table,
    write_posterior_tree_pgls_summary_table,
    write_posterior_tree_pgls_tree_table,
)
from bijux_phylogenetics.comparative.regression_model_selection import (
    compare_comparative_regression_models,
    write_comparative_regression_excluded_taxa_table,
    write_comparative_regression_model_ranking_table,
    write_comparative_regression_pairwise_table,
)
from bijux_phylogenetics.comparative.report_package import (
    build_comparative_report_package,
)
from bijux_phylogenetics.comparative.reporting import (
    build_comparative_method_report,
    build_trait_influence_report,
    compare_comparative_results_across_pruning,
    compare_comparative_results_across_trees,
    write_comparative_method_report,
    write_comparative_methods_summary_text,
)
from bijux_phylogenetics.comparative.trait_imputation import (
    summarize_trait_imputation,
    write_trait_imputation_exclusion_table,
    write_trait_imputation_holdout_table,
    write_trait_imputation_summary_table,
    write_trait_imputation_table,
)
from bijux_phylogenetics.comparative.trait_outliers import (
    summarize_trait_outliers,
    write_trait_outlier_exclusion_table,
    write_trait_outlier_summary_table,
    write_trait_outlier_taxon_table,
)
from bijux_phylogenetics.core.environment import inspect_environment
from bijux_phylogenetics.core.taxonomy import (
    normalize_tree_taxa,
    write_taxon_mapping,
)
from bijux_phylogenetics.engines import (
    list_external_engine_workflows,
    list_mafft_alignment_modes,
    list_trimal_trimming_modes,
    read_engine_version,
    render_inference_workflow_report,
    run_alignment_trimming,
    run_bootstrap_consensus_tree,
    run_bootstrap_support_estimation,
    run_codon_aware_multiple_sequence_alignment,
    run_fasta_to_tree_workflow,
    run_maximum_likelihood_tree_inference,
    run_model_selection,
    run_multiple_sequence_alignment,
    run_sh_alrt_support_estimation,
)
from bijux_phylogenetics.runtime.errors import EngineUnavailableError, PhylogeneticsError
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylogeography import (
    render_coordinate_movement_visualization,
    render_geographic_map_html,
    summarize_continuous_phylogeography,
    summarize_continuous_phylogeography_map,
    summarize_discrete_region_map,
    write_coordinate_estimate_table,
    write_coordinate_movement_branch_table,
    write_coordinate_movement_exclusion_table,
    write_coordinate_movement_outlier_table,
    write_coordinate_movement_summary_table,
    write_geographic_map_exclusion_table,
    write_geographic_map_line_table,
    write_geographic_map_marker_table,
    write_geographic_map_summary_table,
)
from bijux_phylogenetics.provenance.method_tiers import (
    method_tier_metrics,
    method_tier_warnings,
)
from bijux_phylogenetics.runtime.results import build_command_result, build_error_result
_COMMAND_LINE_BENCHMARK_API = {
    "benchmark_alignment_diagnostics": benchmark_alignment_diagnostics,
    "benchmark_large_alignment_scaling": benchmark_large_alignment_scaling,
    "benchmark_large_dataset_stress_suite": benchmark_large_dataset_stress_suite,
    "benchmark_large_tree_model_fitting": benchmark_large_tree_model_fitting,
    "benchmark_large_tree_scaling": benchmark_large_tree_scaling,
    "benchmark_large_tree_set_scaling": benchmark_large_tree_set_scaling,
    "benchmark_real_dataset_macroevolution": benchmark_real_dataset_macroevolution,
    "benchmark_tree_comparison": benchmark_tree_comparison,
    "benchmark_tree_validation": benchmark_tree_validation,
    "benchmark_workflow_practical_limits": benchmark_workflow_practical_limits,
}

_COMMAND_LINE_PARITY_API = {
    "build_generated_geiger_parity_report": build_generated_geiger_parity_report,
}

_COMMAND_LINE_BIOGEOGRAPHY_API = (
    _parse_time_bin_definition,
    build_biogeography_report_package,
    summarize_biogeographic_transition_chronology,
    summarize_constrained_geographic_model,
    summarize_constrained_geographic_report,
    summarize_geographic_migration_event_tree_set,
    summarize_geographic_migration_events,
    summarize_geographic_sampling_bias,
    summarize_geographic_state_model,
    summarize_time_stratified_geographic_transitions,
    write_constrained_geographic_exclusion_table,
    write_constrained_geographic_fit_table,
    write_constrained_geographic_summary_table,
    write_constrained_geographic_transition_table,
    write_dated_biogeography_event_table,
    write_dated_biogeography_exclusion_table,
    write_dated_biogeography_node_table,
    write_dated_biogeography_summary_table,
    write_dated_biogeography_time_bin_table,
    write_geographic_exclusion_table,
    write_geographic_migration_event_summary_table,
    write_geographic_migration_event_table,
    write_geographic_migration_exclusion_table,
    write_geographic_migration_tree_set_event_summary_table,
    write_geographic_migration_tree_set_event_table,
    write_geographic_migration_tree_set_exclusion_table,
    write_geographic_migration_tree_set_summary_table,
    write_geographic_migration_tree_set_tree_table,
    write_geographic_region_probability_table,
    write_geographic_sampling_bias_exclusion_table,
    write_geographic_sampling_bias_node_table,
    write_geographic_sampling_bias_summary_table,
    write_geographic_sampling_bias_transition_table,
    write_geographic_sampling_count_table,
    write_geographic_state_summary_table,
    write_geographic_transition_event_table,
    write_geographic_transition_rate_table,
    write_time_stratified_branch_table,
    write_time_stratified_exclusion_table,
    write_time_stratified_transition_matrix_table,
    write_time_stratified_transition_summary_table,
    write_unsupported_geographic_transition_claim_table,
)

_COMMAND_LINE_PHYLOGEOGRAPHY_API = (
    render_coordinate_movement_visualization,
    render_geographic_map_html,
    summarize_continuous_phylogeography,
    summarize_continuous_phylogeography_map,
    summarize_discrete_region_map,
    write_coordinate_estimate_table,
    write_coordinate_movement_branch_table,
    write_coordinate_movement_exclusion_table,
    write_coordinate_movement_outlier_table,
    write_coordinate_movement_summary_table,
    write_geographic_map_exclusion_table,
    write_geographic_map_line_table,
    write_geographic_map_marker_table,
    write_geographic_map_summary_table,
)

_COMMAND_LINE_ANCESTRAL_API = (
    _parse_assignment_map,
    _parse_transition_pairs,
    _validate_ancestral_discrete_model_arguments,
    build_ancestral_figure_package,
    build_ancestral_report_package,
    build_ancestral_sensitivity_report,
    build_continuous_ancestral_confidence_rows,
    build_continuous_ancestral_tree_set_confidence_rows,
    build_discrete_ancestral_confidence_rows,
    build_discrete_ancestral_tree_set_confidence_rows,
    compare_continuous_ancestral_models,
    compare_discrete_ancestral_reconstructions,
    continuous_ancestral_exclusions,
    discrete_ancestral_exclusions,
    reconstruct_continuous_ancestral_states,
    reconstruct_discrete_ancestral_states,
    render_ancestral_state_report,
    render_ancestral_state_visualization,
    summarize_ancestral_root_sensitivity,
    summarize_ancestral_root_sensitivity_report,
    summarize_ancestral_transition_report,
    summarize_ancestral_transition_tree_set,
    summarize_ancestral_transition_tree_set_report,
    summarize_ancestral_transitions,
    summarize_continuous_ancestral_confidence,
    summarize_continuous_ancestral_report,
    summarize_continuous_ancestral_tree_set,
    summarize_continuous_ancestral_tree_set_confidence,
    summarize_continuous_ancestral_tree_set_report,
    summarize_discrete_ancestral_confidence,
    summarize_discrete_ancestral_report,
    summarize_discrete_ancestral_tree_set,
    summarize_discrete_ancestral_tree_set_confidence,
    summarize_discrete_ancestral_tree_set_report,
    summarize_irreversible_discrete_reconstruction,
    summarize_irreversible_discrete_report,
    summarize_ordered_discrete_reconstruction,
    summarize_ordered_discrete_report,
    validate_discrete_ancestral_reference_examples,
    write_ancestral_confidence_summary_table,
    write_ancestral_root_assumption_table,
    write_ancestral_root_sensitivity_node_table,
    write_ancestral_root_sensitivity_summary_table,
    write_ancestral_state_table,
    write_ancestral_transition_branch_table,
    write_ancestral_transition_count_table,
    write_ancestral_transition_exclusion_table,
    write_ancestral_transition_summary_table,
    write_ancestral_transition_tree_set_branch_table,
    write_ancestral_transition_tree_set_count_table,
    write_ancestral_transition_tree_set_summary_table,
    write_ancestral_transition_tree_set_tree_table,
    write_ancestral_tree_set_exclusion_table,
    write_ancestral_tree_set_tree_table,
    write_continuous_ancestral_confidence_table,
    write_continuous_ancestral_exclusion_table,
    write_continuous_ancestral_summary_table,
    write_continuous_ancestral_tree_set_clade_table,
    write_continuous_ancestral_tree_set_confidence_table,
    write_continuous_ancestral_tree_set_node_table,
    write_continuous_ancestral_tree_set_summary_table,
    write_continuous_ancestral_uncertainty_table,
    write_discrete_ancestral_comparison_table,
    write_discrete_ancestral_confidence_table,
    write_discrete_ancestral_exclusion_table,
    write_discrete_ancestral_fit_table,
    write_discrete_ancestral_probability_table,
    write_discrete_ancestral_summary_table,
    write_discrete_ancestral_transition_table,
    write_discrete_ancestral_tree_set_clade_table,
    write_discrete_ancestral_tree_set_confidence_table,
    write_discrete_ancestral_tree_set_node_table,
    write_discrete_ancestral_tree_set_summary_table,
    write_irreversible_discrete_fit_table,
    write_irreversible_discrete_node_table,
    write_irreversible_discrete_summary_table,
    write_irreversible_discrete_transition_table,
    write_ordered_discrete_fit_table,
    write_ordered_discrete_node_table,
    write_ordered_discrete_summary_table,
    write_ordered_discrete_transition_table,
)


def build_parser() -> argparse.ArgumentParser:
    """Build the repository CLI parser."""
    parser = argparse.ArgumentParser(prog="bijux-phylogenetics")
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    commands = subparsers.add_parser(
        "commands", help="List the registered command taxonomy."
    )
    commands.add_argument("--format", choices=("text", "json"), default="text")

    env = subparsers.add_parser(
        get_command_spec("env").name, help=get_command_spec("env").summary
    )
    env_subparsers = env.add_subparsers(dest="env_command", required=True)
    env_inspect = env_subparsers.add_parser(
        "inspect", help="Inspect runtime dependency availability."
    )
    env_inspect.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(env_inspect)

    phylo = subparsers.add_parser(
        get_command_spec("phylo").name, help=get_command_spec("phylo").summary
    )
    phylo_subparsers = phylo.add_subparsers(dest="phylo_command", required=True)
    phylo_preflight = phylo_subparsers.add_parser(
        "preflight",
        help="Inspect external engine availability, version support, and workflow readiness.",
    )
    phylo_preflight.add_argument(
        "--workflow",
        choices=list_external_engine_workflows(),
        help="Require one selected external-engine workflow to be runnable in the current environment.",
    )
    _add_preflight_executable_arguments(phylo_preflight)
    phylo_preflight.add_argument(
        "--json", action="store_true", help="Emit the preflight report as JSON."
    )
    _add_manifest_argument(phylo_preflight)
    phylo_run = phylo_subparsers.add_parser(
        "run",
        help="Run one governed workflow from one YAML or JSON config file and export a validated result bundle.",
    )
    phylo_run.add_argument("config_path", type=Path)
    phylo_run.add_argument(
        "--json", action="store_true", help="Emit the config-run report as JSON."
    )
    _add_manifest_argument(phylo_run)
    phylo_replay = phylo_subparsers.add_parser(
        "replay",
        help="Rerun one governed phylogenetics workflow from its manifest and compare the replayed outputs.",
    )
    phylo_replay.add_argument("manifest_path", type=Path)
    phylo_replay.add_argument("--out-dir", type=Path)
    _add_preflight_executable_arguments(phylo_replay)
    phylo_replay.add_argument(
        "--json", action="store_true", help="Emit the replay report as JSON."
    )
    _add_manifest_argument(phylo_replay)
    phylo_bundle = phylo_subparsers.add_parser(
        "bundle",
        help="Export one portable workflow-result bundle from a governed workflow manifest.",
    )
    phylo_bundle.add_argument("manifest_path", type=Path)
    phylo_bundle.add_argument("--out-dir", required=True, type=Path)
    phylo_bundle.add_argument(
        "--json", action="store_true", help="Emit the bundle report as JSON."
    )
    _add_manifest_argument(phylo_bundle)
    phylo_validate_bundle = phylo_subparsers.add_parser(
        "validate-bundle",
        help="Validate one workflow-result bundle for checksum integrity and required workflow contents.",
    )
    phylo_validate_bundle.add_argument("bundle_root", type=Path)
    phylo_validate_bundle.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(phylo_validate_bundle)

    add_metadata_commands(subparsers)

    add_traits_commands(subparsers)

    add_prune_command(subparsers)

    add_alignment_commands(subparsers)

    comparative = subparsers.add_parser(
        get_command_spec("comparative").name,
        help=get_command_spec("comparative").summary,
    )
    comparative_subparsers = comparative.add_subparsers(
        dest="comparative_command", required=True
    )
    add_comparative_continuous_commands(comparative_subparsers)
    add_comparative_evolution_commands(comparative_subparsers)
    comparative_clade_traits = comparative_subparsers.add_parser(
        "clade-traits",
        help="Summarize one continuous or categorical trait across internal clades.",
    )
    comparative_clade_traits.add_argument("tree", type=Path)
    comparative_clade_traits.add_argument("table", type=Path)
    comparative_clade_traits.add_argument("--trait", required=True)
    comparative_clade_traits.add_argument("--taxon-column")
    comparative_clade_traits.add_argument(
        "--trait-kind",
        choices=("auto", "continuous", "categorical"),
        default="auto",
        help="Infer trait kind automatically or force continuous/categorical handling.",
    )
    comparative_clade_traits.add_argument(
        "--min-clade-size",
        type=int,
        default=2,
        help="Only summarize internal clades with at least this many analyzed taxa.",
    )
    comparative_clade_traits.add_argument(
        "--summary-out",
        type=Path,
        help="Write one clade-trait summary ledger as TSV or CSV.",
    )
    comparative_clade_traits.add_argument(
        "--clades-out",
        type=Path,
        help="Write one internal clade-trait ledger as TSV or CSV.",
    )
    comparative_clade_traits.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for clade trait summarization as TSV or CSV.",
    )
    comparative_clade_traits.add_argument(
        "--json", action="store_true", help="Emit the clade-trait report as JSON."
    )
    _add_manifest_argument(comparative_clade_traits)
    comparative_phylogenetic_residuals = comparative_subparsers.add_parser(
        "phylogenetic-residuals",
        help="Summarize tree-aware fitted values and residuals for one continuous response and predictor.",
    )
    comparative_phylogenetic_residuals.add_argument("tree", type=Path)
    comparative_phylogenetic_residuals.add_argument("table", type=Path)
    comparative_phylogenetic_residuals.add_argument("--response", required=True)
    comparative_phylogenetic_residuals.add_argument("--predictor", required=True)
    comparative_phylogenetic_residuals.add_argument("--taxon-column")
    comparative_phylogenetic_residuals.add_argument(
        "--method",
        choices=("brownian", "lambda"),
        default="lambda",
        help="Use fixed Brownian covariance or estimate Pagel lambda before computing residuals.",
    )
    comparative_phylogenetic_residuals.add_argument(
        "--summary-out",
        type=Path,
        help="Write one phylogenetic-residual summary ledger as TSV or CSV.",
    )
    comparative_phylogenetic_residuals.add_argument(
        "--residuals-out",
        type=Path,
        help="Write one taxon-level phylogenetic-residual ledger as TSV or CSV.",
    )
    comparative_phylogenetic_residuals.add_argument(
        "--coefficients-out",
        type=Path,
        help="Write one phylogenetic-residual coefficient ledger as TSV or CSV.",
    )
    comparative_phylogenetic_residuals.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for phylogenetic residual review as TSV or CSV.",
    )
    comparative_phylogenetic_residuals.add_argument(
        "--json",
        action="store_true",
        help="Emit the phylogenetic-residual review as JSON.",
    )
    _add_manifest_argument(comparative_phylogenetic_residuals)
    comparative_phylogenetic_anova = comparative_subparsers.add_parser(
        "phylogenetic-anova",
        help="Run a simulation-based phylogenetic ANOVA over one continuous response and one categorical group.",
    )
    comparative_phylogenetic_anova.add_argument("tree", type=Path)
    comparative_phylogenetic_anova.add_argument("table", type=Path)
    comparative_phylogenetic_anova.add_argument("--response", required=True)
    comparative_phylogenetic_anova.add_argument("--group", required=True)
    comparative_phylogenetic_anova.add_argument("--taxon-column")
    comparative_phylogenetic_anova.add_argument(
        "--simulations",
        type=int,
        default=1000,
        help="Number of observed-plus-null F statistics to evaluate.",
    )
    comparative_phylogenetic_anova.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Seed for the Brownian null simulation sequence.",
    )
    comparative_phylogenetic_anova.add_argument(
        "--summary-out",
        type=Path,
        help="Write one phylogenetic-ANOVA summary ledger as TSV or CSV.",
    )
    comparative_phylogenetic_anova.add_argument(
        "--groups-out",
        type=Path,
        help="Write one group-summary ledger for phylogenetic ANOVA as TSV or CSV.",
    )
    comparative_phylogenetic_anova.add_argument(
        "--pairwise-out",
        type=Path,
        help="Write one pairwise-comparison ledger for phylogenetic ANOVA as TSV or CSV.",
    )
    comparative_phylogenetic_anova.add_argument(
        "--simulations-out",
        type=Path,
        help="Write one observed-plus-null F-statistic ledger for phylogenetic ANOVA as TSV or CSV.",
    )
    comparative_phylogenetic_anova.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for phylogenetic ANOVA as TSV or CSV.",
    )
    comparative_phylogenetic_anova.add_argument(
        "--json",
        action="store_true",
        help="Emit the phylogenetic-ANOVA report as JSON.",
    )
    _add_manifest_argument(comparative_phylogenetic_anova)
    comparative_trait_outliers = comparative_subparsers.add_parser(
        "trait-outliers",
        help="Rank continuous-trait taxa by conditional phylogenetic residual size.",
    )
    comparative_trait_outliers.add_argument("tree", type=Path)
    comparative_trait_outliers.add_argument("table", type=Path)
    comparative_trait_outliers.add_argument("--trait", required=True)
    comparative_trait_outliers.add_argument("--taxon-column")
    comparative_trait_outliers.add_argument(
        "--summary-out",
        type=Path,
        help="Write one trait-outlier summary ledger as TSV or CSV.",
    )
    comparative_trait_outliers.add_argument(
        "--outliers-out",
        type=Path,
        help="Write one ranked taxon outlier ledger as TSV or CSV.",
    )
    comparative_trait_outliers.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for trait outlier review as TSV or CSV.",
    )
    comparative_trait_outliers.add_argument(
        "--json",
        action="store_true",
        help="Emit the trait outlier review as JSON.",
    )
    _add_manifest_argument(comparative_trait_outliers)
    comparative_trait_imputation = comparative_subparsers.add_parser(
        "trait-imputation",
        help="Impute missing continuous-trait values under a Brownian phylogenetic model.",
    )
    comparative_trait_imputation.add_argument("tree", type=Path)
    comparative_trait_imputation.add_argument("table", type=Path)
    comparative_trait_imputation.add_argument("--trait", required=True)
    comparative_trait_imputation.add_argument("--taxon-column")
    comparative_trait_imputation.add_argument(
        "--summary-out",
        type=Path,
        help="Write one Brownian trait-imputation summary ledger as TSV or CSV.",
    )
    comparative_trait_imputation.add_argument(
        "--imputations-out",
        type=Path,
        help="Write one imputed-value ledger as TSV or CSV.",
    )
    comparative_trait_imputation.add_argument(
        "--holdout-out",
        type=Path,
        help="Write one leave-one-observed-out validation ledger as TSV or CSV.",
    )
    comparative_trait_imputation.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for trait imputation as TSV or CSV.",
    )
    comparative_trait_imputation.add_argument(
        "--json",
        action="store_true",
        help="Emit the trait-imputation review as JSON.",
    )
    _add_manifest_argument(comparative_trait_imputation)
    comparative_maturity = comparative_subparsers.add_parser(
        "maturity",
        help="Audit comparative residual diagnostics and sensitivity for one response trait workflow.",
    )
    comparative_maturity.add_argument("tree", type=Path)
    comparative_maturity.add_argument("table", type=Path)
    comparative_maturity.add_argument("--response")
    comparative_maturity.add_argument("--predictors", nargs="+")
    comparative_maturity.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass * habitat'.",
    )
    comparative_maturity.add_argument("--taxon-column")
    comparative_maturity.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_maturity.add_argument(
        "--json", action="store_true", help="Emit the maturity audit as JSON."
    )
    _add_manifest_argument(comparative_maturity)
    comparative_covariance_audit = comparative_subparsers.add_parser(
        "covariance-audit",
        help="Audit comparative covariance readiness for PGLS and continuous trait models.",
    )
    comparative_covariance_audit.add_argument("tree", type=Path)
    comparative_covariance_audit.add_argument("table", type=Path)
    comparative_covariance_audit.add_argument(
        "--analysis",
        choices=("pgls", "brownian-trait", "ou-trait"),
        required=True,
    )
    comparative_covariance_audit.add_argument("--trait")
    comparative_covariance_audit.add_argument("--response")
    comparative_covariance_audit.add_argument("--predictors", nargs="+")
    comparative_covariance_audit.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass * habitat'.",
    )
    comparative_covariance_audit.add_argument("--taxon-column")
    comparative_covariance_audit.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_covariance_audit.add_argument(
        "--alpha",
        default="estimate",
        help="Use 'estimate' or a positive OU alpha value.",
    )
    comparative_covariance_audit.add_argument(
        "--summary-out",
        type=Path,
        help="Write one covariance-audit summary row as TSV or CSV.",
    )
    comparative_covariance_audit.add_argument(
        "--candidates-out",
        type=Path,
        help="Write candidate covariance audit rows as TSV or CSV.",
    )
    comparative_covariance_audit.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write excluded covariance-audit taxa as TSV or CSV.",
    )
    comparative_covariance_audit.add_argument(
        "--json", action="store_true", help="Emit the covariance audit as JSON."
    )
    _add_manifest_argument(comparative_covariance_audit)
    comparative_pgls = comparative_subparsers.add_parser(
        "pgls",
        help="Fit a phylogenetic generalized least-squares model.",
    )
    comparative_pgls.add_argument("tree", type=Path)
    comparative_pgls.add_argument("table", type=Path)
    comparative_pgls.add_argument("--response")
    comparative_pgls.add_argument("--predictors", nargs="+")
    comparative_pgls.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass * habitat'.",
    )
    comparative_pgls.add_argument("--taxon-column")
    comparative_pgls.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_pgls.add_argument(
        "--model-matrix-out",
        type=Path,
        help="Write the encoded comparative model matrix as TSV or CSV.",
    )
    comparative_pgls.add_argument(
        "--categorical-contrasts-out",
        type=Path,
        help="Write categorical predictor contrast rows as TSV or CSV.",
    )
    comparative_pgls.add_argument(
        "--interaction-coefficients-out",
        type=Path,
        help="Write interaction coefficient rows as TSV or CSV.",
    )
    comparative_pgls.add_argument(
        "--lambda-profile-out",
        type=Path,
        help="Write the fitted Pagel lambda likelihood profile as TSV or CSV.",
    )
    comparative_pgls.add_argument(
        "--json", action="store_true", help="Emit the model result as JSON."
    )
    _add_manifest_argument(comparative_pgls)
    comparative_brownian_pgls = comparative_subparsers.add_parser(
        "brownian-pgls",
        help="Fit a PGLS model under fixed Brownian shared-path covariance.",
    )
    comparative_brownian_pgls.add_argument("tree", type=Path)
    comparative_brownian_pgls.add_argument("table", type=Path)
    comparative_brownian_pgls.add_argument("--response")
    comparative_brownian_pgls.add_argument("--predictors", nargs="+")
    comparative_brownian_pgls.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass * habitat'.",
    )
    comparative_brownian_pgls.add_argument("--taxon-column")
    comparative_brownian_pgls.add_argument(
        "--covariance-out",
        type=Path,
        help="Write the pairwise Brownian covariance ledger as TSV or CSV.",
    )
    comparative_brownian_pgls.add_argument(
        "--json", action="store_true", help="Emit the Brownian PGLS result as JSON."
    )
    _add_manifest_argument(comparative_brownian_pgls)
    comparative_ou_pgls = comparative_subparsers.add_parser(
        "ou-pgls",
        help="Fit a PGLS model under stationary-root OU covariance.",
    )
    comparative_ou_pgls.add_argument("tree", type=Path)
    comparative_ou_pgls.add_argument("table", type=Path)
    comparative_ou_pgls.add_argument("--response")
    comparative_ou_pgls.add_argument("--predictors", nargs="+")
    comparative_ou_pgls.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass * habitat'.",
    )
    comparative_ou_pgls.add_argument("--taxon-column")
    comparative_ou_pgls.add_argument(
        "--alpha",
        default="estimate",
        help="Use 'estimate' or a positive numeric OU alpha value.",
    )
    comparative_ou_pgls.add_argument(
        "--covariance-out",
        type=Path,
        help="Write the pairwise OU covariance ledger as TSV or CSV.",
    )
    comparative_ou_pgls.add_argument(
        "--alpha-profile-out",
        type=Path,
        help="Write the fitted OU alpha likelihood profile as TSV or CSV.",
    )
    comparative_ou_pgls.add_argument(
        "--json",
        action="store_true",
        help="Emit the OU covariance PGLS result as JSON.",
    )
    _add_manifest_argument(comparative_ou_pgls)
    comparative_multiple_testing = comparative_subparsers.add_parser(
        "multiple-testing",
        help="Adjust PGLS coefficient p-values across many response traits.",
    )
    comparative_multiple_testing.add_argument("tree", type=Path)
    comparative_multiple_testing.add_argument("table", type=Path)
    comparative_multiple_testing.add_argument("--responses", nargs="+", required=True)
    comparative_multiple_testing.add_argument("--predictors", nargs="+", required=True)
    comparative_multiple_testing.add_argument("--taxon-column")
    comparative_multiple_testing.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_multiple_testing.add_argument(
        "--json", action="store_true", help="Emit the correction report as JSON."
    )
    _add_manifest_argument(comparative_multiple_testing)
    comparative_logistic = comparative_subparsers.add_parser(
        "logistic",
        help="Fit a binary phylogenetic logistic approximation with a phylogenetic working correlation.",
    )
    comparative_logistic.add_argument("tree", type=Path)
    comparative_logistic.add_argument("table", type=Path)
    comparative_logistic.add_argument("--response")
    comparative_logistic.add_argument("--predictors", nargs="+")
    comparative_logistic.add_argument(
        "--formula",
        help="Formula-style specification such as 'presence ~ body_mass + habitat'.",
    )
    comparative_logistic.add_argument("--taxon-column")
    comparative_logistic.add_argument(
        "--lambda-value",
        default="1.0",
        help="Use a numeric Pagel lambda value between 0 and 1 for the working correlation.",
    )
    comparative_logistic.add_argument(
        "--coefficients-out",
        type=Path,
        help="Write the fitted logistic coefficient ledger as TSV or CSV.",
    )
    comparative_logistic.add_argument(
        "--fitted-out",
        type=Path,
        help="Write the fitted taxon-level probability ledger as TSV or CSV.",
    )
    comparative_logistic.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write the excluded-taxa ledger as TSV or CSV.",
    )
    comparative_logistic.add_argument(
        "--json", action="store_true", help="Emit the logistic result as JSON."
    )
    _add_manifest_argument(comparative_logistic)
    comparative_clade_residuals = comparative_subparsers.add_parser(
        "clade-residuals",
        help="Aggregate comparative model residuals across internal clades.",
    )
    comparative_clade_residuals.add_argument("tree", type=Path)
    comparative_clade_residuals.add_argument("table", type=Path)
    comparative_clade_residuals.add_argument("--response")
    comparative_clade_residuals.add_argument("--predictors", nargs="+")
    comparative_clade_residuals.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass + habitat'.",
    )
    comparative_clade_residuals.add_argument("--taxon-column")
    comparative_clade_residuals.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1. Binary-response residual aggregation requires a numeric value.",
    )
    comparative_clade_residuals.add_argument(
        "--taxa-out",
        type=Path,
        help="Write the analyzed taxon residual ledger as TSV or CSV.",
    )
    comparative_clade_residuals.add_argument(
        "--clades-out",
        type=Path,
        help="Write the internal clade residual aggregation ledger as TSV or CSV.",
    )
    comparative_clade_residuals.add_argument(
        "--json", action="store_true", help="Emit the clade residual report as JSON."
    )
    _add_manifest_argument(comparative_clade_residuals)
    comparative_clade_stability = comparative_subparsers.add_parser(
        "clade-stability",
        help="Refit one comparative model after removing each major internal clade.",
    )
    comparative_clade_stability.add_argument("tree", type=Path)
    comparative_clade_stability.add_argument("table", type=Path)
    comparative_clade_stability.add_argument("--response")
    comparative_clade_stability.add_argument("--predictors", nargs="+")
    comparative_clade_stability.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass + habitat'.",
    )
    comparative_clade_stability.add_argument("--taxon-column")
    comparative_clade_stability.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1. Binary-response clade stability requires a numeric value.",
    )
    comparative_clade_stability.add_argument(
        "--clades-out",
        type=Path,
        help="Write the leave-one-clade-out stability summary ledger as TSV or CSV.",
    )
    comparative_clade_stability.add_argument(
        "--terms-out",
        type=Path,
        help="Write the coefficient-delta ledger across clade removals as TSV or CSV.",
    )
    comparative_clade_stability.add_argument(
        "--json", action="store_true", help="Emit the clade-stability report as JSON."
    )
    _add_manifest_argument(comparative_clade_stability)
    comparative_posterior_pgls = comparative_subparsers.add_parser(
        "posterior-pgls",
        help="Fit one continuous-trait PGLS model across a posterior or bootstrap tree set.",
    )
    comparative_posterior_pgls.add_argument("tree_set", type=Path)
    comparative_posterior_pgls.add_argument("table", type=Path)
    comparative_posterior_pgls.add_argument("--response")
    comparative_posterior_pgls.add_argument("--predictors", nargs="+")
    comparative_posterior_pgls.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass + habitat'.",
    )
    comparative_posterior_pgls.add_argument("--taxon-column")
    comparative_posterior_pgls.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1 for each retained tree fit.",
    )
    comparative_posterior_pgls.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this leading fraction of the tree set before refitting.",
    )
    comparative_posterior_pgls.add_argument(
        "--significance-threshold",
        type=float,
        default=0.05,
        help="Treat coefficient p-values at or below this threshold as supported.",
    )
    comparative_posterior_pgls.add_argument(
        "--trees-out",
        type=Path,
        help="Write the per-tree posterior PGLS fit ledger as TSV or CSV.",
    )
    comparative_posterior_pgls.add_argument(
        "--coefficients-out",
        type=Path,
        help="Write the per-tree coefficient ledger as TSV or CSV.",
    )
    comparative_posterior_pgls.add_argument(
        "--summary-out",
        type=Path,
        help="Write the coefficient-distribution summary ledger as TSV or CSV.",
    )
    comparative_posterior_pgls.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior-tree PGLS report as JSON.",
    )
    _add_manifest_argument(comparative_posterior_pgls)
    comparative_model_selection = comparative_subparsers.add_parser(
        "model-selection",
        help="Rank competing comparative regression formulas on one shared taxon set.",
    )
    comparative_model_selection.add_argument("tree", type=Path)
    comparative_model_selection.add_argument("table", type=Path)
    comparative_model_selection.add_argument(
        "--formula",
        action="append",
        required=True,
        dest="formulas",
        help="Candidate comparative formula. Repeat this option once per model.",
    )
    comparative_model_selection.add_argument("--taxon-column")
    comparative_model_selection.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1. Binary-response model selection requires a numeric value.",
    )
    comparative_model_selection.add_argument(
        "--ranking-out",
        type=Path,
        help="Write the ranked comparative model table as TSV or CSV.",
    )
    comparative_model_selection.add_argument(
        "--pairwise-out",
        type=Path,
        help="Write nested-versus-non-nested pairwise comparison rows as TSV or CSV.",
    )
    comparative_model_selection.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write the shared-complete-case excluded-taxa ledger as TSV or CSV.",
    )
    comparative_model_selection.add_argument(
        "--json", action="store_true", help="Emit the model-selection report as JSON."
    )
    _add_manifest_argument(comparative_model_selection)
    comparative_multivariate = comparative_subparsers.add_parser(
        "multivariate",
        help="Fit shared-taxon comparative regressions across multiple response traits.",
    )
    comparative_multivariate.add_argument("tree", type=Path)
    comparative_multivariate.add_argument("table", type=Path)
    comparative_multivariate.add_argument("--responses", nargs="+", required=True)
    comparative_multivariate.add_argument("--predictors", nargs="+", required=True)
    comparative_multivariate.add_argument("--taxon-column")
    comparative_multivariate.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_multivariate.add_argument(
        "--covariance-out",
        type=Path,
        help="Write the residual covariance ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--correlation-out",
        type=Path,
        help="Write the residual correlation ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--associations-out",
        type=Path,
        help="Write the residual trait-association ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--coefficients-out",
        type=Path,
        help="Write the per-response coefficient ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--response-models-out",
        type=Path,
        help="Write the per-response model summary ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write the explicit excluded-taxa ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--json",
        action="store_true",
        help="Emit the multivariate regression report as JSON.",
    )
    _add_manifest_argument(comparative_multivariate)
    comparative_report = comparative_subparsers.add_parser(
        "report",
        help="Build an integrated comparative-method report.",
    )
    comparative_report.add_argument("tree", type=Path)
    comparative_report.add_argument("table", type=Path)
    comparative_report.add_argument("--response")
    comparative_report.add_argument("--predictors", nargs="+")
    comparative_report.add_argument("--formula")
    comparative_report.add_argument("--taxon-column")
    comparative_report.add_argument("--out", type=Path)
    comparative_report.add_argument(
        "--methods-summary-out",
        type=Path,
        help="Write reviewer-facing Markdown methods text for the comparative analysis.",
    )
    comparative_report.add_argument(
        "--out-dir",
        type=Path,
        help="Write a full comparative analysis package directory with HTML and reviewer TSV ledgers.",
    )
    comparative_report.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_report.add_argument(
        "--json", action="store_true", help="Emit the comparative report as JSON."
    )
    _add_manifest_argument(comparative_report)
    comparative_influence = comparative_subparsers.add_parser(
        "influence",
        help="Identify predictor terms and taxa driving one comparative result.",
    )
    comparative_influence.add_argument("tree", type=Path)
    comparative_influence.add_argument("table", type=Path)
    comparative_influence.add_argument("--response")
    comparative_influence.add_argument("--predictors", nargs="+")
    comparative_influence.add_argument("--formula")
    comparative_influence.add_argument("--taxon-column")
    comparative_influence.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_influence.add_argument(
        "--json", action="store_true", help="Emit the influence report as JSON."
    )
    _add_manifest_argument(comparative_influence)
    comparative_compare_trees = comparative_subparsers.add_parser(
        "compare-trees",
        help="Compare comparative results across two alternative trees.",
    )
    comparative_compare_trees.add_argument("left_tree", type=Path)
    comparative_compare_trees.add_argument("right_tree", type=Path)
    comparative_compare_trees.add_argument("table", type=Path)
    comparative_compare_trees.add_argument("--response")
    comparative_compare_trees.add_argument("--predictors", nargs="+")
    comparative_compare_trees.add_argument("--formula")
    comparative_compare_trees.add_argument("--taxon-column")
    comparative_compare_trees.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_compare_trees.add_argument(
        "--json", action="store_true", help="Emit the comparison as JSON."
    )
    _add_manifest_argument(comparative_compare_trees)
    comparative_compare_pruning = comparative_subparsers.add_parser(
        "compare-pruning",
        help="Compare comparative results before and after explicit pruning.",
    )
    comparative_compare_pruning.add_argument("tree", type=Path)
    comparative_compare_pruning.add_argument("table", type=Path)
    comparative_compare_pruning.add_argument("--response")
    comparative_compare_pruning.add_argument("--predictors", nargs="+")
    comparative_compare_pruning.add_argument("--formula")
    comparative_compare_pruning.add_argument("--drop-taxa", nargs="+")
    comparative_compare_pruning.add_argument("--keep-taxa", nargs="+")
    comparative_compare_pruning.add_argument("--taxon-column")
    comparative_compare_pruning.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_compare_pruning.add_argument(
        "--json", action="store_true", help="Emit the pruning comparison as JSON."
    )
    _add_manifest_argument(comparative_compare_pruning)

    add_ancestral_commands(subparsers)

    add_biogeography_commands(subparsers)
    add_host_association_commands(subparsers)
    add_ecological_niche_commands(subparsers)
    add_phylogeography_commands(subparsers)
    add_discrete_evolution_commands(subparsers)
    add_diversification_commands(subparsers)

    add_distance_commands(subparsers)
    add_tree_set_commands(subparsers)

    add_simulate_command(subparsers)

    add_benchmark_commands(subparsers)

    add_parity_command(subparsers)

    add_tree_inspection_commands(subparsers)

    normalize = subparsers.add_parser(
        get_command_spec("normalize").name, help=get_command_spec("normalize").summary
    )
    normalize.add_argument("tree", type=Path)
    normalize.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    normalize.add_argument("--out", required=True, type=Path)
    normalize.add_argument(
        "--json", action="store_true", help="Emit the normalization result as JSON."
    )
    _add_manifest_argument(normalize)

    normalize_taxa = subparsers.add_parser(
        get_command_spec("normalize-taxa").name,
        help=get_command_spec("normalize-taxa").summary,
    )
    normalize_taxa.add_argument("tree", type=Path)
    normalize_taxa.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    normalize_taxa.add_argument(
        "--policy", choices=("spaces-to-underscores",), required=True
    )
    normalize_taxa.add_argument("--out", required=True, type=Path)
    normalize_taxa.add_argument("--mapping-out", type=Path)
    normalize_taxa.add_argument(
        "--json", action="store_true", help="Emit the normalization result as JSON."
    )
    _add_manifest_argument(normalize_taxa)

    add_taxonomy_commands(subparsers)

    add_topology_commands(subparsers)

    add_compare_command(subparsers)

    add_annotate_command(subparsers)

    add_diagnose_command(subparsers)

    add_render_command(subparsers)

    add_evidence_command(subparsers)

    add_report_command(subparsers)

    add_demo_command(subparsers)

    adapter = subparsers.add_parser(
        get_command_spec("adapter").name, help=get_command_spec("adapter").summary
    )
    adapter_subparsers = adapter.add_subparsers(dest="adapter_command", required=True)
    adapter_inspect = adapter_subparsers.add_parser(
        "inspect", help="Report external engine version metadata."
    )
    adapter_inspect.add_argument(
        "engine_name", choices=("mafft", "trimal", "iqtree", "FastTree", "MrBayes")
    )
    adapter_inspect.add_argument("--executable", type=str)
    adapter_inspect.add_argument(
        "--json", action="store_true", help="Emit the adapter report as JSON."
    )
    _add_manifest_argument(adapter_inspect)
    adapter_report = adapter_subparsers.add_parser(
        "report", help="Render an HTML report from an engine workflow manifest."
    )
    adapter_report.add_argument("manifest_path", type=Path)
    adapter_report.add_argument("--out", required=True, type=Path)
    adapter_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(adapter_report)
    adapter_align = adapter_subparsers.add_parser(
        "align", help="Run multiple-sequence alignment on unaligned FASTA."
    )
    adapter_align.add_argument("input_path", type=Path)
    adapter_align.add_argument("--out", required=True, type=Path)
    adapter_align.add_argument("--executable", type=str)
    adapter_align.add_argument(
        "--mode",
        choices=list_mafft_alignment_modes(),
        default="auto",
        help="Select the named MAFFT alignment strategy.",
    )
    adapter_align.add_argument(
        "--codon-aware",
        action="store_true",
        help="Translate accepted coding nucleotide sequences to an amino-acid guide, align that guide, then back-translate codon triplets.",
    )
    adapter_align.add_argument(
        "--sequence-type",
        choices=("dna", "rna"),
        help="Declare the coding nucleotide type for codon-aware alignment when explicit forcing is needed.",
    )
    adapter_align.add_argument(
        "--genetic-code",
        default="1",
        help="Use an NCBI genetic code id or codon-table name for codon-aware validation and translation.",
    )
    adapter_align.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_align)
    _add_manifest_argument(adapter_align)
    adapter_trim = adapter_subparsers.add_parser(
        "trim", help="Run external alignment trimming."
    )
    adapter_trim.add_argument("input_path", type=Path)
    adapter_trim.add_argument("--out", required=True, type=Path)
    adapter_trim.add_argument(
        "--mode",
        choices=list_trimal_trimming_modes(),
        default="gap-threshold",
        help="Select the named trimAl trimming strategy.",
    )
    adapter_trim.add_argument("--gap-threshold", type=float, default=0.1)
    adapter_trim.add_argument("--executable", type=str)
    adapter_trim.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_trim)
    _add_manifest_argument(adapter_trim)
    adapter_model = adapter_subparsers.add_parser(
        "model-select", help="Run external sequence-model selection."
    )
    adapter_model.add_argument("input_path", type=Path)
    adapter_model.add_argument("--out-dir", required=True, type=Path)
    adapter_model.add_argument("--prefix", default="model-selection")
    adapter_model.add_argument(
        "--partitions",
        type=Path,
        help="Validate and apply a partition scheme for partitioned model selection.",
    )
    adapter_model.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_model.add_argument("--executable", type=str)
    adapter_model.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_model)
    _add_manifest_argument(adapter_model)
    adapter_ml = adapter_subparsers.add_parser(
        "infer-ml", help="Run maximum-likelihood tree inference."
    )
    adapter_ml.add_argument("input_path", type=Path)
    adapter_ml.add_argument("--out-dir", required=True, type=Path)
    adapter_ml.add_argument("--model", required=True)
    adapter_ml.add_argument("--prefix", default="maximum-likelihood")
    adapter_ml.add_argument(
        "--partitions",
        type=Path,
        help="Validate and apply a partition scheme for partitioned maximum-likelihood inference.",
    )
    adapter_ml.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_ml.add_argument("--executable", type=str)
    adapter_ml.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_ml)
    _add_manifest_argument(adapter_ml)
    adapter_bootstrap = adapter_subparsers.add_parser(
        "bootstrap", help="Run bootstrap support estimation."
    )
    adapter_bootstrap.add_argument("input_path", type=Path)
    adapter_bootstrap.add_argument("--out-dir", required=True, type=Path)
    adapter_bootstrap.add_argument("--model", required=True)
    adapter_bootstrap.add_argument("--replicates", type=int, default=1000)
    adapter_bootstrap.add_argument("--prefix", default="bootstrap-support")
    adapter_bootstrap.add_argument(
        "--partitions",
        type=Path,
        help="Validate and apply a partition scheme for partitioned bootstrap support estimation.",
    )
    adapter_bootstrap.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_bootstrap.add_argument("--executable", type=str)
    adapter_bootstrap.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_bootstrap)
    _add_manifest_argument(adapter_bootstrap)
    adapter_sh_alrt = adapter_subparsers.add_parser(
        "sh-alrt",
        help="Run combined sh-alrt and ultrafast bootstrap support estimation.",
    )
    adapter_sh_alrt.add_argument("input_path", type=Path)
    adapter_sh_alrt.add_argument("--out-dir", required=True, type=Path)
    adapter_sh_alrt.add_argument("--model", required=True)
    adapter_sh_alrt.add_argument("--alrt-replicates", type=int, default=1000)
    adapter_sh_alrt.add_argument("--bootstrap-replicates", type=int, default=1000)
    adapter_sh_alrt.add_argument("--prefix", default="sh-alrt-support")
    adapter_sh_alrt.add_argument(
        "--partitions",
        type=Path,
        help="Validate and apply a partition scheme for combined sh-alrt and ultrafast bootstrap support estimation.",
    )
    adapter_sh_alrt.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_sh_alrt.add_argument("--executable", type=str)
    adapter_sh_alrt.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_sh_alrt)
    _add_manifest_argument(adapter_sh_alrt)
    adapter_fasta_to_tree = adapter_subparsers.add_parser(
        "fasta-to-tree", help="Run alignment-to-tree inference from raw FASTA."
    )
    adapter_fasta_to_tree.add_argument("input_path", type=Path)
    adapter_fasta_to_tree.add_argument("--out-dir", required=True, type=Path)
    adapter_fasta_to_tree.add_argument("--prefix")
    adapter_fasta_to_tree.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_fasta_to_tree.add_argument("--mafft-executable", type=str)
    adapter_fasta_to_tree.add_argument(
        "--alignment-mode",
        choices=list_mafft_alignment_modes(),
        default="auto",
        help="Select the named MAFFT alignment strategy for the raw-input workflow.",
    )
    adapter_fasta_to_tree.add_argument("--trimal-executable", type=str)
    adapter_fasta_to_tree.add_argument(
        "--trimming-mode",
        choices=list_trimal_trimming_modes(),
        default="gap-threshold",
        help="Select the named trimAl trimming strategy for the aligned workflow.",
    )
    adapter_fasta_to_tree.add_argument("--iqtree-executable", type=str)
    adapter_fasta_to_tree.add_argument(
        "--iqtree-seed",
        type=int,
        default=1,
        help="Set the IQ-TREE random seed for deterministic model selection and support estimation.",
    )
    adapter_fasta_to_tree.add_argument(
        "--iqtree-threads",
        type=int,
        default=1,
        help="Set the IQ-TREE thread count used across model selection and inference.",
    )
    adapter_fasta_to_tree.add_argument("--trim-gap-threshold", type=float, default=0.1)
    adapter_fasta_to_tree.add_argument("--bootstrap-replicates", type=int, default=1000)
    adapter_fasta_to_tree.add_argument(
        "--normalize-identifiers",
        action="store_true",
        help="Normalize FASTA identifiers before alignment and resolve any collisions.",
    )
    adapter_fasta_to_tree.add_argument(
        "--remove-invalid-records",
        action="store_true",
        help="Remove empty or illegal sequence records before alignment.",
    )
    adapter_fasta_to_tree.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_fasta_to_tree)
    _add_manifest_argument(adapter_fasta_to_tree)
    adapter_consensus = adapter_subparsers.add_parser(
        "consensus", help="Build a consensus tree from bootstrap trees."
    )
    adapter_consensus.add_argument("input_path", type=Path)
    adapter_consensus.add_argument("--out-dir", required=True, type=Path)
    adapter_consensus.add_argument("--prefix", default="bootstrap-consensus")
    adapter_consensus.add_argument("--minimum-support", type=float, default=0.5)
    adapter_consensus.add_argument("--executable", type=str)
    adapter_consensus.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_consensus)
    _add_manifest_argument(adapter_consensus)
    add_inference_adapter_commands(adapter_subparsers)
    add_mrbayes_adapter_commands(adapter_subparsers)
    add_beast_adapter_commands(adapter_subparsers)
    add_bayesian_adapter_commands(adapter_subparsers)

    return parser


def run_command(args: Any, *, parser: argparse.ArgumentParser) -> int:
    """Run the selected command."""
    try:
        if args.command == "commands":
            _print_commands(output_format=args.format)
            return 0
        if args.command == "env":
            report = inspect_environment()
            outputs = _finalize_outputs(args, command="env", inputs=[])
            _print_result(
                build_command_result(
                    command="env",
                    inputs=[],
                    outputs=outputs,
                    metrics={"dependency_count": len(report.dependencies)},
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "phylo":
            return run_phylo_command(args)
        if args.command == "metadata":
            return run_metadata_command(args)
        if args.command == "traits":
            return run_traits_command(args)
        if args.command == "validate":
            return run_validate_command(args)
        if args.command == "prune":
            return run_prune_command(args)
        if args.command == "alignment":
            alignment_exit_code = run_alignment_command(args)
            if alignment_exit_code is not None:
                return alignment_exit_code
        if args.command == "comparative":
            continuous_exit_code = run_comparative_continuous_command(
                args,
                parser=parser,
            )
            if continuous_exit_code is not None:
                return continuous_exit_code
            evolution_exit_code = run_comparative_evolution_command(
                args,
                parser=parser,
            )
            if evolution_exit_code is not None:
                return evolution_exit_code
            if args.comparative_command == "clade-traits":
                report = summarize_clade_traits(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    minimum_clade_size=args.min_clade_size,
                    trait_kind=args.trait_kind,
                )
                if args.summary_out:
                    write_clade_trait_summary_table(args.summary_out, report)
                if args.clades_out:
                    write_clade_trait_clade_table(args.clades_out, report)
                if args.excluded_taxa_out:
                    write_clade_trait_exclusion_table(
                        args.excluded_taxa_out,
                        report,
                    )
                outputs = _finalize_outputs(
                    args, command="comparative", inputs=[args.tree, args.table]
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "tree_taxon_count": report.tree_taxon_count,
                            "analyzed_taxon_count": report.analyzed_taxon_count,
                            "excluded_taxon_count": len(report.excluded_taxa),
                            "trait_kind": report.trait_kind,
                            "clade_count": len(report.clade_rows),
                            "exceptional_clade_count": len(report.exceptional_clades),
                            "top_exceptional_clade": report.top_exceptional_clade,
                            "top_exceptionality_score": report.top_exceptionality_score,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "phylogenetic-residuals":
                report = summarize_phylogenetic_residuals(
                    args.tree,
                    args.table,
                    response=args.response,
                    predictor=args.predictor,
                    taxon_column=args.taxon_column,
                    method=args.method,
                )
                if args.summary_out:
                    write_phylogenetic_residual_summary_table(args.summary_out, report)
                if args.residuals_out:
                    write_phylogenetic_residual_taxon_table(args.residuals_out, report)
                if args.coefficients_out:
                    write_phylogenetic_residual_coefficient_table(
                        args.coefficients_out,
                        report,
                    )
                if args.excluded_taxa_out:
                    write_phylogenetic_residual_exclusion_table(
                        args.excluded_taxa_out,
                        report,
                    )
                outputs = _finalize_outputs(
                    args, command="comparative", inputs=[args.tree, args.table]
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "tree_taxon_count": report.tree_taxon_count,
                            "analyzed_taxon_count": report.analyzed_taxon_count,
                            "excluded_taxon_count": len(report.excluded_taxa),
                            "method": report.method,
                            "outlier_count": len(report.outlier_taxa),
                            "top_outlier_taxon": (
                                None
                                if not report.taxon_rows
                                else max(
                                    report.taxon_rows,
                                    key=lambda row: row.abs_standardized_residual,
                                ).taxon
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "phylogenetic-anova":
                report = summarize_phylogenetic_anova(
                    args.tree,
                    args.table,
                    response=args.response,
                    group=args.group,
                    taxon_column=args.taxon_column,
                    simulations=args.simulations,
                    seed=args.seed,
                )
                if args.summary_out:
                    write_phylogenetic_anova_summary_table(args.summary_out, report)
                if args.groups_out:
                    write_phylogenetic_anova_group_table(args.groups_out, report)
                if args.pairwise_out:
                    write_phylogenetic_anova_pairwise_table(args.pairwise_out, report)
                if args.simulations_out:
                    write_phylogenetic_anova_simulation_table(
                        args.simulations_out,
                        report,
                    )
                if args.excluded_taxa_out:
                    write_phylogenetic_anova_exclusion_table(
                        args.excluded_taxa_out,
                        report,
                    )
                outputs = _finalize_outputs(
                    args, command="comparative", inputs=[args.tree, args.table]
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "tree_taxon_count": report.tree_taxon_count,
                            "analyzed_taxon_count": report.analyzed_taxon_count,
                            "excluded_taxon_count": len(report.excluded_taxa),
                            "group_count": report.group_count,
                            "simulation_count": report.simulation_count,
                            "p_value": report.p_value,
                            "low_sample_group_count": report.low_sample_group_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "trait-outliers":
                report = summarize_trait_outliers(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                )
                if args.summary_out:
                    write_trait_outlier_summary_table(args.summary_out, report)
                if args.outliers_out:
                    write_trait_outlier_taxon_table(args.outliers_out, report)
                if args.excluded_taxa_out:
                    write_trait_outlier_exclusion_table(
                        args.excluded_taxa_out,
                        report,
                    )
                outputs = _finalize_outputs(
                    args, command="comparative", inputs=[args.tree, args.table]
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "tree_taxon_count": report.tree_taxon_count,
                            "analyzed_taxon_count": report.analyzed_taxon_count,
                            "excluded_taxon_count": len(report.excluded_taxa),
                            "selected_model": report.selected_model,
                            "outlier_count": len(report.outlier_taxa),
                            "top_outlier_taxon": report.top_outlier_taxon,
                            "top_abs_standardized_residual": (
                                report.top_abs_standardized_residual
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "trait-imputation":
                report = summarize_trait_imputation(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                )
                if args.summary_out:
                    write_trait_imputation_summary_table(args.summary_out, report)
                if args.imputations_out:
                    write_trait_imputation_table(args.imputations_out, report)
                if args.holdout_out:
                    write_trait_imputation_holdout_table(args.holdout_out, report)
                if args.excluded_taxa_out:
                    write_trait_imputation_exclusion_table(
                        args.excluded_taxa_out,
                        report,
                    )
                outputs = _finalize_outputs(
                    args, command="comparative", inputs=[args.tree, args.table]
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "tree_taxon_count": report.tree_taxon_count,
                            "observed_taxon_count": report.observed_taxon_count,
                            "imputed_taxon_count": len(report.imputation_rows),
                            "excluded_taxon_count": len(report.excluded_taxa),
                            "holdout_validation_status": report.holdout_validation_status,
                            "holdout_count": len(report.holdout_rows),
                            "holdout_mean_absolute_error": (
                                report.holdout_mean_absolute_error
                            ),
                            "holdout_interval_coverage": (
                                report.holdout_interval_coverage
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            lambda_value: float | str
            if hasattr(args, "lambda_value"):
                if args.lambda_value == "estimate":
                    lambda_value = "estimate"
                else:
                    lambda_value = float(args.lambda_value)
            else:
                lambda_value = "estimate"
            if args.comparative_command == "maturity":
                report = assess_comparative_method_maturity(
                    args.tree,
                    args.table,
                    response=args.response,
                    predictors=list(args.predictors or []),
                    formula=args.formula,
                    taxon_column=args.taxon_column,
                    lambda_value=lambda_value,
                )
                outputs = _finalize_outputs(
                    args, command="comparative", inputs=[args.tree, args.table]
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "selected_model": report.selected_model,
                            "residual_surface_count": len(report.residual_diagnostics),
                            "influential_taxa": len(
                                report.sensitivity.influential_taxa
                            ),
                            "reference_validation_passed": report.reference_validation_passed,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "covariance-audit":
                resolved_lambda: float | str
                if args.lambda_value == "estimate":
                    resolved_lambda = "estimate"
                else:
                    resolved_lambda = float(args.lambda_value)
                resolved_alpha: float | str
                if args.alpha == "estimate":
                    resolved_alpha = "estimate"
                else:
                    resolved_alpha = float(args.alpha)
                report = summarize_comparative_covariance_audit(
                    args.tree,
                    args.table,
                    analysis=args.analysis,
                    trait=args.trait,
                    response=args.response,
                    predictors=list(args.predictors or []),
                    formula=args.formula,
                    taxon_column=args.taxon_column,
                    lambda_value=resolved_lambda,
                    alpha=resolved_alpha,
                )
                outputs: list[Path | str] = []
                if args.summary_out is not None:
                    outputs.append(
                        write_comparative_covariance_audit_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                if args.candidates_out is not None:
                    outputs.append(
                        write_comparative_covariance_audit_candidate_table(
                            args.candidates_out,
                            report,
                        )
                    )
                if args.excluded_taxa_out is not None:
                    outputs.append(
                        write_comparative_covariance_audit_excluded_taxa_table(
                            args.excluded_taxa_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="comparative",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "analysis": report.analysis,
                            "covariance_model": report.covariance_model,
                            "matrix_dimension": report.matrix_dimension,
                            "matrix_rank": report.matrix_rank,
                            "condition_number": report.condition_number,
                            "fit_strategy": report.fit_strategy,
                            "singular": report.singular,
                            "near_singular": report.near_singular,
                            "matched_taxon_count": len(report.matched_taxa),
                            "analysis_taxon_count": len(report.analysis_taxa),
                            "duplicate_tree_taxon_count": len(
                                report.duplicate_tree_taxa
                            ),
                            "duplicate_trait_taxon_count": len(
                                report.duplicate_trait_taxa
                            ),
                            "candidate_row_count": len(report.candidate_rows),
                            "blocker_count": len(report.blockers),
                            "warning_count": len(report.warnings),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "multiple-testing":
                report = run_pgls_multiple_testing(
                    args.tree,
                    args.table,
                    responses=list(args.responses),
                    predictors=list(args.predictors),
                    taxon_column=args.taxon_column,
                    lambda_value=lambda_value,
                )
                outputs = _finalize_outputs(
                    args, command="comparative", inputs=[args.tree, args.table]
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "response_count": len(report.responses),
                            "test_count": len(report.rows),
                            "family_size": report.family_size,
                            "raw_significant_count": report.raw_significant_count,
                            "significant_count": sum(
                                1 for row in report.rows if row.significant
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "logistic":
                report = summarize_phylogenetic_logistic(
                    args.tree,
                    args.table,
                    response=args.response,
                    predictors=list(args.predictors or []),
                    formula=args.formula,
                    taxon_column=args.taxon_column,
                    lambda_value=float(args.lambda_value),
                )
                outputs: list[Path | str] = []
                if args.coefficients_out is not None:
                    outputs.append(
                        write_phylogenetic_logistic_coefficient_table(
                            args.coefficients_out,
                            report,
                        )
                    )
                if args.fitted_out is not None:
                    outputs.append(
                        write_phylogenetic_logistic_fitted_table(
                            args.fitted_out,
                            report,
                        )
                    )
                if args.excluded_taxa_out is not None:
                    outputs.append(
                        write_phylogenetic_logistic_excluded_taxa_table(
                            args.excluded_taxa_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="comparative",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "success_count": report.success_count,
                            "failure_count": report.failure_count,
                            "coefficient_count": len(report.coefficients),
                            "fitted_row_count": len(report.fitted_rows),
                            "lambda_value": report.lambda_value,
                            "approximation_method": report.approximation_method,
                            "converged": report.converged,
                            "iteration_count": report.iteration_count,
                            "binomial_log_likelihood": report.binomial_log_likelihood,
                            "separation_detected": report.separation_detected,
                            "warning_count": len(report.warnings)
                            + len(method_tier_warnings(report.method_tier)),
                            "coefficient_inference_distribution": (
                                report.coefficients[0].inference_distribution
                                if report.coefficients
                                else None
                            ),
                            **method_tier_metrics(report.method_tier),
                        },
                        warnings=method_tier_warnings(report.method_tier)
                        + [warning.message for warning in report.warnings],
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "model-selection":
                report = compare_comparative_regression_models(
                    args.tree,
                    args.table,
                    formulas=list(args.formulas),
                    taxon_column=args.taxon_column,
                    lambda_value=lambda_value,
                )
                outputs: list[Path | str] = []
                if args.ranking_out is not None:
                    outputs.append(
                        write_comparative_regression_model_ranking_table(
                            args.ranking_out,
                            report,
                        )
                    )
                if args.pairwise_out is not None:
                    outputs.append(
                        write_comparative_regression_pairwise_table(
                            args.pairwise_out,
                            report,
                        )
                    )
                if args.excluded_taxa_out is not None:
                    outputs.append(
                        write_comparative_regression_excluded_taxa_table(
                            args.excluded_taxa_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="comparative",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                selected_row = next(row for row in report.rows if row.selected)
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "response": report.response,
                            "model_family": report.model_family,
                            "model_count": len(report.rows),
                            "analysis_taxon_count": len(report.analysis_taxa),
                            "excluded_taxon_count": len(report.excluded_taxa),
                            "pairwise_comparison_count": len(report.pairwise_rows),
                            "best_formula": report.best_formula,
                            "selected_criterion": report.selected_criterion,
                            "selected_log_likelihood": selected_row.log_likelihood,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "clade-residuals":
                report = analyze_comparative_residual_clades(
                    args.tree,
                    args.table,
                    response=args.response,
                    predictors=list(args.predictors or []),
                    formula=args.formula,
                    taxon_column=args.taxon_column,
                    lambda_value=lambda_value,
                )
                outputs: list[Path | str] = []
                if args.taxa_out is not None:
                    outputs.append(
                        write_comparative_residual_taxon_table(
                            args.taxa_out,
                            report,
                        )
                    )
                if args.clades_out is not None:
                    outputs.append(
                        write_comparative_residual_clade_table(
                            args.clades_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="comparative",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                top_clade = (
                    min(report.clade_rows, key=lambda row: row.rank).clade_id
                    if report.clade_rows
                    else None
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "model_family": report.model_family,
                            "taxon_count": len(report.taxon_rows),
                            "clade_count": len(report.clade_rows),
                            "residual_heavy_clade_count": len(
                                report.residual_heavy_clades
                            ),
                            "top_influential_clade": top_clade,
                            "standardized_residual_method": (
                                report.standardized_residual_method
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "clade-stability":
                report = analyze_comparative_clade_stability(
                    args.tree,
                    args.table,
                    response=args.response,
                    predictors=list(args.predictors or []),
                    formula=args.formula,
                    taxon_column=args.taxon_column,
                    lambda_value=lambda_value,
                )
                outputs: list[Path | str] = []
                if args.clades_out is not None:
                    outputs.append(
                        write_comparative_clade_stability_table(
                            args.clades_out,
                            report,
                        )
                    )
                if args.terms_out is not None:
                    outputs.append(
                        write_comparative_clade_coefficient_change_table(
                            args.terms_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="comparative",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                top_clade = (
                    min(
                        (
                            row
                            for row in report.clade_rows
                            if row.fit_status == "fit" and row.rank > 0
                        ),
                        key=lambda row: row.rank,
                    ).clade_id
                    if any(
                        row.fit_status == "fit" and row.rank > 0
                        for row in report.clade_rows
                    )
                    else None
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "model_family": report.model_family,
                            "baseline_taxon_count": len(report.baseline_taxa),
                            "baseline_term_count": report.baseline_term_count,
                            "candidate_clade_count": report.candidate_clade_count,
                            "blocked_clade_count": report.blocked_clade_count,
                            "coefficient_change_row_count": len(
                                report.coefficient_rows
                            ),
                            "top_influential_clade": top_clade,
                            "major_clade_fraction": report.major_clade_fraction,
                            "minimum_major_clade_size": report.minimum_major_clade_size,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "posterior-pgls":
                report = run_posterior_tree_pgls(
                    args.tree_set,
                    args.table,
                    response=args.response,
                    predictors=list(args.predictors or []),
                    formula=args.formula,
                    taxon_column=args.taxon_column,
                    lambda_value=lambda_value,
                    burnin_fraction=args.burnin_fraction,
                    significance_threshold=args.significance_threshold,
                )
                outputs: list[Path | str] = []
                if args.trees_out is not None:
                    outputs.append(
                        write_posterior_tree_pgls_tree_table(
                            args.trees_out,
                            report,
                        )
                    )
                if args.coefficients_out is not None:
                    outputs.append(
                        write_posterior_tree_pgls_coefficient_table(
                            args.coefficients_out,
                            report,
                        )
                    )
                if args.summary_out is not None:
                    outputs.append(
                        write_posterior_tree_pgls_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="comparative",
                    inputs=[args.tree_set, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree_set, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "total_tree_count": report.total_tree_count,
                            "burnin_tree_count": report.burnin_tree_count,
                            "kept_tree_count": report.kept_tree_count,
                            "analysis_taxon_count": len(report.analysis_taxa),
                            "rooted_topology_count": report.rooted_topology_count,
                            "unrooted_topology_count": report.unrooted_topology_count,
                            "tree_fit_row_count": len(report.tree_rows),
                            "coefficient_row_count": len(report.coefficient_rows),
                            "coefficient_summary_count": len(
                                report.coefficient_summaries
                            ),
                            "stable_supported_term_count": sum(
                                row.conclusion_stability == "stable_supported"
                                for row in report.coefficient_summaries
                            ),
                            "direction_conflict_term_count": sum(
                                row.conclusion_stability == "direction_conflict"
                                for row in report.coefficient_summaries
                            ),
                            "lambda_mode": report.lambda_mode,
                            "significance_threshold": report.significance_threshold,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "multivariate":
                report = run_multivariate_comparative_regression(
                    args.tree,
                    args.table,
                    responses=list(args.responses),
                    predictors=list(args.predictors),
                    taxon_column=args.taxon_column,
                    lambda_value=lambda_value,
                )
                if args.covariance_out is not None:
                    write_multivariate_residual_covariance_table(
                        args.covariance_out, report
                    )
                if args.correlation_out is not None:
                    write_multivariate_residual_correlation_table(
                        args.correlation_out, report
                    )
                if args.associations_out is not None:
                    write_multivariate_residual_association_table(
                        args.associations_out, report
                    )
                if args.coefficients_out is not None:
                    write_multivariate_response_coefficient_table(
                        args.coefficients_out, report
                    )
                if args.response_models_out is not None:
                    write_multivariate_response_model_table(
                        args.response_models_out, report
                    )
                if args.excluded_taxa_out is not None:
                    write_multivariate_excluded_taxa_table(
                        args.excluded_taxa_out, report
                    )
                outputs = _finalize_outputs(
                    args, command="comparative", inputs=[args.tree, args.table]
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "response_count": len(report.responses),
                            "predictor_count": len(report.predictors),
                            "analysis_taxa": len(report.analysis_taxa),
                            "excluded_taxa": len(report.excluded_taxa),
                            "residual_covariance_response_count": (
                                report.covariance_diagnostics.response_count
                            ),
                            "residual_covariance_matrix_rank": (
                                report.covariance_diagnostics.matrix_rank
                            ),
                            "residual_covariance_condition_number": (
                                report.covariance_diagnostics.condition_number
                            ),
                            "residual_covariance_singular": (
                                report.covariance_diagnostics.is_singular
                            ),
                            "residual_covariance_near_singular": (
                                report.covariance_diagnostics.is_near_singular
                            ),
                            "residual_covariance_row_count": len(
                                report.covariance_rows
                            ),
                            "residual_correlation_row_count": len(
                                report.correlation_rows
                            ),
                            "residual_association_count": len(report.association_rows),
                            "response_model_count": len(report.response_model_rows),
                            "coefficient_row_count": len(report.coefficient_rows),
                            "warning_count": len(report.warnings),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "report":
                package_result = None
                if args.out_dir is not None:
                    package_result = build_comparative_report_package(
                        args.tree,
                        args.table,
                        out_dir=args.out_dir,
                        response=args.response,
                        predictors=list(args.predictors or []),
                        formula=args.formula,
                        taxon_column=args.taxon_column,
                        lambda_value=lambda_value,
                    )
                    report = package_result.report
                else:
                    report = build_comparative_method_report(
                        args.tree,
                        args.table,
                        response=args.response,
                        predictors=list(args.predictors or []),
                        formula=args.formula,
                        taxon_column=args.taxon_column,
                        lambda_value=lambda_value,
                    )
                if args.out is not None:
                    write_comparative_method_report(args.out, report)
                if args.methods_summary_out is not None:
                    write_comparative_methods_summary_text(
                        args.methods_summary_out, report
                    )
                output_paths: list[Path | str] = [args.out] if args.out else []
                if args.methods_summary_out is not None:
                    output_paths.append(args.methods_summary_out)
                if package_result is not None:
                    output_paths.extend(
                        [
                            package_result.report_path,
                            package_result.methods_summary_path,
                            package_result.reviewer_audit_checklist_path,
                            package_result.summary_table_path,
                            package_result.coefficient_table_path,
                            package_result.residual_table_path,
                            package_result.signal_table_path,
                            package_result.model_comparison_table_path,
                            package_result.interpretation_table_path,
                            package_result.audit_table_path,
                            package_result.contrast_table_path,
                            package_result.manifest_path,
                        ]
                    )
                outputs = _finalize_outputs(
                    args,
                    command="comparative",
                    inputs=[args.tree, args.table],
                    outputs=output_paths,
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "taxon_count": report.snapshot.pgls_model.taxon_count,
                            "selected_model": report.snapshot.model_comparison.better_model,
                            "audit_row_count": len(report.snapshot.audit_rows),
                            "excluded_taxa": len(
                                report.snapshot.pgls_inputs.formula_audit.excluded_taxa
                            ),
                            "limitation_count": len(report.snapshot.limitations),
                            "coefficient_count": len(
                                report.snapshot.pgls_model.coefficients
                            ),
                            "package_output_count": len(outputs),
                        },
                        data=report if package_result is None else package_result,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "influence":
                report = build_trait_influence_report(
                    args.tree,
                    args.table,
                    response=args.response,
                    predictors=list(args.predictors or []),
                    formula=args.formula,
                    taxon_column=args.taxon_column,
                    lambda_value=lambda_value,
                )
                outputs = _finalize_outputs(
                    args, command="comparative", inputs=[args.tree, args.table]
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "predictor_count": len(report.predictor_rows),
                            "taxon_count": len(report.taxon_rows),
                            "top_predictor_terms": len(report.top_predictor_terms),
                            "top_taxa": len(report.top_taxa),
                            "selected_model": report.selected_model,
                        },
                        warnings=report.warnings,
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "compare-trees":
                report = compare_comparative_results_across_trees(
                    args.left_tree,
                    args.right_tree,
                    args.table,
                    response=args.response,
                    predictors=list(args.predictors or []),
                    formula=args.formula,
                    taxon_column=args.taxon_column,
                    lambda_value=lambda_value,
                )
                outputs = _finalize_outputs(
                    args,
                    command="comparative",
                    inputs=[args.left_tree, args.right_tree, args.table],
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.left_tree, args.right_tree, args.table],
                        outputs=outputs,
                        metrics={
                            "coefficient_delta_count": len(report.coefficient_deltas),
                            "sign_changed_terms": len(report.sign_changed_terms),
                            "conclusion_changed": report.conclusion_changed,
                            "left_selected_model": report.left_selected_model,
                            "right_selected_model": report.right_selected_model,
                        },
                        warnings=report.warnings,
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "compare-pruning":
                report = compare_comparative_results_across_pruning(
                    args.tree,
                    args.table,
                    response=args.response,
                    predictors=list(args.predictors or []),
                    formula=args.formula,
                    drop_taxa=list(args.drop_taxa or []),
                    keep_taxa=list(args.keep_taxa or []),
                    taxon_column=args.taxon_column,
                    lambda_value=lambda_value,
                )
                outputs = _finalize_outputs(
                    args, command="comparative", inputs=[args.tree, args.table]
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "baseline_taxa": len(report.baseline_taxa),
                            "pruned_taxa": len(report.pruned_taxa),
                            "dropped_taxa": len(report.dropped_taxa),
                            "sign_changed_terms": len(report.sign_changed_terms),
                            "conclusion_changed": report.conclusion_changed,
                        },
                        warnings=report.warnings,
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "brownian-pgls":
                report = summarize_brownian_covariance_pgls(
                    args.tree,
                    args.table,
                    response=args.response,
                    predictors=list(args.predictors or []),
                    formula=args.formula,
                    taxon_column=args.taxon_column,
                )
                outputs: list[Path | str] = []
                if args.covariance_out is not None:
                    outputs.append(
                        write_brownian_covariance_table(
                            args.covariance_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="comparative",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "predictor_count": len(report.model.predictors),
                            "coefficient_count": len(report.model.coefficients),
                            "covariance_row_count": len(report.rows),
                            "lambda_value": report.model.lambda_value,
                            "covariance_model": "brownian-shared-path",
                            "tree_is_ultrametric": report.tree_is_ultrametric,
                            "minimum_root_to_tip_depth": (
                                report.minimum_root_to_tip_depth
                            ),
                            "maximum_root_to_tip_depth": (
                                report.maximum_root_to_tip_depth
                            ),
                            "raw_log_determinant": report.raw_log_determinant,
                            "positive_definite_before_stabilization": (
                                report.positive_definite_before_stabilization
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.comparative_command == "ou-pgls":
                report = summarize_ou_covariance_pgls(
                    args.tree,
                    args.table,
                    response=args.response,
                    predictors=list(args.predictors or []),
                    formula=args.formula,
                    taxon_column=args.taxon_column,
                    alpha=args.alpha,
                )
                outputs: list[Path | str] = []
                if args.covariance_out is not None:
                    outputs.append(
                        write_ou_covariance_table(args.covariance_out, report)
                    )
                if args.alpha_profile_out is not None:
                    outputs.append(
                        write_ou_alpha_profile_table(args.alpha_profile_out, report)
                    )
                outputs = _finalize_outputs(
                    args,
                    command="comparative",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="comparative",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "predictor_count": len(report.model.predictors),
                            "coefficient_count": len(report.model.coefficients),
                            "covariance_row_count": len(report.rows),
                            "alpha": report.alpha,
                            "alpha_estimation_mode": report.alpha_estimation_mode,
                            "alpha_profile_point_count": len(report.alpha_profile_rows),
                            "alpha_lower_95_confidence_interval": (
                                report.lower_95_confidence_interval
                            ),
                            "alpha_upper_95_confidence_interval": (
                                report.upper_95_confidence_interval
                            ),
                            "covariance_model": "ou-stationary-root",
                            "tree_is_ultrametric": report.tree_is_ultrametric,
                            "raw_log_determinant": report.raw_log_determinant,
                            "positive_definite_before_stabilization": (
                                report.positive_definite_before_stabilization
                            ),
                            "log_likelihood": report.model.log_likelihood,
                            "aic": report.model.aic,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            input_report = inspect_pgls_inputs(
                args.tree,
                args.table,
                response=args.response,
                predictors=list(args.predictors or []),
                formula=args.formula,
                taxon_column=args.taxon_column,
            )
            report = run_pgls(
                args.tree,
                args.table,
                response=args.response,
                predictors=list(args.predictors or []),
                formula=args.formula,
                taxon_column=args.taxon_column,
                lambda_value=lambda_value,
            )
            categorical_contrasts = summarize_pgls_categorical_contrasts(
                args.tree,
                args.table,
                response=args.response,
                predictors=list(args.predictors or []),
                formula=args.formula,
                taxon_column=args.taxon_column,
                lambda_value=lambda_value,
            )
            interaction_coefficients = summarize_pgls_interaction_coefficients(
                args.tree,
                args.table,
                response=args.response,
                predictors=list(args.predictors or []),
                formula=args.formula,
                taxon_column=args.taxon_column,
                lambda_value=lambda_value,
            )
            outputs: list[Path | str] = []
            if args.model_matrix_out is not None:
                outputs.append(
                    write_pgls_model_matrix_table(
                        args.model_matrix_out,
                        input_report.model_matrix,
                    )
                )
            if args.categorical_contrasts_out is not None:
                outputs.append(
                    write_pgls_categorical_contrast_table(
                        args.categorical_contrasts_out,
                        categorical_contrasts,
                    )
                )
            if args.interaction_coefficients_out is not None:
                outputs.append(
                    write_pgls_interaction_coefficient_table(
                        args.interaction_coefficients_out,
                        interaction_coefficients,
                    )
                )
            if args.lambda_profile_out is not None:
                outputs.append(
                    write_pgls_lambda_profile_table(
                        args.lambda_profile_out,
                        report.lambda_fit,
                    )
                )
            outputs = _finalize_outputs(
                args,
                command="comparative",
                inputs=[args.tree, args.table],
                outputs=outputs,
            )
            _print_result(
                build_command_result(
                    command="comparative",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                    warnings=input_report.warnings,
                    metrics={
                        "taxon_count": report.taxon_count,
                        "predictor_count": len(report.predictors),
                        "coefficient_count": len(report.coefficients),
                        "confidence_interval_count": len(report.coefficients),
                        "categorical_contrast_predictor_count": (
                            categorical_contrasts.categorical_predictor_count
                        ),
                        "categorical_contrast_row_count": len(
                            categorical_contrasts.rows
                        ),
                        "interaction_term_count": (
                            interaction_coefficients.interaction_term_count
                        ),
                        "interaction_coefficient_row_count": len(
                            interaction_coefficients.rows
                        ),
                        "intercept_included": input_report.formula.include_intercept,
                        "model_matrix_row_count": input_report.model_matrix.row_count,
                        "model_matrix_column_count": len(
                            input_report.model_matrix.encoded_columns
                        ),
                        "residual_degrees_of_freedom": (
                            report.coefficients[0].degrees_of_freedom
                            if report.coefficients
                            else 0
                        ),
                        "coefficient_inference_distribution": (
                            report.coefficients[0].inference_distribution
                            if report.coefficients
                            else None
                        ),
                        "encoded_predictor_count": len(
                            input_report.model_matrix.encoded_columns
                        )
                        - (1 if input_report.formula.include_intercept else 0),
                        "categorical_predictor_count": len(
                            input_report.categorical_predictors
                        ),
                        "transformed_term_count": len(
                            input_report.formula_audit.transformed_terms
                        ),
                        "lambda_value": report.lambda_value,
                        "lambda_estimation_mode": report.lambda_fit.mode,
                        "lambda_profile_point_count": len(
                            report.lambda_fit.profile_rows
                        ),
                        "lambda_lower_95_confidence_interval": (
                            report.lambda_fit.lower_95_confidence_interval
                        ),
                        "lambda_upper_95_confidence_interval": (
                            report.lambda_fit.upper_95_confidence_interval
                        ),
                        "aic": report.aic,
                        "r_squared": report.r_squared,
                    },
                    data={
                        "inputs": input_report,
                        "model": report,
                        "categorical_contrasts": categorical_contrasts,
                        "interaction_coefficients": interaction_coefficients,
                    },
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "ancestral":
            return run_ancestral_command(args, parser=parser)
        if args.command == "biogeography":
            return run_biogeography_command(args)
        if args.command == "host-association":
            return run_host_association_command(args)
        if args.command == "ecological-niche":
            return run_ecological_niche_command(args)
        if args.command == "phylogeography":
            return run_phylogeography_command(args)
        if args.command == "discrete-evolution":
            return run_discrete_evolution_command(args)
        if args.command == "diversification":
            return run_diversification_command(args)
        if args.command == "distance":
            return run_distance_command(args)
        if args.command == "tree-set":
            return run_tree_set_command(args)
        if args.command == "simulate":
            return run_simulate_command(args, parser=parser)
        if args.command == "benchmark":
            return run_benchmark_command(args)
        if args.command == "parity":
            return run_parity_command(args)
        if args.command == "inspect":
            return run_inspect_command(args)
        if args.command == "normalize":
            tree = load_tree(args.tree, source_format=args.format)
            output_path = write_newick(args.out, tree)
            outputs = _finalize_outputs(
                args, command="normalize", inputs=[args.tree], outputs=[output_path]
            )
            if args.json:
                _print_result(
                    build_command_result(
                        command="normalize",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={"tip_count": tree.tip_count},
                        data={
                            "source_format": tree.source_format,
                            "output_format": "newick",
                        },
                    ),
                    json_output=True,
                )
            else:
                print(output_path)
            return 0
        if args.command == "normalize-taxa":
            tree = load_tree(args.tree, source_format=args.format)
            normalized_tree, report = normalize_tree_taxa(tree, policy=args.policy)
            output_path = write_newick(args.out, normalized_tree)
            mapping_path = args.mapping_out or args.out.with_suffix(
                f"{args.out.suffix}.mapping.tsv"
            )
            write_taxon_mapping(mapping_path, report.renamed_taxa)
            outputs = _finalize_outputs(
                args,
                command="normalize-taxa",
                inputs=[args.tree],
                outputs=[output_path, mapping_path],
            )
            if args.json:
                _print_result(
                    build_command_result(
                        command="normalize-taxa",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={"renamed_taxa": len(report.renamed_taxa)},
                        data=report,
                    ),
                    json_output=True,
                )
            else:
                print(output_path)
            return 0
        if args.command == "taxonomy":
            return run_taxonomy_command(args)
        if args.command == "topology":
            return run_topology_command(args)
        if args.command == "diagnose":
            return run_diagnose_command(args, parser=parser)
        if args.command == "compare":
            return run_compare_command(args, parser=parser)
        if args.command == "annotate":
            return run_annotate_command(args)
        if args.command == "render":
            return run_render_command(args)
        if args.command == "evidence":
            return run_evidence_command(args)
        if args.command == "demo":
            return run_demo_command(args)
        if args.command == "report":
            return run_report_command(args)
        if args.command == "adapter":
            if args.adapter_command == "inspect":
                executable = args.executable or args.engine_name
                report = read_engine_version(
                    args.engine_name,
                    executable,
                    version_args=_adapter_version_args(args.engine_name),
                )
                outputs = _finalize_outputs(
                    args, command="adapter", inputs=[args.engine_name]
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.engine_name],
                        outputs=outputs,
                        metrics={"version_line_count": len(report.text.splitlines())},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "report":
                report = render_inference_workflow_report(
                    manifest_path=args.manifest_path, out_path=args.out
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.manifest_path],
                    outputs=[report.output_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.manifest_path],
                        outputs=outputs,
                        metrics={"warning_count": report.warning_count},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "align":
                if args.codon_aware:
                    report = run_codon_aware_multiple_sequence_alignment(
                        args.input_path,
                        args.out,
                        executable=args.executable or "mafft",
                        mode=args.mode,
                        sequence_type=args.sequence_type,
                        genetic_code=args.genetic_code,
                        resume=args.resume,
                        timeout_seconds=args.timeout_seconds,
                        incomplete_run_policy=args.incomplete_run_policy,
                    )
                    outputs = _finalize_outputs(
                        args,
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=[*report.output_paths.values(), report.manifest_path],
                    )
                    _print_result(
                        build_command_result(
                            command="adapter",
                            inputs=[args.input_path],
                            outputs=outputs,
                            warnings=report.warnings,
                            metrics={
                                "mode": args.mode,
                                "codon_aware": True,
                                "sequence_type": report.sequence_type,
                                "genetic_code_id": report.genetic_code_id,
                                "accepted_sequence_count": report.accepted_sequence_count,
                                "excluded_sequence_count": len(
                                    report.excluded_sequences
                                ),
                                "invalid_codon_sequence_count": report.invalid_codon_sequence_count,
                                "terminal_stop_sequence_count": report.terminal_stop_sequence_count,
                                "resumed": report.resumed,
                                "timeout_seconds": report.run.timeout_seconds,
                            },
                            data=report,
                        ),
                        json_output=args.json,
                    )
                    return 0
                report = run_multiple_sequence_alignment(
                    args.input_path,
                    args.out,
                    executable=args.executable or "mafft",
                    mode=args.mode,
                    resume=args.resume,
                    timeout_seconds=args.timeout_seconds,
                    incomplete_run_policy=args.incomplete_run_policy,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[*report.output_paths.values(), report.manifest_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.run.warning_lines,
                        metrics={
                            "mode": args.mode,
                            "codon_aware": False,
                            "warning_count": len(report.run.warning_lines),
                            "resumed": report.resumed,
                            "timeout_seconds": report.run.timeout_seconds,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "trim":
                report = run_alignment_trimming(
                    args.input_path,
                    args.out,
                    executable=args.executable or "trimal",
                    mode=args.mode,
                    gap_threshold=args.gap_threshold,
                    resume=args.resume,
                    timeout_seconds=args.timeout_seconds,
                    incomplete_run_policy=args.incomplete_run_policy,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[*report.output_paths.values(), report.manifest_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.run.warning_lines,
                        metrics={
                            "mode": args.mode,
                            "warning_count": len(report.run.warning_lines),
                            "retained_site_count": (
                                None
                                if report.trimming_summary is None
                                else report.trimming_summary.retained_site_count
                            ),
                            "removed_site_count": (
                                None
                                if report.trimming_summary is None
                                else report.trimming_summary.removed_site_count
                            ),
                            "input_gap_percentage": (
                                None
                                if report.trimming_summary is None
                                else report.trimming_summary.input_gap_percentage
                            ),
                            "trimmed_gap_percentage": (
                                None
                                if report.trimming_summary is None
                                else report.trimming_summary.trimmed_gap_percentage
                            ),
                            "resumed": report.resumed,
                            "timeout_seconds": report.run.timeout_seconds,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "model-select":
                report = run_model_selection(
                    args.input_path,
                    out_dir=args.out_dir,
                    prefix=args.prefix,
                    executable=args.executable or "iqtree2",
                    sequence_type=args.sequence_type,
                    partition_path=args.partitions,
                    resume=args.resume,
                    timeout_seconds=args.timeout_seconds,
                    incomplete_run_policy=args.incomplete_run_policy,
                )
                adapter_inputs = (
                    [args.input_path]
                    if args.partitions is None
                    else [args.input_path, args.partitions]
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=adapter_inputs,
                    outputs=[*report.output_paths.values(), report.manifest_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=adapter_inputs,
                        outputs=outputs,
                        warnings=report.run.warning_lines,
                        metrics={
                            "selected_model": report.selected_model,
                            "selected_criterion": (
                                None
                                if report.model_selection_summary is None
                                else report.model_selection_summary.selected_criterion
                            ),
                            "candidate_model_count": (
                                0
                                if report.model_selection_summary is None
                                else report.model_selection_summary.candidate_count
                            ),
                            "best_model_aic": (
                                None
                                if report.model_selection_summary is None
                                else report.model_selection_summary.best_model_aic
                            ),
                            "best_model_aicc": (
                                None
                                if report.model_selection_summary is None
                                else report.model_selection_summary.best_model_aicc
                            ),
                            "best_model_bic": (
                                None
                                if report.model_selection_summary is None
                                else report.model_selection_summary.best_model_bic
                            ),
                            "log_likelihood": report.log_likelihood,
                            "partitioned": args.partitions is not None,
                            "resumed": report.resumed,
                            "timeout_seconds": report.run.timeout_seconds,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "infer-ml":
                report = run_maximum_likelihood_tree_inference(
                    args.input_path,
                    out_dir=args.out_dir,
                    model=args.model,
                    prefix=args.prefix,
                    executable=args.executable or "iqtree2",
                    sequence_type=args.sequence_type,
                    partition_path=args.partitions,
                    resume=args.resume,
                    timeout_seconds=args.timeout_seconds,
                    incomplete_run_policy=args.incomplete_run_policy,
                )
                adapter_inputs = (
                    [args.input_path]
                    if args.partitions is None
                    else [args.input_path, args.partitions]
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=adapter_inputs,
                    outputs=[*report.output_paths.values(), report.manifest_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=adapter_inputs,
                        outputs=outputs,
                        warnings=report.run.warning_lines,
                        metrics={
                            "selected_model": report.selected_model,
                            "selected_criterion": (
                                None
                                if report.model_selection_summary is None
                                else report.model_selection_summary.selected_criterion
                            ),
                            "candidate_model_count": (
                                0
                                if report.model_selection_summary is None
                                else report.model_selection_summary.candidate_count
                            ),
                            "log_likelihood": report.log_likelihood,
                            "support_value_count": (
                                0
                                if report.iqtree_summary is None
                                else report.iqtree_summary.support_value_count
                            ),
                            "partitioned": args.partitions is not None,
                            "resumed": report.resumed,
                            "timeout_seconds": report.run.timeout_seconds,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "bootstrap":
                report = run_bootstrap_support_estimation(
                    args.input_path,
                    out_dir=args.out_dir,
                    model=args.model,
                    replicates=args.replicates,
                    prefix=args.prefix,
                    executable=args.executable or "iqtree2",
                    sequence_type=args.sequence_type,
                    partition_path=args.partitions,
                    resume=args.resume,
                    timeout_seconds=args.timeout_seconds,
                    incomplete_run_policy=args.incomplete_run_policy,
                )
                adapter_inputs = (
                    [args.input_path]
                    if args.partitions is None
                    else [args.input_path, args.partitions]
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=adapter_inputs,
                    outputs=[*report.output_paths.values(), report.manifest_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=adapter_inputs,
                        outputs=outputs,
                        warnings=report.run.warning_lines,
                        metrics={
                            "bootstrap_replicates": args.replicates,
                            "selected_model": report.selected_model,
                            "selected_criterion": (
                                None
                                if report.model_selection_summary is None
                                else report.model_selection_summary.selected_criterion
                            ),
                            "candidate_model_count": (
                                0
                                if report.model_selection_summary is None
                                else report.model_selection_summary.candidate_count
                            ),
                            "log_likelihood": report.log_likelihood,
                            "support_value_count": (
                                0
                                if report.iqtree_summary is None
                                else report.iqtree_summary.support_value_count
                            ),
                            "minimum_support": (
                                None
                                if report.bootstrap_support_summary is None
                                else report.bootstrap_support_summary.minimum_support
                            ),
                            "maximum_support": (
                                None
                                if report.bootstrap_support_summary is None
                                else report.bootstrap_support_summary.maximum_support
                            ),
                            "weakly_supported_clade_count": (
                                0
                                if report.bootstrap_support_summary is None
                                else report.bootstrap_support_summary.weakly_supported_clade_count
                            ),
                            "weak_backbone_node_count": (
                                0
                                if report.weak_backbone_report is None
                                else report.weak_backbone_report.weak_backbone_node_count
                            ),
                            "support_histogram": (
                                {}
                                if report.bootstrap_support_summary is None
                                else report.bootstrap_support_summary.support_histogram
                            ),
                            "partitioned": args.partitions is not None,
                            "resumed": report.resumed,
                            "timeout_seconds": report.run.timeout_seconds,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "sh-alrt":
                report = run_sh_alrt_support_estimation(
                    args.input_path,
                    out_dir=args.out_dir,
                    model=args.model,
                    sh_alrt_replicates=args.alrt_replicates,
                    bootstrap_replicates=args.bootstrap_replicates,
                    prefix=args.prefix,
                    executable=args.executable or "iqtree2",
                    sequence_type=args.sequence_type,
                    partition_path=args.partitions,
                    resume=args.resume,
                    timeout_seconds=args.timeout_seconds,
                    incomplete_run_policy=args.incomplete_run_policy,
                )
                adapter_inputs = (
                    [args.input_path]
                    if args.partitions is None
                    else [args.input_path, args.partitions]
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=adapter_inputs,
                    outputs=[*report.output_paths.values(), report.manifest_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=adapter_inputs,
                        outputs=outputs,
                        warnings=report.run.warning_lines,
                        metrics={
                            "sh_alrt_replicates": args.alrt_replicates,
                            "bootstrap_replicates": args.bootstrap_replicates,
                            "selected_model": report.selected_model,
                            "log_likelihood": report.log_likelihood,
                            "support_value_count": (
                                0
                                if report.iqtree_summary is None
                                else report.iqtree_summary.support_value_count
                            ),
                            "sh_alrt_supported_node_count": (
                                0
                                if report.sh_alrt_support_summary is None
                                else report.sh_alrt_support_summary.annotated_node_count
                            ),
                            "conflicting_support_signal_count": (
                                0
                                if report.sh_alrt_support_summary is None
                                else report.sh_alrt_support_summary.conflicting_support_signal_count
                            ),
                            "minimum_sh_alrt_support": (
                                None
                                if report.sh_alrt_support_summary is None
                                else report.sh_alrt_support_summary.minimum_sh_alrt_support
                            ),
                            "maximum_sh_alrt_support": (
                                None
                                if report.sh_alrt_support_summary is None
                                else report.sh_alrt_support_summary.maximum_sh_alrt_support
                            ),
                            "minimum_ufboot_support": (
                                None
                                if report.sh_alrt_support_summary is None
                                else report.sh_alrt_support_summary.minimum_ufboot_support
                            ),
                            "maximum_ufboot_support": (
                                None
                                if report.sh_alrt_support_summary is None
                                else report.sh_alrt_support_summary.maximum_ufboot_support
                            ),
                            "partitioned": args.partitions is not None,
                            "resumed": report.resumed,
                            "timeout_seconds": report.run.timeout_seconds,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "fasta-to-tree":
                report = run_fasta_to_tree_workflow(
                    args.input_path,
                    out_dir=args.out_dir,
                    prefix=args.prefix,
                    sequence_type=args.sequence_type,
                    mafft_executable=args.mafft_executable or "mafft",
                    alignment_mode=args.alignment_mode,
                    trimal_executable=args.trimal_executable or "trimal",
                    trimming_mode=args.trimming_mode,
                    iqtree_executable=args.iqtree_executable or "iqtree2",
                    iqtree_seed=args.iqtree_seed,
                    iqtree_threads=args.iqtree_threads,
                    trim_gap_threshold=args.trim_gap_threshold,
                    bootstrap_replicates=args.bootstrap_replicates,
                    normalize_identifiers=args.normalize_identifiers,
                    remove_invalid_records=args.remove_invalid_records,
                    resume=args.resume,
                    timeout_seconds=args.timeout_seconds,
                    incomplete_run_policy=args.incomplete_run_policy,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[report.engine_artifact_dir, *report.output_paths.values()],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "alignment_mode": args.alignment_mode,
                            "trimming_mode": args.trimming_mode,
                            "iqtree_seed": args.iqtree_seed,
                            "iqtree_threads": args.iqtree_threads,
                            "bootstrap_replicates": args.bootstrap_replicates,
                            "retained_site_count": (
                                None
                                if report.trimming_workflow.trimming_summary is None
                                else report.trimming_workflow.trimming_summary.retained_site_count
                            ),
                            "removed_site_count": (
                                None
                                if report.trimming_workflow.trimming_summary is None
                                else report.trimming_workflow.trimming_summary.removed_site_count
                            ),
                            "selected_model": report.selected_model,
                            "sequence_type": report.sequence_type,
                            "sequence_type_confidence": (
                                report.input_validation.sequence_type_report.confidence
                                if report.repaired_input_validation is None
                                else report.repaired_input_validation.sequence_type_report.confidence
                            ),
                            "normalized_identifier_count": 0
                            if report.input_repair is None
                            else len(report.input_repair.normalized_identifiers),
                            "removed_record_count": 0
                            if report.input_repair is None
                            else len(report.input_repair.removed_records),
                            "resumed": any(
                                workflow.resumed
                                for workflow in (
                                    report.alignment_workflow,
                                    report.trimming_workflow,
                                    report.model_selection_workflow,
                                    report.maximum_likelihood_workflow,
                                    report.bootstrap_workflow,
                                )
                            ),
                            "timeout_seconds": args.timeout_seconds,
                            **method_tier_metrics(report.method_tier),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "consensus":
                report = run_bootstrap_consensus_tree(
                    args.input_path,
                    out_dir=args.out_dir,
                    prefix=args.prefix,
                    executable=args.executable or "iqtree2",
                    minimum_support=args.minimum_support,
                    resume=args.resume,
                    timeout_seconds=args.timeout_seconds,
                    incomplete_run_policy=args.incomplete_run_policy,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[*report.output_paths.values(), report.manifest_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.run.warning_lines,
                        metrics={
                            "minimum_support": args.minimum_support,
                            "resumed": report.resumed,
                            "timeout_seconds": report.run.timeout_seconds,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            mrbayes_exit_code = run_mrbayes_adapter_command(args)
            if mrbayes_exit_code is not None:
                return mrbayes_exit_code
            beast_exit_code = run_beast_adapter_command(args)
            if beast_exit_code is not None:
                return beast_exit_code
            bayesian_exit_code = run_bayesian_adapter_command(args)
            if bayesian_exit_code is not None:
                return bayesian_exit_code
            inference_exit_code = run_inference_adapter_command(args)
            if inference_exit_code is not None:
                return inference_exit_code
            raise EngineUnavailableError(
                f"unsupported adapter command: {args.adapter_command}"
            )
    except PhylogeneticsError as error:
        if _json_requested(args):
            _print_result(
                build_error_result(
                    command=args.command, inputs=_command_inputs(args), error=error
                ),
                json_output=True,
            )
            return 2
        parser.exit(status=2, message=f"{error.code}: {error.message}\n")
    except FileNotFoundError as error:
        parser.exit(status=2, message=f"{error}\n")
    except ValueError as error:
        parser.exit(status=2, message=f"{error}\n")
    except NotImplementedError as error:
        parser.exit(status=2, message=f"{error}\n")
    except Exception as error:  # pragma: no cover - defensive CLI guard
        parser.exit(status=1, message=f"unexpected error: {error}\n")

    parser.print_help(sys.stderr)
    return 2


def main(argv: list[str] | None = None) -> int:
    """Run the canonical phylogenetics command-line interface."""
    parser = build_parser()
    args = parser.parse_args(argv)
    args._argv = list(argv) if argv is not None else list(sys.argv[1:])
    return run_command(args, parser=parser)
