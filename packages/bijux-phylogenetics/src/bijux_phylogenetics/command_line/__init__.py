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
from bijux_phylogenetics.parity import (
    run_ape_parity_cases,
    run_geiger_parity_cases,
    write_ape_parity_observation_table,
    write_ape_parity_summary_table,
    write_geiger_boundary_warning_table,
    write_geiger_model_confidence_table,
    write_geiger_parity_observation_table,
    write_geiger_likelihood_policy_table,
    write_geiger_optimizer_triage_table,
    write_geiger_parity_summary_table,
    write_geiger_parameterization_registry_table,
)
from bijux_phylogenetics.benchmark import (
    benchmark_alignment_diagnostics,
    benchmark_large_alignment_scaling,
    benchmark_real_dataset_macroevolution,
    benchmark_large_tree_model_fitting,
    benchmark_large_tree_set_scaling,
    benchmark_large_tree_scaling,
    benchmark_large_dataset_stress_suite,
    benchmark_workflow_practical_limits,
    benchmark_tree_comparison,
    benchmark_tree_validation,
)
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
from bijux_phylogenetics.command_line.alignment_distance import (
    add_alignment_distance_commands,
    run_alignment_distance_command,
)
from bijux_phylogenetics.command_line.alignment_review import (
    add_alignment_review_commands,
    run_alignment_review_command,
)
from bijux_phylogenetics.command_line.alignment_matrix import (
    add_alignment_matrix_commands,
    run_alignment_matrix_command,
)
from bijux_phylogenetics.command_line.annotate import (
    add_annotate_command,
    run_annotate_command,
)
from bijux_phylogenetics.command_line.alignment_coding import (
    add_alignment_coding_commands,
    run_alignment_coding_command,
)
from bijux_phylogenetics.command_line.alignment_linkage import (
    add_alignment_linkage_commands,
    run_alignment_linkage_command,
)
from bijux_phylogenetics.command_line.comparative_continuous import (
    add_comparative_continuous_commands,
    run_comparative_continuous_command,
)
from bijux_phylogenetics.command_line.comparative_evolution import (
    add_comparative_evolution_commands,
    run_comparative_evolution_command,
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
    _parse_float_csv_row,
    _parse_probability_assignments,
    _parse_rate_rows,
    _parse_time_bin_definition,
    _parse_transition_pairs,
    _split_csv_values,
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
from bijux_phylogenetics.diagnostics.validation import (
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.comparative.discrete_evolution import (
    compare_discrete_state_models,
    count_discrete_stochastic_map_transitions,
    detect_state_imbalance_problems,
    estimate_ancestral_geographic_states,
    load_stochastic_map_collection,
    render_discrete_state_evolution_report,
    render_stochastic_map_density_artifact,
    render_tree_with_geographic_states,
    simulate_discrete_stochastic_maps,
    summarize_discrete_stochastic_map_density,
    summarize_discrete_stochastic_maps,
    validate_discrete_state_coding,
    validate_discrete_transition_reference_examples,
    write_discrete_model_comparison_table,
    write_node_state_probability_table,
    write_stochastic_map_aggregate_transition_matrix,
    write_stochastic_map_branch_occupancy_table,
    write_stochastic_map_branch_probability_table,
    write_stochastic_map_branch_transition_count_table,
    write_stochastic_map_collection,
    write_stochastic_map_density_branch_table,
    write_stochastic_map_density_slice_table,
    write_stochastic_map_event_table,
    write_stochastic_map_segment_table,
    write_stochastic_map_state_time_table,
    write_stochastic_map_summary_table,
    write_stochastic_map_transition_count_matrix,
    write_transition_summary_table,
)
from bijux_phylogenetics.comparative import (
    build_diversification_method_report,
    build_diversification_figure_package,
    compare_diversification_models,
    compute_diversification_gamma_statistic,
    compute_lineage_through_time_curve,
    detect_diversification_outlier_clades,
    detect_incomplete_taxon_sampling_metadata,
    estimate_diversification_rate,
    render_diversification_report,
    run_trait_dependent_diversification_analysis,
    summarize_geiger_birth_death_exclusion,
    summarize_medusa_exclusion,
    write_diversification_methods_summary_text,
    write_clade_diversification_table,
    write_diversification_gamma_statistic_table,
    write_lineage_through_time_table,
    write_trait_dependent_diversification_table,
)
from bijux_phylogenetics.ecological_niche import (
    summarize_niche_transitions,
    write_niche_state_node_table,
    write_niche_transition_branch_table,
    write_niche_transition_clade_table,
    write_niche_transition_count_table,
    write_niche_transition_exclusion_table,
    write_niche_transition_rate_table,
    write_niche_transition_summary_table,
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
from bijux_phylogenetics.runtime.errors import (
    DiversificationAnalysisError,
    EngineUnavailableError,
    PhylogeneticsError,
)
from bijux_phylogenetics.host_association import (
    summarize_host_switching,
    write_host_state_node_table,
    write_host_switch_branch_table,
    write_host_switch_count_table,
    write_host_switch_exclusion_table,
    write_host_switch_fit_table,
    write_host_switch_summary_table,
    write_unsupported_host_switch_claim_table,
)
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
from bijux_phylogenetics.parity import (
    run_phytools_parity_cases,
    write_phytools_parity_observation_table,
    write_phytools_parity_summary_table,
)
from bijux_phylogenetics.provenance.method_tiers import (
    method_tier_metrics,
    method_tier_warnings,
)
from bijux_phylogenetics.parity import (
    validate_reference_parity_examples,
    write_reference_parity_observation_table,
    write_reference_parity_summary_table,
)
from bijux_phylogenetics.runtime.results import build_command_result, build_error_result
from bijux_phylogenetics.simulation import (
    simulate_birth_death_trees,
    simulate_brownian_traits,
    simulate_coalescent_trees,
    simulate_correlated_brownian_trait_collection,
    simulate_discrete_histories,
    simulate_discrete_traits,
    simulate_dna_alignment,
    simulate_early_burst_traits,
    simulate_ou_traits,
    simulate_speciational_traits,
    simulate_protein_alignment,
    simulate_random_trees,
    validate_geiger_sim_char_reference_examples,
    write_continuous_trait_table,
    write_correlated_continuous_trait_collection_summary_table,
    write_correlated_continuous_trait_collection_table,
    write_discrete_history_branch_truth_table,
    write_discrete_history_event_table,
    write_discrete_history_node_truth_table,
    write_discrete_history_segment_table,
    write_discrete_history_summary_table,
    write_discrete_history_tip_truth_table,
    write_discrete_trait_table,
    write_simulated_alignment,
    write_tree_set,
    write_tree_simulation_envelope_table,
    write_tree_simulation_record_table,
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

    alignment = subparsers.add_parser(
        get_command_spec("alignment").name, help=get_command_spec("alignment").summary
    )
    alignment_subparsers = alignment.add_subparsers(
        dest="alignment_command", required=True
    )
    add_alignment_review_commands(alignment_subparsers)
    add_alignment_matrix_commands(alignment_subparsers)
    add_alignment_distance_commands(alignment_subparsers)
    add_alignment_coding_commands(alignment_subparsers)
    add_alignment_linkage_commands(alignment_subparsers)

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

    ancestral = subparsers.add_parser(
        get_command_spec("ancestral").name,
        help=get_command_spec("ancestral").summary,
    )
    ancestral_subparsers = ancestral.add_subparsers(
        dest="ancestral_command", required=True
    )
    ancestral_continuous = ancestral_subparsers.add_parser(
        "continuous",
        help="Reconstruct ancestral states for a continuous trait.",
    )
    ancestral_continuous.add_argument("tree", type=Path)
    ancestral_continuous.add_argument("table", type=Path)
    ancestral_continuous.add_argument("--trait", required=True)
    ancestral_continuous.add_argument("--taxon-column")
    ancestral_continuous.add_argument(
        "--model", choices=("brownian", "ou"), default="brownian"
    )
    ancestral_continuous.add_argument(
        "--estimator",
        choices=("ace-pic", "anc-ml", "fast-anc", "generalized-least-squares"),
        help="Override the continuous ancestral estimator; default follows the selected model.",
    )
    ancestral_continuous.add_argument("--alpha", type=float, default=1.0)
    ancestral_continuous.add_argument("--table-out", type=Path)
    ancestral_continuous.add_argument("--summary-out", type=Path)
    ancestral_continuous.add_argument("--uncertainty-out", type=Path)
    ancestral_continuous.add_argument("--exclusions-out", type=Path)
    ancestral_continuous.add_argument(
        "--json", action="store_true", help="Emit the reconstruction as JSON."
    )
    _add_manifest_argument(ancestral_continuous)
    ancestral_discrete = ancestral_subparsers.add_parser(
        "discrete",
        help="Reconstruct ancestral states for a discrete trait.",
    )
    ancestral_discrete.add_argument("tree", type=Path)
    ancestral_discrete.add_argument("table", type=Path)
    ancestral_discrete.add_argument("--trait", required=True)
    ancestral_discrete.add_argument("--taxon-column")
    ancestral_discrete.add_argument(
        "--model",
        choices=("fitch", "equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="fitch",
    )
    ancestral_discrete.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_discrete.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_discrete.add_argument(
        "--compare-model",
        choices=("fitch", "equal-rates", "symmetric", "all-rates-different", "meristic"),
        help="Optionally compare this reconstruction directly against another discrete model.",
    )
    ancestral_discrete.add_argument(
        "--root-prior-mode",
        choices=("equal", "empirical", "fixed"),
        default="equal",
        help="Likelihood-only root prior policy for discrete ancestral reconstruction.",
    )
    ancestral_discrete.add_argument(
        "--fixed-root-state",
        help="Required when --root-prior-mode fixed; names the forced root state.",
    )
    ancestral_discrete.add_argument("--table-out", type=Path)
    ancestral_discrete.add_argument("--summary-out", type=Path)
    ancestral_discrete.add_argument("--probabilities-out", type=Path)
    ancestral_discrete.add_argument("--transitions-out", type=Path)
    ancestral_discrete.add_argument("--fit-out", type=Path)
    ancestral_discrete.add_argument("--comparison-out", type=Path)
    ancestral_discrete.add_argument("--exclusions-out", type=Path)
    ancestral_discrete.add_argument(
        "--json", action="store_true", help="Emit the reconstruction as JSON."
    )
    _add_manifest_argument(ancestral_discrete)
    ancestral_discrete_reference = ancestral_subparsers.add_parser(
        "discrete-reference",
        help="Validate built-in discrete ancestral reference examples.",
    )
    ancestral_discrete_reference.add_argument(
        "--json",
        action="store_true",
        help="Emit the reference validation report as JSON.",
    )
    ancestral_tree_set = ancestral_subparsers.add_parser(
        "tree-set",
        help="Summarize ancestral reconstruction stability across a tree set.",
    )
    ancestral_tree_set.add_argument("tree_set", type=Path)
    ancestral_tree_set.add_argument("table", type=Path)
    ancestral_tree_set.add_argument("--trait", required=True)
    ancestral_tree_set.add_argument(
        "--kind", choices=("continuous", "discrete"), required=True
    )
    ancestral_tree_set.add_argument("--taxon-column")
    ancestral_tree_set.add_argument(
        "--model",
        choices=(
            "brownian",
            "ou",
            "fitch",
            "equal-rates",
            "symmetric",
            "all-rates-different",
        ),
    )
    ancestral_tree_set.add_argument("--alpha", type=float, default=1.0)
    ancestral_tree_set.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_tree_set.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_tree_set.add_argument("--burnin-fraction", type=float, default=0.0)
    ancestral_tree_set.add_argument("--summary-out", type=Path)
    ancestral_tree_set.add_argument("--trees-out", type=Path)
    ancestral_tree_set.add_argument("--nodes-out", type=Path)
    ancestral_tree_set.add_argument("--clades-out", type=Path)
    ancestral_tree_set.add_argument("--exclusions-out", type=Path)
    ancestral_tree_set.add_argument(
        "--json", action="store_true", help="Emit the tree-set summary as JSON."
    )
    _add_manifest_argument(ancestral_tree_set)
    ancestral_confidence = ancestral_subparsers.add_parser(
        "confidence",
        help="Summarize ancestral state confidence on one tree or tree set.",
    )
    ancestral_confidence.add_argument("tree", type=Path)
    ancestral_confidence.add_argument("table", type=Path)
    ancestral_confidence.add_argument("--trait", required=True)
    ancestral_confidence.add_argument(
        "--kind", choices=("continuous", "discrete"), required=True
    )
    ancestral_confidence.add_argument("--taxon-column")
    ancestral_confidence.add_argument(
        "--model",
        choices=(
            "brownian",
            "ou",
            "fitch",
            "equal-rates",
            "symmetric",
            "all-rates-different",
        ),
    )
    ancestral_confidence.add_argument("--tree-set", action="store_true")
    ancestral_confidence.add_argument("--alpha", type=float, default=1.0)
    ancestral_confidence.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_confidence.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_confidence.add_argument("--burnin-fraction", type=float, default=0.0)
    ancestral_confidence.add_argument("--summary-out", type=Path)
    ancestral_confidence.add_argument("--confidence-out", type=Path)
    ancestral_confidence.add_argument(
        "--json", action="store_true", help="Emit the confidence review as JSON."
    )
    _add_manifest_argument(ancestral_confidence)
    ancestral_root_sensitivity = ancestral_subparsers.add_parser(
        "root-sensitivity",
        help="Summarize how discrete likelihood ancestral reconstructions change under explicit root assumptions.",
    )
    ancestral_root_sensitivity.add_argument("tree", type=Path)
    ancestral_root_sensitivity.add_argument("table", type=Path)
    ancestral_root_sensitivity.add_argument("--trait", required=True)
    ancestral_root_sensitivity.add_argument("--taxon-column")
    ancestral_root_sensitivity.add_argument(
        "--model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
    )
    ancestral_root_sensitivity.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_root_sensitivity.add_argument(
        "--ordered-states",
        help="Comma-delimited explicit ordered state vocabulary.",
    )
    ancestral_root_sensitivity.add_argument("--fixed-root-state")
    ancestral_root_sensitivity.add_argument("--summary-out", type=Path)
    ancestral_root_sensitivity.add_argument("--assumptions-out", type=Path)
    ancestral_root_sensitivity.add_argument("--nodes-out", type=Path)
    ancestral_root_sensitivity.add_argument(
        "--json", action="store_true", help="Emit the root-sensitivity review as JSON."
    )
    _add_manifest_argument(ancestral_root_sensitivity)
    ancestral_ordered_discrete = ancestral_subparsers.add_parser(
        "ordered-discrete",
        help="Compare ordered and unordered discrete likelihood ancestral reconstructions.",
    )
    ancestral_ordered_discrete.add_argument("tree", type=Path)
    ancestral_ordered_discrete.add_argument("table", type=Path)
    ancestral_ordered_discrete.add_argument("--trait", required=True)
    ancestral_ordered_discrete.add_argument("--taxon-column")
    ancestral_ordered_discrete.add_argument(
        "--model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
    )
    ancestral_ordered_discrete.add_argument(
        "--ordered-states",
        required=True,
        help="Comma-delimited explicit ordered state vocabulary.",
    )
    ancestral_ordered_discrete.add_argument("--summary-out", type=Path)
    ancestral_ordered_discrete.add_argument("--fits-out", type=Path)
    ancestral_ordered_discrete.add_argument("--nodes-out", type=Path)
    ancestral_ordered_discrete.add_argument("--transitions-out", type=Path)
    ancestral_ordered_discrete.add_argument(
        "--json", action="store_true", help="Emit the ordered discrete review as JSON."
    )
    _add_manifest_argument(ancestral_ordered_discrete)
    ancestral_irreversible_discrete = ancestral_subparsers.add_parser(
        "irreversible-discrete",
        help="Compare constrained and unconstrained discrete ancestral likelihood reconstructions under an allowed transition graph.",
    )
    ancestral_irreversible_discrete.add_argument("tree", type=Path)
    ancestral_irreversible_discrete.add_argument("table", type=Path)
    ancestral_irreversible_discrete.add_argument("--trait", required=True)
    ancestral_irreversible_discrete.add_argument("--taxon-column")
    ancestral_irreversible_discrete.add_argument(
        "--model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="all-rates-different",
    )
    ancestral_irreversible_discrete.add_argument(
        "--allowed-transitions",
        required=True,
        help="Comma-delimited directed transition graph in SOURCE->TARGET form.",
    )
    ancestral_irreversible_discrete.add_argument("--summary-out", type=Path)
    ancestral_irreversible_discrete.add_argument("--fits-out", type=Path)
    ancestral_irreversible_discrete.add_argument("--nodes-out", type=Path)
    ancestral_irreversible_discrete.add_argument("--transitions-out", type=Path)
    ancestral_irreversible_discrete.add_argument(
        "--json",
        action="store_true",
        help="Emit the irreversible discrete review as JSON.",
    )
    _add_manifest_argument(ancestral_irreversible_discrete)
    ancestral_transitions = ancestral_subparsers.add_parser(
        "transitions",
        help="Count inferred ancestral transitions on one tree or tree set.",
    )
    ancestral_transitions.add_argument("tree", type=Path)
    ancestral_transitions.add_argument("table", type=Path)
    ancestral_transitions.add_argument("--trait", required=True)
    ancestral_transitions.add_argument("--taxon-column")
    ancestral_transitions.add_argument(
        "--model",
        choices=("fitch", "equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="fitch",
    )
    ancestral_transitions.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_transitions.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_transitions.add_argument("--tree-set", action="store_true")
    ancestral_transitions.add_argument("--burnin-fraction", type=float, default=0.0)
    ancestral_transitions.add_argument("--summary-out", type=Path)
    ancestral_transitions.add_argument("--trees-out", type=Path)
    ancestral_transitions.add_argument("--branches-out", type=Path)
    ancestral_transitions.add_argument("--counts-out", type=Path)
    ancestral_transitions.add_argument("--exclusions-out", type=Path)
    ancestral_transitions.add_argument(
        "--json", action="store_true", help="Emit the transition report as JSON."
    )
    _add_manifest_argument(ancestral_transitions)
    ancestral_compare = ancestral_subparsers.add_parser(
        "compare",
        help="Compare two continuous ancestral-state models node by node.",
    )
    ancestral_compare.add_argument("tree", type=Path)
    ancestral_compare.add_argument("table", type=Path)
    ancestral_compare.add_argument("--trait", required=True)
    ancestral_compare.add_argument("--taxon-column")
    ancestral_compare.add_argument(
        "--left-model", choices=("brownian", "ou"), default="brownian"
    )
    ancestral_compare.add_argument(
        "--right-model", choices=("brownian", "ou"), default="ou"
    )
    ancestral_compare.add_argument("--left-alpha", type=float, default=1.0)
    ancestral_compare.add_argument("--right-alpha", type=float, default=1.0)
    ancestral_compare.add_argument(
        "--json", action="store_true", help="Emit the comparison as JSON."
    )
    _add_manifest_argument(ancestral_compare)
    ancestral_sensitivity = ancestral_subparsers.add_parser(
        "sensitivity",
        help="Summarize how ancestral results change across model, tree, pruning, or coding choices.",
    )
    ancestral_sensitivity.add_argument("tree", type=Path)
    ancestral_sensitivity.add_argument("table", type=Path)
    ancestral_sensitivity.add_argument("--trait", required=True)
    ancestral_sensitivity.add_argument(
        "--kind", choices=("continuous", "discrete"), required=True
    )
    ancestral_sensitivity.add_argument("--taxon-column")
    ancestral_sensitivity.add_argument("--model")
    ancestral_sensitivity.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_sensitivity.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_sensitivity.add_argument("--alpha", type=float, default=1.0)
    ancestral_sensitivity.add_argument("--compare-model")
    ancestral_sensitivity.add_argument("--compare-tree", type=Path)
    ancestral_sensitivity.add_argument("--drop-taxa", nargs="+")
    ancestral_sensitivity.add_argument(
        "--coding-map",
        help="Comma-delimited KEY=VALUE recoding map for discrete traits.",
    )
    ancestral_sensitivity.add_argument(
        "--json", action="store_true", help="Emit the sensitivity report as JSON."
    )
    _add_manifest_argument(ancestral_sensitivity)
    ancestral_render = ancestral_subparsers.add_parser(
        "render",
        help="Render a tree annotated with reconstructed ancestral states.",
    )
    ancestral_render.add_argument("tree", type=Path)
    ancestral_render.add_argument("table", type=Path)
    ancestral_render.add_argument("--trait", required=True)
    ancestral_render.add_argument(
        "--kind", choices=("continuous", "discrete"), required=True
    )
    ancestral_render.add_argument("--taxon-column")
    ancestral_render.add_argument("--model")
    ancestral_render.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_render.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_render.add_argument("--alpha", type=float, default=1.0)
    ancestral_render.add_argument(
        "--layout", choices=("cladogram", "phylogram", "circular"), default="phylogram"
    )
    ancestral_render.add_argument(
        "--discrete-node-style", choices=("labels", "pies"), default="labels"
    )
    ancestral_render.add_argument(
        "--branch-coloring", choices=("none", "state", "regime"), default="none"
    )
    ancestral_render.add_argument("--out", required=True, type=Path)
    ancestral_render.add_argument(
        "--json", action="store_true", help="Emit the render result as JSON."
    )
    _add_manifest_argument(ancestral_render)
    ancestral_report = ancestral_subparsers.add_parser(
        "report",
        help="Render an HTML report for ancestral-state reconstruction.",
    )
    ancestral_report.add_argument("tree", type=Path)
    ancestral_report.add_argument("table", type=Path)
    ancestral_report.add_argument("--trait", required=True)
    ancestral_report.add_argument(
        "--kind", choices=("continuous", "discrete"), required=True
    )
    ancestral_report.add_argument("--taxon-column")
    ancestral_report.add_argument("--model")
    ancestral_report.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_report.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_report.add_argument("--alpha", type=float, default=1.0)
    ancestral_report.add_argument("--compare-model")
    ancestral_report.add_argument("--compare-tree", type=Path)
    ancestral_report.add_argument("--drop-taxa", nargs="+")
    ancestral_report.add_argument(
        "--coding-map",
        help="Comma-delimited KEY=VALUE recoding map for discrete traits.",
    )
    ancestral_report.add_argument("--out", type=Path)
    ancestral_report.add_argument(
        "--out-dir",
        type=Path,
        help="Write a full ancestral reconstruction report package directory.",
    )
    ancestral_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(ancestral_report)
    ancestral_package = ancestral_subparsers.add_parser(
        "package",
        help="Write a publication-ready ancestral-state figure package.",
    )
    ancestral_package.add_argument("tree", type=Path)
    ancestral_package.add_argument("table", type=Path)
    ancestral_package.add_argument("--trait", required=True)
    ancestral_package.add_argument(
        "--kind", choices=("continuous", "discrete"), required=True
    )
    ancestral_package.add_argument("--taxon-column")
    ancestral_package.add_argument("--model")
    ancestral_package.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    ancestral_package.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    ancestral_package.add_argument("--alpha", type=float, default=1.0)
    ancestral_package.add_argument(
        "--layout", choices=("cladogram", "phylogram", "circular"), default="phylogram"
    )
    ancestral_package.add_argument("--out-dir", required=True, type=Path)
    ancestral_package.add_argument(
        "--json", action="store_true", help="Emit the package build result as JSON."
    )
    _add_manifest_argument(ancestral_package)

    biogeography = subparsers.add_parser(
        get_command_spec("biogeography").name,
        help=get_command_spec("biogeography").summary,
    )
    biogeography_subparsers = biogeography.add_subparsers(
        dest="biogeography_command",
        required=True,
    )
    biogeography_model = biogeography_subparsers.add_parser(
        "model",
        help="Reconstruct ancestral geographic regions under an ER, SYM, or ARD transition model.",
    )
    biogeography_model.add_argument("tree", type=Path)
    biogeography_model.add_argument("table", type=Path)
    biogeography_model.add_argument("--trait", required=True)
    biogeography_model.add_argument("--taxon-column")
    biogeography_model.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="er",
    )
    biogeography_model.add_argument(
        "--allowed-regions",
        help="Comma-delimited explicit region vocabulary.",
    )
    biogeography_model.add_argument("--summary-out", type=Path)
    biogeography_model.add_argument("--nodes-out", type=Path)
    biogeography_model.add_argument("--rates-out", type=Path)
    biogeography_model.add_argument("--events-out", type=Path)
    biogeography_model.add_argument("--exclusions-out", type=Path)
    biogeography_model.add_argument(
        "--json", action="store_true", help="Emit the biogeography review as JSON."
    )
    _add_manifest_argument(biogeography_model)
    biogeography_constrained = biogeography_subparsers.add_parser(
        "constrained",
        help="Compare constrained and unconstrained geographic fits under an explicit region adjacency matrix.",
    )
    biogeography_constrained.add_argument("tree", type=Path)
    biogeography_constrained.add_argument("table", type=Path)
    biogeography_constrained.add_argument("adjacency", type=Path)
    biogeography_constrained.add_argument("--trait", required=True)
    biogeography_constrained.add_argument("--taxon-column")
    biogeography_constrained.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="ard",
    )
    biogeography_constrained.add_argument("--summary-out", type=Path)
    biogeography_constrained.add_argument("--fits-out", type=Path)
    biogeography_constrained.add_argument("--transitions-out", type=Path)
    biogeography_constrained.add_argument("--unsupported-out", type=Path)
    biogeography_constrained.add_argument("--exclusions-out", type=Path)
    biogeography_constrained.add_argument(
        "--json", action="store_true", help="Emit the biogeography review as JSON."
    )
    _add_manifest_argument(biogeography_constrained)
    biogeography_time_stratified = biogeography_subparsers.add_parser(
        "time-stratified",
        help="Estimate interval-specific geographic transitions across explicit root-depth bins.",
    )
    biogeography_time_stratified.add_argument("tree", type=Path)
    biogeography_time_stratified.add_argument("table", type=Path)
    biogeography_time_stratified.add_argument("--trait", required=True)
    biogeography_time_stratified.add_argument("--taxon-column")
    biogeography_time_stratified.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="er",
    )
    biogeography_time_stratified.add_argument(
        "--allowed-regions",
        help="Comma-delimited explicit region vocabulary.",
    )
    biogeography_time_stratified.add_argument(
        "--time-bin",
        action="append",
        required=True,
        metavar="LABEL:START:END",
        help="Explicit root-depth interval in LABEL:START:END form. Repeat for multiple intervals.",
    )
    biogeography_time_stratified.add_argument("--summary-out", type=Path)
    biogeography_time_stratified.add_argument("--matrix-out", type=Path)
    biogeography_time_stratified.add_argument("--branches-out", type=Path)
    biogeography_time_stratified.add_argument("--exclusions-out", type=Path)
    biogeography_time_stratified.add_argument(
        "--json", action="store_true", help="Emit the biogeography review as JSON."
    )
    _add_manifest_argument(biogeography_time_stratified)
    biogeography_sampling_bias = biogeography_subparsers.add_parser(
        "sampling-bias",
        help="Reweight geographic region sampling and compare weighted versus unweighted ancestral conclusions.",
    )
    biogeography_sampling_bias.add_argument("tree", type=Path)
    biogeography_sampling_bias.add_argument("table", type=Path)
    biogeography_sampling_bias.add_argument("--trait", required=True)
    biogeography_sampling_bias.add_argument("--taxon-column")
    biogeography_sampling_bias.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="er",
    )
    biogeography_sampling_bias.add_argument(
        "--allowed-regions",
        help="Comma-delimited explicit region vocabulary.",
    )
    biogeography_sampling_bias.add_argument(
        "--weights",
        type=Path,
        help="Optional region-weight table with explicit region and weight columns.",
    )
    biogeography_sampling_bias.add_argument(
        "--region-column",
        default="region",
        help="Region column in the optional region-weight table.",
    )
    biogeography_sampling_bias.add_argument(
        "--weight-column",
        default="weight",
        help="Numeric weight column in the optional region-weight table.",
    )
    biogeography_sampling_bias.add_argument("--summary-out", type=Path)
    biogeography_sampling_bias.add_argument("--regions-out", type=Path)
    biogeography_sampling_bias.add_argument("--nodes-out", type=Path)
    biogeography_sampling_bias.add_argument("--transitions-out", type=Path)
    biogeography_sampling_bias.add_argument("--exclusions-out", type=Path)
    biogeography_sampling_bias.add_argument(
        "--json", action="store_true", help="Emit the biogeography review as JSON."
    )
    _add_manifest_argument(biogeography_sampling_bias)
    biogeography_chronology = biogeography_subparsers.add_parser(
        "chronology",
        help="Place inferred geographic transitions into dated-tree age context with automatic age bins.",
    )
    biogeography_chronology.add_argument("tree", type=Path)
    biogeography_chronology.add_argument("table", type=Path)
    biogeography_chronology.add_argument("--trait", required=True)
    biogeography_chronology.add_argument("--taxon-column")
    biogeography_chronology.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="er",
    )
    biogeography_chronology.add_argument(
        "--allowed-regions",
        help="Comma-delimited explicit region vocabulary.",
    )
    biogeography_chronology.add_argument(
        "--time-bin-count",
        type=int,
        default=4,
        help="Number of equal-width age bins between the tips and the root age.",
    )
    biogeography_chronology.add_argument("--summary-out", type=Path)
    biogeography_chronology.add_argument("--nodes-out", type=Path)
    biogeography_chronology.add_argument("--events-out", type=Path)
    biogeography_chronology.add_argument("--bins-out", type=Path)
    biogeography_chronology.add_argument("--exclusions-out", type=Path)
    biogeography_chronology.add_argument(
        "--json", action="store_true", help="Emit the biogeography review as JSON."
    )
    _add_manifest_argument(biogeography_chronology)
    biogeography_events = biogeography_subparsers.add_parser(
        "events",
        help="Extract inferred geographic movement events on one tree or across a retained tree set.",
    )
    biogeography_events.add_argument("tree", type=Path)
    biogeography_events.add_argument("table", type=Path)
    biogeography_events.add_argument("--trait", required=True)
    biogeography_events.add_argument("--taxon-column")
    biogeography_events.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="er",
    )
    biogeography_events.add_argument(
        "--allowed-regions",
        help="Comma-delimited explicit region vocabulary.",
    )
    biogeography_events.add_argument(
        "--tree-set",
        action="store_true",
        help="Interpret the input tree path as a posterior or bootstrap tree set.",
    )
    biogeography_events.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Fraction of leading trees to discard before tree-set event summarization.",
    )
    biogeography_events.add_argument("--summary-out", type=Path)
    biogeography_events.add_argument("--events-out", type=Path)
    biogeography_events.add_argument("--trees-out", type=Path)
    biogeography_events.add_argument("--event-summaries-out", type=Path)
    biogeography_events.add_argument("--exclusions-out", type=Path)
    biogeography_events.add_argument(
        "--json", action="store_true", help="Emit the biogeography review as JSON."
    )
    _add_manifest_argument(biogeography_events)
    biogeography_report = biogeography_subparsers.add_parser(
        "report",
        help="Build a full biogeography report package with region counts, transition evidence, ancestral-region tree, and map output.",
    )
    biogeography_report.add_argument("tree", type=Path)
    biogeography_report.add_argument("table", type=Path)
    biogeography_report.add_argument("--trait", required=True)
    biogeography_report.add_argument("centroids", type=Path)
    biogeography_report.add_argument("--taxon-column")
    biogeography_report.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="er",
    )
    biogeography_report.add_argument(
        "--region-column",
        default="region",
        help="Region key column in the centroid table.",
    )
    biogeography_report.add_argument(
        "--latitude-column",
        default="latitude",
        help="Latitude column in the centroid table.",
    )
    biogeography_report.add_argument(
        "--longitude-column",
        default="longitude",
        help="Longitude column in the centroid table.",
    )
    biogeography_report.add_argument("--out-dir", required=True, type=Path)
    biogeography_report.add_argument(
        "--json", action="store_true", help="Emit the biogeography package as JSON."
    )
    _add_manifest_argument(biogeography_report)
    host_association = subparsers.add_parser(
        get_command_spec("host-association").name,
        help=get_command_spec("host-association").summary,
    )
    host_association_subparsers = host_association.add_subparsers(
        dest="host_association_command",
        required=True,
    )
    host_association_switches = host_association_subparsers.add_parser(
        "switches",
        help="Reconstruct host states, count host switches, and compare constrained host-transition models.",
    )
    host_association_switches.add_argument("tree", type=Path)
    host_association_switches.add_argument("table", type=Path)
    host_association_switches.add_argument("--trait", required=True)
    host_association_switches.add_argument("--taxon-column")
    host_association_switches.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="ard",
    )
    host_association_switches.add_argument(
        "--constraints",
        type=Path,
        help="Optional host-transition constraint ledger with source_host and target_host columns.",
    )
    host_association_switches.add_argument("--summary-out", type=Path)
    host_association_switches.add_argument("--nodes-out", type=Path)
    host_association_switches.add_argument("--branches-out", type=Path)
    host_association_switches.add_argument("--counts-out", type=Path)
    host_association_switches.add_argument("--fits-out", type=Path)
    host_association_switches.add_argument("--unsupported-out", type=Path)
    host_association_switches.add_argument("--exclusions-out", type=Path)
    host_association_switches.add_argument(
        "--json",
        action="store_true",
        help="Emit the host-association review as JSON.",
    )
    _add_manifest_argument(host_association_switches)
    ecological_niche = subparsers.add_parser(
        get_command_spec("ecological-niche").name,
        help=get_command_spec("ecological-niche").summary,
    )
    ecological_niche_subparsers = ecological_niche.add_subparsers(
        dest="ecological_niche_command",
        required=True,
    )
    ecological_niche_transitions = ecological_niche_subparsers.add_parser(
        "transitions",
        help="Fit ecological niche transitions, reconstruct ancestral niches, and rank clade-specific shift burden.",
    )
    ecological_niche_transitions.add_argument("tree", type=Path)
    ecological_niche_transitions.add_argument("table", type=Path)
    ecological_niche_transitions.add_argument("--trait", required=True)
    ecological_niche_transitions.add_argument("--taxon-column")
    ecological_niche_transitions.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="er",
    )
    ecological_niche_transitions.add_argument("--summary-out", type=Path)
    ecological_niche_transitions.add_argument("--nodes-out", type=Path)
    ecological_niche_transitions.add_argument("--rates-out", type=Path)
    ecological_niche_transitions.add_argument("--branches-out", type=Path)
    ecological_niche_transitions.add_argument("--counts-out", type=Path)
    ecological_niche_transitions.add_argument("--clades-out", type=Path)
    ecological_niche_transitions.add_argument("--exclusions-out", type=Path)
    ecological_niche_transitions.add_argument(
        "--json",
        action="store_true",
        help="Emit the ecological-niche review as JSON.",
    )
    _add_manifest_argument(ecological_niche_transitions)
    phylogeography = subparsers.add_parser(
        get_command_spec("phylogeography").name,
        help=get_command_spec("phylogeography").summary,
    )
    phylogeography_subparsers = phylogeography.add_subparsers(
        dest="phylogeography_command",
        required=True,
    )
    phylogeography_coordinates = phylogeography_subparsers.add_parser(
        "coordinates",
        help="Reconstruct continuous geographic coordinates, review branch movement, and render coordinate-space movement output.",
    )
    phylogeography_coordinates.add_argument("tree", type=Path)
    phylogeography_coordinates.add_argument("table", type=Path)
    phylogeography_coordinates.add_argument("--latitude-column", required=True)
    phylogeography_coordinates.add_argument("--longitude-column", required=True)
    phylogeography_coordinates.add_argument("--taxon-column")
    phylogeography_coordinates.add_argument(
        "--model",
        choices=("brownian", "ou"),
        default="brownian",
    )
    phylogeography_coordinates.add_argument("--alpha", type=float, default=1.0)
    phylogeography_coordinates.add_argument("--summary-out", type=Path)
    phylogeography_coordinates.add_argument("--estimates-out", type=Path)
    phylogeography_coordinates.add_argument("--branches-out", type=Path)
    phylogeography_coordinates.add_argument("--outliers-out", type=Path)
    phylogeography_coordinates.add_argument("--exclusions-out", type=Path)
    phylogeography_coordinates.add_argument("--visualization-out", type=Path)
    phylogeography_coordinates.add_argument(
        "--json",
        action="store_true",
        help="Emit the phylogeography review as JSON.",
    )
    _add_manifest_argument(phylogeography_coordinates)
    phylogeography_coordinates_map = phylogeography_subparsers.add_parser(
        "coordinates-map",
        help="Render one HTML world map from continuous geographic coordinate reconstruction.",
    )
    phylogeography_coordinates_map.add_argument("tree", type=Path)
    phylogeography_coordinates_map.add_argument("table", type=Path)
    phylogeography_coordinates_map.add_argument("--latitude-column", required=True)
    phylogeography_coordinates_map.add_argument("--longitude-column", required=True)
    phylogeography_coordinates_map.add_argument("--taxon-column")
    phylogeography_coordinates_map.add_argument(
        "--model",
        choices=("brownian", "ou"),
        default="brownian",
    )
    phylogeography_coordinates_map.add_argument("--alpha", type=float, default=1.0)
    phylogeography_coordinates_map.add_argument(
        "--minimum-midpoint-depth",
        type=float,
    )
    phylogeography_coordinates_map.add_argument(
        "--maximum-midpoint-depth",
        type=float,
    )
    phylogeography_coordinates_map.add_argument("--summary-out", type=Path)
    phylogeography_coordinates_map.add_argument("--markers-out", type=Path)
    phylogeography_coordinates_map.add_argument("--lines-out", type=Path)
    phylogeography_coordinates_map.add_argument("--exclusions-out", type=Path)
    phylogeography_coordinates_map.add_argument("--html-out", type=Path)
    phylogeography_coordinates_map.add_argument(
        "--json",
        action="store_true",
        help="Emit the mapped phylogeography review as JSON.",
    )
    _add_manifest_argument(phylogeography_coordinates_map)
    phylogeography_regions_map = phylogeography_subparsers.add_parser(
        "regions-map",
        help="Render one HTML world map from discrete ancestral geographic region reconstruction.",
    )
    phylogeography_regions_map.add_argument("tree", type=Path)
    phylogeography_regions_map.add_argument("table", type=Path)
    phylogeography_regions_map.add_argument("--trait", required=True)
    phylogeography_regions_map.add_argument("--centroids", type=Path, required=True)
    phylogeography_regions_map.add_argument("--taxon-column")
    phylogeography_regions_map.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="er",
    )
    phylogeography_regions_map.add_argument(
        "--region-column",
        default="region",
    )
    phylogeography_regions_map.add_argument(
        "--latitude-column",
        default="latitude",
    )
    phylogeography_regions_map.add_argument(
        "--longitude-column",
        default="longitude",
    )
    phylogeography_regions_map.add_argument(
        "--minimum-midpoint-depth",
        type=float,
    )
    phylogeography_regions_map.add_argument(
        "--maximum-midpoint-depth",
        type=float,
    )
    phylogeography_regions_map.add_argument("--summary-out", type=Path)
    phylogeography_regions_map.add_argument("--markers-out", type=Path)
    phylogeography_regions_map.add_argument("--lines-out", type=Path)
    phylogeography_regions_map.add_argument("--exclusions-out", type=Path)
    phylogeography_regions_map.add_argument("--html-out", type=Path)
    phylogeography_regions_map.add_argument(
        "--json",
        action="store_true",
        help="Emit the mapped regional reconstruction review as JSON.",
    )
    _add_manifest_argument(phylogeography_regions_map)

    discrete_evolution = subparsers.add_parser(
        get_command_spec("discrete-evolution").name,
        help=get_command_spec("discrete-evolution").summary,
    )
    discrete_evolution_subparsers = discrete_evolution.add_subparsers(
        dest="discrete_evolution_command",
        required=True,
    )
    discrete_validate = discrete_evolution_subparsers.add_parser(
        "validate-coding",
        help="Validate discrete-state labels against tree-overlapping taxa.",
    )
    discrete_validate.add_argument("tree", type=Path)
    discrete_validate.add_argument("table", type=Path)
    discrete_validate.add_argument("--trait", required=True)
    discrete_validate.add_argument("--taxon-column")
    discrete_validate.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, accept any single token state label.",
    )
    discrete_validate.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    discrete_validate.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    discrete_validate.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(discrete_validate)
    discrete_imbalance = discrete_evolution_subparsers.add_parser(
        "imbalance",
        help="Detect rare, dominant, or degenerate state balance problems.",
    )
    discrete_imbalance.add_argument("tree", type=Path)
    discrete_imbalance.add_argument("table", type=Path)
    discrete_imbalance.add_argument("--trait", required=True)
    discrete_imbalance.add_argument("--taxon-column")
    discrete_imbalance.add_argument(
        "--json", action="store_true", help="Emit the imbalance report as JSON."
    )
    _add_manifest_argument(discrete_imbalance)
    discrete_reference = discrete_evolution_subparsers.add_parser(
        "reference",
        help="Validate deterministic discrete-state transition examples against built-in reference expectations.",
    )
    discrete_reference.add_argument(
        "--json",
        action="store_true",
        help="Emit the reference-validation report as JSON.",
    )
    _add_manifest_argument(discrete_reference)
    discrete_model = discrete_evolution_subparsers.add_parser(
        "model",
        help="Run one discrete-state transition model and export node or branch summaries.",
    )
    discrete_model.add_argument("tree", type=Path)
    discrete_model.add_argument("table", type=Path)
    discrete_model.add_argument("--trait", required=True)
    discrete_model.add_argument("--taxon-column")
    discrete_model.add_argument(
        "--model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
    )
    discrete_model.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    discrete_model.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    discrete_model.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, infer observed states from the table.",
    )
    discrete_model.add_argument(
        "--node-table-out", type=Path, help="Write node-state probabilities as TSV."
    )
    discrete_model.add_argument(
        "--transitions-out", type=Path, help="Write branch transition summaries as TSV."
    )
    discrete_model.add_argument(
        "--json", action="store_true", help="Emit the model report as JSON."
    )
    _add_manifest_argument(discrete_model)
    discrete_compare = discrete_evolution_subparsers.add_parser(
        "compare-models",
        help="Compare two supported discrete-state evolution models node by node.",
    )
    discrete_compare.add_argument("tree", type=Path)
    discrete_compare.add_argument("table", type=Path)
    discrete_compare.add_argument("--trait", required=True)
    discrete_compare.add_argument("--taxon-column")
    discrete_compare.add_argument(
        "--left-model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
    )
    discrete_compare.add_argument(
        "--right-model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="all-rates-different",
    )
    discrete_compare.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    discrete_compare.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    discrete_compare.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, infer observed states from the table.",
    )
    discrete_compare.add_argument(
        "--table-out", type=Path, help="Write node-wise model differences as TSV."
    )
    discrete_compare.add_argument(
        "--json", action="store_true", help="Emit the comparison report as JSON."
    )
    _add_manifest_argument(discrete_compare)
    discrete_stochastic = discrete_evolution_subparsers.add_parser(
        "stochastic-map",
        help="Generate seeded stochastic character maps from a fitted discrete-state CTMC.",
    )
    discrete_stochastic.add_argument("tree", type=Path)
    discrete_stochastic.add_argument("table", type=Path)
    discrete_stochastic.add_argument("--trait", required=True)
    discrete_stochastic.add_argument("--taxon-column")
    discrete_stochastic.add_argument(
        "--model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
    )
    discrete_stochastic.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    discrete_stochastic.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    discrete_stochastic.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, infer observed states from the table.",
    )
    discrete_stochastic.add_argument("--replicates", type=int, default=100)
    discrete_stochastic.add_argument("--seed", type=int, default=0)
    discrete_stochastic.add_argument(
        "--collection-out", type=Path, help="Write stochastic maps as JSON."
    )
    discrete_stochastic.add_argument(
        "--summary-out", type=Path, help="Write stochastic-map summary as TSV."
    )
    discrete_stochastic.add_argument(
        "--state-times-out",
        type=Path,
        help="Write per-state time-in-state summaries as TSV.",
    )
    discrete_stochastic.add_argument(
        "--branch-occupancy-out",
        type=Path,
        help="Write per-branch state-occupancy summaries as TSV.",
    )
    discrete_stochastic.add_argument(
        "--count-matrix-out",
        type=Path,
        help="Write one per-replicate transition count matrix as TSV.",
    )
    discrete_stochastic.add_argument(
        "--aggregate-matrix-out",
        type=Path,
        help="Write one aggregate mean transition matrix as TSV.",
    )
    discrete_stochastic.add_argument(
        "--branch-transition-out",
        type=Path,
        help="Write per-branch transition-count summaries as TSV.",
    )
    discrete_stochastic.add_argument(
        "--segments-out",
        type=Path,
        help="Write flat branch-state segment rows as TSV.",
    )
    discrete_stochastic.add_argument(
        "--events-out",
        type=Path,
        help="Write flat stochastic transition-event rows as TSV.",
    )
    discrete_stochastic.add_argument(
        "--focal-state",
        help="Resolve one focal state for density slices and density rendering. When omitted, binary collections default to the second state in the fitted state order.",
    )
    discrete_stochastic.add_argument(
        "--density-resolution",
        type=int,
        default=100,
        help="Set the branch-slice resolution used for density summaries.",
    )
    discrete_stochastic.add_argument(
        "--branch-probabilities-out",
        type=Path,
        help="Write per-branch state-probability summaries as TSV.",
    )
    discrete_stochastic.add_argument(
        "--density-branches-out",
        type=Path,
        help="Write per-branch focal-state density summaries as TSV.",
    )
    discrete_stochastic.add_argument(
        "--density-slices-out",
        type=Path,
        help="Write flat branch-slice density rows as TSV.",
    )
    discrete_stochastic.add_argument(
        "--density-figure-out",
        type=Path,
        help="Write one branch-colored density artifact as .svg or .html.",
    )
    discrete_stochastic.add_argument(
        "--layout",
        choices=("phylogram", "cladogram", "circular"),
        default="phylogram",
        help="Choose the layout for any density artifact written from this command.",
    )
    discrete_stochastic.add_argument(
        "--json",
        action="store_true",
        help="Emit the stochastic-map collection as JSON.",
    )
    _add_manifest_argument(discrete_stochastic)
    discrete_summarize_maps = discrete_evolution_subparsers.add_parser(
        "summarize-maps",
        help="Summarize a previously written stochastic-map collection.",
    )
    discrete_summarize_maps.add_argument("input_path", type=Path)
    discrete_summarize_maps.add_argument(
        "--summary-out", type=Path, help="Write stochastic-map summary as TSV."
    )
    discrete_summarize_maps.add_argument(
        "--state-times-out",
        type=Path,
        help="Write per-state time-in-state summaries as TSV.",
    )
    discrete_summarize_maps.add_argument(
        "--branch-occupancy-out",
        type=Path,
        help="Write per-branch state-occupancy summaries as TSV.",
    )
    discrete_summarize_maps.add_argument(
        "--json", action="store_true", help="Emit the stochastic-map summary as JSON."
    )
    _add_manifest_argument(discrete_summarize_maps)
    discrete_count_maps = discrete_evolution_subparsers.add_parser(
        "count-maps",
        help="Count directional transitions in a previously written stochastic-map collection.",
    )
    discrete_count_maps.add_argument("input_path", type=Path)
    discrete_count_maps.add_argument(
        "--count-matrix-out",
        type=Path,
        help="Write one per-replicate transition count matrix as TSV.",
    )
    discrete_count_maps.add_argument(
        "--aggregate-matrix-out",
        type=Path,
        help="Write one aggregate mean transition matrix as TSV.",
    )
    discrete_count_maps.add_argument(
        "--branch-transition-out",
        type=Path,
        help="Write per-branch transition-count summaries as TSV.",
    )
    discrete_count_maps.add_argument(
        "--events-out",
        type=Path,
        help="Write flat stochastic transition-event rows as TSV.",
    )
    discrete_count_maps.add_argument(
        "--json",
        action="store_true",
        help="Emit the stochastic-map count report as JSON.",
    )
    _add_manifest_argument(discrete_count_maps)
    discrete_density_maps = discrete_evolution_subparsers.add_parser(
        "density-maps",
        help="Summarize posterior density over a previously written stochastic-map collection.",
    )
    discrete_density_maps.add_argument("input_path", type=Path)
    discrete_density_maps.add_argument(
        "--focal-state",
        help="Resolve one focal state for density slices and density rendering. When omitted, binary collections default to the second state in the fitted state order.",
    )
    discrete_density_maps.add_argument(
        "--resolution",
        type=int,
        default=100,
        help="Set the branch-slice resolution used for density summaries.",
    )
    discrete_density_maps.add_argument(
        "--branch-probabilities-out",
        type=Path,
        help="Write per-branch state-probability summaries as TSV.",
    )
    discrete_density_maps.add_argument(
        "--density-branches-out",
        type=Path,
        help="Write per-branch focal-state density summaries as TSV.",
    )
    discrete_density_maps.add_argument(
        "--density-slices-out",
        type=Path,
        help="Write flat branch-slice density rows as TSV.",
    )
    discrete_density_maps.add_argument(
        "--out",
        type=Path,
        help="Write one branch-colored density artifact as .svg or .html.",
    )
    discrete_density_maps.add_argument(
        "--layout",
        choices=("phylogram", "cladogram", "circular"),
        default="phylogram",
    )
    discrete_density_maps.add_argument(
        "--json",
        action="store_true",
        help="Emit the stochastic-map density report as JSON.",
    )
    _add_manifest_argument(discrete_density_maps)
    discrete_render = discrete_evolution_subparsers.add_parser(
        "render",
        help="Render a tree annotated with reconstructed geographic or other discrete states.",
    )
    discrete_render.add_argument("tree", type=Path)
    discrete_render.add_argument("table", type=Path)
    discrete_render.add_argument("--trait", required=True)
    discrete_render.add_argument("--taxon-column")
    discrete_render.add_argument(
        "--model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
    )
    discrete_render.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    discrete_render.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    discrete_render.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, infer observed states from the table.",
    )
    discrete_render.add_argument(
        "--layout", choices=("cladogram", "phylogram", "circular"), default="phylogram"
    )
    discrete_render.add_argument("--out", required=True, type=Path)
    discrete_render.add_argument(
        "--json", action="store_true", help="Emit the render result as JSON."
    )
    _add_manifest_argument(discrete_render)
    discrete_report = discrete_evolution_subparsers.add_parser(
        "report",
        help="Render an HTML report for one discrete-state evolution analysis.",
    )
    discrete_report.add_argument("tree", type=Path)
    discrete_report.add_argument("table", type=Path)
    discrete_report.add_argument("--trait", required=True)
    discrete_report.add_argument("--taxon-column")
    discrete_report.add_argument(
        "--model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
    )
    discrete_report.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    discrete_report.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    discrete_report.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, infer observed states from the table.",
    )
    discrete_report.add_argument(
        "--compare-model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
    )
    discrete_report.add_argument("--out", required=True, type=Path)
    discrete_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(discrete_report)

    diversification = subparsers.add_parser(
        get_command_spec("diversification").name,
        help=get_command_spec("diversification").summary,
    )
    diversification_subparsers = diversification.add_subparsers(
        dest="diversification_command", required=True
    )
    diversification_ltt = diversification_subparsers.add_parser(
        "ltt",
        help="Compute a lineage-through-time curve for one rooted ultrametric tree.",
    )
    diversification_ltt.add_argument("tree", type=Path)
    diversification_ltt.add_argument(
        "--out", type=Path, help="Write the lineage-through-time table as TSV."
    )
    diversification_ltt.add_argument(
        "--json", action="store_true", help="Emit the LTT report as JSON."
    )
    _add_manifest_argument(diversification_ltt)
    diversification_sampling = diversification_subparsers.add_parser(
        "sampling",
        help="Inspect taxon sampling-fraction metadata against the tree tips.",
    )
    diversification_sampling.add_argument("tree", type=Path)
    diversification_sampling.add_argument("table", type=Path)
    diversification_sampling.add_argument("--taxon-column")
    diversification_sampling.add_argument("--sampling-column")
    diversification_sampling.add_argument(
        "--json", action="store_true", help="Emit the sampling report as JSON."
    )
    _add_manifest_argument(diversification_sampling)
    diversification_estimate = diversification_subparsers.add_parser(
        "estimate",
        help="Estimate a simple Yule or birth-death diversification model.",
    )
    diversification_estimate.add_argument("tree", type=Path)
    diversification_estimate.add_argument("--metadata", type=Path)
    diversification_estimate.add_argument("--taxon-column")
    diversification_estimate.add_argument("--sampling-column")
    diversification_estimate.add_argument(
        "--model", choices=("yule", "birth-death"), default="birth-death"
    )
    diversification_estimate.add_argument(
        "--json", action="store_true", help="Emit the diversification estimate as JSON."
    )
    _add_manifest_argument(diversification_estimate)
    diversification_gamma = diversification_subparsers.add_parser(
        "gamma-stat",
        help="Compute the Pybus-Harvey diversification gamma statistic.",
    )
    diversification_gamma.add_argument("tree", type=Path)
    diversification_gamma.add_argument("--metadata", type=Path)
    diversification_gamma.add_argument("--taxon-column")
    diversification_gamma.add_argument("--sampling-column")
    diversification_gamma.add_argument(
        "--out",
        type=Path,
        help="Write the diversification gamma-statistic table as TSV.",
    )
    diversification_gamma.add_argument(
        "--json",
        action="store_true",
        help="Emit the diversification gamma-statistic report as JSON.",
    )
    _add_manifest_argument(diversification_gamma)
    diversification_compare = diversification_subparsers.add_parser(
        "compare-models",
        help="Compare Yule and birth-death diversification fits.",
    )
    diversification_compare.add_argument("tree", type=Path)
    diversification_compare.add_argument("--metadata", type=Path)
    diversification_compare.add_argument("--taxon-column")
    diversification_compare.add_argument("--sampling-column")
    diversification_compare.add_argument(
        "--json", action="store_true", help="Emit the model comparison as JSON."
    )
    _add_manifest_argument(diversification_compare)
    diversification_clades = diversification_subparsers.add_parser(
        "clades",
        help="Detect clades with unusually high or low diversification.",
    )
    diversification_clades.add_argument("tree", type=Path)
    diversification_clades.add_argument(
        "--model", choices=("yule", "birth-death"), default="birth-death"
    )
    diversification_clades.add_argument("--min-tip-count", type=int, default=2)
    diversification_clades.add_argument(
        "--out", type=Path, help="Write the clade diversification table as TSV."
    )
    diversification_clades.add_argument(
        "--json", action="store_true", help="Emit the clade scan report as JSON."
    )
    _add_manifest_argument(diversification_clades)
    diversification_trait = diversification_subparsers.add_parser(
        "trait-dependent",
        help="Summarize simple trait-linked diversification rates when states form interpretable clades.",
    )
    diversification_trait.add_argument("tree", type=Path)
    diversification_trait.add_argument("table", type=Path)
    diversification_trait.add_argument("--trait", required=True)
    diversification_trait.add_argument("--taxon-column")
    diversification_trait.add_argument(
        "--out",
        type=Path,
        help="Write the trait-dependent diversification table as TSV.",
    )
    diversification_trait.add_argument(
        "--json", action="store_true", help="Emit the trait-dependent report as JSON."
    )
    _add_manifest_argument(diversification_trait)
    diversification_package = diversification_subparsers.add_parser(
        "package",
        help="Build a publication-oriented diversification figure package.",
    )
    diversification_package.add_argument("tree", type=Path)
    diversification_package.add_argument("--metadata", type=Path)
    diversification_package.add_argument("--taxon-column")
    diversification_package.add_argument("--sampling-column")
    diversification_package.add_argument(
        "--model", choices=("yule", "birth-death"), default="birth-death"
    )
    diversification_package.add_argument("--min-tip-count", type=int, default=2)
    diversification_package.add_argument("--out-dir", required=True, type=Path)
    diversification_package.add_argument(
        "--json", action="store_true", help="Emit the package build result as JSON."
    )
    _add_manifest_argument(diversification_package)
    diversification_report = diversification_subparsers.add_parser(
        "report",
        help="Render an HTML diversification and macroevolution report.",
    )
    diversification_report.add_argument("tree", type=Path)
    diversification_report.add_argument("--metadata", type=Path)
    diversification_report.add_argument("--taxon-column")
    diversification_report.add_argument("--sampling-column")
    diversification_report.add_argument("--traits", type=Path)
    diversification_report.add_argument("--trait")
    diversification_report.add_argument("--out", required=True, type=Path)
    diversification_report.add_argument(
        "--methods-summary-out",
        type=Path,
        help="Write reviewer-facing Markdown methods text for the diversification analysis.",
    )
    diversification_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(diversification_report)
    diversification_methods_summary = diversification_subparsers.add_parser(
        "methods-summary",
        help="Write reviewer-facing Markdown methods text for one diversification analysis.",
    )
    diversification_methods_summary.add_argument("tree", type=Path)
    diversification_methods_summary.add_argument("--metadata", type=Path)
    diversification_methods_summary.add_argument("--taxon-column")
    diversification_methods_summary.add_argument("--sampling-column")
    diversification_methods_summary.add_argument("--traits", type=Path)
    diversification_methods_summary.add_argument("--trait")
    diversification_methods_summary.add_argument(
        "--estimate-model",
        choices=("yule", "birth-death"),
        default="birth-death",
    )
    diversification_methods_summary.add_argument(
        "--clade-model",
        choices=("yule", "birth-death"),
        default="birth-death",
    )
    diversification_methods_summary.add_argument(
        "--min-tip-count",
        type=int,
        default=2,
        help="Minimum clade size included in the clade outlier review.",
    )
    diversification_methods_summary.add_argument("--out", required=True, type=Path)
    diversification_methods_summary.add_argument(
        "--json", action="store_true", help="Emit the methods summary metrics as JSON."
    )
    _add_manifest_argument(diversification_methods_summary)
    diversification_medusa = diversification_subparsers.add_parser(
        "medusa",
        help="Explain the explicit exclusion boundary for geiger::medusa parity.",
    )
    diversification_medusa.add_argument("tree", type=Path)
    diversification_medusa.add_argument("--metadata", type=Path)
    diversification_medusa.add_argument("--taxon-column")
    diversification_medusa.add_argument("--sampling-column")
    diversification_medusa.add_argument(
        "--json", action="store_true", help="Emit the MEDUSA exclusion as JSON."
    )
    _add_manifest_argument(diversification_medusa)
    diversification_bd_ms = diversification_subparsers.add_parser(
        "bd-ms",
        help="Explain the explicit exclusion boundary for geiger::bd.ms birth-death parity.",
    )
    diversification_bd_ms.add_argument("tree", type=Path)
    diversification_bd_ms.add_argument("--metadata", type=Path)
    diversification_bd_ms.add_argument("--taxon-column")
    diversification_bd_ms.add_argument("--sampling-column")
    diversification_bd_ms.add_argument(
        "--json",
        action="store_true",
        help="Emit the birth-death exclusion as JSON.",
    )
    _add_manifest_argument(diversification_bd_ms)

    add_distance_commands(subparsers)
    add_tree_set_commands(subparsers)

    simulate = subparsers.add_parser(
        get_command_spec("simulate").name, help=get_command_spec("simulate").summary
    )
    simulate_subparsers = simulate.add_subparsers(
        dest="simulate_command", required=True
    )
    simulate_birth_death = simulate_subparsers.add_parser(
        "tree-birth-death",
        help="Simulate one or more trees under a birth-death process.",
    )
    simulate_birth_death.add_argument("--tree-count", type=int, default=1)
    simulate_birth_death.add_argument("--tip-count", type=int, required=True)
    simulate_birth_death.add_argument("--birth-rate", type=float, default=1.0)
    simulate_birth_death.add_argument("--death-rate", type=float, default=0.25)
    simulate_birth_death.add_argument("--seed", type=int, default=1)
    simulate_birth_death.add_argument("--out", required=True, type=Path)
    simulate_birth_death.add_argument("--record-table-out", type=Path)
    simulate_birth_death.add_argument("--envelope-table-out", type=Path)
    simulate_birth_death.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_birth_death)
    simulate_random_tree = simulate_subparsers.add_parser(
        "tree-random",
        help="Simulate one or more rooted random trees with uniform branch lengths.",
    )
    simulate_random_tree.add_argument("--tree-count", type=int, default=1)
    simulate_random_tree.add_argument("--tip-count", type=int, required=True)
    simulate_random_tree.add_argument("--seed", type=int, default=1)
    simulate_random_tree.add_argument("--out", required=True, type=Path)
    simulate_random_tree.add_argument("--record-table-out", type=Path)
    simulate_random_tree.add_argument("--envelope-table-out", type=Path)
    simulate_random_tree.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_random_tree)
    simulate_coalescent = simulate_subparsers.add_parser(
        "tree-coalescent",
        help="Simulate one or more trees under a coalescent model.",
    )
    simulate_coalescent.add_argument("--tree-count", type=int, default=1)
    simulate_coalescent.add_argument("--tip-count", type=int, required=True)
    simulate_coalescent.add_argument("--population-size", type=float, default=1.0)
    simulate_coalescent.add_argument("--seed", type=int, default=1)
    simulate_coalescent.add_argument("--out", required=True, type=Path)
    simulate_coalescent.add_argument("--record-table-out", type=Path)
    simulate_coalescent.add_argument("--envelope-table-out", type=Path)
    simulate_coalescent.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_coalescent)
    simulate_brownian = simulate_subparsers.add_parser(
        "traits-brownian",
        help="Simulate a continuous tip trait under Brownian motion.",
    )
    simulate_brownian.add_argument("tree", type=Path)
    simulate_brownian.add_argument("--root-state", type=float, default=0.0)
    simulate_brownian.add_argument("--sigma", type=float)
    simulate_brownian.add_argument("--sigma-squared", type=float)
    simulate_brownian.add_argument("--seed", type=int, default=1)
    simulate_brownian.add_argument("--out", required=True, type=Path)
    simulate_brownian.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_brownian)
    simulate_speciational = simulate_subparsers.add_parser(
        "traits-speciational",
        help="Simulate a continuous tip trait under geiger-style speciational Brownian motion.",
    )
    simulate_speciational.add_argument("tree", type=Path)
    simulate_speciational.add_argument("--root-state", type=float, default=0.0)
    simulate_speciational.add_argument("--sigma", type=float)
    simulate_speciational.add_argument("--sigma-squared", type=float)
    simulate_speciational.add_argument("--seed", type=int, default=1)
    simulate_speciational.add_argument("--out", required=True, type=Path)
    simulate_speciational.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_speciational)
    simulate_brownian_correlated = simulate_subparsers.add_parser(
        "traits-brownian-correlated",
        help="Simulate correlated continuous tip traits under multivariate Brownian motion.",
    )
    simulate_brownian_correlated.add_argument("tree", type=Path)
    simulate_brownian_correlated.add_argument(
        "--trait",
        action="append",
        required=True,
        help="One trait name. Repeat to define the multivariate trait order.",
    )
    simulate_brownian_correlated.add_argument(
        "--root-state",
        action="append",
        type=float,
        default=[],
        help="One root state value aligned to the declared trait order.",
    )
    covariance_group = simulate_brownian_correlated.add_mutually_exclusive_group(
        required=True
    )
    covariance_group.add_argument(
        "--covariance-row",
        action="append",
        dest="covariance_rows",
        help="One comma-separated covariance-matrix row. Repeat to build the full matrix.",
    )
    covariance_group.add_argument(
        "--correlation-row",
        action="append",
        dest="correlation_rows",
        help="One comma-separated correlation-matrix row. Repeat to build the full matrix.",
    )
    simulate_brownian_correlated.add_argument(
        "--trait-standard-deviation",
        action="append",
        type=float,
        default=[],
        help="One trait standard deviation aligned to the declared trait order when using --correlation-row.",
    )
    simulate_brownian_correlated.add_argument("--replicates", type=int, default=128)
    simulate_brownian_correlated.add_argument("--seed", type=int, default=1)
    simulate_brownian_correlated.add_argument("--out", required=True, type=Path)
    simulate_brownian_correlated.add_argument("--summary-out", type=Path)
    simulate_brownian_correlated.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_brownian_correlated)
    simulate_ou = simulate_subparsers.add_parser(
        "traits-ou",
        help="Simulate a continuous tip trait under an OU process.",
    )
    simulate_ou.add_argument("tree", type=Path)
    simulate_ou.add_argument("--root-state", type=float, default=0.0)
    simulate_ou.add_argument("--sigma", type=float, default=1.0)
    simulate_ou.add_argument("--alpha", type=float, default=1.0)
    simulate_ou.add_argument("--theta", type=float, default=0.0)
    simulate_ou.add_argument("--seed", type=int, default=1)
    simulate_ou.add_argument("--out", required=True, type=Path)
    simulate_ou.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_ou)
    simulate_early_burst = simulate_subparsers.add_parser(
        "traits-early-burst",
        help="Simulate a continuous tip trait under an early-burst branch-rate process.",
    )
    simulate_early_burst.add_argument("tree", type=Path)
    simulate_early_burst.add_argument("--root-state", type=float, default=0.0)
    simulate_early_burst.add_argument("--sigma", type=float, default=1.0)
    simulate_early_burst.add_argument("--rate-change", type=float, default=1.0)
    simulate_early_burst.add_argument("--seed", type=int, default=1)
    simulate_early_burst.add_argument("--out", required=True, type=Path)
    simulate_early_burst.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_early_burst)
    simulate_discrete = simulate_subparsers.add_parser(
        "traits-discrete",
        help="Simulate a discrete tip trait under a symmetric jump process.",
    )
    simulate_discrete.add_argument("tree", type=Path)
    simulate_discrete.add_argument("--states", nargs="+", required=True)
    simulate_discrete.add_argument("--transition-rate", type=float, default=1.0)
    simulate_discrete.add_argument("--root-state")
    simulate_discrete.add_argument("--seed", type=int, default=1)
    simulate_discrete.add_argument("--out", required=True, type=Path)
    simulate_discrete.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_discrete)
    simulate_history_discrete = simulate_subparsers.add_parser(
        "history-discrete",
        help="Simulate true discrete histories on a fixed tree from one explicit rate matrix.",
    )
    simulate_history_discrete.add_argument("tree", type=Path)
    simulate_history_discrete.add_argument("--states", nargs="+", required=True)
    simulate_history_discrete.add_argument(
        "--rate",
        action="append",
        required=True,
        help="One SOURCE->TARGET=RATE entry. Repeat to build the rate matrix.",
    )
    simulate_history_discrete.add_argument("--root-state")
    simulate_history_discrete.add_argument(
        "--root-probability",
        action="append",
        default=[],
        help="One STATE=PROBABILITY entry. Repeat to define the root prior.",
    )
    simulate_history_discrete.add_argument("--replicates", type=int, default=1)
    simulate_history_discrete.add_argument("--seed", type=int, default=1)
    simulate_history_discrete.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Tip-state truth table output path.",
    )
    simulate_history_discrete.add_argument("--nodes-out", type=Path)
    simulate_history_discrete.add_argument("--branches-out", type=Path)
    simulate_history_discrete.add_argument("--events-out", type=Path)
    simulate_history_discrete.add_argument("--segments-out", type=Path)
    simulate_history_discrete.add_argument("--summary-out", type=Path)
    simulate_history_discrete.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_history_discrete)
    simulate_dna = simulate_subparsers.add_parser(
        "alignment-dna",
        help="Simulate a DNA alignment along a rooted tree.",
    )
    simulate_dna.add_argument("tree", type=Path)
    simulate_dna.add_argument("--sequence-length", type=int, required=True)
    simulate_dna.add_argument("--substitution-rate", type=float, default=1.0)
    simulate_dna.add_argument("--seed", type=int, default=1)
    simulate_dna.add_argument("--out", required=True, type=Path)
    simulate_dna.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_dna)
    simulate_protein = simulate_subparsers.add_parser(
        "alignment-protein",
        help="Simulate a protein alignment along a rooted tree.",
    )
    simulate_protein.add_argument("tree", type=Path)
    simulate_protein.add_argument("--sequence-length", type=int, required=True)
    simulate_protein.add_argument("--substitution-rate", type=float, default=1.0)
    simulate_protein.add_argument("--seed", type=int, default=1)
    simulate_protein.add_argument("--out", required=True, type=Path)
    simulate_protein.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_protein)
    simulate_validate_sim_char_reference = simulate_subparsers.add_parser(
        "validate-sim-char-reference",
        help="Validate governed geiger::sim.char summary envelopes.",
    )
    simulate_validate_sim_char_reference.add_argument(
        "--json",
        action="store_true",
        help="Emit the governed validation report as JSON.",
    )
    _add_manifest_argument(simulate_validate_sim_char_reference)

    benchmark = subparsers.add_parser(
        get_command_spec("benchmark").name, help=get_command_spec("benchmark").summary
    )
    benchmark_subparsers = benchmark.add_subparsers(
        dest="benchmark_command", required=True
    )
    benchmark_validate = benchmark_subparsers.add_parser(
        "tree-validation",
        help="Benchmark tree validation across size classes.",
    )
    benchmark_validate.add_argument("--replicates", type=int, default=3)
    benchmark_validate.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_validate)
    benchmark_compare = benchmark_subparsers.add_parser(
        "tree-comparison",
        help="Benchmark tree comparison across increasing taxon counts.",
    )
    benchmark_compare.add_argument("--replicates", type=int, default=3)
    benchmark_compare.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_compare)
    benchmark_large_tree = benchmark_subparsers.add_parser(
        "large-tree-scaling",
        help="Benchmark large-tree validation, comparison, rendering, and reporting.",
    )
    benchmark_large_tree.add_argument("--replicates", type=int, default=1)
    benchmark_large_tree.add_argument(
        "--tip-count",
        action="append",
        dest="tip_counts",
        type=int,
        help="Add one governed tree size to benchmark. Repeat to benchmark multiple sizes.",
    )
    benchmark_large_tree.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_large_tree)
    benchmark_large_alignment = benchmark_subparsers.add_parser(
        "large-alignment-scaling",
        help="Benchmark large-alignment diagnostics, trimming, distance, and readiness.",
    )
    benchmark_large_alignment.add_argument("--replicates", type=int, default=1)
    benchmark_large_alignment.add_argument(
        "--sequence-count",
        action="append",
        dest="sequence_counts",
        type=int,
        help="Add one sequence count to benchmark. Repeat to benchmark multiple size classes.",
    )
    benchmark_large_alignment.add_argument(
        "--alignment-length",
        action="append",
        dest="alignment_lengths",
        type=int,
        help="Add one alignment length to benchmark. Repeat to benchmark multiple size classes.",
    )
    benchmark_large_alignment.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_large_alignment)
    benchmark_large_tree_set = benchmark_subparsers.add_parser(
        "large-tree-set-scaling",
        help="Benchmark large-tree-set consensus, RF diversity, clustering, and uncertainty summaries.",
    )
    benchmark_large_tree_set.add_argument("--replicates", type=int, default=1)
    benchmark_large_tree_set.add_argument(
        "--tree-count",
        action="append",
        dest="tree_counts",
        type=int,
        help="Add one posterior tree count to benchmark. Repeat to benchmark multiple size classes.",
    )
    benchmark_large_tree_set.add_argument(
        "--tip-count",
        action="append",
        dest="tip_counts",
        type=int,
        help="Add one taxon count to benchmark. Repeat to benchmark multiple size classes.",
    )
    benchmark_large_tree_set.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_large_tree_set)
    benchmark_practical_limits = benchmark_subparsers.add_parser(
        "workflow-practical-limits",
        help="Report the largest governed workflow classes currently exercised in benchmark and stress lanes.",
    )
    benchmark_practical_limits.add_argument("--replicates", type=int, default=1)
    benchmark_practical_limits.add_argument(
        "--tree-tip-count",
        action="append",
        dest="tree_tip_counts",
        type=int,
        help="Add one large-tree taxon count. Repeat to override the governed tree-size classes.",
    )
    benchmark_practical_limits.add_argument(
        "--sequence-count",
        action="append",
        dest="sequence_counts",
        type=int,
        help="Add one sequence count for the large-alignment classes. Repeat alongside --alignment-length.",
    )
    benchmark_practical_limits.add_argument(
        "--alignment-length",
        action="append",
        dest="alignment_lengths",
        type=int,
        help="Add one aligned-site count for the large-alignment classes. Repeat alongside --sequence-count.",
    )
    benchmark_practical_limits.add_argument(
        "--posterior-tree-count",
        action="append",
        dest="posterior_tree_counts",
        type=int,
        help="Add one posterior tree count for the tree-set classes. Repeat alongside --tree-set-tip-count.",
    )
    benchmark_practical_limits.add_argument(
        "--tree-set-tip-count",
        action="append",
        dest="tree_set_tip_counts",
        type=int,
        help="Add one taxon count for the tree-set classes. Repeat alongside --posterior-tree-count.",
    )
    benchmark_practical_limits.add_argument(
        "--stress-tier",
        action="append",
        dest="stress_tiers",
        choices=("small", "heavy"),
        help="Include one governed stress tier. Repeat to aggregate multiple tiers.",
    )
    benchmark_practical_limits.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_practical_limits)
    benchmark_alignment = benchmark_subparsers.add_parser(
        "alignment-diagnostics",
        help="Benchmark alignment diagnostics across increasing sequence counts.",
    )
    benchmark_alignment.add_argument("--replicates", type=int, default=3)
    benchmark_alignment.add_argument("--sequence-length", type=int, default=128)
    benchmark_alignment.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_alignment)
    benchmark_stress = benchmark_subparsers.add_parser(
        "stress-suite",
        help="Benchmark large-dataset stress workloads across governed tiers.",
    )
    benchmark_stress.add_argument(
        "--tier",
        choices=("small", "heavy"),
        default="small",
        help="Select the governed stress tier to execute.",
    )
    benchmark_stress.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_stress)
    benchmark_large_tree_model = benchmark_subparsers.add_parser(
        "large-tree-model-fitting",
        help="Benchmark 100+ taxon continuous and discrete model fitting with governed geiger comparison and heavy-tier review.",
    )
    benchmark_large_tree_model.add_argument(
        "--tier",
        choices=("small", "heavy"),
        default="small",
        help="Select the governed model-fitting tier to execute.",
    )
    benchmark_large_tree_model.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_large_tree_model)
    benchmark_real_dataset = benchmark_subparsers.add_parser(
        "real-dataset-macroevolution",
        help="Benchmark continuous and discrete macroevolution model fitting on the published Central European seashore flora dataset against stored local geiger references.",
    )
    benchmark_real_dataset.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_real_dataset)

    parity = subparsers.add_parser(
        get_command_spec("parity").name, help=get_command_spec("parity").summary
    )
    parity.add_argument(
        "--reference-source",
        choices=("checked-fixture", "ape-live", "phytools-live", "geiger-live"),
        default="checked-fixture",
        help="Choose the checked fixture parity suite or one of the live external execution harnesses.",
    )
    parity.add_argument(
        "--extended",
        action="store_true",
        help="Include the optional larger posterior tree-set parity fixtures.",
    )
    parity.add_argument(
        "--ape-case",
        action="append",
        dest="ape_cases",
        help="Restrict the live ape parity harness to one or more governed case ids.",
    )
    parity.add_argument(
        "--ape-rscript-executable",
        default="Rscript",
        help="Executable used to launch the live ape parity runner.",
    )
    parity.add_argument(
        "--ape-failure-root",
        type=Path,
        help="Directory for reproducible live ape mismatch and skip artifacts.",
    )
    parity.add_argument(
        "--phytools-case",
        action="append",
        dest="phytools_cases",
        help="Restrict the live phytools parity harness to one or more governed case ids.",
    )
    parity.add_argument(
        "--phytools-rscript-executable",
        default="Rscript",
        help="Executable used to launch the live phytools parity runner.",
    )
    parity.add_argument(
        "--phytools-failure-root",
        type=Path,
        help="Directory for reproducible live phytools mismatch and skip artifacts.",
    )
    parity.add_argument(
        "--geiger-case",
        action="append",
        dest="geiger_cases",
        help="Restrict the live geiger parity harness to one or more governed case ids.",
    )
    parity.add_argument(
        "--geiger-rscript-executable",
        default="Rscript",
        help="Executable used to launch the live geiger parity runner.",
    )
    parity.add_argument(
        "--geiger-failure-root",
        type=Path,
        help="Directory for reproducible live geiger mismatch and skip artifacts.",
    )
    parity.add_argument("--summary-out", type=Path)
    parity.add_argument("--observations-out", type=Path)
    parity.add_argument("--optimizer-triage-out", type=Path)
    parity.add_argument("--boundary-warning-out", type=Path)
    parity.add_argument("--likelihood-policy-out", type=Path)
    parity.add_argument("--model-confidence-out", type=Path)
    parity.add_argument("--parameterization-registry-out", type=Path)
    parity.add_argument(
        "--json", action="store_true", help="Emit the parity report as JSON."
    )
    _add_manifest_argument(parity)

    validate = subparsers.add_parser(
        get_command_spec("validate").name, help=get_command_spec("validate").summary
    )
    validate.add_argument("tree", type=Path)
    validate.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    validate.add_argument("--allow-duplicates", action="store_true")
    validate.add_argument("--allow-negative-branches", action="store_true")
    validate.add_argument("--require-rooted", action="store_true")
    validate.add_argument("--require-ultrametric", action="store_true")
    validate.add_argument("--strict", action="store_true")
    validate.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(validate)

    inspect = subparsers.add_parser(
        get_command_spec("inspect").name, help=get_command_spec("inspect").summary
    )
    inspect.add_argument("tree", type=Path)
    inspect.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    inspect.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(inspect)

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
            report = validate_tree_path(
                args.tree,
                source_format=args.format,
                allow_duplicates=args.allow_duplicates,
                strict=args.strict,
                allow_negative_branch_lengths=args.allow_negative_branches,
                require_rooted=args.require_rooted,
                require_ultrametric=args.require_ultrametric,
            )
            outputs = _finalize_outputs(args, command="validate", inputs=[args.tree])
            _print_result(
                build_command_result(
                    command="validate",
                    inputs=[args.tree],
                    outputs=outputs,
                    warnings=report.warnings,
                    metrics={
                        "tip_count": report.tip_count,
                        "internal_node_count": report.internal_node_count,
                        "validity_decision": report.validity_decision,
                        "syntax_valid": report.syntax_valid,
                        "biologically_safe": report.biologically_safe,
                        "polytomy_count": report.polytomy_count,
                        "missing_internal_branch_count": len(
                            report.missing_internal_branch_nodes
                        ),
                        "missing_terminal_branch_count": len(
                            report.missing_terminal_branch_taxa
                        ),
                        "singleton_internal_node_count": len(
                            report.singleton_internal_nodes
                        ),
                        "integrity_issue_count": len(report.integrity_issues),
                        "unsafe_external_label_count": len(
                            report.unsafe_external_labels
                        ),
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "prune":
            return run_prune_command(args)
        if args.command == "alignment":
            review_exit_code = run_alignment_review_command(args)
            if review_exit_code is not None:
                return review_exit_code
            matrix_exit_code = run_alignment_matrix_command(args)
            if matrix_exit_code is not None:
                return matrix_exit_code
            distance_exit_code = run_alignment_distance_command(args)
            if distance_exit_code is not None:
                return distance_exit_code
            coding_exit_code = run_alignment_coding_command(args)
            if coding_exit_code is not None:
                return coding_exit_code
            linkage_exit_code = run_alignment_linkage_command(args)
            if linkage_exit_code is not None:
                return linkage_exit_code
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
            if args.ancestral_command == "continuous":
                if (
                    args.model == "brownian"
                    and args.estimator == "generalized-least-squares"
                ):
                    parser.error(
                        "continuous ancestral estimator generalized-least-squares requires model ou"
                    )
                if args.model == "ou" and args.estimator in {
                    "ace-pic",
                    "anc-ml",
                    "fast-anc",
                }:
                    parser.error(
                        "continuous ancestral estimators ace-pic, anc-ml, and fast-anc require model brownian"
                    )
                report = reconstruct_continuous_ancestral_states(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    estimator=args.estimator,
                    alpha=args.alpha,
                )
                summary = summarize_continuous_ancestral_report(report)
                exclusions = continuous_ancestral_exclusions(report)
                outputs: list[Path | str] = []
                if args.table_out is not None:
                    outputs.append(write_ancestral_state_table(args.table_out, report))
                if args.summary_out is not None:
                    outputs.append(
                        write_continuous_ancestral_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                if args.uncertainty_out is not None:
                    outputs.append(
                        write_continuous_ancestral_uncertainty_table(
                            args.uncertainty_out,
                            report,
                        )
                    )
                if args.exclusions_out is not None:
                    outputs.append(
                        write_continuous_ancestral_exclusion_table(
                            args.exclusions_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "estimate_count": len(report.estimates),
                            "internal_node_count": summary.internal_node_count,
                            "excluded_taxon_count": len(exclusions),
                            "unstable_node_count": summary.unstable_node_count,
                            "model": report.model,
                            "estimator": report.estimator,
                            "tree_is_ultrametric": summary.tree_is_ultrametric,
                            "covariance_near_singular": (
                                summary.covariance_near_singular
                            ),
                            "covariance_condition_number": (
                                summary.covariance_condition_number
                            ),
                            "log_likelihood": summary.log_likelihood,
                            "residual_sigma_squared": (summary.residual_sigma_squared),
                            "optimizer_name": summary.optimizer_name,
                            "optimizer_converged": summary.optimizer_converged,
                            "optimizer_iteration_count": (
                                summary.optimizer_iteration_count
                            ),
                            "optimizer_function_evaluation_count": (
                                summary.optimizer_function_evaluation_count
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.ancestral_command == "discrete":
                if args.state_ordering == "ordered" and args.model == "fitch":
                    parser.error(
                        "ordered ancestral discrete reconstruction requires a likelihood model"
                    )
                if args.compare_model == args.model:
                    parser.error(
                        "discrete ancestral compare-model must differ from the primary model"
                    )
                report = reconstruct_discrete_ancestral_states(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                    root_prior_mode=args.root_prior_mode,
                    fixed_root_state=args.fixed_root_state,
                )
                summary = summarize_discrete_ancestral_report(report)
                exclusions = discrete_ancestral_exclusions(report)
                comparison = (
                    compare_discrete_ancestral_reconstructions(
                        args.tree,
                        args.table,
                        trait=args.trait,
                        taxon_column=args.taxon_column,
                        left_model=args.model,
                        right_model=args.compare_model,
                        state_ordering=args.state_ordering,
                        ordered_states=_split_csv_values(args.ordered_states) or None,
                        root_prior_mode=args.root_prior_mode,
                        fixed_root_state=args.fixed_root_state,
                    )
                    if args.compare_model is not None
                    else None
                )
                outputs: list[Path | str] = []
                if args.table_out is not None:
                    outputs.append(write_ancestral_state_table(args.table_out, report))
                if args.summary_out is not None:
                    outputs.append(
                        write_discrete_ancestral_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                if args.probabilities_out is not None:
                    outputs.append(
                        write_discrete_ancestral_probability_table(
                            args.probabilities_out,
                            report,
                        )
                    )
                if args.transitions_out is not None:
                    outputs.append(
                        write_discrete_ancestral_transition_table(
                            args.transitions_out,
                            report,
                        )
                    )
                if args.fit_out is not None:
                    outputs.append(
                        write_discrete_ancestral_fit_table(
                            args.fit_out,
                            report,
                        )
                    )
                if args.comparison_out is not None and comparison is not None:
                    outputs.append(
                        write_discrete_ancestral_comparison_table(
                            args.comparison_out,
                            comparison,
                        )
                    )
                if args.exclusions_out is not None:
                    outputs.append(
                        write_discrete_ancestral_exclusion_table(
                            args.exclusions_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "estimate_count": len(report.estimates),
                            "internal_node_count": summary.internal_node_count,
                            "ambiguous_internal_node_count": (
                                summary.ambiguous_internal_node_count
                            ),
                            "excluded_taxon_count": len(exclusions),
                            "state_count": len(report.observed_states),
                            "minimal_change_count": report.minimal_change_count,
                            "parsimonious_root_state_count": (
                                report.parsimonious_root_state_count
                            ),
                            "unstable_node_count": summary.unstable_node_count,
                            "log_likelihood": report.log_likelihood,
                            "parameter_count": report.parameter_count,
                            "aic": report.aic,
                            "root_prior_mode": report.root_prior_mode,
                            "fixed_root_state": report.fixed_root_state,
                            "phytools_rerooting_method_comparable": (
                                report.rerooting_method_compatibility.comparable
                            ),
                            "optimizer_converged": (
                                None
                                if report.optimizer_diagnostics is None
                                else report.optimizer_diagnostics.converged
                            ),
                            "optimizer_iteration_count": (
                                None
                                if report.optimizer_diagnostics is None
                                else report.optimizer_diagnostics.iteration_count
                            ),
                            "optimizer_function_evaluation_count": (
                                None
                                if report.optimizer_diagnostics is None
                                else report.optimizer_diagnostics.function_evaluation_count
                            ),
                            "overparameterized": report.overparameterized,
                            "transition_rate_count": len(report.transition_rate_rows),
                            "baseline_model": (
                                None
                                if report.baseline_comparison is None
                                else report.baseline_comparison.baseline_model
                            ),
                            "baseline_delta_aic": (
                                None
                                if report.baseline_comparison is None
                                else report.baseline_comparison.delta_aic
                            ),
                            "preferred_model_by_aic": (
                                None
                                if report.baseline_comparison is None
                                else report.baseline_comparison.preferred_model_by_aic
                            ),
                            "comparison_node_count": (
                                len(comparison.rows) if comparison is not None else 0
                            ),
                            "comparison_differing_node_count": (
                                comparison.differing_node_count
                                if comparison is not None
                                else 0
                            ),
                            "model": report.model,
                        },
                        data={
                            "report": report,
                            "comparison": comparison,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            if args.ancestral_command == "discrete-reference":
                report = validate_discrete_ancestral_reference_examples()
                outputs = _finalize_outputs(args, command="ancestral", inputs=[])
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[],
                        outputs=outputs,
                        metrics={
                            "case_count": report.case_count,
                            "external_case_count": report.external_case_count,
                            "all_passed": report.all_passed,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.ancestral_command == "tree-set":
                if args.kind == "continuous":
                    resolved_model = args.model or "brownian"
                    if resolved_model not in {"brownian", "ou"}:
                        parser.error(
                            "continuous ancestral tree-set reconstruction requires model brownian or ou"
                        )
                    report = summarize_continuous_ancestral_tree_set(
                        args.tree_set,
                        args.table,
                        trait=args.trait,
                        taxon_column=args.taxon_column,
                        model=resolved_model,
                        alpha=args.alpha,
                        burnin_fraction=args.burnin_fraction,
                    )
                    summary = summarize_continuous_ancestral_tree_set_report(report)
                    outputs: list[Path | str] = []
                    if args.summary_out is not None:
                        outputs.append(
                            write_continuous_ancestral_tree_set_summary_table(
                                args.summary_out,
                                report,
                            )
                        )
                    if args.trees_out is not None:
                        outputs.append(
                            write_ancestral_tree_set_tree_table(args.trees_out, report)
                        )
                    if args.nodes_out is not None:
                        outputs.append(
                            write_continuous_ancestral_tree_set_node_table(
                                args.nodes_out,
                                report,
                            )
                        )
                    if args.clades_out is not None:
                        outputs.append(
                            write_continuous_ancestral_tree_set_clade_table(
                                args.clades_out,
                                report,
                            )
                        )
                    if args.exclusions_out is not None:
                        outputs.append(
                            write_ancestral_tree_set_exclusion_table(
                                args.exclusions_out,
                                report,
                            )
                        )
                    outputs = _finalize_outputs(
                        args,
                        command="ancestral",
                        inputs=[args.tree_set, args.table],
                        outputs=outputs,
                    )
                    _print_result(
                        build_command_result(
                            command="ancestral",
                            inputs=[args.tree_set, args.table],
                            outputs=outputs,
                            warnings=report.warnings,
                            metrics={
                                "kind": "continuous",
                                "model": report.model,
                                "total_tree_count": report.total_tree_count,
                                "kept_tree_count": report.kept_tree_count,
                                "rooted_topology_count": report.rooted_topology_count,
                                "unrooted_topology_count": report.unrooted_topology_count,
                                "node_row_count": len(report.node_rows),
                                "clade_summary_count": len(report.clade_summaries),
                                "excluded_taxon_count": len(report.exclusions),
                                "unstable_clade_count": summary.unstable_clade_count,
                            },
                            data=report,
                        ),
                        json_output=args.json,
                    )
                    return 0
                resolved_model = args.model or "fitch"
                if resolved_model not in {
                    "fitch",
                    "equal-rates",
                    "symmetric",
                    "all-rates-different",
                }:
                    parser.error(
                        "discrete ancestral tree-set reconstruction requires a discrete model"
                    )
                if args.state_ordering == "ordered" and resolved_model == "fitch":
                    parser.error(
                        "ordered ancestral tree-set discrete reconstruction requires a likelihood model"
                    )
                report = summarize_discrete_ancestral_tree_set(
                    args.tree_set,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=resolved_model,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                    burnin_fraction=args.burnin_fraction,
                )
                summary = summarize_discrete_ancestral_tree_set_report(report)
                outputs = []
                if args.summary_out is not None:
                    outputs.append(
                        write_discrete_ancestral_tree_set_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                if args.trees_out is not None:
                    outputs.append(
                        write_ancestral_tree_set_tree_table(args.trees_out, report)
                    )
                if args.nodes_out is not None:
                    outputs.append(
                        write_discrete_ancestral_tree_set_node_table(
                            args.nodes_out,
                            report,
                        )
                    )
                if args.clades_out is not None:
                    outputs.append(
                        write_discrete_ancestral_tree_set_clade_table(
                            args.clades_out,
                            report,
                        )
                    )
                if args.exclusions_out is not None:
                    outputs.append(
                        write_ancestral_tree_set_exclusion_table(
                            args.exclusions_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="ancestral",
                    inputs=[args.tree_set, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree_set, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "kind": "discrete",
                            "model": report.model,
                            "total_tree_count": report.total_tree_count,
                            "kept_tree_count": report.kept_tree_count,
                            "rooted_topology_count": report.rooted_topology_count,
                            "unrooted_topology_count": report.unrooted_topology_count,
                            "node_row_count": len(report.node_rows),
                            "clade_summary_count": len(report.clade_summaries),
                            "excluded_taxon_count": len(report.exclusions),
                            "unstable_clade_count": summary.unstable_clade_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.ancestral_command == "confidence":
                if not args.tree_set and args.burnin_fraction != 0.0:
                    parser.error("--burnin-fraction requires --tree-set")
                if args.kind == "continuous":
                    resolved_model = args.model or "brownian"
                    if resolved_model not in {"brownian", "ou"}:
                        parser.error(
                            "continuous ancestral confidence review requires model brownian or ou"
                        )
                    if args.tree_set:
                        report = summarize_continuous_ancestral_tree_set(
                            args.tree,
                            args.table,
                            trait=args.trait,
                            taxon_column=args.taxon_column,
                            model=resolved_model,
                            alpha=args.alpha,
                            burnin_fraction=args.burnin_fraction,
                        )
                        confidence_rows = (
                            build_continuous_ancestral_tree_set_confidence_rows(report)
                        )
                        confidence_summary = (
                            summarize_continuous_ancestral_tree_set_confidence(report)
                        )
                        outputs: list[Path | str] = []
                        if args.summary_out is not None:
                            outputs.append(
                                write_ancestral_confidence_summary_table(
                                    args.summary_out,
                                    confidence_summary,
                                )
                            )
                        if args.confidence_out is not None:
                            outputs.append(
                                write_continuous_ancestral_tree_set_confidence_table(
                                    args.confidence_out,
                                    report,
                                )
                            )
                        outputs = _finalize_outputs(
                            args,
                            command="ancestral",
                            inputs=[args.tree, args.table],
                            outputs=outputs,
                        )
                        _print_result(
                            build_command_result(
                                command="ancestral",
                                inputs=[args.tree, args.table],
                                outputs=outputs,
                                warnings=report.warnings,
                                metrics={
                                    "kind": "continuous",
                                    "source_kind": "tree_set",
                                    "model": report.model,
                                    "kept_tree_count": report.kept_tree_count,
                                    "confidence_row_count": len(confidence_rows),
                                    "low_confidence_count": (
                                        confidence_summary.low_confidence_count
                                    ),
                                    "unstable_count": confidence_summary.unstable_count,
                                    "high_entropy_count": (
                                        confidence_summary.high_entropy_count
                                    ),
                                    "top_uncertain_id": (
                                        confidence_summary.top_uncertain_id
                                    ),
                                },
                                data={
                                    "report": report,
                                    "confidence_summary": confidence_summary,
                                },
                            ),
                            json_output=args.json,
                        )
                        return 0
                    report = reconstruct_continuous_ancestral_states(
                        args.tree,
                        args.table,
                        trait=args.trait,
                        taxon_column=args.taxon_column,
                        model=resolved_model,
                        alpha=args.alpha,
                    )
                    confidence_rows = build_continuous_ancestral_confidence_rows(report)
                    confidence_summary = summarize_continuous_ancestral_confidence(
                        report
                    )
                    outputs = []
                    if args.summary_out is not None:
                        outputs.append(
                            write_ancestral_confidence_summary_table(
                                args.summary_out,
                                confidence_summary,
                            )
                        )
                    if args.confidence_out is not None:
                        outputs.append(
                            write_continuous_ancestral_confidence_table(
                                args.confidence_out,
                                report,
                            )
                        )
                    outputs = _finalize_outputs(
                        args,
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                    )
                    _print_result(
                        build_command_result(
                            command="ancestral",
                            inputs=[args.tree, args.table],
                            outputs=outputs,
                            warnings=report.warnings,
                            metrics={
                                "kind": "continuous",
                                "source_kind": "tree",
                                "model": report.model,
                                "confidence_row_count": len(confidence_rows),
                                "low_confidence_count": (
                                    confidence_summary.low_confidence_count
                                ),
                                "unstable_count": confidence_summary.unstable_count,
                                "high_entropy_count": (
                                    confidence_summary.high_entropy_count
                                ),
                                "top_uncertain_id": confidence_summary.top_uncertain_id,
                            },
                            data={
                                "report": report,
                                "confidence_summary": confidence_summary,
                            },
                        ),
                        json_output=args.json,
                    )
                    return 0
                resolved_model = args.model or "fitch"
                if resolved_model not in {
                    "fitch",
                    "equal-rates",
                    "symmetric",
                    "all-rates-different",
                }:
                    parser.error(
                        "discrete ancestral confidence review requires a discrete model"
                    )
                if args.state_ordering == "ordered" and resolved_model == "fitch":
                    parser.error(
                        "ordered ancestral confidence review requires a likelihood model"
                    )
                if args.tree_set:
                    report = summarize_discrete_ancestral_tree_set(
                        args.tree,
                        args.table,
                        trait=args.trait,
                        taxon_column=args.taxon_column,
                        model=resolved_model,
                        state_ordering=args.state_ordering,
                        ordered_states=_split_csv_values(args.ordered_states) or None,
                        burnin_fraction=args.burnin_fraction,
                    )
                    confidence_rows = build_discrete_ancestral_tree_set_confidence_rows(
                        report
                    )
                    confidence_summary = (
                        summarize_discrete_ancestral_tree_set_confidence(report)
                    )
                    outputs = []
                    if args.summary_out is not None:
                        outputs.append(
                            write_ancestral_confidence_summary_table(
                                args.summary_out,
                                confidence_summary,
                            )
                        )
                    if args.confidence_out is not None:
                        outputs.append(
                            write_discrete_ancestral_tree_set_confidence_table(
                                args.confidence_out,
                                report,
                            )
                        )
                    outputs = _finalize_outputs(
                        args,
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                    )
                    _print_result(
                        build_command_result(
                            command="ancestral",
                            inputs=[args.tree, args.table],
                            outputs=outputs,
                            warnings=report.warnings,
                            metrics={
                                "kind": "discrete",
                                "source_kind": "tree_set",
                                "model": report.model,
                                "kept_tree_count": report.kept_tree_count,
                                "confidence_row_count": len(confidence_rows),
                                "low_confidence_count": (
                                    confidence_summary.low_confidence_count
                                ),
                                "unstable_count": confidence_summary.unstable_count,
                                "high_entropy_count": (
                                    confidence_summary.high_entropy_count
                                ),
                                "top_uncertain_id": (
                                    confidence_summary.top_uncertain_id
                                ),
                            },
                            data={
                                "report": report,
                                "confidence_summary": confidence_summary,
                            },
                        ),
                        json_output=args.json,
                    )
                    return 0
                report = reconstruct_discrete_ancestral_states(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=resolved_model,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                )
                confidence_rows = build_discrete_ancestral_confidence_rows(report)
                confidence_summary = summarize_discrete_ancestral_confidence(report)
                outputs = []
                if args.summary_out is not None:
                    outputs.append(
                        write_ancestral_confidence_summary_table(
                            args.summary_out,
                            confidence_summary,
                        )
                    )
                if args.confidence_out is not None:
                    outputs.append(
                        write_discrete_ancestral_confidence_table(
                            args.confidence_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "kind": "discrete",
                            "source_kind": "tree",
                            "model": report.model,
                            "confidence_row_count": len(confidence_rows),
                            "low_confidence_count": (
                                confidence_summary.low_confidence_count
                            ),
                            "unstable_count": confidence_summary.unstable_count,
                            "high_entropy_count": (
                                confidence_summary.high_entropy_count
                            ),
                            "top_uncertain_id": confidence_summary.top_uncertain_id,
                        },
                        data={
                            "report": report,
                            "confidence_summary": confidence_summary,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            if args.ancestral_command == "root-sensitivity":
                report = summarize_ancestral_root_sensitivity(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                    fixed_root_state=args.fixed_root_state,
                )
                summary = summarize_ancestral_root_sensitivity_report(report)
                outputs = []
                if args.summary_out is not None:
                    outputs.append(
                        write_ancestral_root_sensitivity_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                if args.assumptions_out is not None:
                    outputs.append(
                        write_ancestral_root_assumption_table(
                            args.assumptions_out,
                            report,
                        )
                    )
                if args.nodes_out is not None:
                    outputs.append(
                        write_ancestral_root_sensitivity_node_table(
                            args.nodes_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "model": report.model,
                            "state_ordering": report.state_ordering,
                            "analyzed_taxon_count": report.analyzed_taxon_count,
                            "assumption_count": summary.assumption_count,
                            "compared_node_count": summary.compared_node_count,
                            "state_changed_node_count": (
                                summary.state_changed_node_count
                            ),
                            "support_changed_node_count": (
                                summary.support_changed_node_count
                            ),
                            "top_sensitive_node": summary.top_sensitive_node,
                            "fixed_root_state": report.fixed_root_state,
                        },
                        data={
                            "report": report,
                            "summary": summary,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            if args.ancestral_command == "ordered-discrete":
                report = summarize_ordered_discrete_reconstruction(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    ordered_states=_split_csv_values(args.ordered_states) or [],
                )
                summary = summarize_ordered_discrete_report(report)
                outputs = []
                if args.summary_out is not None:
                    outputs.append(
                        write_ordered_discrete_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                if args.fits_out is not None:
                    outputs.append(
                        write_ordered_discrete_fit_table(
                            args.fits_out,
                            report,
                        )
                    )
                if args.nodes_out is not None:
                    outputs.append(
                        write_ordered_discrete_node_table(
                            args.nodes_out,
                            report,
                        )
                    )
                if args.transitions_out is not None:
                    outputs.append(
                        write_ordered_discrete_transition_table(
                            args.transitions_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "model": report.model,
                            "ordered_state_count": len(report.ordered_states),
                            "fit_count": len(report.fit_rows),
                            "differing_node_count": summary.differing_node_count,
                            "ambiguity_change_count": (summary.ambiguity_change_count),
                            "restricted_transition_count": (
                                summary.restricted_transition_count
                            ),
                            "preferred_ordering": summary.preferred_ordering,
                        },
                        data={
                            "report": report,
                            "summary": summary,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            if args.ancestral_command == "irreversible-discrete":
                report = summarize_irreversible_discrete_reconstruction(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    allowed_transition_pairs=_parse_transition_pairs(
                        args.allowed_transitions
                    ),
                )
                summary = summarize_irreversible_discrete_report(report)
                outputs = []
                if args.summary_out is not None:
                    outputs.append(
                        write_irreversible_discrete_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                if args.fits_out is not None:
                    outputs.append(
                        write_irreversible_discrete_fit_table(
                            args.fits_out,
                            report,
                        )
                    )
                if args.nodes_out is not None:
                    outputs.append(
                        write_irreversible_discrete_node_table(
                            args.nodes_out,
                            report,
                        )
                    )
                if args.transitions_out is not None:
                    outputs.append(
                        write_irreversible_discrete_transition_table(
                            args.transitions_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "model": report.model,
                            "allowed_transition_count": len(
                                report.allowed_transition_pairs
                            ),
                            "fit_count": len(report.fit_rows),
                            "differing_node_count": summary.differing_node_count,
                            "ambiguity_change_count": (summary.ambiguity_change_count),
                            "forbidden_transition_count": (
                                summary.forbidden_transition_count
                            ),
                            "preferred_constraint": summary.preferred_constraint,
                        },
                        data={
                            "report": report,
                            "summary": summary,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            if args.ancestral_command == "transitions":
                if args.state_ordering == "ordered" and args.model == "fitch":
                    parser.error(
                        "ordered ancestral transition counting requires a likelihood model"
                    )
                if not args.tree_set and args.trees_out is not None:
                    parser.error("--trees-out requires --tree-set")
                if not args.tree_set and args.burnin_fraction != 0.0:
                    parser.error("--burnin-fraction requires --tree-set")
                if args.tree_set:
                    report = summarize_ancestral_transition_tree_set(
                        args.tree,
                        args.table,
                        trait=args.trait,
                        taxon_column=args.taxon_column,
                        model=args.model,
                        state_ordering=args.state_ordering,
                        ordered_states=_split_csv_values(args.ordered_states) or None,
                        burnin_fraction=args.burnin_fraction,
                    )
                    summary = summarize_ancestral_transition_tree_set_report(report)
                    outputs = []
                    if args.summary_out is not None:
                        outputs.append(
                            write_ancestral_transition_tree_set_summary_table(
                                args.summary_out,
                                report,
                            )
                        )
                    if args.trees_out is not None:
                        outputs.append(
                            write_ancestral_transition_tree_set_tree_table(
                                args.trees_out,
                                report,
                            )
                        )
                    if args.branches_out is not None:
                        outputs.append(
                            write_ancestral_transition_tree_set_branch_table(
                                args.branches_out,
                                report,
                            )
                        )
                    if args.counts_out is not None:
                        outputs.append(
                            write_ancestral_transition_tree_set_count_table(
                                args.counts_out,
                                report,
                            )
                        )
                    if args.exclusions_out is not None:
                        outputs.append(
                            write_ancestral_transition_exclusion_table(
                                args.exclusions_out,
                                report,
                            )
                        )
                    outputs = _finalize_outputs(
                        args,
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                    )
                    _print_result(
                        build_command_result(
                            command="ancestral",
                            inputs=[args.tree, args.table],
                            outputs=outputs,
                            warnings=report.warnings,
                            metrics={
                                "tree_set": True,
                                "model": report.model,
                                "total_tree_count": report.total_tree_count,
                                "kept_tree_count": report.kept_tree_count,
                                "rooted_topology_count": report.rooted_topology_count,
                                "unrooted_topology_count": (
                                    report.unrooted_topology_count
                                ),
                                "transition_pair_count": len(report.transition_rows),
                                "topology_sensitive_transition_pair_count": (
                                    summary.topology_sensitive_transition_pair_count
                                ),
                                "uncertainty_sensitive_transition_pair_count": (
                                    summary.uncertainty_sensitive_transition_pair_count
                                ),
                                "excluded_taxon_count": len(report.exclusions),
                            },
                            data=report,
                        ),
                        json_output=args.json,
                    )
                    return 0
                report = summarize_ancestral_transitions(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                )
                summary = summarize_ancestral_transition_report(report)
                outputs = []
                if args.summary_out is not None:
                    outputs.append(
                        write_ancestral_transition_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                if args.branches_out is not None:
                    outputs.append(
                        write_ancestral_transition_branch_table(
                            args.branches_out,
                            report,
                        )
                    )
                if args.counts_out is not None:
                    outputs.append(
                        write_ancestral_transition_count_table(
                            args.counts_out,
                            report,
                        )
                    )
                if args.exclusions_out is not None:
                    outputs.append(
                        write_ancestral_transition_exclusion_table(
                            args.exclusions_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "tree_set": False,
                            "model": report.model,
                            "total_branch_count": summary.total_branch_count,
                            "changed_branch_count": summary.changed_branch_count,
                            "certain_change_count": summary.certain_change_count,
                            "uncertain_change_count": summary.uncertain_change_count,
                            "transition_pair_count": len(report.transition_rows),
                            "excluded_taxon_count": len(report.exclusions),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.ancestral_command == "compare":
                report = compare_continuous_ancestral_models(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    left_model=args.left_model,
                    right_model=args.right_model,
                    left_alpha=args.left_alpha,
                    right_alpha=args.right_alpha,
                )
                outputs = _finalize_outputs(
                    args, command="ancestral", inputs=[args.tree, args.table]
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "compared_node_count": len(report.rows),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.ancestral_command == "sensitivity":
                _validate_ancestral_discrete_model_arguments(args, parser)
                resolved_model = args.model or (
                    "brownian" if args.kind == "continuous" else "fitch"
                )
                report = build_ancestral_sensitivity_report(
                    tree_path=args.tree,
                    traits_path=args.table,
                    trait=args.trait,
                    reconstruction_kind=args.kind,
                    model=resolved_model,
                    taxon_column=args.taxon_column,
                    alpha=args.alpha,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                    compare_tree_path=args.compare_tree,
                    compare_model=args.compare_model,
                    drop_taxa=args.drop_taxa,
                    coding_map=_parse_assignment_map(args.coding_map) or None,
                )
                outputs = _finalize_outputs(
                    args, command="ancestral", inputs=[args.tree, args.table]
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "baseline_node_count": report.baseline_node_count,
                            "has_model_sensitivity": report.model_sensitivity
                            is not None,
                            "has_tree_sensitivity": report.tree_sensitivity is not None,
                            "has_pruning_sensitivity": report.pruning_sensitivity
                            is not None,
                            "has_trait_coding_sensitivity": report.trait_coding_sensitivity
                            is not None,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.ancestral_command == "render":
                _validate_ancestral_discrete_model_arguments(args, parser)
                if args.kind == "continuous" and args.branch_coloring == "state":
                    parser.error(
                        "continuous ancestral rendering does not support branch coloring 'state'"
                    )
                if args.kind == "discrete" and args.branch_coloring == "regime":
                    parser.error(
                        "discrete ancestral rendering does not support branch coloring 'regime'"
                    )
                if args.kind == "continuous":
                    resolved_model = args.model or "brownian"
                    reconstruction = reconstruct_continuous_ancestral_states(
                        args.tree,
                        args.table,
                        trait=args.trait,
                        taxon_column=args.taxon_column,
                        model=resolved_model,
                        alpha=args.alpha,
                    )
                else:
                    resolved_model = args.model or "fitch"
                    reconstruction = reconstruct_discrete_ancestral_states(
                        args.tree,
                        args.table,
                        trait=args.trait,
                        taxon_column=args.taxon_column,
                        model=resolved_model,
                        state_ordering=args.state_ordering,
                        ordered_states=_split_csv_values(args.ordered_states) or None,
                    )
                result = render_ancestral_state_visualization(
                    args.tree,
                    reconstruction,
                    out_path=args.out,
                    layout=args.layout,
                    discrete_node_style=args.discrete_node_style,
                    branch_coloring=args.branch_coloring,
                )
                rendered_outputs = (
                    [result.output_path]
                    if result.format == "svg"
                    else [result.output_path, result.svg_path]
                )
                outputs = _finalize_outputs(
                    args,
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=rendered_outputs,
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=getattr(reconstruction, "warnings", []),
                        metrics={
                            "tip_count": result.tree_render.tip_count,
                            "format": result.format,
                            "layout": result.layout,
                            "rendered_internal_annotation_count": (
                                result.tree_render.rendered_internal_annotation_count
                            ),
                            "rendered_internal_pie_count": (
                                result.tree_render.rendered_internal_pie_count
                            ),
                            "rendered_branch_color_count": (
                                result.tree_render.rendered_branch_color_count
                            ),
                        },
                        data={
                            "reconstruction": reconstruction,
                            "visualization": result,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            resolved_model = args.model or (
                "brownian" if args.kind == "continuous" else "fitch"
            )
            if args.ancestral_command == "package":
                _validate_ancestral_discrete_model_arguments(args, parser)
                result = build_ancestral_figure_package(
                    tree_path=args.tree,
                    traits_path=args.table,
                    trait=args.trait,
                    reconstruction_kind=args.kind,
                    out_dir=args.out_dir,
                    taxon_column=args.taxon_column,
                    model=resolved_model,
                    alpha=args.alpha,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                    layout=args.layout,
                )
                outputs = _finalize_outputs(
                    args,
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=[
                        result.figure_path,
                        result.figure_png_path,
                        result.figure_html_path,
                        result.review_path,
                        result.node_table_path,
                        result.uncertainty_table_path,
                        result.node_review_path,
                        result.legend_path,
                        result.model_description_path,
                        result.caption_path,
                        result.manifest_path,
                        result.reproducibility_manifest_path,
                    ],
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "output_dir": str(result.output_dir),
                            "artifact_count": 12,
                            "publication_ready": result.audit.publication_ready,
                            "internal_state_visible": result.audit.internal_state_visible,
                            "uncertainty_visible": result.audit.uncertainty_visible,
                            "ambiguous_internal_node_count": (
                                result.audit.ambiguous_internal_node_count
                            ),
                            "unstable_internal_node_count": (
                                result.audit.unstable_internal_node_count
                            ),
                            "rendered_internal_annotation_count": (
                                result.audit.rendered_internal_annotation_count
                            ),
                            "rendered_internal_pie_count": (
                                result.audit.rendered_internal_pie_count
                            ),
                        },
                        data=result,
                    ),
                    json_output=args.json,
                )
                return 0
            _validate_ancestral_discrete_model_arguments(args, parser)
            if args.out is None and args.out_dir is None:
                parser.error("ancestral report requires --out or --out-dir")
            if args.out_dir is not None:
                result = build_ancestral_report_package(
                    tree_path=args.tree,
                    traits_path=args.table,
                    trait=args.trait,
                    reconstruction_kind=args.kind,
                    out_dir=args.out_dir,
                    taxon_column=args.taxon_column,
                    model=resolved_model,
                    alpha=args.alpha,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                    compare_model=args.compare_model,
                    compare_tree_path=args.compare_tree,
                    drop_taxa=args.drop_taxa,
                    coding_map=_parse_assignment_map(args.coding_map) or None,
                )
                output_paths: list[Path | str] = [
                    result.report_path,
                    result.methods_summary_path,
                    result.reviewer_audit_checklist_path,
                    result.figure_path,
                    result.figure_png_path,
                    result.figure_html_path,
                    result.summary_table_path,
                    result.node_table_path,
                    result.uncertainty_table_path,
                    result.transition_count_table_path,
                    result.transition_branch_table_path,
                    result.exclusion_table_path,
                    result.manifest_path,
                ]
                if args.out is not None:
                    args.out.parent.mkdir(parents=True, exist_ok=True)
                    args.out.write_text(
                        result.report_path.read_text(encoding="utf-8"),
                        encoding="utf-8",
                    )
                    svg_out = args.out.with_suffix(".svg")
                    svg_out.write_text(
                        result.figure_path.read_text(encoding="utf-8"),
                        encoding="utf-8",
                    )
                    output_paths.extend([args.out, svg_out])
                outputs = _finalize_outputs(
                    args,
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=output_paths,
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "report_kind": "ancestral-report-package",
                            "reconstruction_kind": result.reconstruction_kind,
                            "output_dir": str(result.output_dir),
                            "artifact_count": 13,
                            "methods_summary_warning_count": (
                                result.methods_summary.warning_count
                            ),
                            "transition_count_row_count": result.machine_manifest[
                                "metrics"
                            ]["transition_count_row_count"],
                        },
                        data=result,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.out is None:
                raise ValueError(
                    "ancestral report rendering requires an explicit output path"
                )
            result = render_ancestral_state_report(
                tree_path=args.tree,
                traits_path=args.table,
                trait=args.trait,
                reconstruction_kind=args.kind,
                out_path=args.out,
                taxon_column=args.taxon_column,
                model=resolved_model,
                alpha=args.alpha,
                state_ordering=args.state_ordering,
                ordered_states=_split_csv_values(args.ordered_states) or None,
                compare_model=args.compare_model,
                compare_tree_path=args.compare_tree,
                drop_taxa=args.drop_taxa,
                coding_map=_parse_assignment_map(args.coding_map) or None,
            )
            outputs = _finalize_outputs(
                args,
                command="ancestral",
                inputs=[args.tree, args.table],
                outputs=[result.output_path, args.out.with_suffix(".svg")],
            )
            _print_result(
                build_command_result(
                    command="ancestral",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                    metrics={
                        "report_kind": result.report_kind,
                        "reconstruction_kind": result.reconstruction_kind,
                    },
                    data=result,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "biogeography":
            if args.biogeography_command == "model":
                report = summarize_geographic_state_model(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    allowed_regions=_split_csv_values(args.allowed_regions) or None,
                )
                outputs: list[Path | str] = []
                if args.summary_out is not None:
                    outputs.append(
                        write_geographic_state_summary_table(args.summary_out, report)
                    )
                if args.nodes_out is not None:
                    outputs.append(
                        write_geographic_region_probability_table(
                            args.nodes_out, report
                        )
                    )
                if args.rates_out is not None:
                    outputs.append(
                        write_geographic_transition_rate_table(args.rates_out, report)
                    )
                if args.events_out is not None:
                    outputs.append(
                        write_geographic_transition_event_table(args.events_out, report)
                    )
                if args.exclusions_out is not None:
                    outputs.append(
                        write_geographic_exclusion_table(args.exclusions_out, report)
                    )
                outputs = _finalize_outputs(
                    args,
                    command="biogeography",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="biogeography",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "model": report.model,
                            "observed_region_count": (
                                report.summary.observed_region_count
                            ),
                            "internal_node_count": report.summary.internal_node_count,
                            "transition_rate_row_count": (
                                report.summary.transition_rate_row_count
                            ),
                            "changed_branch_count": (
                                report.summary.changed_branch_count
                            ),
                            "strongly_supported_transition_count": (
                                report.summary.strongly_supported_transition_count
                            ),
                            "excluded_taxon_count": (
                                report.summary.excluded_taxon_count
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.biogeography_command == "constrained":
                report = summarize_constrained_geographic_model(
                    args.tree,
                    args.table,
                    args.adjacency,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                )
                summary = summarize_constrained_geographic_report(report)
                outputs: list[Path | str] = []
                if args.summary_out is not None:
                    outputs.append(
                        write_constrained_geographic_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                if args.fits_out is not None:
                    outputs.append(
                        write_constrained_geographic_fit_table(
                            args.fits_out,
                            report,
                        )
                    )
                if args.transitions_out is not None:
                    outputs.append(
                        write_constrained_geographic_transition_table(
                            args.transitions_out,
                            report,
                        )
                    )
                if args.unsupported_out is not None:
                    outputs.append(
                        write_unsupported_geographic_transition_claim_table(
                            args.unsupported_out,
                            report,
                        )
                    )
                if args.exclusions_out is not None:
                    outputs.append(
                        write_constrained_geographic_exclusion_table(
                            args.exclusions_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="biogeography",
                    inputs=[args.tree, args.table, args.adjacency],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="biogeography",
                        inputs=[args.tree, args.table, args.adjacency],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "model": report.model,
                            "allowed_transition_count": (
                                summary.allowed_transition_count
                            ),
                            "forbidden_transition_count": (
                                summary.forbidden_transition_count
                            ),
                            "unsupported_transition_claim_count": (
                                summary.unsupported_transition_claim_count
                            ),
                            "preferred_constraint": summary.preferred_constraint,
                            "excluded_taxon_count": summary.excluded_taxon_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.biogeography_command == "time-stratified":
                report = summarize_time_stratified_geographic_transitions(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    allowed_regions=_split_csv_values(args.allowed_regions) or None,
                    time_bins=[
                        _parse_time_bin_definition(raw_time_bin)
                        for raw_time_bin in args.time_bin
                    ],
                )
                outputs: list[Path | str] = []
                if args.summary_out is not None:
                    outputs.append(
                        write_time_stratified_transition_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                if args.matrix_out is not None:
                    outputs.append(
                        write_time_stratified_transition_matrix_table(
                            args.matrix_out,
                            report,
                        )
                    )
                if args.branches_out is not None:
                    outputs.append(
                        write_time_stratified_branch_table(
                            args.branches_out,
                            report,
                        )
                    )
                if args.exclusions_out is not None:
                    outputs.append(
                        write_time_stratified_exclusion_table(
                            args.exclusions_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="biogeography",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="biogeography",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "model": report.model,
                            "time_bin_count": report.summary.time_bin_count,
                            "matrix_row_count": report.summary.matrix_row_count,
                            "changed_branch_count": (
                                report.summary.changed_branch_count
                            ),
                            "allocated_transition_weight_total": (
                                report.summary.allocated_transition_weight_total
                            ),
                            "excluded_taxon_count": (
                                report.summary.excluded_taxon_count
                            ),
                            "warning_count": report.summary.warning_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.biogeography_command == "sampling-bias":
                report = summarize_geographic_sampling_bias(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    allowed_regions=_split_csv_values(args.allowed_regions) or None,
                    weights_path=args.weights,
                    region_column=args.region_column,
                    weight_column=args.weight_column,
                )
                outputs: list[Path | str] = []
                if args.summary_out is not None:
                    outputs.append(
                        write_geographic_sampling_bias_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                if args.regions_out is not None:
                    outputs.append(
                        write_geographic_sampling_count_table(
                            args.regions_out,
                            report,
                        )
                    )
                if args.nodes_out is not None:
                    outputs.append(
                        write_geographic_sampling_bias_node_table(
                            args.nodes_out,
                            report,
                        )
                    )
                if args.transitions_out is not None:
                    outputs.append(
                        write_geographic_sampling_bias_transition_table(
                            args.transitions_out,
                            report,
                        )
                    )
                if args.exclusions_out is not None:
                    outputs.append(
                        write_geographic_sampling_bias_exclusion_table(
                            args.exclusions_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="biogeography",
                    inputs=(
                        [args.tree, args.table]
                        if args.weights is None
                        else [args.tree, args.table, args.weights]
                    ),
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="biogeography",
                        inputs=(
                            [args.tree, args.table]
                            if args.weights is None
                            else [args.tree, args.table, args.weights]
                        ),
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "model": report.summary.model,
                            "weighting_mode": report.summary.weighting_mode,
                            "region_dominated": report.summary.region_dominated,
                            "dominant_region": report.summary.dominant_region,
                            "dominant_region_fraction": (
                                report.summary.dominant_region_fraction
                            ),
                            "root_region_changed": report.summary.root_region_changed,
                            "changed_internal_node_count": (
                                report.summary.changed_internal_node_count
                            ),
                            "changed_transition_count": (
                                report.summary.changed_transition_count
                            ),
                            "excluded_taxon_count": (
                                report.summary.excluded_taxon_count
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.biogeography_command == "chronology":
                report = summarize_biogeographic_transition_chronology(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    allowed_regions=_split_csv_values(args.allowed_regions) or None,
                    time_bin_count=args.time_bin_count,
                )
                outputs: list[Path | str] = []
                if args.summary_out is not None:
                    outputs.append(
                        write_dated_biogeography_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                if args.nodes_out is not None:
                    outputs.append(
                        write_dated_biogeography_node_table(
                            args.nodes_out,
                            report,
                        )
                    )
                if args.events_out is not None:
                    outputs.append(
                        write_dated_biogeography_event_table(
                            args.events_out,
                            report,
                        )
                    )
                if args.bins_out is not None:
                    outputs.append(
                        write_dated_biogeography_time_bin_table(
                            args.bins_out,
                            report,
                        )
                    )
                if args.exclusions_out is not None:
                    outputs.append(
                        write_dated_biogeography_exclusion_table(
                            args.exclusions_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="biogeography",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="biogeography",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "model": report.summary.model,
                            "tree_is_time_scaled": (report.summary.tree_is_time_scaled),
                            "root_age": report.summary.root_age,
                            "event_count": report.summary.event_count,
                            "time_bin_count": report.summary.time_bin_count,
                            "high_uncertainty_bin_count": (
                                report.summary.high_uncertainty_bin_count
                            ),
                            "excluded_taxon_count": (
                                report.summary.excluded_taxon_count
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.biogeography_command == "events":
                outputs: list[Path | str] = []
                if args.tree_set:
                    report = summarize_geographic_migration_event_tree_set(
                        args.tree,
                        args.table,
                        trait=args.trait,
                        taxon_column=args.taxon_column,
                        model=args.model,
                        allowed_regions=_split_csv_values(args.allowed_regions) or None,
                        burnin_fraction=args.burnin_fraction,
                    )
                    if args.summary_out is not None:
                        outputs.append(
                            write_geographic_migration_tree_set_summary_table(
                                args.summary_out,
                                report,
                            )
                        )
                    if args.events_out is not None:
                        outputs.append(
                            write_geographic_migration_tree_set_event_table(
                                args.events_out,
                                report,
                            )
                        )
                    if args.trees_out is not None:
                        outputs.append(
                            write_geographic_migration_tree_set_tree_table(
                                args.trees_out,
                                report,
                            )
                        )
                    if args.event_summaries_out is not None:
                        outputs.append(
                            write_geographic_migration_tree_set_event_summary_table(
                                args.event_summaries_out,
                                report,
                            )
                        )
                    if args.exclusions_out is not None:
                        outputs.append(
                            write_geographic_migration_tree_set_exclusion_table(
                                args.exclusions_out,
                                report,
                            )
                        )
                    outputs = _finalize_outputs(
                        args,
                        command="biogeography",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                    )
                    _print_result(
                        build_command_result(
                            command="biogeography",
                            inputs=[args.tree, args.table],
                            outputs=outputs,
                            warnings=report.warnings,
                            metrics={
                                "report_mode": "tree_set",
                                "model": report.model,
                                "kept_tree_count": report.summary.kept_tree_count,
                                "event_row_count": report.summary.event_row_count,
                                "event_summary_count": (
                                    report.summary.event_summary_count
                                ),
                                "topology_sensitive_event_count": (
                                    report.summary.topology_sensitive_event_count
                                ),
                                "excluded_taxon_count": (
                                    report.summary.excluded_taxon_count
                                ),
                                "warning_count": report.summary.warning_count,
                            },
                            data=report,
                        ),
                        json_output=args.json,
                    )
                    return 0
                report = summarize_geographic_migration_events(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    allowed_regions=_split_csv_values(args.allowed_regions) or None,
                )
                if args.summary_out is not None:
                    outputs.append(
                        write_geographic_migration_event_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                if args.events_out is not None:
                    outputs.append(
                        write_geographic_migration_event_table(
                            args.events_out,
                            report,
                        )
                    )
                if args.exclusions_out is not None:
                    outputs.append(
                        write_geographic_migration_exclusion_table(
                            args.exclusions_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="biogeography",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="biogeography",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "report_mode": "single_tree",
                            "model": report.model,
                            "event_count": report.summary.event_count,
                            "strongly_supported_event_count": (
                                report.summary.strongly_supported_event_count
                            ),
                            "mean_event_support": (report.summary.mean_event_support),
                            "excluded_taxon_count": (
                                report.summary.excluded_taxon_count
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.biogeography_command == "report":
                result = build_biogeography_report_package(
                    tree_path=args.tree,
                    traits_path=args.table,
                    centroids_path=args.centroids,
                    trait=args.trait,
                    out_dir=args.out_dir,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    region_column=args.region_column,
                    latitude_column=args.latitude_column,
                    longitude_column=args.longitude_column,
                )
                outputs = _finalize_outputs(
                    args,
                    command="biogeography",
                    inputs=[args.tree, args.table, args.centroids],
                    outputs=[
                        result.report_path,
                        result.tree_figure_path,
                        result.map_path,
                        result.legend_path,
                        result.caption_path,
                        result.summary_table_path,
                        result.region_count_table_path,
                        result.node_table_path,
                        result.transition_matrix_path,
                        result.event_table_path,
                        result.map_marker_table_path,
                        result.map_line_table_path,
                        result.exclusion_table_path,
                        result.manifest_path,
                        result.reproducibility_manifest_path,
                    ],
                )
                _print_result(
                    build_command_result(
                        command="biogeography",
                        inputs=[args.tree, args.table, args.centroids],
                        outputs=outputs,
                        warnings=result.warnings,
                        metrics={
                            "report_kind": "biogeography-report-package",
                            "model": result.state_report.model,
                            "output_dir": str(result.output_dir),
                            "artifact_count": 15,
                            "observed_region_count": (
                                result.state_report.summary.observed_region_count
                            ),
                            "transition_rate_row_count": (
                                result.state_report.summary.transition_rate_row_count
                            ),
                            "event_count": result.event_report.summary.event_count,
                            "visible_map_line_count": (
                                result.map_report.summary.visible_line_count
                            ),
                            "publication_ready": result.audit.publication_ready,
                            "legend_entry_count": result.audit.legend_entry_count,
                            "caption_ready": result.audit.caption_ready,
                            "rendered_internal_pie_count": (
                                result.audit.rendered_internal_pie_count
                            ),
                        },
                        data=result,
                    ),
                    json_output=args.json,
                )
                return 0
        if (
            args.command == "host-association"
            and args.host_association_command == "switches"
        ):
            report = summarize_host_switching(
                args.tree,
                args.table,
                trait=args.trait,
                taxon_column=args.taxon_column,
                model=args.model,
                constraint_path=args.constraints,
            )
            outputs: list[Path | str] = []
            if args.summary_out is not None:
                outputs.append(
                    write_host_switch_summary_table(
                        args.summary_out,
                        report,
                    )
                )
            if args.nodes_out is not None:
                outputs.append(
                    write_host_state_node_table(
                        args.nodes_out,
                        report,
                    )
                )
            if args.branches_out is not None:
                outputs.append(
                    write_host_switch_branch_table(
                        args.branches_out,
                        report,
                    )
                )
            if args.counts_out is not None:
                outputs.append(
                    write_host_switch_count_table(
                        args.counts_out,
                        report,
                    )
                )
            if args.fits_out is not None:
                outputs.append(
                    write_host_switch_fit_table(
                        args.fits_out,
                        report,
                    )
                )
            if args.unsupported_out is not None:
                outputs.append(
                    write_unsupported_host_switch_claim_table(
                        args.unsupported_out,
                        report,
                    )
                )
            if args.exclusions_out is not None:
                outputs.append(
                    write_host_switch_exclusion_table(
                        args.exclusions_out,
                        report,
                    )
                )
            inputs = [args.tree, args.table]
            if args.constraints is not None:
                inputs.append(args.constraints)
            outputs = _finalize_outputs(
                args,
                command="host-association",
                inputs=inputs,
                outputs=outputs,
            )
            _print_result(
                build_command_result(
                    command="host-association",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=report.warnings,
                    metrics={
                        "model": report.model,
                        "analysis_constraint_mode": (
                            report.summary.analysis_constraint_mode
                        ),
                        "observed_host_count": report.summary.observed_host_count,
                        "host_switch_count": report.summary.host_switch_count,
                        "certain_host_switch_count": (
                            report.summary.certain_host_switch_count
                        ),
                        "uncertain_host_switch_count": (
                            report.summary.uncertain_host_switch_count
                        ),
                        "preferred_constraint": (report.summary.preferred_constraint),
                        "unsupported_switch_claim_count": (
                            report.summary.unsupported_switch_claim_count
                        ),
                        "excluded_taxon_count": (report.summary.excluded_taxon_count),
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if (
            args.command == "ecological-niche"
            and args.ecological_niche_command == "transitions"
        ):
            report = summarize_niche_transitions(
                args.tree,
                args.table,
                trait=args.trait,
                taxon_column=args.taxon_column,
                model=args.model,
            )
            outputs: list[Path | str] = []
            if args.summary_out is not None:
                outputs.append(
                    write_niche_transition_summary_table(
                        args.summary_out,
                        report,
                    )
                )
            if args.nodes_out is not None:
                outputs.append(
                    write_niche_state_node_table(
                        args.nodes_out,
                        report,
                    )
                )
            if args.rates_out is not None:
                outputs.append(
                    write_niche_transition_rate_table(
                        args.rates_out,
                        report,
                    )
                )
            if args.branches_out is not None:
                outputs.append(
                    write_niche_transition_branch_table(
                        args.branches_out,
                        report,
                    )
                )
            if args.counts_out is not None:
                outputs.append(
                    write_niche_transition_count_table(
                        args.counts_out,
                        report,
                    )
                )
            if args.clades_out is not None:
                outputs.append(
                    write_niche_transition_clade_table(
                        args.clades_out,
                        report,
                    )
                )
            if args.exclusions_out is not None:
                outputs.append(
                    write_niche_transition_exclusion_table(
                        args.exclusions_out,
                        report,
                    )
                )
            outputs = _finalize_outputs(
                args,
                command="ecological-niche",
                inputs=[args.tree, args.table],
                outputs=outputs,
            )
            _print_result(
                build_command_result(
                    command="ecological-niche",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                    warnings=report.warnings,
                    metrics={
                        "model": report.model,
                        "observed_niche_count": report.summary.observed_niche_count,
                        "transition_rate_row_count": (
                            report.summary.transition_rate_row_count
                        ),
                        "changed_branch_count": (report.summary.changed_branch_count),
                        "certain_transition_count": (
                            report.summary.certain_transition_count
                        ),
                        "uncertain_transition_count": (
                            report.summary.uncertain_transition_count
                        ),
                        "repeated_shift_clade_count": (
                            report.summary.repeated_shift_clade_count
                        ),
                        "excluded_taxon_count": (report.summary.excluded_taxon_count),
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "phylogeography":
            if args.phylogeography_command == "coordinates":
                report = summarize_continuous_phylogeography(
                    args.tree,
                    args.table,
                    latitude_column=args.latitude_column,
                    longitude_column=args.longitude_column,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    alpha=args.alpha,
                )
                outputs: list[Path | str] = []
                if args.summary_out is not None:
                    outputs.append(
                        write_coordinate_movement_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                if args.estimates_out is not None:
                    outputs.append(
                        write_coordinate_estimate_table(
                            args.estimates_out,
                            report,
                        )
                    )
                if args.branches_out is not None:
                    outputs.append(
                        write_coordinate_movement_branch_table(
                            args.branches_out,
                            report,
                        )
                    )
                if args.outliers_out is not None:
                    outputs.append(
                        write_coordinate_movement_outlier_table(
                            args.outliers_out,
                            report,
                        )
                    )
                if args.exclusions_out is not None:
                    outputs.append(
                        write_coordinate_movement_exclusion_table(
                            args.exclusions_out,
                            report,
                        )
                    )
                if args.visualization_out is not None:
                    outputs.append(
                        render_coordinate_movement_visualization(
                            report,
                            out_path=args.visualization_out,
                        ).output_path
                    )
                outputs = _finalize_outputs(
                    args,
                    command="phylogeography",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="phylogeography",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "model": report.model,
                            "analyzed_taxon_count": report.summary.analyzed_taxon_count,
                            "outlier_jump_count": report.summary.outlier_jump_count,
                            "impossible_jump_count": (
                                report.summary.impossible_jump_count
                            ),
                            "flagged_branch_count": (
                                report.summary.flagged_branch_count
                            ),
                            "maximum_jump_km": report.summary.maximum_jump_km,
                            "excluded_taxon_count": (
                                report.summary.excluded_taxon_count
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.phylogeography_command == "coordinates-map":
                report = summarize_continuous_phylogeography_map(
                    args.tree,
                    args.table,
                    latitude_column=args.latitude_column,
                    longitude_column=args.longitude_column,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    alpha=args.alpha,
                    minimum_midpoint_depth=args.minimum_midpoint_depth,
                    maximum_midpoint_depth=args.maximum_midpoint_depth,
                )
                outputs: list[Path | str] = []
                if args.summary_out is not None:
                    outputs.append(
                        write_geographic_map_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                if args.markers_out is not None:
                    outputs.append(
                        write_geographic_map_marker_table(
                            args.markers_out,
                            report,
                        )
                    )
                if args.lines_out is not None:
                    outputs.append(
                        write_geographic_map_line_table(
                            args.lines_out,
                            report,
                        )
                    )
                if args.exclusions_out is not None:
                    outputs.append(
                        write_geographic_map_exclusion_table(
                            args.exclusions_out,
                            report,
                        )
                    )
                if args.html_out is not None:
                    outputs.append(
                        render_geographic_map_html(
                            report,
                            out_path=args.html_out,
                        ).output_path
                    )
                outputs = _finalize_outputs(
                    args,
                    command="phylogeography",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="phylogeography",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "map_mode": report.summary.mode,
                            "model": report.summary.model,
                            "tip_marker_count": report.summary.tip_marker_count,
                            "internal_marker_count": (
                                report.summary.internal_marker_count
                                + report.summary.root_marker_count
                            ),
                            "line_count": report.summary.line_count,
                            "visible_line_count": report.summary.visible_line_count,
                            "time_filter_applied": (report.summary.time_filter_applied),
                            "excluded_record_count": (
                                report.summary.excluded_record_count
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.phylogeography_command == "regions-map":
                report = summarize_discrete_region_map(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    centroids_path=args.centroids,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    region_column=args.region_column,
                    latitude_column=args.latitude_column,
                    longitude_column=args.longitude_column,
                    minimum_midpoint_depth=args.minimum_midpoint_depth,
                    maximum_midpoint_depth=args.maximum_midpoint_depth,
                )
                outputs: list[Path | str] = []
                if args.summary_out is not None:
                    outputs.append(
                        write_geographic_map_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                if args.markers_out is not None:
                    outputs.append(
                        write_geographic_map_marker_table(
                            args.markers_out,
                            report,
                        )
                    )
                if args.lines_out is not None:
                    outputs.append(
                        write_geographic_map_line_table(
                            args.lines_out,
                            report,
                        )
                    )
                if args.exclusions_out is not None:
                    outputs.append(
                        write_geographic_map_exclusion_table(
                            args.exclusions_out,
                            report,
                        )
                    )
                if args.html_out is not None:
                    outputs.append(
                        render_geographic_map_html(
                            report,
                            out_path=args.html_out,
                        ).output_path
                    )
                inputs = [args.tree, args.table, args.centroids]
                outputs = _finalize_outputs(
                    args,
                    command="phylogeography",
                    inputs=inputs,
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="phylogeography",
                        inputs=inputs,
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "map_mode": report.summary.mode,
                            "model": report.summary.model,
                            "tip_marker_count": report.summary.tip_marker_count,
                            "internal_marker_count": (
                                report.summary.internal_marker_count
                                + report.summary.root_marker_count
                            ),
                            "line_count": report.summary.line_count,
                            "visible_line_count": report.summary.visible_line_count,
                            "time_filter_applied": (report.summary.time_filter_applied),
                            "excluded_record_count": (
                                report.summary.excluded_record_count
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
        if args.command == "discrete-evolution":
            allowed_states = (
                _split_csv_values(args.allowed_states)
                if hasattr(args, "allowed_states")
                else []
            )
            if args.discrete_evolution_command == "validate-coding":
                report = validate_discrete_state_coding(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    allowed_states=allowed_states or None,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                )
                outputs = _finalize_outputs(
                    args, command="discrete-evolution", inputs=[args.tree, args.table]
                )
                _print_result(
                    build_command_result(
                        command="discrete-evolution",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "valid": report.valid,
                            "issue_count": len(report.issues),
                            "observed_state_count": len(report.observed_states),
                            "state_ordering": report.state_ordering,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.discrete_evolution_command == "imbalance":
                report = detect_state_imbalance_problems(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                )
                outputs = _finalize_outputs(
                    args, command="discrete-evolution", inputs=[args.tree, args.table]
                )
                _print_result(
                    build_command_result(
                        command="discrete-evolution",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=[warning.message for warning in report.warnings],
                        metrics={
                            "taxon_count": report.taxon_count,
                            "observed_state_count": len(report.observed_states),
                            "warning_count": len(report.warnings),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.discrete_evolution_command == "reference":
                report = validate_discrete_transition_reference_examples()
                outputs = _finalize_outputs(
                    args, command="discrete-evolution", inputs=[]
                )
                _print_result(
                    build_command_result(
                        command="discrete-evolution",
                        inputs=[],
                        outputs=outputs,
                        metrics={
                            "case_count": report.case_count,
                            "all_passed": report.all_passed,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.discrete_evolution_command == "model":
                report = estimate_ancestral_geographic_states(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    allowed_states=allowed_states or None,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                )
                outputs: list[Path | str] = []
                if args.node_table_out is not None:
                    outputs.append(
                        write_node_state_probability_table(args.node_table_out, report)
                    )
                if args.transitions_out is not None:
                    outputs.append(
                        write_transition_summary_table(args.transitions_out, report)
                    )
                outputs = _finalize_outputs(
                    args,
                    command="discrete-evolution",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="discrete-evolution",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "observed_state_count": len(report.observed_states),
                            "transition_count": report.transition_summary.transition_count,
                            "strongly_supported_transition_count": report.transition_summary.strongly_supported_transition_count,
                            "model": report.model,
                            "state_ordering": report.state_ordering,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.discrete_evolution_command == "stochastic-map":
                report = simulate_discrete_stochastic_maps(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    allowed_states=allowed_states or None,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                    replicates=args.replicates,
                    seed=args.seed,
                )
                count_report = count_discrete_stochastic_map_transitions(report)
                density_report = None
                density_result = None
                if any(
                    output is not None
                    for output in (
                        args.branch_probabilities_out,
                        args.density_branches_out,
                        args.density_slices_out,
                        args.density_figure_out,
                    )
                ):
                    density_report = summarize_discrete_stochastic_map_density(
                        report,
                        resolution=args.density_resolution,
                        focal_state=args.focal_state,
                    )
                outputs: list[Path | str] = []
                if args.collection_out is not None:
                    outputs.append(
                        write_stochastic_map_collection(args.collection_out, report)
                    )
                if args.summary_out is not None:
                    outputs.append(
                        write_stochastic_map_summary_table(
                            args.summary_out, report.summary
                        )
                    )
                if args.state_times_out is not None:
                    outputs.append(
                        write_stochastic_map_state_time_table(
                            args.state_times_out,
                            report.summary,
                        )
                    )
                if args.branch_occupancy_out is not None:
                    outputs.append(
                        write_stochastic_map_branch_occupancy_table(
                            args.branch_occupancy_out,
                            report.summary,
                        )
                    )
                if args.count_matrix_out is not None:
                    outputs.append(
                        write_stochastic_map_transition_count_matrix(
                            args.count_matrix_out,
                            count_report,
                        )
                    )
                if args.aggregate_matrix_out is not None:
                    outputs.append(
                        write_stochastic_map_aggregate_transition_matrix(
                            args.aggregate_matrix_out,
                            count_report,
                        )
                    )
                if args.branch_transition_out is not None:
                    outputs.append(
                        write_stochastic_map_branch_transition_count_table(
                            args.branch_transition_out,
                            count_report,
                        )
                    )
                if args.segments_out is not None:
                    outputs.append(
                        write_stochastic_map_segment_table(
                            args.segments_out,
                            report,
                        )
                    )
                if args.events_out is not None:
                    outputs.append(
                        write_stochastic_map_event_table(
                            args.events_out,
                            report,
                        )
                    )
                if (
                    density_report is not None
                    and args.branch_probabilities_out is not None
                ):
                    outputs.append(
                        write_stochastic_map_branch_probability_table(
                            args.branch_probabilities_out,
                            density_report,
                        )
                    )
                if density_report is not None and args.density_branches_out is not None:
                    outputs.append(
                        write_stochastic_map_density_branch_table(
                            args.density_branches_out,
                            density_report,
                        )
                    )
                if density_report is not None and args.density_slices_out is not None:
                    outputs.append(
                        write_stochastic_map_density_slice_table(
                            args.density_slices_out,
                            density_report,
                        )
                    )
                if density_report is not None and args.density_figure_out is not None:
                    density_result = render_stochastic_map_density_artifact(
                        density_report,
                        tree_path=report.tree_path,
                        out_path=args.density_figure_out,
                        layout=args.layout,
                    )
                    outputs.extend(
                        [density_result.output_path, density_result.svg_path]
                        if density_result.format == "html"
                        else [density_result.output_path]
                    )
                outputs = _finalize_outputs(
                    args,
                    command="discrete-evolution",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="discrete-evolution",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=list(
                            dict.fromkeys(
                                [
                                    *report.warnings,
                                    *(
                                        density_report.warnings
                                        if density_report is not None
                                        else []
                                    ),
                                ]
                            )
                        ),
                        metrics={
                            "requested_replicate_count": report.replicates,
                            "successful_replicate_count": report.summary.replicate_count,
                            "simulation_failure_count": report.summary.simulation_failure_count,
                            "mean_total_transition_count": report.summary.mean_total_transition_count,
                            "branch_state_row_count": len(
                                report.summary.branch_occupancy_rows
                            ),
                            "count_matrix_row_count": len(count_report.matrix_rows),
                            "branch_transition_row_count": len(
                                count_report.branch_rows
                            ),
                            "model": report.model,
                            "state_ordering": report.state_ordering,
                            "conditioned_on_node_estimates": report.conditioned_on_node_estimates,
                            "parameter_count": report.fit_audit.parameter_count,
                            "log_likelihood": report.fit_audit.log_likelihood,
                            "aic": report.fit_audit.aic,
                            "overparameterized": report.fit_audit.overparameterized,
                            "optimizer_converged": report.fit_audit.optimizer_converged,
                            "optimizer_hit_lower_parameter_bound": report.fit_audit.optimizer_hit_lower_parameter_bound,
                            "optimizer_hit_upper_parameter_bound": report.fit_audit.optimizer_hit_upper_parameter_bound,
                            "fit_warning_count": len(report.fit_audit.warnings),
                            "baseline_model": report.fit_audit.baseline_model,
                            "preferred_model_by_aic": report.fit_audit.preferred_model_by_aic,
                            "branch_probability_row_count": (
                                len(density_report.branch_state_rows)
                                if density_report is not None
                                else 0
                            ),
                            "density_branch_row_count": (
                                len(density_report.branch_rows)
                                if density_report is not None
                                else 0
                            ),
                            "density_slice_row_count": (
                                len(density_report.density_rows)
                                if density_report is not None
                                else 0
                            ),
                            "density_focal_state": (
                                density_report.focal_state
                                if density_report is not None
                                else None
                            ),
                            "density_resolution": (
                                density_report.resolution
                                if density_report is not None
                                else None
                            ),
                            "density_rendered_branch_color_count": (
                                density_result.rendered_branch_color_count
                                if density_result is not None
                                else 0
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.discrete_evolution_command == "summarize-maps":
                collection = load_stochastic_map_collection(args.input_path)
                report = summarize_discrete_stochastic_maps(collection)
                outputs: list[Path | str] = []
                if args.summary_out is not None:
                    outputs.append(
                        write_stochastic_map_summary_table(args.summary_out, report)
                    )
                if args.state_times_out is not None:
                    outputs.append(
                        write_stochastic_map_state_time_table(
                            args.state_times_out,
                            report,
                        )
                    )
                if args.branch_occupancy_out is not None:
                    outputs.append(
                        write_stochastic_map_branch_occupancy_table(
                            args.branch_occupancy_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="discrete-evolution",
                    inputs=[args.input_path],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="discrete-evolution",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "replicate_count": report.replicate_count,
                            "mean_total_transition_count": report.mean_total_transition_count,
                            "simulation_failure_count": report.simulation_failure_count,
                            "branch_state_row_count": len(report.branch_occupancy_rows),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.discrete_evolution_command == "count-maps":
                collection = load_stochastic_map_collection(args.input_path)
                report = count_discrete_stochastic_map_transitions(collection)
                outputs: list[Path | str] = []
                if args.count_matrix_out is not None:
                    outputs.append(
                        write_stochastic_map_transition_count_matrix(
                            args.count_matrix_out,
                            report,
                        )
                    )
                if args.aggregate_matrix_out is not None:
                    outputs.append(
                        write_stochastic_map_aggregate_transition_matrix(
                            args.aggregate_matrix_out,
                            report,
                        )
                    )
                if args.branch_transition_out is not None:
                    outputs.append(
                        write_stochastic_map_branch_transition_count_table(
                            args.branch_transition_out,
                            report,
                        )
                    )
                if args.events_out is not None:
                    outputs.append(
                        write_stochastic_map_event_table(
                            args.events_out,
                            collection,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="discrete-evolution",
                    inputs=[args.input_path],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="discrete-evolution",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "replicate_count": report.replicate_count,
                            "mean_total_transition_count": report.mean_total_transition_count,
                            "count_matrix_row_count": len(report.matrix_rows),
                            "branch_transition_row_count": len(report.branch_rows),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.discrete_evolution_command == "density-maps":
                collection = load_stochastic_map_collection(args.input_path)
                report = summarize_discrete_stochastic_map_density(
                    collection,
                    resolution=args.resolution,
                    focal_state=args.focal_state,
                )
                outputs: list[Path | str] = []
                density_result = None
                if args.branch_probabilities_out is not None:
                    outputs.append(
                        write_stochastic_map_branch_probability_table(
                            args.branch_probabilities_out,
                            report,
                        )
                    )
                if args.density_branches_out is not None:
                    outputs.append(
                        write_stochastic_map_density_branch_table(
                            args.density_branches_out,
                            report,
                        )
                    )
                if args.density_slices_out is not None:
                    outputs.append(
                        write_stochastic_map_density_slice_table(
                            args.density_slices_out,
                            report,
                        )
                    )
                if args.out is not None:
                    density_result = render_stochastic_map_density_artifact(
                        report,
                        tree_path=collection.tree_path,
                        out_path=args.out,
                        layout=args.layout,
                    )
                    outputs.extend(
                        [density_result.output_path, density_result.svg_path]
                        if density_result.format == "html"
                        else [density_result.output_path]
                    )
                outputs = _finalize_outputs(
                    args,
                    command="discrete-evolution",
                    inputs=[args.input_path],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="discrete-evolution",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "replicate_count": report.replicate_count,
                            "resolution": report.resolution,
                            "branch_probability_row_count": len(
                                report.branch_state_rows
                            ),
                            "density_branch_row_count": len(report.branch_rows),
                            "density_slice_row_count": len(report.density_rows),
                            "focal_state": report.focal_state,
                            "baseline_state": report.baseline_state,
                            "rendered_branch_color_count": (
                                density_result.rendered_branch_color_count
                                if density_result is not None
                                else 0
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.discrete_evolution_command == "render":
                report = estimate_ancestral_geographic_states(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    allowed_states=allowed_states or None,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                )
                result = render_tree_with_geographic_states(
                    args.tree,
                    report,
                    out_path=args.out,
                    layout=args.layout,
                )
                outputs = _finalize_outputs(
                    args,
                    command="discrete-evolution",
                    inputs=[args.tree, args.table],
                    outputs=[result.output_path],
                )
                _print_result(
                    build_command_result(
                        command="discrete-evolution",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "tip_count": result.tip_count,
                            "rendered_internal_annotation_count": result.rendered_internal_annotation_count,
                            "layout": result.layout,
                            "model": report.model,
                            "state_ordering": report.state_ordering,
                        },
                        data={
                            "reconstruction": report,
                            "render": result,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            if args.discrete_evolution_command == "report":
                result = render_discrete_state_evolution_report(
                    tree_path=args.tree,
                    traits_path=args.table,
                    trait=args.trait,
                    out_path=args.out,
                    taxon_column=args.taxon_column,
                    model=args.model,
                    allowed_states=allowed_states or None,
                    state_ordering=args.state_ordering,
                    ordered_states=_split_csv_values(args.ordered_states) or None,
                    compare_model=args.compare_model,
                )
                outputs = _finalize_outputs(
                    args,
                    command="discrete-evolution",
                    inputs=[args.tree, args.table],
                    outputs=[result.output_path, args.out.with_suffix(".svg")],
                )
                _print_result(
                    build_command_result(
                        command="discrete-evolution",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "report_kind": result.report_kind,
                            "model": result.model,
                            "state_ordering": args.state_ordering,
                        },
                        data=result,
                    ),
                    json_output=args.json,
                )
                return 0
            comparison = compare_discrete_state_models(
                args.tree,
                args.table,
                trait=args.trait,
                taxon_column=args.taxon_column,
                left_model=args.left_model,
                right_model=args.right_model,
                allowed_states=allowed_states or None,
                state_ordering=args.state_ordering,
                ordered_states=_split_csv_values(args.ordered_states) or None,
            )
            outputs: list[Path | str] = []
            if args.table_out is not None:
                outputs.append(
                    write_discrete_model_comparison_table(args.table_out, comparison)
                )
            outputs = _finalize_outputs(
                args,
                command="discrete-evolution",
                inputs=[args.tree, args.table],
                outputs=outputs,
            )
            _print_result(
                build_command_result(
                    command="discrete-evolution",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                    metrics={
                        "better_model": comparison.better_model,
                        "model_count": len(comparison.rows),
                        "differing_node_count": sum(
                            1 for row in comparison.node_differences if row.differs
                        ),
                        "state_ordering": args.state_ordering,
                    },
                    data=comparison,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "diversification":
            if args.diversification_command == "ltt":
                report = compute_lineage_through_time_curve(args.tree)
                outputs: list[Path | str] = []
                if args.out is not None:
                    outputs.append(write_lineage_through_time_table(args.out, report))
                outputs = _finalize_outputs(
                    args,
                    command="diversification",
                    inputs=[args.tree],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="diversification",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={
                            "tip_count": report.tip_count,
                            "root_age": report.root_age,
                            "point_count": len(report.points),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.diversification_command == "sampling":
                report = detect_incomplete_taxon_sampling_metadata(
                    args.tree,
                    args.table,
                    taxon_column=args.taxon_column,
                    sampling_column=args.sampling_column,
                )
                outputs = _finalize_outputs(
                    args, command="diversification", inputs=[args.tree, args.table]
                )
                _print_result(
                    build_command_result(
                        command="diversification",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "complete": report.complete,
                            "matched_taxon_count": len(report.matched_taxa),
                            "missing_taxon_count": len(report.missing_taxa),
                            "invalid_row_count": len(report.invalid_rows),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.diversification_command == "estimate":
                inputs = [args.tree]
                if args.metadata is not None:
                    inputs.append(args.metadata)
                report = estimate_diversification_rate(
                    args.tree,
                    metadata_path=args.metadata,
                    taxon_column=args.taxon_column,
                    sampling_column=args.sampling_column,
                    model=args.model,
                )
                outputs = _finalize_outputs(
                    args, command="diversification", inputs=inputs
                )
                _print_result(
                    build_command_result(
                        command="diversification",
                        inputs=inputs,
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "model": report.model,
                            "sampling_fraction": report.sampling_fraction,
                            "net_diversification_rate": report.net_diversification_rate,
                            "aic": report.aic,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.diversification_command == "gamma-stat":
                inputs = [args.tree]
                if args.metadata is not None:
                    inputs.append(args.metadata)
                report = compute_diversification_gamma_statistic(
                    args.tree,
                    metadata_path=args.metadata,
                    taxon_column=args.taxon_column,
                    sampling_column=args.sampling_column,
                )
                outputs: list[Path | str] = []
                if args.out is not None:
                    outputs.append(
                        write_diversification_gamma_statistic_table(args.out, report)
                    )
                outputs = _finalize_outputs(
                    args, command="diversification", inputs=inputs, outputs=outputs
                )
                _print_result(
                    build_command_result(
                        command="diversification",
                        inputs=inputs,
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "tip_count": report.tip_count,
                            "branching_time_count": report.branching_time_count,
                            "gamma_statistic": report.gamma_statistic,
                            "sampling_fraction": report.sampling_fraction,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.diversification_command == "compare-models":
                inputs = [args.tree]
                if args.metadata is not None:
                    inputs.append(args.metadata)
                report = compare_diversification_models(
                    args.tree,
                    metadata_path=args.metadata,
                    taxon_column=args.taxon_column,
                    sampling_column=args.sampling_column,
                )
                outputs = _finalize_outputs(
                    args, command="diversification", inputs=inputs
                )
                _print_result(
                    build_command_result(
                        command="diversification",
                        inputs=inputs,
                        outputs=outputs,
                        metrics={
                            "better_model": report.better_model,
                            "model_count": len(report.rows),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.diversification_command == "clades":
                report = detect_diversification_outlier_clades(
                    args.tree,
                    min_tip_count=args.min_tip_count,
                    model=args.model,
                )
                outputs: list[Path | str] = []
                if args.out is not None:
                    outputs.append(write_clade_diversification_table(args.out, report))
                outputs = _finalize_outputs(
                    args,
                    command="diversification",
                    inputs=[args.tree],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="diversification",
                        inputs=[args.tree],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "global_rate": report.global_rate,
                            "high_clade_count": len(report.high_diversification_clades),
                            "low_clade_count": len(report.low_diversification_clades),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.diversification_command == "trait-dependent":
                report = run_trait_dependent_diversification_analysis(
                    args.tree,
                    args.table,
                    trait=args.trait,
                    taxon_column=args.taxon_column,
                )
                outputs: list[Path | str] = []
                if args.out is not None:
                    outputs.append(
                        write_trait_dependent_diversification_table(args.out, report)
                    )
                outputs = _finalize_outputs(
                    args,
                    command="diversification",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="diversification",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "state_count": len(report.states),
                            "monophyletic_state_count": sum(
                                1 for row in report.states if row.monophyletic
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.diversification_command == "package":
                inputs = [args.tree]
                if args.metadata is not None:
                    inputs.append(args.metadata)
                result = build_diversification_figure_package(
                    args.tree,
                    out_dir=args.out_dir,
                    metadata_path=args.metadata,
                    taxon_column=args.taxon_column,
                    sampling_column=args.sampling_column,
                    min_tip_count=args.min_tip_count,
                    model=args.model,
                )
                outputs = _finalize_outputs(
                    args,
                    command="diversification",
                    inputs=inputs,
                    outputs=[
                        result.lineage_figure_path,
                        result.clade_figure_path,
                        result.model_figure_path,
                        result.lineage_table_path,
                        result.clade_table_path,
                        result.model_table_path,
                        result.legend_path,
                        result.caption_path,
                        result.methods_summary_path,
                        result.review_path,
                        result.manifest_path,
                        result.reproducibility_manifest_path,
                    ],
                )
                warnings = (
                    []
                    if result.sampling_report is None
                    else list(result.sampling_report.warnings)
                )
                _print_result(
                    build_command_result(
                        command="diversification",
                        inputs=inputs,
                        outputs=outputs,
                        warnings=warnings,
                        metrics={
                            "publication_ready": result.audit.publication_ready,
                            "sampling_metadata_complete": (
                                result.audit.sampling_metadata_complete
                            ),
                            "plotted_ltt_point_count": (
                                result.audit.plotted_ltt_point_count
                            ),
                            "plotted_clade_count": result.audit.plotted_clade_count,
                            "highlighted_outlier_count": (
                                result.audit.highlighted_outlier_count
                            ),
                            "plotted_model_count": result.audit.plotted_model_count,
                            "better_model": result.audit.better_model,
                            "methods_summary_warning_count": (
                                result.methods_summary.warning_count
                            ),
                        },
                        data=result.machine_manifest,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.diversification_command == "methods-summary":
                inputs = [args.tree]
                if args.metadata is not None:
                    inputs.append(args.metadata)
                if args.traits is not None:
                    inputs.append(args.traits)
                report = build_diversification_method_report(
                    args.tree,
                    metadata_path=args.metadata,
                    taxon_column=args.taxon_column,
                    sampling_column=args.sampling_column,
                    traits_path=args.traits,
                    trait=args.trait,
                    estimate_model=args.estimate_model,
                    clade_model=args.clade_model,
                    clade_min_tip_count=args.min_tip_count,
                )
                result = write_diversification_methods_summary_text(args.out, report)
                outputs = _finalize_outputs(
                    args,
                    command="diversification",
                    inputs=inputs,
                    outputs=[result.output_path],
                )
                warnings = (
                    []
                    if result.report.sampling_report is None
                    else list(result.report.sampling_report.warnings)
                )
                _print_result(
                    build_command_result(
                        command="diversification",
                        inputs=inputs,
                        outputs=outputs,
                        warnings=warnings,
                        metrics={
                            "warning_count": result.warning_count,
                            "better_model": result.better_model,
                            "sampling_metadata_complete": (
                                result.sampling_metadata_complete
                            ),
                            "clade_observation_count": (
                                result.clade_observation_count
                            ),
                            "trait_state_count": (
                                0
                                if result.report.trait_report is None
                                else len(result.report.trait_report.states)
                            ),
                        },
                        data=result,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.diversification_command == "medusa":
                report = summarize_medusa_exclusion(
                    args.tree,
                    metadata_path=args.metadata,
                    taxon_column=args.taxon_column,
                    sampling_column=args.sampling_column,
                )
                raise DiversificationAnalysisError(
                    report.exclusion_reason,
                    code="diversification_medusa_explicitly_excluded",
                    details={
                        "failure_reason": report.exclusion_code,
                        "supported_surfaces": report.supported_surfaces,
                        "missing_surfaces": report.missing_surfaces,
                        "tip_count": report.validation.tip_count,
                        "rooted": report.validation.rooted,
                        "ultrametric": report.validation.ultrametric,
                        "sampling_metadata_complete": (
                            None
                            if report.sampling_report is None
                            else report.sampling_report.complete
                        ),
                    },
                )
            if args.diversification_command == "bd-ms":
                report = summarize_geiger_birth_death_exclusion(
                    args.tree,
                    metadata_path=args.metadata,
                    taxon_column=args.taxon_column,
                    sampling_column=args.sampling_column,
                )
                raise DiversificationAnalysisError(
                    report.exclusion_reason,
                    code="diversification_birth_death_explicitly_excluded",
                    details={
                        "failure_reason": report.exclusion_code,
                        "geiger_reference_surface": report.geiger_reference_surface,
                        "geiger_reference_arguments": (
                            report.geiger_reference_arguments
                        ),
                        "owned_surface": report.owned_surface,
                        "tip_count": report.validation.tip_count,
                        "rooted": report.validation.rooted,
                        "ultrametric": report.validation.ultrametric,
                        "sampling_metadata_complete": (
                            None
                            if report.sampling_report is None
                            else report.sampling_report.complete
                        ),
                    },
                )
            inputs = [args.tree]
            if args.metadata is not None:
                inputs.append(args.metadata)
            if args.traits is not None:
                inputs.append(args.traits)
            result = render_diversification_report(
                tree_path=args.tree,
                out_path=args.out,
                metadata_path=args.metadata,
                taxon_column=args.taxon_column,
                sampling_column=args.sampling_column,
                traits_path=args.traits,
                trait=args.trait,
                methods_summary_path=args.methods_summary_out,
            )
            output_paths: list[Path | str] = [result.output_path]
            if result.methods_summary_path is not None:
                output_paths.append(result.methods_summary_path)
            outputs = _finalize_outputs(
                args,
                command="diversification",
                inputs=inputs,
                outputs=output_paths,
            )
            _print_result(
                build_command_result(
                    command="diversification",
                    inputs=inputs,
                    outputs=outputs,
                    metrics={
                        "report_kind": result.report_kind,
                        "methods_summary_warning_count": (
                            result.methods_summary_warning_count
                        ),
                        "better_model": result.report.model_comparison.better_model,
                    },
                    data=result,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "distance":
            return run_distance_command(args)
        if args.command == "tree-set":
            return run_tree_set_command(args)
        if args.command == "simulate":
            if args.simulate_command == "tree-birth-death":
                trees, report = simulate_birth_death_trees(
                    tree_count=args.tree_count,
                    tip_count=args.tip_count,
                    birth_rate=args.birth_rate,
                    death_rate=args.death_rate,
                    seed=args.seed,
                )
                output_path = write_tree_set(args.out, trees)
                outputs_to_finalize = [output_path]
                if args.record_table_out is not None:
                    outputs_to_finalize.append(
                        write_tree_simulation_record_table(
                            args.record_table_out, report
                        )
                    )
                if args.envelope_table_out is not None:
                    outputs_to_finalize.append(
                        write_tree_simulation_envelope_table(
                            args.envelope_table_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="simulate",
                    inputs=[],
                    outputs=outputs_to_finalize,
                )
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[],
                        outputs=outputs,
                        metrics={
                            "tree_count": report.tree_count,
                            "tip_count": report.tip_count,
                            "pooled_branch_count": report.pooled_branch_count,
                            "envelope_metric_count": len(report.envelope_metrics),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "tree-random":
                trees, report = simulate_random_trees(
                    tree_count=args.tree_count,
                    tip_count=args.tip_count,
                    seed=args.seed,
                )
                output_path = write_tree_set(args.out, trees)
                outputs_to_finalize = [output_path]
                if args.record_table_out is not None:
                    outputs_to_finalize.append(
                        write_tree_simulation_record_table(
                            args.record_table_out, report
                        )
                    )
                if args.envelope_table_out is not None:
                    outputs_to_finalize.append(
                        write_tree_simulation_envelope_table(
                            args.envelope_table_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="simulate",
                    inputs=[],
                    outputs=outputs_to_finalize,
                )
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[],
                        outputs=outputs,
                        metrics={
                            "tree_count": report.tree_count,
                            "tip_count": report.tip_count,
                            "pooled_branch_count": report.pooled_branch_count,
                            "envelope_metric_count": len(report.envelope_metrics),
                            "branch_length_model": report.branch_length_model,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "tree-coalescent":
                trees, report = simulate_coalescent_trees(
                    tree_count=args.tree_count,
                    tip_count=args.tip_count,
                    population_size=args.population_size,
                    seed=args.seed,
                )
                output_path = write_tree_set(args.out, trees)
                outputs_to_finalize = [output_path]
                if args.record_table_out is not None:
                    outputs_to_finalize.append(
                        write_tree_simulation_record_table(
                            args.record_table_out, report
                        )
                    )
                if args.envelope_table_out is not None:
                    outputs_to_finalize.append(
                        write_tree_simulation_envelope_table(
                            args.envelope_table_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="simulate",
                    inputs=[],
                    outputs=outputs_to_finalize,
                )
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[],
                        outputs=outputs,
                        metrics={
                            "tree_count": report.tree_count,
                            "tip_count": report.tip_count,
                            "pooled_branch_count": report.pooled_branch_count,
                            "envelope_metric_count": len(report.envelope_metrics),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "traits-brownian":
                report = simulate_brownian_traits(
                    args.tree,
                    root_state=args.root_state,
                    sigma=args.sigma,
                    sigma_squared=args.sigma_squared,
                    seed=args.seed,
                )
                output_path = write_continuous_trait_table(args.out, report)
                outputs = _finalize_outputs(
                    args, command="simulate", inputs=[args.tree], outputs=[output_path]
                )
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={
                            "tip_count": report.tip_count,
                            "trait_count": len(report.traits),
                            "sigma_squared": report.sigma_squared,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "traits-speciational":
                report = simulate_speciational_traits(
                    args.tree,
                    root_state=args.root_state,
                    sigma=args.sigma,
                    sigma_squared=args.sigma_squared,
                    seed=args.seed,
                )
                output_path = write_continuous_trait_table(args.out, report)
                outputs = _finalize_outputs(
                    args, command="simulate", inputs=[args.tree], outputs=[output_path]
                )
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={
                            "tip_count": report.tip_count,
                            "trait_count": len(report.traits),
                            "sigma_squared": report.sigma_squared,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "traits-brownian-correlated":
                if args.covariance_rows is None and args.correlation_rows is None:
                    parser.error(
                        "correlated Brownian simulation requires either --covariance-row or --correlation-row"
                    )
                if args.covariance_rows is not None and args.trait_standard_deviation:
                    parser.error(
                        "trait standard deviations can only be used with --correlation-row"
                    )
                if (
                    args.correlation_rows is not None
                    and not args.trait_standard_deviation
                ):
                    parser.error(
                        "correlated Brownian simulation requires --trait-standard-deviation with --correlation-row"
                    )
                report = simulate_correlated_brownian_trait_collection(
                    args.tree,
                    trait_names=args.trait,
                    evolutionary_covariance_matrix=(
                        None
                        if args.covariance_rows is None
                        else [
                            _parse_float_csv_row(raw_row)
                            for raw_row in args.covariance_rows
                        ]
                    ),
                    evolutionary_correlation_matrix=(
                        None
                        if args.correlation_rows is None
                        else [
                            _parse_float_csv_row(raw_row)
                            for raw_row in args.correlation_rows
                        ]
                    ),
                    trait_standard_deviations=(
                        None
                        if not args.trait_standard_deviation
                        else args.trait_standard_deviation
                    ),
                    root_states=None if not args.root_state else args.root_state,
                    replicates=args.replicates,
                    seed=args.seed,
                )
                outputs: list[Path | str] = [
                    write_correlated_continuous_trait_collection_table(args.out, report)
                ]
                if args.summary_out is not None:
                    outputs.append(
                        write_correlated_continuous_trait_collection_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="simulate",
                    inputs=[args.tree],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={
                            "tip_count": report.tip_count,
                            "trait_count": len(report.trait_names),
                            "replicate_count": report.replicate_count,
                            "summary_row_count": len(report.rows),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "traits-ou":
                report = simulate_ou_traits(
                    args.tree,
                    root_state=args.root_state,
                    sigma=args.sigma,
                    alpha=args.alpha,
                    theta=args.theta,
                    seed=args.seed,
                )
                output_path = write_continuous_trait_table(args.out, report)
                outputs = _finalize_outputs(
                    args, command="simulate", inputs=[args.tree], outputs=[output_path]
                )
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={
                            "tip_count": report.tip_count,
                            "trait_count": len(report.traits),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "traits-early-burst":
                report = simulate_early_burst_traits(
                    args.tree,
                    root_state=args.root_state,
                    sigma=args.sigma,
                    rate_change=args.rate_change,
                    seed=args.seed,
                )
                output_path = write_continuous_trait_table(args.out, report)
                outputs = _finalize_outputs(
                    args, command="simulate", inputs=[args.tree], outputs=[output_path]
                )
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={
                            "tip_count": report.tip_count,
                            "trait_count": len(report.traits),
                            "rate_change": report.rate_change,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "traits-discrete":
                report = simulate_discrete_traits(
                    args.tree,
                    states=args.states,
                    transition_rate=args.transition_rate,
                    root_state=args.root_state,
                    seed=args.seed,
                )
                output_path = write_discrete_trait_table(args.out, report)
                outputs = _finalize_outputs(
                    args, command="simulate", inputs=[args.tree], outputs=[output_path]
                )
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={
                            "tip_count": report.tip_count,
                            "trait_count": len(report.traits),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "history-discrete":
                root_probability_rows = _parse_probability_assignments(
                    args.root_probability
                )
                report = simulate_discrete_histories(
                    args.tree,
                    states=args.states,
                    rate_rows=_parse_rate_rows(args.rate),
                    root_state=args.root_state,
                    root_state_probabilities=(root_probability_rows or None),
                    replicates=args.replicates,
                    seed=args.seed,
                )
                outputs_to_finalize = [
                    write_discrete_history_tip_truth_table(args.out, report)
                ]
                if args.nodes_out is not None:
                    outputs_to_finalize.append(
                        write_discrete_history_node_truth_table(args.nodes_out, report)
                    )
                if args.branches_out is not None:
                    outputs_to_finalize.append(
                        write_discrete_history_branch_truth_table(
                            args.branches_out,
                            report,
                        )
                    )
                if args.events_out is not None:
                    outputs_to_finalize.append(
                        write_discrete_history_event_table(args.events_out, report)
                    )
                if args.segments_out is not None:
                    outputs_to_finalize.append(
                        write_discrete_history_segment_table(args.segments_out, report)
                    )
                if args.summary_out is not None:
                    outputs_to_finalize.append(
                        write_discrete_history_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="simulate",
                    inputs=[args.tree],
                    outputs=outputs_to_finalize,
                )
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={
                            "tip_count": report.tip_count,
                            "branch_count": report.branch_count,
                            "replicate_count": report.replicate_count,
                            "state_count": len(report.states),
                            "mean_total_transition_count": (
                                report.mean_total_transition_count
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "alignment-dna":
                report = simulate_dna_alignment(
                    args.tree,
                    sequence_length=args.sequence_length,
                    substitution_rate=args.substitution_rate,
                    seed=args.seed,
                )
                output_path = write_simulated_alignment(args.out, report)
                outputs = _finalize_outputs(
                    args, command="simulate", inputs=[args.tree], outputs=[output_path]
                )
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={
                            "tip_count": report.tip_count,
                            "sequence_length": report.sequence_length,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "validate-sim-char-reference":
                report = validate_geiger_sim_char_reference_examples()
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[],
                        outputs=[],
                        metrics={
                            "case_count": report.case_count,
                            "all_passed": report.all_passed,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.simulate_command == "alignment-protein":
                report = simulate_protein_alignment(
                    args.tree,
                    sequence_length=args.sequence_length,
                    substitution_rate=args.substitution_rate,
                    seed=args.seed,
                )
                output_path = write_simulated_alignment(args.out, report)
                outputs = _finalize_outputs(
                    args, command="simulate", inputs=[args.tree], outputs=[output_path]
                )
                _print_result(
                    build_command_result(
                        command="simulate",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={
                            "tip_count": report.tip_count,
                            "sequence_length": report.sequence_length,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
        if args.command == "benchmark":
            if args.benchmark_command == "tree-validation":
                report = benchmark_tree_validation(replicates=args.replicates)
            elif args.benchmark_command == "tree-comparison":
                report = benchmark_tree_comparison(replicates=args.replicates)
            elif args.benchmark_command == "large-tree-scaling":
                report = benchmark_large_tree_scaling(
                    replicates=args.replicates,
                    tip_counts=args.tip_counts,
                )
            elif args.benchmark_command == "large-alignment-scaling":
                classes = None
                if args.sequence_counts is not None or args.alignment_lengths is not None:
                    sequence_counts = args.sequence_counts or []
                    alignment_lengths = args.alignment_lengths or []
                    if len(sequence_counts) != len(alignment_lengths):
                        raise ValueError(
                            "large-alignment-scaling requires the same number of --sequence-count and --alignment-length values"
                        )
                    classes = [
                        (
                            f"sequences-{sequence_count}-sites-{alignment_length}",
                            sequence_count,
                            alignment_length,
                        )
                        for sequence_count, alignment_length in zip(
                            sequence_counts,
                            alignment_lengths,
                            strict=True,
                        )
                    ]
                report = benchmark_large_alignment_scaling(
                    replicates=args.replicates,
                    size_classes=classes,
                )
            elif args.benchmark_command == "large-tree-set-scaling":
                classes = None
                if args.tree_counts is not None or args.tip_counts is not None:
                    tree_counts = args.tree_counts or []
                    tip_counts = args.tip_counts or []
                    if len(tree_counts) != len(tip_counts):
                        raise ValueError(
                            "large-tree-set-scaling requires the same number of --tree-count and --tip-count values"
                        )
                    classes = [
                        (f"trees-{tree_count}-taxa-{tip_count}", tree_count, tip_count)
                        for tree_count, tip_count in zip(
                            tree_counts,
                            tip_counts,
                            strict=True,
                        )
                    ]
                report = benchmark_large_tree_set_scaling(
                    replicates=args.replicates,
                    size_classes=classes,
                )
            elif args.benchmark_command == "workflow-practical-limits":
                alignment_classes = None
                if (
                    args.sequence_counts is not None
                    or args.alignment_lengths is not None
                ):
                    sequence_counts = args.sequence_counts or []
                    alignment_lengths = args.alignment_lengths or []
                    if len(sequence_counts) != len(alignment_lengths):
                        raise ValueError(
                            "workflow-practical-limits requires the same number of --sequence-count and --alignment-length values"
                        )
                    alignment_classes = [
                        (
                            f"sequences-{sequence_count}-sites-{alignment_length}",
                            sequence_count,
                            alignment_length,
                        )
                        for sequence_count, alignment_length in zip(
                            sequence_counts,
                            alignment_lengths,
                            strict=True,
                        )
                    ]
                tree_set_classes = None
                if (
                    args.posterior_tree_counts is not None
                    or args.tree_set_tip_counts is not None
                ):
                    posterior_tree_counts = args.posterior_tree_counts or []
                    tree_set_tip_counts = args.tree_set_tip_counts or []
                    if len(posterior_tree_counts) != len(tree_set_tip_counts):
                        raise ValueError(
                            "workflow-practical-limits requires the same number of --posterior-tree-count and --tree-set-tip-count values"
                        )
                    tree_set_classes = [
                        (
                            f"trees-{tree_count}-taxa-{tip_count}",
                            tree_count,
                            tip_count,
                        )
                        for tree_count, tip_count in zip(
                            posterior_tree_counts,
                            tree_set_tip_counts,
                            strict=True,
                        )
                    ]
                report = benchmark_workflow_practical_limits(
                    replicates=args.replicates,
                    tree_tip_counts=args.tree_tip_counts,
                    alignment_size_classes=alignment_classes,
                    tree_set_size_classes=tree_set_classes,
                    stress_tiers=args.stress_tiers,
                )
            elif args.benchmark_command == "stress-suite":
                report = benchmark_large_dataset_stress_suite(tier=args.tier)
            elif args.benchmark_command == "large-tree-model-fitting":
                report = benchmark_large_tree_model_fitting(tier=args.tier)
            elif args.benchmark_command == "real-dataset-macroevolution":
                report = benchmark_real_dataset_macroevolution()
            else:
                report = benchmark_alignment_diagnostics(
                    replicates=args.replicates,
                    sequence_length=args.sequence_length,
                )
            outputs = _finalize_outputs(args, command="benchmark", inputs=[])
            if hasattr(report, "entries"):
                metrics = {
                    "entry_count": len(report.entries),
                }
            elif hasattr(report, "summary_rows"):
                metrics = {
                    "summary_row_count": len(report.summary_rows),
                    "model_row_count": len(report.model_rows),
                    "alignment_review_row_count": len(report.alignment_review_rows),
                    "parity_row_count": len(report.parity_rows),
                }
            else:
                metrics = {
                    "observation_count": (
                        len(report.observations)
                        if hasattr(report, "observations")
                        else sum(len(row.observations) for row in report.workflows)
                    ),
                }
            if hasattr(report, "replicates"):
                metrics["replicates"] = report.replicates
            if hasattr(report, "tier"):
                metrics["tier"] = report.tier
            if hasattr(report, "case_count"):
                metrics["case_count"] = report.case_count
            if hasattr(report, "geiger_match_case_count"):
                metrics["geiger_match_case_count"] = report.geiger_match_case_count
            if hasattr(report, "threshold_pass_case_count"):
                metrics["threshold_pass_case_count"] = (
                    report.threshold_pass_case_count
                )
            if hasattr(report, "too_slow_case_count"):
                metrics["too_slow_case_count"] = report.too_slow_case_count
            if hasattr(report, "unstable_case_count"):
                metrics["unstable_case_count"] = report.unstable_case_count
            if hasattr(report, "stress_tiers"):
                metrics["stress_tier_count"] = len(report.stress_tiers)
            if hasattr(report, "workflows"):
                metrics["workflow_count"] = len(report.workflows)
                if hasattr(report, "tip_counts"):
                    metrics["max_tip_count"] = max(report.tip_counts)
                if hasattr(report, "tree_counts"):
                    metrics["max_tree_count"] = max(report.tree_counts)
                if hasattr(report, "alignment_lengths"):
                    metrics["max_alignment_length"] = max(report.alignment_lengths)
                    metrics["max_sequence_count"] = max(report.sequence_counts)
            if hasattr(report, "entries"):
                metrics["workflow_count"] = len(report.entries)
                taxon_limits = [
                    row.tested_taxon_limit
                    for row in report.entries
                    if row.tested_taxon_limit is not None
                ]
                site_limits = [
                    row.tested_site_limit
                    for row in report.entries
                    if row.tested_site_limit is not None
                ]
                tree_limits = [
                    row.tested_tree_limit
                    for row in report.entries
                    if row.tested_tree_limit is not None
                ]
                posterior_limits = [
                    row.tested_posterior_size
                    for row in report.entries
                    if row.tested_posterior_size is not None
                ]
                if taxon_limits:
                    metrics["max_taxon_limit"] = max(taxon_limits)
                if site_limits:
                    metrics["max_site_limit"] = max(site_limits)
                if tree_limits:
                    metrics["max_tree_limit"] = max(tree_limits)
                if posterior_limits:
                    metrics["max_posterior_size"] = max(posterior_limits)
            if hasattr(report, "observations") and report.observations:
                taxon_counts = [
                    row.taxon_count
                    for row in report.observations
                    if hasattr(row, "taxon_count") and row.taxon_count is not None
                ]
                if taxon_counts:
                    metrics["max_taxon_count"] = max(taxon_counts)
            _print_result(
                build_command_result(
                    command="benchmark",
                    inputs=[],
                    outputs=outputs,
                    metrics=metrics,
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "parity":
            if (
                args.reference_source in {"ape-live", "phytools-live", "geiger-live"}
                and args.extended
            ):
                raise ValueError(
                    "--extended is only supported for the checked fixture parity suite"
                )
            if args.reference_source == "ape-live":
                report = run_ape_parity_cases(
                    case_ids=args.ape_cases,
                    rscript_executable=args.ape_rscript_executable,
                    failure_root=args.ape_failure_root,
                )
                output_paths: list[Path | str] = []
                summary_path = None
                observation_path = None
                if args.summary_out is not None:
                    summary_path = write_ape_parity_summary_table(
                        args.summary_out,
                        report,
                    )
                    output_paths.append(summary_path)
                if args.observations_out is not None:
                    observation_path = write_ape_parity_observation_table(
                        args.observations_out,
                        report,
                    )
                    output_paths.append(observation_path)
                outputs = _finalize_outputs(
                    args,
                    command="parity",
                    inputs=[],
                    outputs=output_paths,
                )
                _print_result(
                    build_command_result(
                        command="parity",
                        inputs=[],
                        outputs=outputs,
                        metrics={
                            "all_passed": report.all_passed,
                            "case_count": report.case_count,
                            "function_count": len(report.summary_rows),
                            "failed_case_count": report.failed_case_count,
                            "skipped_case_count": report.skipped_case_count,
                            "reference_source": args.reference_source,
                        },
                        data={
                            "report": report,
                            "summary_table": summary_path,
                            "observation_table": observation_path,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            if args.reference_source == "phytools-live":
                report = run_phytools_parity_cases(
                    case_ids=args.phytools_cases,
                    rscript_executable=args.phytools_rscript_executable,
                    failure_root=args.phytools_failure_root,
                )
                output_paths: list[Path | str] = []
                summary_path = None
                observation_path = None
                if args.summary_out is not None:
                    summary_path = write_phytools_parity_summary_table(
                        args.summary_out,
                        report,
                    )
                    output_paths.append(summary_path)
                if args.observations_out is not None:
                    observation_path = write_phytools_parity_observation_table(
                        args.observations_out,
                        report,
                    )
                    output_paths.append(observation_path)
                outputs = _finalize_outputs(
                    args,
                    command="parity",
                    inputs=[],
                    outputs=output_paths,
                )
                _print_result(
                    build_command_result(
                        command="parity",
                        inputs=[],
                        outputs=outputs,
                        metrics={
                            "all_passed": report.all_passed,
                            "case_count": report.case_count,
                            "function_count": len(report.summary_rows),
                            "failed_case_count": report.failed_case_count,
                            "skipped_case_count": report.skipped_case_count,
                            "reference_source": args.reference_source,
                        },
                        data={
                            "report": report,
                            "summary_table": summary_path,
                            "observation_table": observation_path,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            if args.reference_source == "geiger-live":
                report = run_geiger_parity_cases(
                    case_ids=args.geiger_cases,
                    rscript_executable=args.geiger_rscript_executable,
                    failure_root=args.geiger_failure_root,
                )
                output_paths: list[Path | str] = []
                summary_path = None
                observation_path = None
                optimizer_triage_path = None
                boundary_warning_path = None
                likelihood_policy_path = None
                model_confidence_path = None
                parameterization_registry_path = None
                if args.summary_out is not None:
                    summary_path = write_geiger_parity_summary_table(
                        args.summary_out,
                        report,
                    )
                    output_paths.append(summary_path)
                if args.observations_out is not None:
                    observation_path = write_geiger_parity_observation_table(
                        args.observations_out,
                        report,
                    )
                    output_paths.append(observation_path)
                if args.optimizer_triage_out is not None:
                    optimizer_triage_path = write_geiger_optimizer_triage_table(
                        args.optimizer_triage_out,
                        report,
                    )
                    output_paths.append(optimizer_triage_path)
                if args.boundary_warning_out is not None:
                    boundary_warning_path = write_geiger_boundary_warning_table(
                        args.boundary_warning_out,
                        report,
                    )
                    output_paths.append(boundary_warning_path)
                if args.likelihood_policy_out is not None:
                    likelihood_policy_path = write_geiger_likelihood_policy_table(
                        args.likelihood_policy_out,
                        report,
                    )
                    output_paths.append(likelihood_policy_path)
                if args.model_confidence_out is not None:
                    model_confidence_path = write_geiger_model_confidence_table(
                        args.model_confidence_out,
                        report,
                    )
                    output_paths.append(model_confidence_path)
                if args.parameterization_registry_out is not None:
                    parameterization_registry_path = (
                        write_geiger_parameterization_registry_table(
                            args.parameterization_registry_out,
                            report,
                        )
                    )
                    output_paths.append(parameterization_registry_path)
                outputs = _finalize_outputs(
                    args,
                    command="parity",
                    inputs=[],
                    outputs=output_paths,
                )
                _print_result(
                    build_command_result(
                        command="parity",
                        inputs=[],
                        outputs=outputs,
                        metrics={
                            "all_passed": report.all_passed,
                            "case_count": report.case_count,
                            "function_count": len(report.summary_rows),
                            "failed_case_count": report.failed_case_count,
                            "skipped_case_count": report.skipped_case_count,
                            "boundary_warning_row_count": len(
                                report.boundary_warning_rows
                            ),
                            "model_confidence_row_count": len(
                                report.model_confidence_rows
                            ),
                            "reference_source": args.reference_source,
                        },
                        data={
                            "report": report,
                            "summary_table": summary_path,
                            "observation_table": observation_path,
                            "optimizer_triage_table": optimizer_triage_path,
                            "boundary_warning_table": boundary_warning_path,
                            "likelihood_policy_table": likelihood_policy_path,
                            "model_confidence_table": model_confidence_path,
                            "parameterization_registry_table": (
                                parameterization_registry_path
                            ),
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            report = validate_reference_parity_examples(include_extended=args.extended)
            output_paths: list[Path | str] = []
            summary_path = None
            observation_path = None
            if args.summary_out is not None:
                summary_path = write_reference_parity_summary_table(
                    args.summary_out,
                    report,
                )
                output_paths.append(summary_path)
            if args.observations_out is not None:
                observation_path = write_reference_parity_observation_table(
                    args.observations_out,
                    report,
                )
                output_paths.append(observation_path)
            outputs = _finalize_outputs(
                args,
                command="parity",
                inputs=[],
                outputs=output_paths,
            )
            _print_result(
                build_command_result(
                    command="parity",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "all_passed": report.all_passed,
                        "case_count": report.case_count,
                        "method_count": len(report.covered_methods),
                        "failed_case_count": report.failed_case_count,
                        "reference_source": args.reference_source,
                        "extended": args.extended,
                    },
                    data={
                        "report": report,
                        "summary_table": summary_path,
                        "observation_table": observation_path,
                    },
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "inspect":
            report = inspect_tree_path(args.tree, source_format=args.format)
            outputs = _finalize_outputs(args, command="inspect", inputs=[args.tree])
            _print_result(
                build_command_result(
                    command="inspect",
                    inputs=[args.tree],
                    outputs=outputs,
                    warnings=report.warnings,
                    metrics={
                        "tip_count": report.tip_count,
                        "node_count": report.node_count,
                        "internal_node_count": report.internal_node_count,
                        "edge_count": report.edge_count,
                        "clade_count": report.clade_count,
                        "is_binary": report.is_binary,
                        "polytomy_count": report.polytomy_count,
                        "branch_length_status": report.branch_length_status,
                        "is_ultrametric": report.is_ultrametric,
                        "tree_diameter": report.tree_diameter,
                        "colless_imbalance_index": report.colless_imbalance_index,
                        "sackin_imbalance_index": report.sackin_imbalance_index,
                        "tree_quality_score": report.tree_quality_score,
                        "zero_length_branch_count": report.zero_length_branch_count,
                        "cherry_count": report.cherry_count,
                        "missing_internal_branch_count": len(
                            report.missing_internal_branch_nodes
                        ),
                        "missing_terminal_branch_count": len(
                            report.missing_terminal_branch_taxa
                        ),
                        "singleton_internal_node_count": len(
                            report.singleton_internal_nodes
                        ),
                        "long_branch_outlier_count": len(report.long_branch_outliers),
                        "short_branch_outlier_count": len(report.short_branch_outliers),
                        "likely_support_label_count": len(report.likely_support_labels),
                        "likely_named_internal_label_count": len(
                            report.likely_named_internal_labels
                        ),
                        "suspicious_support_range_count": len(
                            report.suspicious_support_value_ranges
                        ),
                        "root_classification": report.root_state_confidence.classification,
                        "internal_label_conflict_count": len(
                            report.internal_label_conflicts
                        ),
                        "unsafe_external_label_count": len(
                            report.unsafe_external_labels
                        ),
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
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
