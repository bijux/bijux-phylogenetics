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
    write_ape_parity_observation_table,
    write_ape_parity_summary_table,
)
from bijux_phylogenetics.bayesian import (
    assess_beast_burnin_sensitivity,
    assess_beast_convergence,
    assess_mrbayes_burnin_sensitivity,
    assess_mrbayes_convergence,
    build_bayesian_evidence_package,
    compute_mrbayes_effective_sample_sizes,
    parse_beast_log,
    parse_beast_posterior_tree_samples,
    parse_mrbayes_consensus_tree,
    parse_mrbayes_mcmc_diagnostics,
    parse_mrbayes_parameter_traces,
    parse_mrbayes_posterior_tree_samples,
    prepare_beast_time_tree_analysis,
    prepare_mrbayes_analysis,
    render_bayesian_posterior_report,
    render_calibration_audit_report,
    run_beast_posterior_inference,
    run_mrbayes_posterior_inference,
    subsample_beast_posterior_tree_set,
    subsample_mrbayes_posterior_tree_set,
    summarize_beast_analysis_xml,
    summarize_beast_log,
    summarize_beast_posterior_topology_diversity,
    summarize_beast_posterior_trees,
    summarize_mrbayes_parameter_diagnostics,
    summarize_mrbayes_posterior_trees,
    validate_fossil_calibration_table,
    validate_tip_dating_metadata,
    write_bayesian_methods_summary_text,
    write_beast_burnin_sensitivity_slice_table,
    write_beast_log_summary_table,
    write_beast_posterior_tree_set,
    write_mrbayes_burnin_sensitivity_slice_table,
    write_mrbayes_parameter_summary_table,
    write_posterior_tree_subsample,
    write_posterior_tree_subsample_table,
    write_supplementary_bayesian_diagnostics_table,
)
from bijux_phylogenetics.bayesian.burnin import (
    DEFAULT_BURNIN_FRACTIONS,
    write_burnin_clade_shift_table,
    write_burnin_parameter_shift_table,
)
from bijux_phylogenetics.benchmark import (
    benchmark_alignment_diagnostics,
    benchmark_large_dataset_stress_suite,
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
    _evidence_book_metrics,
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
from bijux_phylogenetics.command_line.distance import (
    add_distance_commands,
    run_distance_command,
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
    _build_annotation_strips,
    _build_numeric_trait_map,
    _build_string_trait_map,
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
from bijux_phylogenetics.core.demo import run_capability_demo
from bijux_phylogenetics.core.environment import inspect_environment
from bijux_phylogenetics.core.metadata import load_taxon_table
from bijux_phylogenetics.core.taxonomy import (
    normalize_tree_taxa,
    write_taxon_mapping,
)
from bijux_phylogenetics.datasets import (
    run_avian_reproductive_trait_demo,
    run_catarrhine_mitogenome_five_locus_panel_demo,
    run_central_european_seashore_flora_demo,
    run_gnathostome_ortholog_protein_benchmark_demo,
    run_influenza_a_ha_reference_demo,
    run_pleistocene_bear_cytb_fragment_demo,
    run_primate_comparative_demo,
    run_rabies_cross_host_geography_panel_demo,
    run_rabies_cross_host_panel_demo,
    run_rabies_geographic_transition_panel_demo,
    run_rabies_method_sensitivity_panel_demo,
)
from bijux_phylogenetics.datasets.continuous_mode_recovery import (
    run_continuous_mode_recovery_panel_demo,
)
from bijux_phylogenetics.datasets.data_quality_stress import (
    run_catarrhine_data_quality_stress_panel_demo,
)
from bijux_phylogenetics.datasets.known_answer_reference import (
    run_known_answer_reference_demo,
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
    compare_diversification_models,
    compute_diversification_gamma_statistic,
    compute_lineage_through_time_curve,
    detect_diversification_outlier_clades,
    detect_incomplete_taxon_sampling_metadata,
    estimate_diversification_rate,
    render_diversification_report,
    run_trait_dependent_diversification_analysis,
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
    compare_fast_and_ml_trees,
    list_external_engine_workflows,
    list_mafft_alignment_modes,
    list_trimal_trimming_modes,
    read_engine_version,
    render_inference_workflow_report,
    run_alignment_trimming,
    run_bootstrap_consensus_tree,
    run_bootstrap_support_estimation,
    run_codon_aware_multiple_sequence_alignment,
    run_fast_tree_inference,
    run_fasta_to_tree_workflow,
    run_maximum_likelihood_tree_inference,
    run_model_selection,
    run_multiple_sequence_alignment,
    run_sh_alrt_support_estimation,
    run_tree_inference_comparison,
)
from bijux_phylogenetics.engines.inference_reproducibility import (
    run_inference_reproducibility_check,
)
from bijux_phylogenetics.engines.large_alignment_inference import (
    run_large_alignment_inference,
)
from bijux_phylogenetics.runtime.errors import (
    EngineUnavailableError,
    EvidenceContractError,
    MetadataJoinError,
    PhylogeneticsError,
)
from bijux_phylogenetics.evidence.book import validate_evidence_book
from bijux_phylogenetics.evidence.bundles import bundle_directory, validate_bundle
from bijux_phylogenetics.evidence.workbench import (
    DOCS_EVIDENCE_OVERVIEW,
    build_evidence_book_selection,
    build_evidence_book_study,
    list_registered_evidence_studies,
    refresh_evidence_book,
    rerun_evidence_book_selection,
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
from bijux_phylogenetics.render.package import build_tree_figure_package
from bijux_phylogenetics.render.svg import audit_support_label_rendering, render_tree_svg
from bijux_phylogenetics.reports.service import (
    render_alignment_report,
    render_dataset_report,
    render_level_one_release_gate_report,
    render_phylo_inputs_report,
    render_release_truth_report,
    render_taxon_report,
    render_tree_report,
    render_workflow_validation_report,
)
from bijux_phylogenetics.reports.tree_package import build_tree_report_package
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
    simulate_protein_alignment,
    simulate_random_trees,
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
from bijux_phylogenetics.trees import (
    cluster_trees_by_topology,
    compute_clade_frequency_table,
    compute_tree_distance_matrix,
    detect_unstable_clades,
    write_clade_frequency_table,
    write_topology_cluster_table,
    write_tree_distance_matrix,
    write_unstable_clade_table,
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
        choices=("fitch", "equal-rates", "symmetric", "all-rates-different"),
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
        choices=("fitch", "equal-rates", "symmetric", "all-rates-different"),
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
        choices=("equal-rates", "symmetric", "all-rates-different"),
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
        choices=("equal-rates", "symmetric", "all-rates-different"),
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
        choices=("equal-rates", "symmetric", "all-rates-different"),
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
        choices=("fitch", "equal-rates", "symmetric", "all-rates-different"),
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
        choices=("equal-rates", "symmetric", "all-rates-different"),
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
        choices=("equal-rates", "symmetric", "all-rates-different"),
        default="equal-rates",
    )
    discrete_compare.add_argument(
        "--right-model",
        choices=("equal-rates", "symmetric", "all-rates-different"),
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
        choices=("equal-rates", "symmetric", "all-rates-different"),
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
        choices=("equal-rates", "symmetric", "all-rates-different"),
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
        choices=("equal-rates", "symmetric", "all-rates-different"),
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
        "--compare-model", choices=("equal-rates", "symmetric", "all-rates-different")
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
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(diversification_report)

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

    parity = subparsers.add_parser(
        get_command_spec("parity").name, help=get_command_spec("parity").summary
    )
    parity.add_argument(
        "--reference-source",
        choices=("checked-fixture", "ape-live", "phytools-live"),
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
    parity.add_argument("--summary-out", type=Path)
    parity.add_argument("--observations-out", type=Path)
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

    render = subparsers.add_parser(
        get_command_spec("render").name, help=get_command_spec("render").summary
    )
    render.add_argument("tree", type=Path)
    render.add_argument("--metadata", type=Path)
    render.add_argument("--traits", type=Path)
    render.add_argument("--taxon-column")
    render.add_argument("--label-column")
    render.add_argument(
        "--layout", choices=["cladogram", "phylogram", "circular"], default="cladogram"
    )
    render.add_argument("--support-labels", action="store_true")
    render.add_argument("--categorical-column")
    render.add_argument("--continuous-column")
    render.add_argument("--metadata-strip-columns")
    render.add_argument("--heatmap-columns")
    render.add_argument("--collapse-clades")
    render.add_argument("--package-dir", type=Path)
    render.add_argument("--out", required=True, type=Path)
    render.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(render)

    evidence = subparsers.add_parser(
        get_command_spec("evidence").name, help=get_command_spec("evidence").summary
    )
    evidence_subparsers = evidence.add_subparsers(
        dest="evidence_command", required=True
    )
    evidence_bundle = evidence_subparsers.add_parser(
        "bundle",
        help="Bundle explicit phylogenetics inputs and outputs as evidence.",
    )
    evidence_bundle.add_argument("--inputs", nargs="+", required=True, type=Path)
    evidence_bundle.add_argument("--outputs", nargs="+", required=True, type=Path)
    evidence_bundle.add_argument("--out", required=True, type=Path)
    evidence_bundle.add_argument(
        "--json", action="store_true", help="Emit the bundle report as JSON."
    )
    _add_manifest_argument(evidence_bundle)
    evidence_validate = evidence_subparsers.add_parser(
        "validate", help="Validate an existing evidence bundle."
    )
    evidence_validate.add_argument("bundle_root", type=Path)
    evidence_validate.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(evidence_validate)
    evidence_book = evidence_subparsers.add_parser(
        "book",
        help="Govern evidence-book generation, validation, and partial reruns.",
    )
    evidence_book_subparsers = evidence_book.add_subparsers(
        dest="evidence_book_command",
        required=True,
    )
    evidence_book_studies = evidence_book_subparsers.add_parser(
        "studies",
        help="List governed evidence-book studies and partial rerun capabilities.",
    )
    evidence_book_studies.add_argument(
        "--json", action="store_true", help="Emit the study registry as JSON."
    )
    _add_manifest_argument(evidence_book_studies)
    evidence_book_build = evidence_book_subparsers.add_parser(
        "build",
        help="Refresh governed evidence-book outputs or rebuild one registered study.",
    )
    evidence_book_build.add_argument(
        "study_id",
        nargs="?",
        help="Optional registered study identifier to rebuild before refreshing the evidence-book.",
    )
    evidence_book_build.add_argument(
        "--evidence-id",
        dest="evidence_ids",
        action="append",
        default=[],
        help="Optional Evidence ID to rebuild within the selected study. May be repeated.",
    )
    evidence_book_build.add_argument(
        "--json", action="store_true", help="Emit the build report as JSON."
    )
    _add_manifest_argument(evidence_book_build)
    evidence_book_validate = evidence_book_subparsers.add_parser(
        "validate",
        help="Validate the governed evidence-book surface and summarize coverage gaps.",
    )
    evidence_book_validate.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(evidence_book_validate)
    evidence_book_rerun = evidence_book_subparsers.add_parser(
        "rerun",
        help="Regenerate selected Evidence IDs for a study and refresh governed outputs.",
    )
    evidence_book_rerun.add_argument("study_id")
    evidence_book_rerun.add_argument("evidence_ids", nargs="+")
    evidence_book_rerun.add_argument(
        "--json", action="store_true", help="Emit the rerun report as JSON."
    )
    _add_manifest_argument(evidence_book_rerun)

    report = subparsers.add_parser(
        get_command_spec("report").name, help=get_command_spec("report").summary
    )
    report_subparsers = report.add_subparsers(dest="report_command", required=True)
    report_tree = report_subparsers.add_parser(
        "tree", help="Render a deterministic single-tree HTML report."
    )
    report_tree.add_argument("tree", type=Path)
    report_tree.add_argument("--out", required=True, type=Path)
    report_tree.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_tree)
    report_tree_package = report_subparsers.add_parser(
        "tree-package",
        help="Build a full tree report package with figure and TSV ledgers.",
    )
    report_tree_package.add_argument("tree", type=Path)
    report_tree_package.add_argument("--out-dir", required=True, type=Path)
    report_tree_package.add_argument(
        "--json", action="store_true", help="Emit the package build result as JSON."
    )
    _add_manifest_argument(report_tree_package)
    report_alignment = report_subparsers.add_parser(
        "alignment", help="Render an alignment-only HTML diagnostic report."
    )
    report_alignment.add_argument("--alignment", required=True, type=Path)
    report_alignment.add_argument("--out", required=True, type=Path)
    report_alignment.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_alignment)
    report_dataset = report_subparsers.add_parser(
        "dataset", help="Render a tree plus table dataset HTML report."
    )
    report_dataset.add_argument("--tree", required=True, type=Path)
    report_dataset.add_argument("--metadata", required=True, type=Path)
    report_dataset.add_argument("--traits", type=Path)
    report_dataset.add_argument("--alignment", type=Path)
    report_dataset.add_argument("--tip-dates", type=Path)
    report_dataset.add_argument("--calibrations", type=Path)
    report_dataset.add_argument("--out", required=True, type=Path)
    report_dataset.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_dataset)
    report_phylo_inputs = report_subparsers.add_parser(
        "phylo-inputs",
        help="Render a tree plus alignment HTML input report.",
    )
    report_phylo_inputs.add_argument("--tree", required=True, type=Path)
    report_phylo_inputs.add_argument("--alignment", required=True, type=Path)
    report_phylo_inputs.add_argument("--out", required=True, type=Path)
    report_phylo_inputs.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_phylo_inputs)
    report_taxonomy = report_subparsers.add_parser(
        "taxonomy", help="Render a reviewer-facing taxon audit HTML report."
    )
    report_taxonomy.add_argument("--tree", required=True, type=Path)
    report_taxonomy.add_argument("--synonym-table", type=Path)
    report_taxonomy.add_argument("--metadata", type=Path)
    report_taxonomy.add_argument("--traits", type=Path)
    report_taxonomy.add_argument("--alignment", type=Path)
    report_taxonomy.add_argument("--filtered-alignment", type=Path)
    report_taxonomy.add_argument("--inference-tree", type=Path)
    report_taxonomy.add_argument("--reported-taxa", type=Path)
    report_taxonomy.add_argument("--out", required=True, type=Path)
    report_taxonomy.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_taxonomy)
    report_workflow_validation = report_subparsers.add_parser(
        "workflow-validation",
        help="Render the Level 1 workflow validation fixture report.",
    )
    report_workflow_validation.add_argument("--fixtures-root", type=Path)
    report_workflow_validation.add_argument("--out", required=True, type=Path)
    report_workflow_validation.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_workflow_validation)
    report_release_gate = report_subparsers.add_parser(
        "release-gate",
        help="Render the Level 1 release gate for the checked-in workflow fixtures.",
    )
    report_release_gate.add_argument("--fixtures-root", type=Path)
    report_release_gate.add_argument("--out", required=True, type=Path)
    report_release_gate.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_release_gate)
    report_release_truth = report_subparsers.add_parser(
        "release-truth",
        help="Render one machine-produced release truth report from pytest and workflow evidence.",
    )
    report_release_truth.add_argument(
        "--test-report",
        type=Path,
        action="append",
        required=True,
        help="Path to one pytest JUnit XML report for the full test surface. Repeat to aggregate multiple sessions.",
    )
    report_release_truth.add_argument(
        "--real-engine-test-report",
        type=Path,
        action="append",
        required=True,
        help="Path to one pytest JUnit XML report for real-engine tests. Repeat to aggregate multiple sessions.",
    )
    report_release_truth.add_argument("--fixtures-root", type=Path)
    report_release_truth.add_argument(
        "--stress-tier",
        choices=("small", "heavy"),
        default="small",
        help="Governed stress tier to benchmark during release truth generation.",
    )
    report_release_truth.add_argument(
        "--parity-extended",
        action="store_true",
        help="Include the governed extended reference-parity suite in the release truth report.",
    )
    report_release_truth.add_argument("--out", required=True, type=Path)
    report_release_truth.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_release_truth)

    demo = subparsers.add_parser(
        get_command_spec("demo").name, help=get_command_spec("demo").summary
    )
    demo_subparsers = demo.add_subparsers(dest="demo_command", required=True)
    demo_run = demo_subparsers.add_parser(
        "run", help="Run the repository capability demo workflow."
    )
    demo_run.add_argument("--out", required=True, type=Path)
    demo_run.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_run)
    demo_primate = demo_subparsers.add_parser(
        "primate-comparative",
        help="Materialize the packaged primate dataset and comparative workflow outputs.",
    )
    demo_primate.add_argument("--out", required=True, type=Path)
    demo_primate.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_primate)
    demo_birds = demo_subparsers.add_parser(
        "avian-reproductive-traits",
        help="Materialize the packaged avian reproductive dataset and workflow outputs.",
    )
    demo_birds.add_argument("--out", required=True, type=Path)
    demo_birds.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_birds)
    demo_plants = demo_subparsers.add_parser(
        "central-european-seashore-flora",
        help="Materialize the packaged Central European plant dataset and workflow outputs.",
    )
    demo_plants.add_argument("--out", required=True, type=Path)
    demo_plants.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_plants)
    demo_viruses = demo_subparsers.add_parser(
        "influenza-a-ha-reference-panel",
        help="Materialize the packaged influenza A HA dataset and rerun the sequence-to-tree workflow outputs.",
    )
    demo_viruses.add_argument("--out", required=True, type=Path)
    demo_viruses.add_argument("--mafft-executable", type=str)
    demo_viruses.add_argument("--trimal-executable", type=str)
    demo_viruses.add_argument("--iqtree-executable", type=str)
    demo_viruses.add_argument("--iqtree-seed", type=int, default=1)
    demo_viruses.add_argument("--iqtree-threads", type=int, default=1)
    demo_viruses.add_argument("--bootstrap-replicates", type=int, default=1000)
    demo_viruses.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_viruses)
    demo_protein_benchmark = demo_subparsers.add_parser(
        "gnathostome-ortholog-protein-benchmark",
        help="Materialize the packaged gnathostome protein benchmark and rerun the governed amino-acid sequence-to-tree outputs.",
    )
    demo_protein_benchmark.add_argument("--out", required=True, type=Path)
    demo_protein_benchmark.add_argument("--mafft-executable", type=str)
    demo_protein_benchmark.add_argument("--trimal-executable", type=str)
    demo_protein_benchmark.add_argument("--iqtree-executable", type=str)
    demo_protein_benchmark.add_argument("--iqtree-seed", type=int, default=1)
    demo_protein_benchmark.add_argument("--iqtree-threads", type=int, default=1)
    demo_protein_benchmark.add_argument(
        "--bootstrap-replicates", type=int, default=1000
    )
    demo_protein_benchmark.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_protein_benchmark)
    demo_ancient_dna = demo_subparsers.add_parser(
        "pleistocene-bear-cytb-fragments",
        help="Materialize the packaged ancient-DNA-style bear dataset and rerun the sequence-to-tree workflow outputs.",
    )
    demo_ancient_dna.add_argument("--out", required=True, type=Path)
    demo_ancient_dna.add_argument("--mafft-executable", type=str)
    demo_ancient_dna.add_argument("--trimal-executable", type=str)
    demo_ancient_dna.add_argument("--iqtree-executable", type=str)
    demo_ancient_dna.add_argument("--iqtree-seed", type=int, default=1)
    demo_ancient_dna.add_argument("--iqtree-threads", type=int, default=1)
    demo_ancient_dna.add_argument("--bootstrap-replicates", type=int, default=1000)
    demo_ancient_dna.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_ancient_dna)
    demo_rabies = demo_subparsers.add_parser(
        "rabies-cross-host-panel",
        help="Materialize the packaged rabies host-switching dataset and rerun the governed host-transition review outputs.",
    )
    demo_rabies.add_argument("--out", required=True, type=Path)
    demo_rabies.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_rabies)
    demo_rabies_geography = demo_subparsers.add_parser(
        "rabies-geographic-transition-panel",
        help="Materialize the packaged rabies geography dataset and rerun the governed geographic transition review outputs.",
    )
    demo_rabies_geography.add_argument("--out", required=True, type=Path)
    demo_rabies_geography.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_rabies_geography)
    demo_rabies_host_geography = demo_subparsers.add_parser(
        "rabies-cross-host-geography-panel",
        help="Materialize the packaged rabies integrated dataset and rerun the full sequence-to-tree, host, and geography workflow outputs.",
    )
    demo_rabies_host_geography.add_argument("--out", required=True, type=Path)
    demo_rabies_host_geography.add_argument(
        "--config",
        type=Path,
        help="Optional workflow config JSON. Defaults to the packaged dataset config.",
    )
    demo_rabies_host_geography.add_argument("--mafft-executable", type=str)
    demo_rabies_host_geography.add_argument("--trimal-executable", type=str)
    demo_rabies_host_geography.add_argument("--iqtree-executable", type=str)
    demo_rabies_host_geography.add_argument("--fasttree-executable", type=str)
    demo_rabies_host_geography.add_argument("--iqtree-seed", type=int, default=1)
    demo_rabies_host_geography.add_argument("--iqtree-threads", type=int, default=1)
    demo_rabies_host_geography.add_argument(
        "--bootstrap-replicates", type=int, default=1000
    )
    demo_rabies_host_geography.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_rabies_host_geography)
    demo_rabies_method_sensitivity = demo_subparsers.add_parser(
        "rabies-method-sensitivity-panel",
        help="Materialize the packaged rabies method-sensitivity dataset and rerun the governed preprocessing and engine-comparison workflow outputs.",
    )
    demo_rabies_method_sensitivity.add_argument("--out", required=True, type=Path)
    demo_rabies_method_sensitivity.add_argument("--mafft-executable", type=str)
    demo_rabies_method_sensitivity.add_argument("--trimal-executable", type=str)
    demo_rabies_method_sensitivity.add_argument("--iqtree-executable", type=str)
    demo_rabies_method_sensitivity.add_argument("--fasttree-executable", type=str)
    demo_rabies_method_sensitivity.add_argument("--iqtree-seed", type=int, default=1)
    demo_rabies_method_sensitivity.add_argument("--iqtree-threads", type=int, default=1)
    demo_rabies_method_sensitivity.add_argument(
        "--bootstrap-replicates", type=int, default=1000
    )
    demo_rabies_method_sensitivity.add_argument(
        "--parallel-workers",
        type=int,
        default=None,
        help="Number of isolated variant workers to run in parallel. Defaults to the packaged workflow config.",
    )
    demo_rabies_method_sensitivity.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_rabies_method_sensitivity)
    demo_catarrhine_mitogenome = demo_subparsers.add_parser(
        "catarrhine-mitogenome-five-locus-panel",
        help="Materialize the packaged catarrhine multi-locus dataset and rerun the governed concatenation and partitioned inference outputs.",
    )
    demo_catarrhine_mitogenome.add_argument("--out", required=True, type=Path)
    demo_catarrhine_mitogenome.add_argument("--iqtree-executable", type=str)
    demo_catarrhine_mitogenome.add_argument("--iqtree-seed", type=int, default=1)
    demo_catarrhine_mitogenome.add_argument("--iqtree-threads", type=int, default=1)
    demo_catarrhine_mitogenome.add_argument(
        "--bootstrap-replicates", type=int, default=1000
    )
    demo_catarrhine_mitogenome.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_catarrhine_mitogenome)
    demo_catarrhine_stress = demo_subparsers.add_parser(
        "catarrhine-data-quality-stress-panel",
        help="Materialize the packaged catarrhine dirty-data stress dataset and rerun the governed audit and cleanup outputs.",
    )
    demo_catarrhine_stress.add_argument("--out", required=True, type=Path)
    demo_catarrhine_stress.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_catarrhine_stress)
    demo_continuous_mode_recovery = demo_subparsers.add_parser(
        "continuous-mode-recovery-panel",
        help="Materialize the packaged continuous-trait recovery dataset and rerun the governed simulation-recovery outputs.",
    )
    demo_continuous_mode_recovery.add_argument("--out", required=True, type=Path)
    demo_continuous_mode_recovery.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_continuous_mode_recovery)
    demo_known_answer = demo_subparsers.add_parser(
        "known-answer-reference-panel",
        help="Materialize the packaged known-answer simulation dataset and rerun the governed recovery outputs.",
    )
    demo_known_answer.add_argument("--out", required=True, type=Path)
    demo_known_answer.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_known_answer)

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
    adapter_fast = adapter_subparsers.add_parser(
        "infer-fast", help="Run fast approximate tree inference."
    )
    adapter_fast.add_argument("input_path", type=Path)
    adapter_fast.add_argument("--out", required=True, type=Path)
    adapter_fast.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_fast.add_argument("--executable", type=str)
    adapter_fast.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_fast)
    _add_manifest_argument(adapter_fast)
    adapter_large = adapter_subparsers.add_parser(
        "infer-large",
        help="Run large-alignment FastTree inference with streamed preflight and resource reporting.",
    )
    adapter_large.add_argument("input_path", type=Path)
    adapter_large.add_argument("--out-dir", required=True, type=Path)
    adapter_large.add_argument("--prefix", default="large-alignment-inference")
    adapter_large.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_large.add_argument("--executable", type=str)
    adapter_large.add_argument("--resume", action="store_true")
    adapter_large.add_argument(
        "--timeout-seconds",
        type=float,
        help="Stop the FastTree inference step if it exceeds this wall-clock budget.",
    )
    adapter_large.add_argument(
        "--incomplete-run-policy",
        choices=("reject", "clean"),
        default="reject",
        help="Reject or clean incomplete FastTree outputs before starting a new large-alignment run.",
    )
    adapter_large.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_manifest_argument(adapter_large)
    adapter_compare = adapter_subparsers.add_parser(
        "compare", help="Compare fast approximate and ML trees."
    )
    adapter_compare.add_argument("--fast-tree", required=True, type=Path)
    adapter_compare.add_argument("--ml-tree", required=True, type=Path)
    adapter_compare.add_argument("--out", required=True, type=Path)
    adapter_compare.add_argument(
        "--json", action="store_true", help="Emit the comparison report as JSON."
    )
    _add_manifest_argument(adapter_compare)
    adapter_compare_engines = adapter_subparsers.add_parser(
        "compare-engines",
        help="Run IQ-TREE and FastTree on one alignment and compare the inferred trees.",
    )
    adapter_compare_engines.add_argument("input_path", type=Path)
    adapter_compare_engines.add_argument("--out-dir", required=True, type=Path)
    adapter_compare_engines.add_argument("--prefix", default="engine-comparison")
    adapter_compare_engines.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_compare_engines.add_argument("--iqtree-executable", type=str)
    adapter_compare_engines.add_argument("--fasttree-executable", type=str)
    adapter_compare_engines.add_argument(
        "--iqtree-seed",
        type=int,
        default=1,
        help="Set the IQ-TREE random seed for deterministic comparison runs.",
    )
    adapter_compare_engines.add_argument(
        "--iqtree-threads",
        type=int,
        default=1,
        help="Set the IQ-TREE thread count used during the comparison run.",
    )
    adapter_compare_engines.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=1000,
        help="Set the ultrafast bootstrap replicate count used for the IQ-TREE support workflow.",
    )
    adapter_compare_engines.add_argument(
        "--json",
        action="store_true",
        help="Emit the comparison workflow report as JSON.",
    )
    _add_external_adapter_execution_arguments(adapter_compare_engines)
    _add_manifest_argument(adapter_compare_engines)
    adapter_reproducibility = adapter_subparsers.add_parser(
        "reproducibility",
        help="Rerun supported IQ-TREE inference and classify deterministic versus unstable outputs.",
    )
    adapter_reproducibility.add_argument("input_path", type=Path)
    adapter_reproducibility.add_argument("--out-dir", required=True, type=Path)
    adapter_reproducibility.add_argument(
        "--prefix", default="inference-reproducibility"
    )
    adapter_reproducibility.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_reproducibility.add_argument("--iqtree-executable", type=str)
    adapter_reproducibility.add_argument(
        "--iqtree-seed",
        type=int,
        default=1,
        help="Set the IQ-TREE random seed used for every rerun.",
    )
    adapter_reproducibility.add_argument(
        "--iqtree-threads",
        type=int,
        default=1,
        help="Set the IQ-TREE thread count used for every rerun.",
    )
    adapter_reproducibility.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=1000,
        help="Set the ultrafast bootstrap replicate count used for every rerun.",
    )
    adapter_reproducibility.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="Set how many repeated supported-inference runs to compare.",
    )
    adapter_reproducibility.add_argument(
        "--json",
        action="store_true",
        help="Emit the reproducibility workflow report as JSON.",
    )
    _add_manifest_argument(adapter_reproducibility)
    adapter_mrbayes_prepare = adapter_subparsers.add_parser(
        "mrbayes-prepare",
        help="Prepare a MrBayes NEXUS analysis from an aligned FASTA file.",
    )
    adapter_mrbayes_prepare.add_argument("input_path", type=Path)
    adapter_mrbayes_prepare.add_argument("--out", required=True, type=Path)
    adapter_mrbayes_prepare.add_argument("--partitions", type=Path)
    adapter_mrbayes_prepare.add_argument("--model", default="gtr")
    adapter_mrbayes_prepare.add_argument("--rates", default="gamma")
    adapter_mrbayes_prepare.add_argument("--ngen", type=int, default=10000)
    adapter_mrbayes_prepare.add_argument("--nchains", type=int, default=4)
    adapter_mrbayes_prepare.add_argument("--samplefreq", type=int, default=100)
    adapter_mrbayes_prepare.add_argument("--printfreq", type=int, default=100)
    adapter_mrbayes_prepare.add_argument("--burnin-fraction", type=float, default=0.25)
    adapter_mrbayes_prepare.add_argument(
        "--json", action="store_true", help="Emit the preparation report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_prepare)
    adapter_mrbayes_run = adapter_subparsers.add_parser(
        "mrbayes-run",
        help="Run a prepared MrBayes posterior inference workflow.",
    )
    adapter_mrbayes_run.add_argument("input_path", type=Path)
    adapter_mrbayes_run.add_argument("--executable", type=str)
    adapter_mrbayes_run.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_mrbayes_run)
    _add_manifest_argument(adapter_mrbayes_run)
    adapter_mrbayes_summarize = adapter_subparsers.add_parser(
        "mrbayes-summarize",
        help="Summarize MrBayes posterior trees after burn-in removal.",
    )
    adapter_mrbayes_summarize.add_argument("input_path", type=Path)
    adapter_mrbayes_summarize.add_argument(
        "--burnin-fraction", type=float, default=0.25
    )
    adapter_mrbayes_summarize.add_argument(
        "--json", action="store_true", help="Emit the posterior summary as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_summarize)
    adapter_mrbayes_traces = adapter_subparsers.add_parser(
        "mrbayes-traces",
        help="Parse a MrBayes parameter trace table.",
    )
    adapter_mrbayes_traces.add_argument("input_path", type=Path)
    adapter_mrbayes_traces.add_argument(
        "--json", action="store_true", help="Emit the trace report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_traces)
    adapter_mrbayes_trees = adapter_subparsers.add_parser(
        "mrbayes-trees",
        help="Parse a MrBayes posterior tree set into sampled trees.",
    )
    adapter_mrbayes_trees.add_argument("input_path", type=Path)
    adapter_mrbayes_trees.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior tree set report as JSON.",
    )
    _add_manifest_argument(adapter_mrbayes_trees)
    adapter_mrbayes_subsample = adapter_subparsers.add_parser(
        "mrbayes-subsample",
        help="Subsample MrBayes posterior trees while preserving generation metadata.",
    )
    adapter_mrbayes_subsample.add_argument("input_path", type=Path)
    adapter_mrbayes_subsample.add_argument(
        "--method",
        required=True,
        choices=("evenly-spaced", "random"),
        help="Select evenly spaced thinning or a seeded random retained subset.",
    )
    adapter_mrbayes_subsample.add_argument("--burnin-fraction", type=float, default=0.0)
    adapter_mrbayes_subsample.add_argument("--thinning-interval", type=int)
    adapter_mrbayes_subsample.add_argument("--sample-count", type=int)
    adapter_mrbayes_subsample.add_argument("--seed", type=int)
    adapter_mrbayes_subsample.add_argument(
        "--tree-set-out",
        type=Path,
        help="Write the retained posterior tree set as normalized Newick.",
    )
    adapter_mrbayes_subsample.add_argument(
        "--sample-table-out",
        type=Path,
        help="Write a TSV ledger of retained posterior-tree metadata.",
    )
    adapter_mrbayes_subsample.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior subsampling report as JSON.",
    )
    _add_manifest_argument(adapter_mrbayes_subsample)
    adapter_mrbayes_mcmc = adapter_subparsers.add_parser(
        "mrbayes-mcmc",
        help="Parse a MrBayes MCMC diagnostics table.",
    )
    adapter_mrbayes_mcmc.add_argument("input_path", type=Path)
    adapter_mrbayes_mcmc.add_argument(
        "--json", action="store_true", help="Emit the MCMC diagnostics report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_mcmc)
    adapter_mrbayes_consensus = adapter_subparsers.add_parser(
        "mrbayes-consensus",
        help="Parse a MrBayes consensus tree with posterior-probability annotations.",
    )
    adapter_mrbayes_consensus.add_argument("input_path", type=Path)
    adapter_mrbayes_consensus.add_argument(
        "--json", action="store_true", help="Emit the consensus tree report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_consensus)
    adapter_mrbayes_ess = adapter_subparsers.add_parser(
        "mrbayes-ess",
        help="Compute effective sample sizes from a MrBayes trace table.",
    )
    adapter_mrbayes_ess.add_argument("input_path", type=Path)
    adapter_mrbayes_ess.add_argument(
        "--json", action="store_true", help="Emit the ESS report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_ess)
    adapter_mrbayes_parameters = adapter_subparsers.add_parser(
        "mrbayes-parameters",
        help="Summarize burn-in-aware posterior parameter diagnostics from a MrBayes trace table.",
    )
    adapter_mrbayes_parameters.add_argument("input_path", type=Path)
    adapter_mrbayes_parameters.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this fraction of early samples before reporting posterior summaries.",
    )
    adapter_mrbayes_parameters.add_argument(
        "--summary-out",
        type=Path,
        help="Write a TSV parameter-summary table for the retained trace samples.",
    )
    adapter_mrbayes_parameters.add_argument(
        "--json", action="store_true", help="Emit the parameter diagnostics as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_parameters)
    adapter_mrbayes_burnin = adapter_subparsers.add_parser(
        "mrbayes-burnin-sensitivity",
        help="Compare MrBayes posterior summaries across multiple burn-in fractions.",
    )
    adapter_mrbayes_burnin.add_argument("posterior_trees", type=Path)
    adapter_mrbayes_burnin.add_argument("--traces", type=Path)
    adapter_mrbayes_burnin.add_argument(
        "--burnin-fractions",
        nargs="+",
        type=float,
        default=list(DEFAULT_BURNIN_FRACTIONS),
    )
    adapter_mrbayes_burnin.add_argument("--slice-out", type=Path)
    adapter_mrbayes_burnin.add_argument("--parameter-out", type=Path)
    adapter_mrbayes_burnin.add_argument("--clade-out", type=Path)
    adapter_mrbayes_burnin.add_argument(
        "--json",
        action="store_true",
        help="Emit the burn-in sensitivity report as JSON.",
    )
    _add_manifest_argument(adapter_mrbayes_burnin)
    adapter_mrbayes_convergence = adapter_subparsers.add_parser(
        "mrbayes-convergence",
        help="Assess MrBayes trace convergence from ESS and trace drift.",
    )
    adapter_mrbayes_convergence.add_argument("input_path", type=Path)
    adapter_mrbayes_convergence.add_argument(
        "--ess-threshold", type=float, default=200.0
    )
    adapter_mrbayes_convergence.add_argument(
        "--mean-shift-threshold", type=float, default=0.5
    )
    adapter_mrbayes_convergence.add_argument(
        "--json", action="store_true", help="Emit the convergence report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_convergence)
    adapter_mrbayes_report = adapter_subparsers.add_parser(
        "mrbayes-report",
        help="Render an HTML Bayesian posterior report from posterior trees and traces.",
    )
    adapter_mrbayes_report.add_argument("posterior_trees", type=Path)
    adapter_mrbayes_report.add_argument("--traces", required=True, type=Path)
    adapter_mrbayes_report.add_argument("--out", required=True, type=Path)
    adapter_mrbayes_report.add_argument("--burnin-fraction", type=float, default=0.25)
    adapter_mrbayes_report.add_argument("--ess-threshold", type=float, default=200.0)
    adapter_mrbayes_report.add_argument(
        "--mean-shift-threshold", type=float, default=0.5
    )
    adapter_mrbayes_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_report)
    adapter_beast_prepare = adapter_subparsers.add_parser(
        "beast-prepare",
        help="Prepare a BEAST-style time-tree XML analysis from aligned sequences and dating inputs.",
    )
    adapter_beast_prepare.add_argument("input_path", type=Path)
    adapter_beast_prepare.add_argument("--out", required=True, type=Path)
    adapter_beast_prepare.add_argument("--tree", type=Path)
    adapter_beast_prepare.add_argument("--calibrations", type=Path)
    adapter_beast_prepare.add_argument("--tip-dates", type=Path)
    adapter_beast_prepare.add_argument("--clock-model", default="strict")
    adapter_beast_prepare.add_argument("--tree-prior", default="yule")
    adapter_beast_prepare.add_argument("--chain-length", type=int, default=1000000)
    adapter_beast_prepare.add_argument("--log-every", type=int, default=1000)
    adapter_beast_prepare.add_argument(
        "--json", action="store_true", help="Emit the preparation report as JSON."
    )
    _add_manifest_argument(adapter_beast_prepare)
    adapter_beast_xml = adapter_subparsers.add_parser(
        "beast-xml",
        help="Summarize and validate one prepared BEAST analysis XML.",
    )
    adapter_beast_xml.add_argument("input_path", type=Path)
    adapter_beast_xml.add_argument(
        "--json", action="store_true", help="Emit the XML summary report as JSON."
    )
    _add_manifest_argument(adapter_beast_xml)
    adapter_beast_run = adapter_subparsers.add_parser(
        "beast-run",
        help="Run a prepared BEAST posterior inference workflow.",
    )
    adapter_beast_run.add_argument("input_path", type=Path)
    adapter_beast_run.add_argument("--executable", type=str)
    adapter_beast_run.add_argument("--threads", type=int, default=1)
    adapter_beast_run.add_argument("--seed", type=int, default=1)
    adapter_beast_run.add_argument(
        "--no-overwrite",
        action="store_true",
        help="Keep any existing posterior outputs instead of passing the BEAST overwrite flag.",
    )
    adapter_beast_run.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_beast_run)
    _add_manifest_argument(adapter_beast_run)
    adapter_beast_calibrations = adapter_subparsers.add_parser(
        "beast-calibrations",
        help="Validate a fossil calibration table against a tree.",
    )
    adapter_beast_calibrations.add_argument("tree_path", type=Path)
    adapter_beast_calibrations.add_argument("calibration_path", type=Path)
    adapter_beast_calibrations.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(adapter_beast_calibrations)
    adapter_beast_tip_dates = adapter_subparsers.add_parser(
        "beast-tip-dates",
        help="Validate tip-dating metadata against a tree and optional alignment.",
    )
    adapter_beast_tip_dates.add_argument("tree_path", type=Path)
    adapter_beast_tip_dates.add_argument("tip_dates_path", type=Path)
    adapter_beast_tip_dates.add_argument("--alignment", type=Path)
    adapter_beast_tip_dates.add_argument("--date-column", default="date")
    adapter_beast_tip_dates.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(adapter_beast_tip_dates)
    adapter_beast_log = adapter_subparsers.add_parser(
        "beast-log",
        help="Parse a BEAST log file into a deterministic numeric trace table.",
    )
    adapter_beast_log.add_argument("input_path", type=Path)
    adapter_beast_log.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this fraction of early samples before reporting summaries.",
    )
    adapter_beast_log.add_argument(
        "--summary-out",
        type=Path,
        help="Write a TSV parameter-summary table for the parsed log.",
    )
    adapter_beast_log.add_argument(
        "--json", action="store_true", help="Emit the parsed log report as JSON."
    )
    _add_manifest_argument(adapter_beast_log)
    adapter_beast_burnin = adapter_subparsers.add_parser(
        "beast-burnin-sensitivity",
        help="Compare BEAST posterior summaries across multiple burn-in fractions.",
    )
    adapter_beast_burnin.add_argument("posterior_trees", type=Path)
    adapter_beast_burnin.add_argument("--log", type=Path)
    adapter_beast_burnin.add_argument(
        "--burnin-fractions",
        nargs="+",
        type=float,
        default=list(DEFAULT_BURNIN_FRACTIONS),
    )
    adapter_beast_burnin.add_argument("--slice-out", type=Path)
    adapter_beast_burnin.add_argument("--parameter-out", type=Path)
    adapter_beast_burnin.add_argument("--clade-out", type=Path)
    adapter_beast_burnin.add_argument(
        "--json",
        action="store_true",
        help="Emit the burn-in sensitivity report as JSON.",
    )
    _add_manifest_argument(adapter_beast_burnin)
    adapter_beast_parameters = adapter_subparsers.add_parser(
        "beast-parameters",
        help="Summarize burn-in-aware posterior parameter diagnostics from a BEAST log.",
    )
    adapter_beast_parameters.add_argument("input_path", type=Path)
    adapter_beast_parameters.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this fraction of early samples before reporting posterior summaries.",
    )
    adapter_beast_parameters.add_argument(
        "--summary-out",
        type=Path,
        help="Write a TSV parameter-summary table for the retained log samples.",
    )
    adapter_beast_parameters.add_argument(
        "--json", action="store_true", help="Emit the parameter diagnostics as JSON."
    )
    _add_manifest_argument(adapter_beast_parameters)
    adapter_beast_trees = adapter_subparsers.add_parser(
        "beast-trees",
        help="Parse a BEAST posterior tree set into state-tagged normalized trees.",
    )
    adapter_beast_trees.add_argument("input_path", type=Path)
    adapter_beast_trees.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this fraction of early sampled trees before reporting summaries.",
    )
    adapter_beast_trees.add_argument(
        "--tree-set-out",
        type=Path,
        help="Write the retained posterior tree set as normalized Newick.",
    )
    adapter_beast_trees.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior tree set report as JSON.",
    )
    _add_manifest_argument(adapter_beast_trees)
    adapter_beast_subsample = adapter_subparsers.add_parser(
        "beast-subsample",
        help="Subsample BEAST posterior trees while preserving state metadata.",
    )
    adapter_beast_subsample.add_argument("input_path", type=Path)
    adapter_beast_subsample.add_argument(
        "--method",
        required=True,
        choices=("evenly-spaced", "random"),
        help="Select evenly spaced thinning or a seeded random retained subset.",
    )
    adapter_beast_subsample.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this fraction of early sampled trees before subsampling.",
    )
    adapter_beast_subsample.add_argument("--thinning-interval", type=int)
    adapter_beast_subsample.add_argument("--sample-count", type=int)
    adapter_beast_subsample.add_argument("--seed", type=int)
    adapter_beast_subsample.add_argument(
        "--tree-set-out",
        type=Path,
        help="Write the retained posterior tree set as normalized Newick.",
    )
    adapter_beast_subsample.add_argument(
        "--sample-table-out",
        type=Path,
        help="Write a TSV ledger of retained posterior-tree metadata.",
    )
    adapter_beast_subsample.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior subsampling report as JSON.",
    )
    _add_manifest_argument(adapter_beast_subsample)
    adapter_beast_consensus = adapter_subparsers.add_parser(
        "beast-consensus",
        help="Build a majority-rule consensus tree from BEAST posterior tree samples.",
    )
    adapter_beast_consensus.add_argument("input_path", type=Path)
    adapter_beast_consensus.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.25,
        help="Discard this fraction of early sampled trees before consensus building.",
    )
    adapter_beast_consensus.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Write the posterior-probability-annotated consensus tree as Newick.",
    )
    adapter_beast_consensus.add_argument(
        "--tree-set-out",
        type=Path,
        help="Write the retained posterior tree set as normalized Newick.",
    )
    adapter_beast_consensus.add_argument(
        "--clade-table-out",
        type=Path,
        help="Write the retained clade-frequency ledger as TSV.",
    )
    adapter_beast_consensus.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior consensus report as JSON.",
    )
    _add_manifest_argument(adapter_beast_consensus)
    adapter_beast_diversity = adapter_subparsers.add_parser(
        "beast-diversity",
        help="Summarize topology diversity across BEAST posterior tree samples.",
    )
    adapter_beast_diversity.add_argument("input_path", type=Path)
    adapter_beast_diversity.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.25,
        help="Discard this fraction of early sampled trees before topology review.",
    )
    adapter_beast_diversity.add_argument(
        "--tree-set-out",
        type=Path,
        help="Write the retained posterior tree set as normalized Newick.",
    )
    adapter_beast_diversity.add_argument(
        "--distance-out",
        type=Path,
        help="Write the pairwise RF distance ledger as TSV.",
    )
    adapter_beast_diversity.add_argument(
        "--topology-out",
        type=Path,
        help="Write the rooted topology cluster ledger as TSV.",
    )
    adapter_beast_diversity.add_argument(
        "--unstable-clade-out",
        type=Path,
        help="Write the unstable-clade ledger as TSV.",
    )
    adapter_beast_diversity.add_argument(
        "--json",
        action="store_true",
        help="Emit the topology diversity report as JSON.",
    )
    _add_manifest_argument(adapter_beast_diversity)
    adapter_beast_convergence = adapter_subparsers.add_parser(
        "beast-convergence",
        help="Assess BEAST log convergence from ESS and trace drift.",
    )
    adapter_beast_convergence.add_argument("input_path", type=Path)
    adapter_beast_convergence.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this fraction of early samples before assessing convergence.",
    )
    adapter_beast_convergence.add_argument("--ess-threshold", type=float, default=200.0)
    adapter_beast_convergence.add_argument(
        "--mean-shift-threshold", type=float, default=0.5
    )
    adapter_beast_convergence.add_argument(
        "--json", action="store_true", help="Emit the convergence report as JSON."
    )
    _add_manifest_argument(adapter_beast_convergence)
    adapter_beast_calibration_report = adapter_subparsers.add_parser(
        "beast-calibration-report",
        help="Render an HTML calibration audit report.",
    )
    adapter_beast_calibration_report.add_argument("tree_path", type=Path)
    adapter_beast_calibration_report.add_argument("calibration_path", type=Path)
    adapter_beast_calibration_report.add_argument("--out", required=True, type=Path)
    adapter_beast_calibration_report.add_argument("--tip-dates", type=Path)
    adapter_beast_calibration_report.add_argument("--alignment", type=Path)
    adapter_beast_calibration_report.add_argument("--date-column", default="date")
    adapter_beast_calibration_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(adapter_beast_calibration_report)
    adapter_bayesian_evidence = adapter_subparsers.add_parser(
        "bayesian-evidence",
        help="Bundle Bayesian configs, trees, logs, diagnostics, and reports into one evidence package.",
    )
    adapter_bayesian_evidence.add_argument("--out-dir", required=True, type=Path)
    adapter_bayesian_evidence.add_argument(
        "--inputs", nargs="+", required=True, type=Path
    )
    adapter_bayesian_evidence.add_argument(
        "--configs", nargs="+", required=True, type=Path
    )
    adapter_bayesian_evidence.add_argument(
        "--trees", nargs="+", required=True, type=Path
    )
    adapter_bayesian_evidence.add_argument(
        "--logs", nargs="+", required=True, type=Path
    )
    adapter_bayesian_evidence.add_argument(
        "--diagnostics", nargs="+", required=True, type=Path
    )
    adapter_bayesian_evidence.add_argument(
        "--reports", nargs="+", required=True, type=Path
    )
    adapter_bayesian_evidence.add_argument(
        "--json", action="store_true", help="Emit the evidence-package report as JSON."
    )
    _add_manifest_argument(adapter_bayesian_evidence)
    adapter_bayesian_table = adapter_subparsers.add_parser(
        "bayesian-diagnostics-table",
        help="Write a supplementary Bayesian diagnostics table from posterior logs.",
    )
    adapter_bayesian_table.add_argument("posterior_trees", type=Path)
    adapter_bayesian_table.add_argument("--log", required=True, type=Path)
    adapter_bayesian_table.add_argument("--additional-logs", nargs="*", type=Path)
    adapter_bayesian_table.add_argument("--out", required=True, type=Path)
    adapter_bayesian_table.add_argument(
        "--burnin-fractions", nargs="+", type=float, default=[0.1, 0.25, 0.5]
    )
    adapter_bayesian_table.add_argument("--ess-threshold", type=float, default=200.0)
    adapter_bayesian_table.add_argument(
        "--mean-shift-threshold", type=float, default=0.5
    )
    adapter_bayesian_table.add_argument(
        "--cross-chain-mean-shift-threshold", type=float, default=0.75
    )
    adapter_bayesian_table.add_argument(
        "--json", action="store_true", help="Emit the diagnostics-table result as JSON."
    )
    _add_manifest_argument(adapter_bayesian_table)
    adapter_bayesian_methods = adapter_subparsers.add_parser(
        "bayesian-methods",
        help="Write reviewer-facing Bayesian methods summary text.",
    )
    adapter_bayesian_methods.add_argument("posterior_trees", type=Path)
    adapter_bayesian_methods.add_argument("--log", required=True, type=Path)
    adapter_bayesian_methods.add_argument("--additional-logs", nargs="*", type=Path)
    adapter_bayesian_methods.add_argument("--analysis-xml", type=Path)
    adapter_bayesian_methods.add_argument("--out", required=True, type=Path)
    adapter_bayesian_methods.add_argument("--tree-prior", default="unspecified")
    adapter_bayesian_methods.add_argument("--clock-model", default="unspecified")
    adapter_bayesian_methods.add_argument("--calibration-path", type=Path)
    adapter_bayesian_methods.add_argument("--tip-dates-path", type=Path)
    adapter_bayesian_methods.add_argument(
        "--burnin-fractions", nargs="+", type=float, default=[0.1, 0.25, 0.5]
    )
    adapter_bayesian_methods.add_argument("--ess-threshold", type=float, default=200.0)
    adapter_bayesian_methods.add_argument(
        "--mean-shift-threshold", type=float, default=0.5
    )
    adapter_bayesian_methods.add_argument(
        "--cross-chain-mean-shift-threshold", type=float, default=0.75
    )
    adapter_bayesian_methods.add_argument(
        "--json", action="store_true", help="Emit the methods-summary result as JSON."
    )
    _add_manifest_argument(adapter_bayesian_methods)

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
                output_paths: list[Path | str] = [args.out] if args.out else []
                if package_result is not None:
                    output_paths.extend(
                        [
                            package_result.report_path,
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
                            "package_output_count": 0 if package_result is None else 10,
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
                        result.node_table_path,
                        result.uncertainty_table_path,
                        result.legend_path,
                        result.model_description_path,
                        result.caption_path,
                        result.manifest_path,
                    ],
                )
                _print_result(
                    build_command_result(
                        command="ancestral",
                        inputs=[args.tree, args.table],
                        outputs=outputs,
                        metrics={
                            "output_dir": str(result.output_dir),
                            "artifact_count": 9,
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
                            "artifact_count": 11,
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
                        result.summary_table_path,
                        result.region_count_table_path,
                        result.node_table_path,
                        result.transition_matrix_path,
                        result.event_table_path,
                        result.map_marker_table_path,
                        result.map_line_table_path,
                        result.exclusion_table_path,
                        result.manifest_path,
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
                            "artifact_count": 12,
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
            )
            outputs = _finalize_outputs(
                args,
                command="diversification",
                inputs=inputs,
                outputs=[result.output_path],
            )
            _print_result(
                build_command_result(
                    command="diversification",
                    inputs=inputs,
                    outputs=outputs,
                    metrics={
                        "report_kind": result.report_kind,
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
            elif args.benchmark_command == "stress-suite":
                report = benchmark_large_dataset_stress_suite(tier=args.tier)
            else:
                report = benchmark_alignment_diagnostics(
                    replicates=args.replicates,
                    sequence_length=args.sequence_length,
                )
            outputs = _finalize_outputs(args, command="benchmark", inputs=[])
            metrics = {
                "observation_count": len(report.observations),
            }
            if hasattr(report, "replicates"):
                metrics["replicates"] = report.replicates
            if hasattr(report, "tier"):
                metrics["tier"] = report.tier
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
            if args.reference_source in {"ape-live", "phytools-live"} and args.extended:
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
            metadata_table = (
                load_taxon_table(args.metadata, taxon_column=args.taxon_column)
                if args.metadata is not None
                else None
            )
            traits_table = (
                load_taxon_table(args.traits, taxon_column=args.taxon_column)
                if args.traits is not None
                else None
            )
            labels: dict[str, str] | None = None
            if metadata_table is not None and args.label_column is not None:
                if args.label_column not in metadata_table.columns:
                    raise MetadataJoinError(
                        f"metadata table does not contain label column '{args.label_column}'"
                    )
                labels = {
                    row[metadata_table.taxon_column]: row[args.label_column]
                    for row in metadata_table.rows
                    if row[args.label_column]
                }
            categorical_traits = (
                _build_string_trait_map(traits_table, args.categorical_column)
                if traits_table is not None and args.categorical_column is not None
                else None
            )
            continuous_traits = (
                _build_numeric_trait_map(traits_table, args.continuous_column)
                if traits_table is not None and args.continuous_column is not None
                else None
            )
            metadata_strips = (
                _build_annotation_strips(
                    metadata_table, _split_csv_values(args.metadata_strip_columns)
                )
                if metadata_table is not None
                else []
            )
            heatmap_columns = (
                _build_annotation_strips(
                    traits_table, _split_csv_values(args.heatmap_columns)
                )
                if traits_table is not None
                else []
            )
            collapsed_clades = _split_csv_values(args.collapse_clades)
            support_audit = (
                audit_support_label_rendering(args.tree)
                if args.support_labels
                else None
            )
            result = render_tree_svg(
                args.tree,
                out_path=args.out,
                labels=labels,
                layout=args.layout,
                show_support_values=args.support_labels
                and (support_audit.validated if support_audit is not None else False),
                categorical_traits=categorical_traits,
                continuous_traits=continuous_traits,
                metadata_strips=metadata_strips,
                heatmap_columns=heatmap_columns,
                collapsed_clades=collapsed_clades,
                validated_support_labels={}
                if support_audit is None
                else support_audit.labels_by_node,
                support_validation_warnings=[]
                if support_audit is None
                else support_audit.warnings,
            )
            inputs = [args.tree]
            if args.metadata is not None:
                inputs.append(args.metadata)
            if args.traits is not None:
                inputs.append(args.traits)
            outputs = [result.output_path]
            package_result = None
            if args.package_dir is not None:
                package_result = build_tree_figure_package(
                    args.tree,
                    out_dir=args.package_dir,
                    labels=labels,
                    layout=args.layout,
                    show_support_values=args.support_labels,
                    categorical_traits=categorical_traits,
                    continuous_traits=continuous_traits,
                    metadata_strips=metadata_strips,
                    heatmap_columns=heatmap_columns,
                    collapsed_clades=collapsed_clades,
                )
                outputs.append(package_result.output_dir)
            outputs = _finalize_outputs(
                args, command="render", inputs=inputs, outputs=outputs
            )
            if args.json:
                _print_result(
                    build_command_result(
                        command="render",
                        inputs=inputs,
                        outputs=outputs,
                        warnings=result.missing_metadata_labels
                        + ([] if support_audit is None else support_audit.warnings),
                        metrics={
                            "tip_count": result.tip_count,
                            "visible_tip_count": result.visible_tip_count,
                            "label_count": result.label_count,
                            "rendered_support_count": result.rendered_support_count,
                            "rendered_categorical_trait_count": result.rendered_categorical_trait_count,
                            "rendered_continuous_trait_count": result.rendered_continuous_trait_count,
                            "rendered_metadata_strip_count": result.rendered_metadata_strip_count,
                            "rendered_heatmap_column_count": result.rendered_heatmap_column_count,
                            "collapsed_clade_count": result.collapsed_clade_count,
                        },
                        data={
                            "render": result,
                            "figure_package_dir": package_result.output_dir
                            if package_result is not None
                            else None,
                            "figure_package_audit": None
                            if package_result is None
                            else package_result.audit,
                            "support_audit": support_audit,
                        },
                    ),
                    json_output=True,
                )
                return 0
            print(result.output_path)
            return 0
        if args.command == "evidence":
            if args.evidence_command == "bundle":
                report = bundle_directory(args.inputs, args.outputs, args.out)
                inputs = [*args.inputs, *args.outputs]
                outputs = _finalize_outputs(
                    args, command="evidence", inputs=inputs, outputs=[args.out]
                )
                _print_result(
                    build_command_result(
                        command="evidence",
                        inputs=inputs,
                        outputs=outputs,
                        metrics={
                            "file_count": report.file_count,
                            "input_file_count": report.input_file_count,
                            "output_file_count": report.output_file_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.evidence_command == "validate":
                report = validate_bundle(args.bundle_root)
                if not report.valid:
                    raise EvidenceContractError(
                        f"evidence bundle validation failed with {len(report.mismatches)} mismatch(es)"
                    )
                outputs = _finalize_outputs(
                    args, command="evidence", inputs=[args.bundle_root]
                )
                _print_result(
                    build_command_result(
                        command="evidence",
                        inputs=[args.bundle_root],
                        outputs=outputs,
                        metrics={"file_count": report.file_count},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            repo_root = Path.cwd()
            if args.evidence_book_command == "studies":
                studies = list_registered_evidence_studies(repo_root)
                outputs = _finalize_outputs(args, command="evidence", inputs=[])
                _print_result(
                    build_command_result(
                        command="evidence",
                        inputs=[],
                        outputs=outputs,
                        metrics={
                            "study_count": len(studies),
                            "partial_rerun_capable_count": sum(
                                1 for study in studies if study.supports_partial_rerun
                            ),
                        },
                        data={"studies": studies},
                    ),
                    json_output=args.json,
                )
                return 0
            if args.evidence_book_command == "build":
                if args.study_id is None and args.evidence_ids:
                    raise EvidenceContractError(
                        "--evidence-id requires a study_id for evidence book build"
                    )
                if args.study_id is None:
                    refresh_report = refresh_evidence_book(repo_root)
                    outputs = _finalize_outputs(
                        args,
                        command="evidence",
                        inputs=[],
                        outputs=refresh_report.updated_paths,
                    )
                    metrics: dict[str, object] = {
                        "reviewer_summary_count": refresh_report.reviewer_summary_count,
                        "updated_path_count": len(refresh_report.updated_paths),
                        **_evidence_book_metrics(repo_root),
                    }
                    _print_result(
                        build_command_result(
                            command="evidence",
                            inputs=[],
                            outputs=outputs,
                            metrics=metrics,
                            data=refresh_report,
                        ),
                        json_output=args.json,
                    )
                    return 0
                if args.evidence_ids:
                    report = build_evidence_book_selection(
                        repo_root,
                        args.study_id,
                        args.evidence_ids,
                    )
                    outputs = _finalize_outputs(
                        args,
                        command="evidence",
                        inputs=[],
                        outputs=report.refresh_report.updated_paths,
                    )
                    metrics: dict[str, object] = {
                        "selected_study_count": 1,
                        "selected_evidence_count": len(report.selected_evidence_ids),
                        "updated_path_count": len(report.refresh_report.updated_paths),
                        "reviewer_summary_count": report.refresh_report.reviewer_summary_count,
                        **_evidence_book_metrics(repo_root),
                    }
                    _print_result(
                        build_command_result(
                            command="evidence",
                            inputs=[],
                            outputs=outputs,
                            metrics=metrics,
                            data=report,
                        ),
                        json_output=args.json,
                    )
                    return 0
                report = build_evidence_book_study(repo_root, args.study_id)
                build_inputs = (
                    []
                    if report.study_report.build_script_path is None
                    else [Path(report.study_report.build_script_path)]
                )
                outputs = _finalize_outputs(
                    args,
                    command="evidence",
                    inputs=build_inputs,
                    outputs=report.refresh_report.updated_paths,
                )
                metrics: dict[str, object] = {
                    "selected_study_count": 1,
                    "updated_path_count": len(report.refresh_report.updated_paths),
                    "reviewer_summary_count": report.refresh_report.reviewer_summary_count,
                    **_evidence_book_metrics(repo_root),
                }
                _print_result(
                    build_command_result(
                        command="evidence",
                        inputs=build_inputs,
                        outputs=outputs,
                        metrics=metrics,
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.evidence_book_command == "validate":
                report = validate_evidence_book(repo_root)
                if not report.valid:
                    raise EvidenceContractError(
                        f"evidence-book validation failed with {len(report.issues)} issue(s)"
                    )
                outputs = _finalize_outputs(
                    args,
                    command="evidence",
                    inputs=[repo_root / "evidence-book"],
                    outputs=[
                        repo_root / "evidence-book" / "index" / "coverage-gaps.json",
                        repo_root / "evidence-book" / "index" / "freshness-report.json",
                        repo_root / "evidence-book" / "index" / "integrity-report.json",
                        repo_root / DOCS_EVIDENCE_OVERVIEW,
                    ],
                )
                _print_result(
                    build_command_result(
                        command="evidence",
                        inputs=[repo_root / "evidence-book"],
                        outputs=outputs,
                        metrics={
                            "issue_count": len(report.issues),
                            **_evidence_book_metrics(repo_root),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            report = rerun_evidence_book_selection(
                repo_root, args.study_id, args.evidence_ids
            )
            outputs = _finalize_outputs(
                args,
                command="evidence",
                inputs=[],
                outputs=report.refresh_report.updated_paths,
            )
            _print_result(
                build_command_result(
                    command="evidence",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "selected_evidence_count": len(
                            report.rerun_report.selected_evidence_ids
                        ),
                        "updated_path_count": len(report.refresh_report.updated_paths),
                        **_evidence_book_metrics(repo_root),
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "demo":
            if args.demo_command == "run":
                result = run_capability_demo(args.out)
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.tree_report,
                        result.dataset_report,
                        result.phylo_inputs_report,
                        result.comparison_report,
                        result.capability_summary,
                    ],
                )
                if args.json:
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={"artifact_count": 5},
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            if args.demo_command == "primate-comparative":
                result = run_primate_comparative_demo(args.out)
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.dataset_export.readme_path,
                        result.dataset_export.tree_path,
                        result.dataset_export.traits_path,
                        result.workflow_bundle.summary_path,
                        result.workflow_bundle.pgls_lambda_profile_path,
                        result.workflow_bundle.brownian_summary_path,
                        result.workflow_bundle.ou_summary_path,
                        result.workflow_bundle.signal_summary_path,
                        result.workflow_bundle.signal_permutations_path,
                        result.workflow_bundle.continuous_ancestral_summary_path,
                        result.workflow_bundle.continuous_ancestral_uncertainty_path,
                        result.workflow_bundle.discrete_ancestral_summary_path,
                        result.workflow_bundle.discrete_ancestral_probability_path,
                        result.overview_path,
                    ],
                )
                if args.json:
                    expected_output_count = len(
                        list(result.dataset_export.expected_output_root.glob("*"))
                    )
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={
                                "artifact_count": len(outputs),
                                "dataset_taxon_count": result.dataset.taxon_count,
                                "reference_output_count": expected_output_count,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            if args.demo_command == "avian-reproductive-traits":
                result = run_avian_reproductive_trait_demo(args.out)
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.dataset_export.readme_path,
                        result.dataset_export.tree_path,
                        result.dataset_export.traits_path,
                        result.workflow_bundle.summary_path,
                        result.workflow_bundle.pgls_lambda_profile_path,
                        result.workflow_bundle.brownian_summary_path,
                        result.workflow_bundle.ou_summary_path,
                        result.workflow_bundle.signal_summary_path,
                        result.workflow_bundle.signal_permutations_path,
                        result.workflow_bundle.continuous_ancestral_summary_path,
                        result.workflow_bundle.continuous_ancestral_uncertainty_path,
                        result.workflow_bundle.discrete_ancestral_summary_path,
                        result.workflow_bundle.discrete_ancestral_probability_path,
                        result.workflow_bundle.clade_summary_path,
                        result.workflow_bundle.clade_rows_path,
                        result.overview_path,
                    ],
                )
                if args.json:
                    expected_output_count = len(
                        list(result.dataset_export.expected_output_root.glob("*"))
                    )
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={
                                "artifact_count": len(outputs),
                                "dataset_taxon_count": result.dataset.taxon_count,
                                "reference_output_count": expected_output_count,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            if args.demo_command == "central-european-seashore-flora":
                result = run_central_european_seashore_flora_demo(args.out)
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.dataset_export.readme_path,
                        result.dataset_export.tree_path,
                        result.dataset_export.traits_path,
                        result.workflow_bundle.summary_path,
                        result.workflow_bundle.pgls_lambda_profile_path,
                        result.workflow_bundle.brownian_summary_path,
                        result.workflow_bundle.ou_summary_path,
                        result.workflow_bundle.signal_summary_path,
                        result.workflow_bundle.signal_permutations_path,
                        result.workflow_bundle.continuous_ancestral_summary_path,
                        result.workflow_bundle.continuous_ancestral_uncertainty_path,
                        result.workflow_bundle.discrete_ancestral_summary_path,
                        result.workflow_bundle.discrete_ancestral_probability_path,
                        result.workflow_bundle.clade_summary_path,
                        result.workflow_bundle.clade_rows_path,
                        result.overview_path,
                    ],
                )
                if args.json:
                    expected_output_count = len(
                        list(result.dataset_export.expected_output_root.glob("*"))
                    )
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={
                                "artifact_count": len(outputs),
                                "dataset_taxon_count": result.dataset.taxon_count,
                                "reference_output_count": expected_output_count,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            if args.demo_command == "influenza-a-ha-reference-panel":
                result = run_influenza_a_ha_reference_demo(
                    args.out,
                    mafft_executable=args.mafft_executable or "mafft",
                    trimal_executable=args.trimal_executable or "trimal",
                    iqtree_executable=args.iqtree_executable or "iqtree2",
                    iqtree_seed=args.iqtree_seed,
                    iqtree_threads=args.iqtree_threads,
                    bootstrap_replicates=args.bootstrap_replicates,
                )
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.dataset_export.readme_path,
                        result.dataset_export.sequences_path,
                        result.workflow_bundle.summary_path,
                        result.workflow_bundle.alignment_path,
                        result.workflow_bundle.trimmed_alignment_path,
                        result.workflow_bundle.tree_path,
                        result.workflow_bundle.model_table_path,
                        result.workflow_bundle.support_table_path,
                        result.workflow_bundle.log_path,
                        result.workflow_bundle.manifest_path,
                        result.overview_path,
                    ],
                )
                if args.json:
                    expected_output_count = len(
                        list(result.dataset_export.expected_output_root.glob("*"))
                    )
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={
                                "artifact_count": len(outputs),
                                "sequence_count": result.dataset.sequence_count,
                                "sequence_type": result.dataset.sequence_type,
                                "selected_model": result.workflow_bundle.selected_model,
                                "minimum_support": result.workflow_bundle.minimum_support,
                                "maximum_support": result.workflow_bundle.maximum_support,
                                "weakly_supported_clade_count": result.workflow_bundle.weakly_supported_clade_count,
                                "reference_output_count": expected_output_count,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            if args.demo_command == "gnathostome-ortholog-protein-benchmark":
                result = run_gnathostome_ortholog_protein_benchmark_demo(
                    args.out,
                    mafft_executable=args.mafft_executable or "mafft",
                    trimal_executable=args.trimal_executable or "trimal",
                    iqtree_executable=args.iqtree_executable or "iqtree2",
                    iqtree_seed=args.iqtree_seed,
                    iqtree_threads=args.iqtree_threads,
                    bootstrap_replicates=args.bootstrap_replicates,
                )
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.dataset_export.readme_path,
                        result.dataset_export.sequences_path,
                        result.workflow_bundle.summary_path,
                        result.workflow_bundle.assumptions_path,
                        result.workflow_bundle.alignment_path,
                        result.workflow_bundle.trimmed_alignment_path,
                        result.workflow_bundle.tree_path,
                        result.workflow_bundle.model_table_path,
                        result.workflow_bundle.support_table_path,
                        result.workflow_bundle.log_path,
                        result.workflow_bundle.manifest_path,
                        result.overview_path,
                    ],
                )
                if args.json:
                    expected_output_count = len(
                        list(result.dataset_export.expected_output_root.glob("*"))
                    )
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={
                                "artifact_count": len(outputs),
                                "sequence_count": result.dataset.sequence_count,
                                "sequence_type": result.dataset.sequence_type,
                                "selected_model": result.workflow_bundle.selected_model,
                                "alignment_length": result.workflow_bundle.alignment_length,
                                "trimmed_alignment_length": result.workflow_bundle.trimmed_alignment_length,
                                "minimum_support": result.workflow_bundle.minimum_support,
                                "maximum_support": result.workflow_bundle.maximum_support,
                                "weakly_supported_clade_count": result.workflow_bundle.weakly_supported_clade_count,
                                "state_space": "amino-acid",
                                "model_selection_scope": "protein-models-only",
                                "reference_output_count": expected_output_count,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            if args.demo_command == "pleistocene-bear-cytb-fragments":
                result = run_pleistocene_bear_cytb_fragment_demo(
                    args.out,
                    mafft_executable=args.mafft_executable or "mafft",
                    trimal_executable=args.trimal_executable or "trimal",
                    iqtree_executable=args.iqtree_executable or "iqtree2",
                    iqtree_seed=args.iqtree_seed,
                    iqtree_threads=args.iqtree_threads,
                    bootstrap_replicates=args.bootstrap_replicates,
                )
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.dataset_export.readme_path,
                        result.dataset_export.sequences_path,
                        result.workflow_bundle.summary_path,
                        result.workflow_bundle.missingness_effects_path,
                        result.workflow_bundle.alignment_path,
                        result.workflow_bundle.trimmed_alignment_path,
                        result.workflow_bundle.cleaned_alignment_path,
                        result.workflow_bundle.tree_path,
                        result.workflow_bundle.model_table_path,
                        result.workflow_bundle.support_table_path,
                        result.overview_path,
                    ],
                )
                if args.json:
                    expected_output_count = len(
                        list(result.dataset_export.expected_output_root.glob("*"))
                    )
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={
                                "artifact_count": len(outputs),
                                "sequence_count": result.dataset.sequence_count,
                                "degraded_sequence_count": len(
                                    result.dataset.degraded_sequence_ids
                                ),
                                "selected_model": result.workflow_bundle.selected_model,
                                "minimum_support": result.workflow_bundle.minimum_support,
                                "maximum_support": result.workflow_bundle.maximum_support,
                                "removed_column_count": result.workflow_bundle.removed_column_count,
                                "cleaned_missing_data_fraction": result.workflow_bundle.cleaned_missing_data_fraction,
                                "reference_output_count": expected_output_count,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            if args.demo_command == "rabies-cross-host-panel":
                result = run_rabies_cross_host_panel_demo(args.out)
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.dataset_export.readme_path,
                        result.dataset_export.sequences_path,
                        result.dataset_export.tree_path,
                        result.dataset_export.hosts_path,
                        result.workflow_bundle.workflow_summary_path,
                        result.workflow_bundle.host_switch_summary_path,
                        result.workflow_bundle.host_state_nodes_path,
                        result.workflow_bundle.host_switch_branches_path,
                        result.workflow_bundle.host_switch_counts_path,
                        result.workflow_bundle.host_switch_fits_path,
                        result.workflow_bundle.host_switch_unsupported_path,
                        result.workflow_bundle.host_switch_exclusions_path,
                        result.overview_path,
                    ],
                )
                if args.json:
                    expected_output_count = len(
                        list(result.dataset_export.expected_output_root.glob("*"))
                    )
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={
                                "artifact_count": len(outputs),
                                "taxon_count": result.dataset.taxon_count,
                                "workflow_trait": result.dataset.workflow_trait,
                                "observed_host_group_count": (
                                    result.dataset.observed_host_group_count
                                ),
                                "analysis_constraint_mode": (
                                    result.workflow_bundle.analysis_constraint_mode
                                ),
                                "root_host": result.workflow_bundle.root_host,
                                "root_confidence": (
                                    result.workflow_bundle.root_confidence
                                ),
                                "host_switch_count": (
                                    result.workflow_bundle.host_switch_count
                                ),
                                "certain_host_switch_count": (
                                    result.workflow_bundle.certain_host_switch_count
                                ),
                                "uncertain_host_switch_count": (
                                    result.workflow_bundle.uncertain_host_switch_count
                                ),
                                "reference_output_count": expected_output_count,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            if args.demo_command == "rabies-geographic-transition-panel":
                result = run_rabies_geographic_transition_panel_demo(args.out)
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.dataset_export.readme_path,
                        result.dataset_export.sequences_path,
                        result.dataset_export.tree_path,
                        result.dataset_export.regions_path,
                        result.workflow_bundle.workflow_summary_path,
                        result.workflow_bundle.geographic_state_summary_path,
                        result.workflow_bundle.geographic_region_probability_path,
                        result.workflow_bundle.geographic_transition_rate_path,
                        result.workflow_bundle.geographic_transition_event_path,
                        result.workflow_bundle.geographic_state_exclusion_path,
                        result.workflow_bundle.geographic_migration_summary_path,
                        result.workflow_bundle.geographic_migration_event_path,
                        result.workflow_bundle.geographic_migration_exclusion_path,
                        result.overview_path,
                    ],
                )
                if args.json:
                    expected_output_count = len(
                        list(result.dataset_export.expected_output_root.glob("*"))
                    )
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={
                                "artifact_count": len(outputs),
                                "taxon_count": result.dataset.taxon_count,
                                "workflow_trait": result.dataset.workflow_trait,
                                "observed_region_group_count": (
                                    result.dataset.observed_region_group_count
                                ),
                                "root_region": result.workflow_bundle.root_region,
                                "root_region_probability": (
                                    result.workflow_bundle.root_region_probability
                                ),
                                "changed_branch_count": (
                                    result.workflow_bundle.changed_branch_count
                                ),
                                "strongly_supported_transition_count": (
                                    result.workflow_bundle.strongly_supported_transition_count
                                ),
                                "migration_event_count": (
                                    result.workflow_bundle.migration_event_count
                                ),
                                "strongly_supported_migration_event_count": (
                                    result.workflow_bundle.strongly_supported_migration_event_count
                                ),
                                "reference_output_count": expected_output_count,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            if args.demo_command == "rabies-cross-host-geography-panel":
                result = run_rabies_cross_host_geography_panel_demo(
                    args.out,
                    config_path=args.config,
                    mafft_executable=args.mafft_executable or "mafft",
                    trimal_executable=args.trimal_executable or "trimal",
                    iqtree_executable=args.iqtree_executable or "iqtree2",
                    fasttree_executable=args.fasttree_executable or "FastTree",
                    iqtree_seed=args.iqtree_seed,
                    iqtree_threads=args.iqtree_threads,
                    bootstrap_replicates=args.bootstrap_replicates,
                )
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.dataset_export.readme_path,
                        result.dataset_export.workflow_config_path,
                        result.dataset_export.sequences_path,
                        result.dataset_export.metadata_path,
                        result.dataset_export.centroids_path,
                        result.dataset_export.accession_table_path,
                        result.workflow_bundle.workflow_summary_path,
                        result.workflow_bundle.resource_observations_path,
                        result.workflow_bundle.config_audit_path,
                        result.workflow_bundle.resolved_config_path,
                        result.workflow_bundle.input_validation_path,
                        result.workflow_bundle.alignment_quality_path,
                        result.workflow_bundle.alignment_sequence_ranking_path,
                        result.workflow_bundle.alignment_path,
                        result.workflow_bundle.trimmed_alignment_path,
                        result.workflow_bundle.tree_path,
                        result.workflow_bundle.rooting_report_path,
                        result.workflow_bundle.model_table_path,
                        result.workflow_bundle.support_table_path,
                        result.workflow_bundle.clade_table_path,
                        result.workflow_bundle.bootstrap_summary_path,
                        result.workflow_bundle.bootstrap_tree_comparison_summary_path,
                        result.workflow_bundle.host_switch_summary_path,
                        result.workflow_bundle.host_switch_counts_path,
                        result.workflow_bundle.biogeography_report_path,
                        result.workflow_bundle.biogeography_tree_figure_path,
                        result.workflow_bundle.biogeography_map_path,
                        result.workflow_bundle.comparative_report_path,
                        result.workflow_bundle.comparative_summary_path,
                        result.workflow_bundle.conclusion_stability_summary_path,
                        result.workflow_bundle.key_clade_stability_path,
                        result.workflow_bundle.support_value_stability_path,
                        result.workflow_bundle.ancestral_state_stability_path,
                        result.workflow_bundle.comparative_coefficient_stability_path,
                        result.workflow_bundle.conclusion_stability_report_path,
                        result.workflow_bundle.scientific_findings_path,
                        result.workflow_bundle.final_report_path,
                        result.workflow_bundle.final_manifest_path,
                        result.overview_path,
                        result.overview_html_path,
                        result.package_manifest_path,
                    ],
                )
                if args.json:
                    expected_output_count = len(
                        [
                            path
                            for path in result.dataset_export.expected_output_root.rglob(
                                "*"
                            )
                            if path.is_file()
                        ]
                    )
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={
                                "artifact_count": len(outputs),
                                "sequence_count": result.dataset.sequence_count,
                                "config_path": str(
                                    result.dataset_export.workflow_config_path
                                ),
                                "biological_question": (
                                    "Do the host-associated rabies lineages in this compact panel occupy one distinct geographic regime while retaining one coherent phylogenetic signal?"
                                ),
                                "short_answer": (
                                    "The rooted panel remains anchored in bat and north_asia, and `host_group[canid]` shows a nominally supported positive longitude association under the selected comparative model, but the inference remains cautionary because the panel is intentionally compact."
                                ),
                                "host_trait": result.dataset.host_trait,
                                "geography_trait": result.dataset.geography_trait,
                                "selected_model": result.workflow_bundle.selected_model,
                                "aligned_quality_score": (
                                    result.workflow_bundle.aligned_quality_score
                                ),
                                "trimmed_quality_score": (
                                    result.workflow_bundle.trimmed_quality_score
                                ),
                                "minimum_support": result.workflow_bundle.minimum_support,
                                "maximum_support": result.workflow_bundle.maximum_support,
                                "root_host": result.workflow_bundle.root_host,
                                "root_region": result.workflow_bundle.root_region,
                                "host_switch_count": result.workflow_bundle.host_switch_count,
                                "migration_event_count": result.workflow_bundle.migration_event_count,
                                "clade_row_count": (
                                    result.workflow_bundle.clade_row_count
                                ),
                                "bootstrap_tree_count": (
                                    result.workflow_bundle.bootstrap_tree_count
                                ),
                                "timeout_seconds": (
                                    result.workflow_bundle.timeout_seconds
                                ),
                                "max_bootstrap_tree_count": (
                                    result.workflow_bundle.max_bootstrap_tree_count
                                ),
                                "max_report_table_rows": (
                                    result.workflow_bundle.max_report_table_rows
                                ),
                                "budget_warning_count": (
                                    result.workflow_bundle.budget_warning_count
                                ),
                                "bootstrap_review_runtime_seconds": (
                                    result.workflow_bundle.bootstrap_review_runtime_seconds
                                ),
                                "bootstrap_review_peak_memory_bytes": (
                                    result.workflow_bundle.bootstrap_review_peak_memory_bytes
                                ),
                                "bootstrap_consensus_rooted_rf_distance": (
                                    result.workflow_bundle.bootstrap_consensus_rooted_rf_distance
                                ),
                                "comparative_formula": (
                                    result.workflow_bundle.comparative_formula
                                ),
                                "comparative_selected_model": (
                                    result.workflow_bundle.comparative_selected_model
                                ),
                                "conclusion_stable_count": (
                                    result.workflow_bundle.conclusion_stable_count
                                ),
                                "conclusion_weak_count": (
                                    result.workflow_bundle.conclusion_weak_count
                                ),
                                "conclusion_unstable_count": (
                                    result.workflow_bundle.conclusion_unstable_count
                                ),
                                "config_check_count": (
                                    result.workflow_bundle.config_check_count
                                ),
                                "scientific_finding_count": (
                                    result.workflow_bundle.scientific_finding_count
                                ),
                                "reference_output_count": expected_output_count,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            if args.demo_command == "rabies-method-sensitivity-panel":
                result = run_rabies_method_sensitivity_panel_demo(
                    args.out,
                    mafft_executable=args.mafft_executable or "mafft",
                    trimal_executable=args.trimal_executable or "trimal",
                    iqtree_executable=args.iqtree_executable or "iqtree2",
                    fasttree_executable=args.fasttree_executable or "FastTree",
                    iqtree_seed=args.iqtree_seed,
                    iqtree_threads=args.iqtree_threads,
                    bootstrap_replicates=args.bootstrap_replicates,
                    parallel_workers=args.parallel_workers,
                )
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.dataset_export.readme_path,
                        result.dataset_export.config_path,
                        result.dataset_export.sequences_path,
                        result.dataset_export.metadata_path,
                        result.workflow_bundle.workflow_summary_path,
                        result.workflow_bundle.variant_summary_path,
                        result.workflow_bundle.parallel_summary_path,
                        result.workflow_bundle.preprocessing_comparison_path,
                        result.workflow_bundle.stable_clades_path,
                        result.workflow_bundle.changed_clades_path,
                        result.workflow_bundle.conclusion_summary_path,
                        result.workflow_bundle.config_path,
                        result.workflow_bundle.manifest_path,
                        result.workflow_bundle.report_manifest_path,
                        result.workflow_bundle.report_path,
                        result.overview_path,
                    ],
                )
                if args.json:
                    expected_output_count = len(
                        [
                            path
                            for path in result.dataset_export.expected_output_root.rglob(
                                "*"
                            )
                            if path.is_file()
                        ]
                    )
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={
                                "artifact_count": len(outputs),
                                "taxon_count": result.dataset.taxon_count,
                                "variant_count": result.workflow_bundle.variant_count,
                                "parallel_workers": (
                                    result.workflow_bundle.parallel_workers
                                ),
                                "execution_mode": (
                                    result.workflow_bundle.execution_mode
                                ),
                                "stable_clade_count": (
                                    result.workflow_bundle.stable_clade_count
                                ),
                                "changed_clade_count": (
                                    result.workflow_bundle.changed_clade_count
                                ),
                                "preprocessing_change_pair_count": (
                                    result.workflow_bundle.preprocessing_change_pair_count
                                ),
                                "rooted_engine_change_variant_count": (
                                    result.workflow_bundle.rooted_engine_change_variant_count
                                ),
                                "serious_conflict_variant_count": (
                                    result.workflow_bundle.serious_conflict_variant_count
                                ),
                                "report_linked_artifact_count": (
                                    result.workflow_bundle.report_linked_artifact_count
                                ),
                                "report_html_size_bytes": (
                                    result.workflow_bundle.report_html_size_bytes
                                ),
                                "report_linked_artifact_bytes": (
                                    result.workflow_bundle.report_linked_artifact_bytes
                                ),
                                "report_total_output_bytes": (
                                    result.workflow_bundle.report_total_output_bytes
                                ),
                                "reference_output_count": expected_output_count,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            if args.demo_command == "catarrhine-mitogenome-five-locus-panel":
                result = run_catarrhine_mitogenome_five_locus_panel_demo(
                    args.out,
                    iqtree_executable=args.iqtree_executable or "iqtree2",
                    iqtree_seed=args.iqtree_seed,
                    iqtree_threads=args.iqtree_threads,
                    bootstrap_replicates=args.bootstrap_replicates,
                )
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.dataset_export.readme_path,
                        result.dataset_export.taxa_path,
                        *sorted(
                            result.dataset_export.locus_alignment_root.glob("*.fasta")
                        ),
                        result.workflow_bundle.workflow_summary_path,
                        result.workflow_bundle.supermatrix_path,
                        result.workflow_bundle.partitions_path,
                        result.workflow_bundle.occupancy_taxa_path,
                        result.workflow_bundle.occupancy_loci_path,
                        result.workflow_bundle.occupancy_matrix_path,
                        result.workflow_bundle.partition_summary_path,
                        result.workflow_bundle.model_candidates_path,
                        result.workflow_bundle.support_tree_path,
                        result.workflow_bundle.support_table_path,
                        result.overview_path,
                    ],
                )
                if args.json:
                    expected_output_count = len(
                        list(result.dataset_export.expected_output_root.glob("*"))
                    )
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={
                                "artifact_count": len(outputs),
                                "taxon_count": result.dataset.taxon_count,
                                "locus_count": result.dataset.locus_count,
                                "alignment_length": (
                                    result.workflow_bundle.alignment_length
                                ),
                                "partition_count": (
                                    result.workflow_bundle.partition_count
                                ),
                                "selected_model": result.workflow_bundle.selected_model,
                                "minimum_support": (
                                    result.workflow_bundle.minimum_support
                                ),
                                "maximum_support": (
                                    result.workflow_bundle.maximum_support
                                ),
                                "weakly_supported_clade_count": (
                                    result.workflow_bundle.weakly_supported_clade_count
                                ),
                                "reference_output_count": expected_output_count,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            if args.demo_command == "catarrhine-data-quality-stress-panel":
                result = run_catarrhine_data_quality_stress_panel_demo(args.out)
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.dataset_export.readme_path,
                        result.dataset_export.raw_alignment_path,
                        result.dataset_export.raw_sequence_input_path,
                        result.dataset_export.raw_coding_sequences_path,
                        result.dataset_export.raw_tree_path,
                        result.dataset_export.raw_traits_path,
                        result.dataset_export.raw_trait_mismatch_path,
                        result.workflow_bundle.workflow_summary_path,
                        result.workflow_bundle.raw_sequence_findings_path,
                        result.workflow_bundle.raw_sequence_repair_path,
                        result.workflow_bundle.repaired_sequence_input_path,
                        result.workflow_bundle.repaired_sequence_validation_path,
                        result.workflow_bundle.coding_sequence_exclusions_path,
                        result.workflow_bundle.prepared_coding_sequences_path,
                        result.workflow_bundle.raw_trait_linkage_path,
                        result.workflow_bundle.trait_duplicates_path,
                        result.workflow_bundle.trait_missing_values_path,
                        result.workflow_bundle.sequence_outliers_path,
                        result.workflow_bundle.tree_issues_path,
                        result.workflow_bundle.repair_actions_path,
                        result.workflow_bundle.cleaned_traits_path,
                        result.workflow_bundle.cleaned_alignment_path,
                        result.workflow_bundle.cleaned_tree_path,
                        result.workflow_bundle.cleaned_linkage_path,
                        result.workflow_bundle.cleaned_validation_path,
                        result.overview_path,
                    ],
                )
                if args.json:
                    expected_output_count = len(
                        list(result.dataset_export.expected_output_root.glob("*"))
                    )
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={
                                "artifact_count": len(outputs),
                                "raw_taxon_count": result.workflow_bundle.raw_taxon_count,
                                "cleaned_taxon_count": result.workflow_bundle.cleaned_taxon_count,
                                "duplicate_sequence_identifier_count": (
                                    result.workflow_bundle.duplicate_sequence_identifier_count
                                ),
                                "illegal_character_count": (
                                    result.workflow_bundle.illegal_character_count
                                ),
                                "empty_sequence_count": (
                                    result.workflow_bundle.empty_sequence_count
                                ),
                                "raw_sequence_length_outlier_count": (
                                    result.workflow_bundle.raw_sequence_length_outlier_count
                                ),
                                "duplicate_trait_taxon_count": (
                                    result.workflow_bundle.duplicate_trait_taxon_count
                                ),
                                "missing_trait_value_count": (
                                    result.workflow_bundle.missing_trait_value_count
                                ),
                                "sequence_outlier_count": (
                                    result.workflow_bundle.sequence_outlier_count
                                ),
                                "tree_zero_length_branch_count": (
                                    result.workflow_bundle.tree_zero_length_branch_count
                                ),
                                "tree_negative_branch_count": (
                                    result.workflow_bundle.tree_negative_branch_count
                                ),
                                "tree_long_branch_outlier_count": (
                                    result.workflow_bundle.tree_long_branch_outlier_count
                                ),
                                "coding_frame_error_count": (
                                    result.workflow_bundle.coding_frame_error_count
                                ),
                                "coding_internal_stop_count": (
                                    result.workflow_bundle.coding_internal_stop_count
                                ),
                                "raw_trait_missing_from_traits_count": (
                                    result.workflow_bundle.raw_trait_missing_from_traits_count
                                ),
                                "raw_trait_extra_taxon_count": (
                                    result.workflow_bundle.raw_trait_extra_taxon_count
                                ),
                                "dropped_taxon_count": (
                                    result.workflow_bundle.dropped_taxon_count
                                ),
                                "repaired_branch_count": (
                                    result.workflow_bundle.repaired_branch_count
                                ),
                                "reference_output_count": expected_output_count,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            if args.demo_command == "continuous-mode-recovery-panel":
                result = run_continuous_mode_recovery_panel_demo(args.out)
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.dataset_export.readme_path,
                        result.dataset_export.reference_tree_path,
                        result.dataset_export.simulation_cases_path,
                        result.workflow_bundle.workflow_summary_path,
                        result.workflow_bundle.recovery_summary_path,
                        result.workflow_bundle.parameter_recovery_path,
                        result.workflow_bundle.model_choice_path,
                        result.workflow_bundle.warning_review_path,
                        result.overview_path,
                    ],
                )
                if args.json:
                    expected_output_count = len(
                        [
                            path
                            for path in result.dataset_export.expected_output_root.rglob(
                                "*"
                            )
                            if path.is_file()
                        ]
                    )
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={
                                "artifact_count": len(outputs),
                                "taxon_count": result.dataset.taxon_count,
                                "case_count": result.dataset.case_count,
                                "selection_match_count": (
                                    result.workflow_bundle.selection_match_count
                                ),
                                "parameter_pass_count": (
                                    result.workflow_bundle.parameter_pass_count
                                ),
                                "parameter_row_count": (
                                    result.workflow_bundle.parameter_row_count
                                ),
                                "expected_warning_case_count": (
                                    result.workflow_bundle.expected_warning_case_count
                                ),
                                "expected_warning_present_count": (
                                    result.workflow_bundle.expected_warning_present_count
                                ),
                                "reference_output_count": expected_output_count,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            if args.demo_command == "known-answer-reference-panel":
                result = run_known_answer_reference_demo(args.out)
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.dataset_export.readme_path,
                        result.dataset_export.true_tree_path,
                        result.dataset_export.alignment_path,
                        result.dataset_export.continuous_traits_path,
                        result.dataset_export.ou_traits_path,
                        result.dataset_export.discrete_traits_path,
                        result.dataset_export.host_traits_path,
                        result.dataset_export.geographic_traits_path,
                        result.dataset_export.true_parameters_path,
                        result.dataset_export.true_continuous_nodes_path,
                        result.dataset_export.true_ou_nodes_path,
                        result.dataset_export.true_discrete_nodes_path,
                        result.dataset_export.true_host_nodes_path,
                        result.dataset_export.true_geographic_nodes_path,
                        result.dataset_export.true_host_switch_events_path,
                        result.dataset_export.true_geographic_transition_events_path,
                        result.dataset_export.recovery_thresholds_path,
                        result.workflow_bundle.workflow_summary_path,
                        result.workflow_bundle.distance_tree_path,
                        result.workflow_bundle.tree_recovery_path,
                        result.workflow_bundle.parameter_recovery_path,
                        result.workflow_bundle.brownian_fit_summary_path,
                        result.workflow_bundle.ou_fit_summary_path,
                        result.workflow_bundle.continuous_ancestral_summary_path,
                        result.workflow_bundle.continuous_ancestral_uncertainty_path,
                        result.workflow_bundle.continuous_node_recovery_path,
                        result.workflow_bundle.discrete_ancestral_summary_path,
                        result.workflow_bundle.discrete_ancestral_probability_path,
                        result.workflow_bundle.discrete_node_recovery_path,
                        result.workflow_bundle.host_switch_summary_path,
                        result.workflow_bundle.host_state_nodes_path,
                        result.workflow_bundle.host_switch_branches_path,
                        result.workflow_bundle.host_node_recovery_path,
                        result.workflow_bundle.host_event_recovery_path,
                        result.workflow_bundle.geographic_ancestral_summary_path,
                        result.workflow_bundle.geographic_state_probability_path,
                        result.workflow_bundle.geographic_transition_summary_path,
                        result.workflow_bundle.geographic_node_recovery_path,
                        result.workflow_bundle.geographic_event_recovery_path,
                        result.workflow_bundle.threshold_evaluation_path,
                        result.overview_path,
                    ],
                )
                if args.json:
                    expected_output_count = len(
                        list(result.dataset_export.expected_output_root.glob("*"))
                    )
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={
                                "artifact_count": len(outputs),
                                "taxon_count": result.dataset.taxon_count,
                                "sequence_length": result.dataset.sequence_length,
                                "distance_method": result.dataset.distance_method,
                                "distance_model": result.dataset.distance_model,
                                "rooted_topology_equal": (
                                    result.workflow_bundle.rooted_topology_equal
                                ),
                                "same_unrooted_topology": (
                                    result.workflow_bundle.same_unrooted_topology
                                ),
                                "same_taxa_different_rooting": (
                                    result.workflow_bundle.same_taxa_different_rooting
                                ),
                                "robinson_foulds_distance": (
                                    result.workflow_bundle.robinson_foulds_distance
                                ),
                                "parameter_row_count": (
                                    result.workflow_bundle.parameter_row_count
                                ),
                                "threshold_pass_count": (
                                    result.workflow_bundle.threshold_pass_count
                                ),
                                "threshold_row_count": (
                                    result.workflow_bundle.threshold_row_count
                                ),
                                "continuous_internal_node_mean_absolute_error": (
                                    result.workflow_bundle.continuous_internal_node_mean_absolute_error
                                ),
                                "discrete_internal_node_accuracy": (
                                    result.workflow_bundle.discrete_internal_node_accuracy
                                ),
                                "host_internal_node_accuracy": (
                                    result.workflow_bundle.host_internal_node_accuracy
                                ),
                                "host_event_accuracy": (
                                    result.workflow_bundle.host_event_accuracy
                                ),
                                "geographic_internal_node_accuracy": (
                                    result.workflow_bundle.geographic_internal_node_accuracy
                                ),
                                "geographic_event_accuracy": (
                                    result.workflow_bundle.geographic_event_accuracy
                                ),
                                "reference_output_count": expected_output_count,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            raise NotImplementedError(f"unsupported demo command: {args.demo_command}")
        if args.command == "report":
            if args.report_command == "tree":
                result = render_tree_report(tree_path=args.tree, out_path=args.out)
                outputs = _finalize_outputs(
                    args,
                    command="report",
                    inputs=[args.tree],
                    outputs=[result.output_path, result.machine_manifest_path],
                )
                if args.json:
                    _print_result(
                        build_command_result(
                            command="report",
                            inputs=[args.tree],
                            outputs=outputs,
                            warnings=result.validation.warnings
                            + result.inspection.warnings,
                            metrics={"tip_count": result.inspection.tip_count},
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_path)
                return 0
            if args.report_command == "tree-package":
                result = build_tree_report_package(args.tree, out_dir=args.out_dir)
                outputs = _finalize_outputs(
                    args,
                    command="report",
                    inputs=[args.tree],
                    outputs=[
                        result.report_path,
                        result.figure_path,
                        result.support_table_path,
                        result.clade_table_path,
                        result.branch_stats_path,
                        result.manifest_path,
                    ],
                )
                if args.json:
                    _print_result(
                        build_command_result(
                            command="report",
                            inputs=[args.tree],
                            outputs=outputs,
                            warnings=result.validation.warnings
                            + result.inspection.warnings,
                            metrics={
                                "tip_count": result.inspection.tip_count,
                                "supported_branch_count": sum(
                                    1
                                    for row in result.support_rows
                                    if row.support is not None
                                ),
                                "rendered_support_count": (
                                    result.figure.rendered_support_count
                                ),
                                "long_outlier_count": (
                                    result.branch_stats.long_outlier_count
                                ),
                                **method_tier_metrics(result.method_tier),
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_dir)
                return 0
            if args.report_command == "alignment":
                result = render_alignment_report(
                    alignment_path=args.alignment, out_path=args.out
                )
                outputs = _finalize_outputs(
                    args,
                    command="report",
                    inputs=[args.alignment],
                    outputs=[result.output_path, result.machine_manifest_path],
                )
                if args.json:
                    _print_result(
                        build_command_result(
                            command="report",
                            inputs=[args.alignment],
                            outputs=outputs,
                            warnings=result.alignment_forensic.warnings,
                            metrics={
                                "sequence_count": result.alignment.sequence_count,
                                "alignment_length": result.alignment.alignment_length,
                                "quality_score": result.alignment_quality.quality_score,
                                "warning_count": len(
                                    result.alignment_forensic.warnings
                                ),
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_path)
                return 0
            if args.report_command == "dataset":
                result = render_dataset_report(
                    tree_path=args.tree,
                    metadata_path=args.metadata,
                    traits_path=args.traits,
                    alignment_path=args.alignment,
                    tip_dates_path=args.tip_dates,
                    calibration_path=args.calibrations,
                    out_path=args.out,
                )
                inputs = [args.tree, args.metadata]
                if args.traits is not None:
                    inputs.append(args.traits)
                if args.alignment is not None:
                    inputs.append(args.alignment)
                if args.tip_dates is not None:
                    inputs.append(args.tip_dates)
                if args.calibrations is not None:
                    inputs.append(args.calibrations)
                outputs = _finalize_outputs(
                    args,
                    command="report",
                    inputs=inputs,
                    outputs=[result.output_path, result.machine_manifest_path],
                )
                if args.json:
                    _print_result(
                        build_command_result(
                            command="report",
                            inputs=inputs,
                            outputs=outputs,
                            warnings=result.validation.warnings
                            + result.inspection.warnings,
                            metrics={
                                "tip_count": result.inspection.tip_count,
                                "linked_taxa": result.metadata_linkage.linked_taxa,
                                "readiness_decision": None
                                if result.dataset_audit is None
                                else result.dataset_audit.readiness_decision,
                                "excluded_taxa": 0
                                if result.dataset_audit is None
                                else len(result.dataset_audit.exclusion_table.rows),
                                "blocked_analysis_count": 0
                                if result.dataset_audit is None
                                else len(result.dataset_audit.blocked_analyses),
                                "risky_analysis_count": 0
                                if result.dataset_audit is None
                                else sum(
                                    1
                                    for row in result.dataset_audit.analysis_decisions
                                    if row.decision == "risky"
                                ),
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_path)
                return 0
            if args.report_command == "phylo-inputs":
                result = render_phylo_inputs_report(
                    tree_path=args.tree,
                    alignment_path=args.alignment,
                    out_path=args.out,
                )
                inputs = [args.tree, args.alignment]
                outputs = _finalize_outputs(
                    args,
                    command="report",
                    inputs=inputs,
                    outputs=[result.output_path, result.machine_manifest_path],
                )
                if args.json:
                    _print_result(
                        build_command_result(
                            command="report",
                            inputs=inputs,
                            outputs=outputs,
                            warnings=result.validation.warnings
                            + result.inspection.warnings,
                            metrics={
                                "tip_count": result.inspection.tip_count,
                                "alignment_length": result.alignment.alignment_length,
                                "linked_taxa": result.alignment_linkage.linked_taxa,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_path)
                return 0
            if args.report_command == "taxonomy":
                result = render_taxon_report(
                    tree_path=args.tree,
                    synonym_table_path=args.synonym_table,
                    metadata_path=args.metadata,
                    traits_path=args.traits,
                    alignment_path=args.alignment,
                    filtered_alignment_path=args.filtered_alignment,
                    inference_tree_path=args.inference_tree,
                    reported_taxa_path=args.reported_taxa,
                    out_path=args.out,
                )
                inputs = [
                    args.tree,
                    *([args.synonym_table] if args.synonym_table is not None else []),
                ]
                if args.metadata is not None:
                    inputs.append(args.metadata)
                if args.traits is not None:
                    inputs.append(args.traits)
                if args.alignment is not None:
                    inputs.append(args.alignment)
                if args.filtered_alignment is not None:
                    inputs.append(args.filtered_alignment)
                if args.inference_tree is not None:
                    inputs.append(args.inference_tree)
                if args.reported_taxa is not None:
                    inputs.append(args.reported_taxa)
                outputs = _finalize_outputs(
                    args,
                    command="report",
                    inputs=inputs,
                    outputs=[result.output_path, result.machine_manifest_path],
                )
                if args.json:
                    _print_result(
                        build_command_result(
                            command="report",
                            inputs=inputs,
                            outputs=outputs,
                            warnings=result.taxon_audit.warnings,
                            metrics={
                                "tree_tip_count": result.taxon_audit.tree_tip_count,
                                "status": result.taxon_audit.status,
                                "conflict_count": len(
                                    result.taxon_audit.mapping_conflicts.rows
                                ),
                                "crosswalk_rows": 0
                                if result.taxon_crosswalk is None
                                else len(result.taxon_crosswalk.rows),
                                "excluded_taxa": 0
                                if result.taxon_exclusions is None
                                else len(result.taxon_exclusions.rows),
                                "loss_stage_count": 0
                                if result.taxon_workflow_loss is None
                                else len(result.taxon_workflow_loss.loss_stage_counts),
                                "unstable_taxa": 0
                                if result.taxon_stability is None
                                else len(result.taxon_stability.unstable_taxa),
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_path)
                return 0
            if args.report_command == "workflow-validation":
                result = render_workflow_validation_report(
                    out_path=args.out,
                    fixtures_root=args.fixtures_root,
                )
                inputs = [] if args.fixtures_root is None else [args.fixtures_root]
                outputs = _finalize_outputs(
                    args,
                    command="report",
                    inputs=inputs,
                    outputs=[result.output_path, result.machine_manifest_path],
                )
                if args.json:
                    _print_result(
                        build_command_result(
                            command="report",
                            inputs=inputs,
                            outputs=outputs,
                            metrics={
                                "total_fixture_count": result.validation.total_fixture_count,
                                "passed_fixture_count": result.validation.passed_fixture_count,
                                "workflow_count": len(result.validation.workflows),
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_path)
                return 0
            if args.report_command == "release-gate":
                result = render_level_one_release_gate_report(
                    out_path=args.out,
                    fixtures_root=args.fixtures_root,
                )
                inputs = [] if args.fixtures_root is None else [args.fixtures_root]
                outputs = _finalize_outputs(
                    args,
                    command="report",
                    inputs=inputs,
                    outputs=[result.output_path, result.machine_manifest_path],
                )
                if args.json:
                    _print_result(
                        build_command_result(
                            command="report",
                            inputs=inputs,
                            outputs=outputs,
                            warnings=result.release_gate.dataset_warnings,
                            metrics={
                                "decision": result.release_gate.gate.decision,
                                "retained_taxa": len(
                                    result.release_gate.gate.retained_taxa
                                ),
                                "excluded_taxa": len(
                                    result.release_gate.gate.excluded_taxa
                                ),
                                "blocked_analysis_count": len(
                                    result.release_gate.gate.blocked_analyses
                                ),
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_path)
                return 0
            if args.report_command == "release-truth":
                result = render_release_truth_report(
                    out_path=args.out,
                    test_report_paths=args.test_report,
                    real_engine_test_report_paths=args.real_engine_test_report,
                    fixtures_root=args.fixtures_root,
                    include_extended_parity=args.parity_extended,
                    stress_tier=args.stress_tier,
                )
                inputs = [
                    *args.test_report,
                    *args.real_engine_test_report,
                    *([args.fixtures_root] if args.fixtures_root is not None else []),
                ]
                outputs = _finalize_outputs(
                    args,
                    command="report",
                    inputs=inputs,
                    outputs=[result.output_path, result.machine_manifest_path],
                )
                if args.json:
                    _print_result(
                        build_command_result(
                            command="report",
                            inputs=inputs,
                            outputs=outputs,
                            warnings=result.release_truth.known_limitations,
                            metrics={
                                "total_tests": result.release_truth.total_tests.total_tests,
                                "total_tests_passed": result.release_truth.total_tests.passed_tests,
                                "total_tests_failed": result.release_truth.total_tests.failed_tests,
                                "total_tests_skipped": result.release_truth.total_tests.skipped_tests,
                                "real_engine_tests": result.release_truth.real_engine_tests.total_tests,
                                "real_engine_tests_passed": result.release_truth.real_engine_tests.passed_tests,
                                "real_engine_tests_failed": result.release_truth.real_engine_tests.failed_tests,
                                "real_engine_tests_skipped": result.release_truth.real_engine_tests.skipped_tests,
                                "supported_workflow_count": len(
                                    result.release_truth.supported_workflows
                                ),
                                "experimental_workflow_count": len(
                                    result.release_truth.experimental_workflows
                                ),
                                "flagship_dataset_count": len(
                                    result.release_truth.flagship_datasets
                                ),
                                "reference_parity_case_count": result.release_truth.reference_parity.case_count,
                                "stress_workload_count": len(
                                    result.release_truth.stress_suite.observations
                                ),
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_path)
                return 0
            raise NotImplementedError(
                f"unsupported report command: {args.report_command}"
            )
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
            if args.adapter_command == "infer-fast":
                report = run_fast_tree_inference(
                    args.input_path,
                    args.out,
                    executable=args.executable or "FastTree",
                    sequence_type=args.sequence_type,
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
                            "warning_count": len(report.run.warning_lines),
                            "approximate_method": (
                                None
                                if report.fasttree_support_summary is None
                                else report.fasttree_support_summary.approximate_method
                            ),
                            "support_label_kind": (
                                None
                                if report.fasttree_support_summary is None
                                else report.fasttree_support_summary.support_label_kind
                            ),
                            "support_scale": (
                                None
                                if report.fasttree_support_summary is None
                                else report.fasttree_support_summary.support_scale
                            ),
                            "annotated_node_count": (
                                0
                                if report.fasttree_support_summary is None
                                else report.fasttree_support_summary.annotated_node_count
                            ),
                            "minimum_local_support": (
                                None
                                if report.fasttree_support_summary is None
                                else report.fasttree_support_summary.minimum_local_support
                            ),
                            "maximum_local_support": (
                                None
                                if report.fasttree_support_summary is None
                                else report.fasttree_support_summary.maximum_local_support
                            ),
                            "weakly_supported_clade_count": (
                                0
                                if report.fasttree_support_summary is None
                                else report.fasttree_support_summary.weakly_supported_clade_count
                            ),
                            "resumed": report.resumed,
                            "timeout_seconds": report.run.timeout_seconds,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-prepare":
                report = prepare_mrbayes_analysis(
                    args.input_path,
                    args.out,
                    partition_path=args.partitions,
                    model=args.model,
                    rates=args.rates,
                    ngen=args.ngen,
                    nchains=args.nchains,
                    samplefreq=args.samplefreq,
                    printfreq=args.printfreq,
                    burnin_fraction=args.burnin_fraction,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[args.out],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "character_count": report.character_count,
                            "partitioned": report.partition_path is not None,
                            "partition_count": report.partition_count,
                            "partition_warning_count": len(report.partition_warnings),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-run":
                report = run_mrbayes_posterior_inference(
                    args.input_path,
                    executable=args.executable or "mb",
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
                            "warning_count": len(report.run.warning_lines),
                            "resumed": report.resumed,
                            "timeout_seconds": report.run.timeout_seconds,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-summarize":
                consensus_tree, report = summarize_mrbayes_posterior_trees(
                    args.input_path,
                    burnin_fraction=args.burnin_fraction,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[report.filtered_tree_set_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "kept_tree_count": report.kept_tree_count,
                            "rooted_topology_count": report.rooted_topology_count,
                            "tip_count": consensus_tree.tip_count,
                        },
                        data={
                            "summary": report,
                            "consensus_newick": report.consensus_newick,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-traces":
                report = parse_mrbayes_parameter_traces(args.input_path)
                outputs = _finalize_outputs(
                    args, command="adapter", inputs=[args.input_path]
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "row_count": report.row_count,
                            "column_count": len(report.columns),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-trees":
                report = parse_mrbayes_posterior_tree_samples(args.input_path)
                outputs = _finalize_outputs(
                    args, command="adapter", inputs=[args.input_path]
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "tree_count": report.tree_count,
                            "rooted_tree_count": report.rooted_tree_count,
                            "sampled_generation_count": len(report.sampled_generations),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-subsample":
                report = subsample_mrbayes_posterior_tree_set(
                    args.input_path,
                    method=args.method,
                    thinning_interval=args.thinning_interval,
                    sample_count=args.sample_count,
                    burnin_fraction=args.burnin_fraction,
                    random_seed=args.seed,
                )
                outputs: list[Path | str] = []
                if args.tree_set_out is not None:
                    outputs.append(
                        write_posterior_tree_subsample(args.tree_set_out, report)
                    )
                if args.sample_table_out is not None:
                    outputs.append(
                        write_posterior_tree_subsample_table(
                            args.sample_table_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "total_tree_count": report.total_tree_count,
                            "burnin_tree_count": report.burnin_tree_count,
                            "pre_subsampling_tree_count": (
                                report.pre_subsampling_tree_count
                            ),
                            "retained_tree_count": report.retained_tree_count,
                            "selection_method": report.selection_method,
                            "retained_generation_count": len(
                                [
                                    tree
                                    for tree in report.trees
                                    if tree.generation is not None
                                ]
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-mcmc":
                report = parse_mrbayes_mcmc_diagnostics(args.input_path)
                outputs = _finalize_outputs(
                    args, command="adapter", inputs=[args.input_path]
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "row_count": report.row_count,
                            "column_count": len(report.columns),
                            "comment_count": len(report.comment_lines),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-consensus":
                tree, report = parse_mrbayes_consensus_tree(args.input_path)
                outputs = _finalize_outputs(
                    args, command="adapter", inputs=[args.input_path]
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "tip_count": tree.tip_count,
                            "annotated_node_count": report.annotated_node_count,
                            "maximum_posterior_probability": (
                                report.maximum_posterior_probability
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-ess":
                report = compute_mrbayes_effective_sample_sizes(args.input_path)
                outputs = _finalize_outputs(
                    args, command="adapter", inputs=[args.input_path]
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={"parameter_count": len(report.effective_sample_sizes)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-parameters":
                report = summarize_mrbayes_parameter_diagnostics(
                    args.input_path,
                    burnin_fraction=args.burnin_fraction,
                )
                outputs: list[Path | str] = []
                if args.summary_out is not None:
                    outputs.append(
                        write_mrbayes_parameter_summary_table(
                            args.summary_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "burnin_fraction": report.burnin_fraction,
                            "kept_row_count": report.kept_row_count,
                            "parameter_count": len(report.parameter_summaries),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-burnin-sensitivity":
                report = assess_mrbayes_burnin_sensitivity(
                    args.posterior_trees,
                    trace_path=args.traces,
                    burnin_fractions=tuple(args.burnin_fractions),
                )
                inputs = [
                    args.posterior_trees,
                    *([args.traces] if args.traces is not None else []),
                ]
                outputs: list[Path | str] = []
                if args.slice_out is not None:
                    outputs.append(
                        write_mrbayes_burnin_sensitivity_slice_table(
                            args.slice_out,
                            report,
                        )
                    )
                if args.parameter_out is not None:
                    outputs.append(
                        write_burnin_parameter_shift_table(
                            args.parameter_out,
                            report.parameter_shifts,
                        )
                    )
                if args.clade_out is not None:
                    outputs.append(
                        write_burnin_clade_shift_table(
                            args.clade_out,
                            report.clade_shifts,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=inputs,
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=inputs,
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "slice_count": len(report.slices),
                            "parameter_shift_count": len(report.parameter_shifts),
                            "unstable_parameter_count": report.unstable_parameter_count,
                            "clade_shift_count": len(report.clade_shifts),
                            "unstable_clade_count": report.unstable_clade_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-convergence":
                report = assess_mrbayes_convergence(
                    args.input_path,
                    ess_threshold=args.ess_threshold,
                    mean_shift_threshold=args.mean_shift_threshold,
                )
                outputs = _finalize_outputs(
                    args, command="adapter", inputs=[args.input_path]
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=[warning["message"] for warning in report.warnings],
                        metrics={
                            "warning_count": len(report.warnings),
                            "converged": report.converged,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "mrbayes-report":
                report = render_bayesian_posterior_report(
                    posterior_tree_path=args.posterior_trees,
                    trace_path=args.traces,
                    out_path=args.out,
                    burnin_fraction=args.burnin_fraction,
                    ess_threshold=args.ess_threshold,
                    mean_shift_threshold=args.mean_shift_threshold,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.posterior_trees, args.traces],
                    outputs=[report.output_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.posterior_trees, args.traces],
                        outputs=outputs,
                        warnings=method_tier_warnings(report.method_tier),
                        metrics={
                            "kept_tree_count": report.kept_tree_count,
                            "warning_count": report.warning_count
                            + len(method_tier_warnings(report.method_tier)),
                            **method_tier_metrics(report.method_tier),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-prepare":
                report = prepare_beast_time_tree_analysis(
                    args.input_path,
                    args.out,
                    tree_path=args.tree,
                    calibration_path=args.calibrations,
                    tip_dates_path=args.tip_dates,
                    clock_model=args.clock_model,
                    tree_prior=args.tree_prior,
                    chain_length=args.chain_length,
                    log_every=args.log_every,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[args.out],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "taxon_count": report.taxon_count,
                            "character_count": report.character_count,
                            "calibration_count": report.calibration_count,
                            "tip_date_count": report.tip_date_count,
                            "warning_count": report.warning_count,
                            "starting_tree_source": report.starting_tree_source,
                            "beast_data_type": report.beast_data_type,
                            "substitution_model": report.substitution_model,
                            "clock_model": report.clock_model,
                            "tree_prior": report.tree_prior,
                            "chain_length": report.chain_length,
                            "log_every": report.log_every,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-xml":
                report = summarize_beast_analysis_xml(args.input_path)
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "valid": report.valid,
                            "issue_count": len(report.issues),
                            "taxon_count": report.taxon_count,
                            "character_count": report.character_count,
                            "calibration_count": report.calibration_count,
                            "tip_date_count": report.tip_date_count,
                            "chain_length": report.chain_length,
                            "logger_count": report.logger_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-run":
                report = run_beast_posterior_inference(
                    args.input_path,
                    executable=args.executable or "beast",
                    overwrite=not args.no_overwrite,
                    threads=args.threads,
                    seed=args.seed,
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
                            "warning_count": len(report.run.warning_lines),
                            "threads": args.threads,
                            "seed": args.seed,
                            "overwrite": not args.no_overwrite,
                            "resumed": report.resumed,
                            "timeout_seconds": report.run.timeout_seconds,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-calibrations":
                report = validate_fossil_calibration_table(
                    args.tree_path, args.calibration_path
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.tree_path, args.calibration_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.tree_path, args.calibration_path],
                        outputs=outputs,
                        metrics={
                            "calibration_count": report.calibration_count,
                            "invalid_calibration_count": report.invalid_calibration_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-tip-dates":
                report = validate_tip_dating_metadata(
                    args.tree_path,
                    args.tip_dates_path,
                    alignment_path=args.alignment,
                    date_column=args.date_column,
                )
                inputs = [
                    args.tree_path,
                    args.tip_dates_path,
                    *([args.alignment] if args.alignment is not None else []),
                ]
                outputs = _finalize_outputs(args, command="adapter", inputs=inputs)
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=inputs,
                        outputs=outputs,
                        metrics={
                            "valid_tip_count": report.valid_tip_count,
                            "invalid_tip_count": report.invalid_tip_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-log":
                report = parse_beast_log(args.input_path)
                summary = summarize_beast_log(
                    args.input_path, burnin_fraction=args.burnin_fraction
                )
                outputs: list[Path | str] = []
                if args.summary_out is not None:
                    outputs.append(
                        write_beast_log_summary_table(args.summary_out, summary)
                    )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "row_count": report.row_count,
                            "column_count": len(report.columns),
                            "burnin_fraction": summary.burnin_fraction,
                            "kept_row_count": summary.kept_row_count,
                            "posterior_parameter_count": len(
                                summary.posterior_parameters
                            ),
                            "likelihood_parameter_count": len(
                                summary.likelihood_parameters
                            ),
                            "prior_parameter_count": len(summary.prior_parameters),
                            "clock_parameter_count": len(summary.clock_parameters),
                            "tree_parameter_count": len(summary.tree_parameters),
                        },
                        data={"log": report, "summary": summary},
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-burnin-sensitivity":
                report = assess_beast_burnin_sensitivity(
                    args.posterior_trees,
                    log_path=args.log,
                    burnin_fractions=tuple(args.burnin_fractions),
                )
                inputs = [
                    args.posterior_trees,
                    *([args.log] if args.log is not None else []),
                ]
                outputs: list[Path | str] = []
                if args.slice_out is not None:
                    outputs.append(
                        write_beast_burnin_sensitivity_slice_table(
                            args.slice_out,
                            report,
                        )
                    )
                if args.parameter_out is not None:
                    outputs.append(
                        write_burnin_parameter_shift_table(
                            args.parameter_out,
                            report.parameter_shifts,
                        )
                    )
                if args.clade_out is not None:
                    outputs.append(
                        write_burnin_clade_shift_table(
                            args.clade_out,
                            report.clade_shifts,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=inputs,
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=inputs,
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "slice_count": len(report.slices),
                            "parameter_shift_count": len(report.parameter_shifts),
                            "unstable_parameter_count": report.unstable_parameter_count,
                            "clade_shift_count": len(report.clade_shifts),
                            "unstable_clade_count": report.unstable_clade_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-parameters":
                report = summarize_beast_log(
                    args.input_path,
                    burnin_fraction=args.burnin_fraction,
                )
                outputs: list[Path | str] = []
                if args.summary_out is not None:
                    outputs.append(
                        write_beast_log_summary_table(args.summary_out, report)
                    )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "burnin_fraction": report.burnin_fraction,
                            "kept_row_count": report.kept_row_count,
                            "parameter_count": len(report.parameter_summaries),
                            "posterior_parameter_count": len(
                                report.posterior_parameters
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-trees":
                report = parse_beast_posterior_tree_samples(
                    args.input_path,
                    burnin_fraction=args.burnin_fraction,
                )
                outputs: list[Path | str] = []
                if args.tree_set_out is not None:
                    outputs.append(
                        write_beast_posterior_tree_set(args.tree_set_out, report)
                    )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "total_tree_count": report.total_tree_count,
                            "kept_tree_count": report.kept_tree_count,
                            "rooted_tree_count": report.rooted_tree_count,
                            "burnin_fraction": report.burnin_fraction,
                            "clade_count": len(report.clades),
                            "sampled_state_count": len(report.sampled_states),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-subsample":
                report = subsample_beast_posterior_tree_set(
                    args.input_path,
                    method=args.method,
                    thinning_interval=args.thinning_interval,
                    sample_count=args.sample_count,
                    burnin_fraction=args.burnin_fraction,
                    random_seed=args.seed,
                )
                outputs: list[Path | str] = []
                if args.tree_set_out is not None:
                    outputs.append(
                        write_posterior_tree_subsample(args.tree_set_out, report)
                    )
                if args.sample_table_out is not None:
                    outputs.append(
                        write_posterior_tree_subsample_table(
                            args.sample_table_out,
                            report,
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "total_tree_count": report.total_tree_count,
                            "burnin_tree_count": report.burnin_tree_count,
                            "pre_subsampling_tree_count": (
                                report.pre_subsampling_tree_count
                            ),
                            "retained_tree_count": report.retained_tree_count,
                            "selection_method": report.selection_method,
                            "retained_state_count": len(
                                [
                                    tree
                                    for tree in report.trees
                                    if tree.state is not None
                                ]
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-consensus":
                consensus_tree, report = summarize_beast_posterior_trees(
                    args.input_path,
                    burnin_fraction=args.burnin_fraction,
                )
                outputs: list[Path | str] = [write_newick(args.out, consensus_tree)]
                if args.tree_set_out is not None:
                    args.tree_set_out.parent.mkdir(parents=True, exist_ok=True)
                    args.tree_set_out.write_text(
                        report.retained_tree_set_path.read_text(encoding="utf-8"),
                        encoding="utf-8",
                    )
                    outputs.append(args.tree_set_out)
                if args.clade_table_out is not None:
                    outputs.append(
                        write_clade_frequency_table(
                            args.clade_table_out,
                            compute_clade_frequency_table(
                                report.retained_tree_set_path
                            ),
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "total_tree_count": report.total_tree_count,
                            "kept_tree_count": report.kept_tree_count,
                            "annotated_node_count": report.annotated_node_count,
                            "clade_frequency_count": report.clade_frequency_count,
                            "burnin_fraction": report.burnin_fraction,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-diversity":
                report = summarize_beast_posterior_topology_diversity(
                    args.input_path,
                    burnin_fraction=args.burnin_fraction,
                )
                outputs: list[Path | str] = []
                if args.tree_set_out is not None:
                    args.tree_set_out.parent.mkdir(parents=True, exist_ok=True)
                    args.tree_set_out.write_text(
                        Path(report.retained_tree_set_path).read_text(encoding="utf-8"),
                        encoding="utf-8",
                    )
                    outputs.append(args.tree_set_out)
                if args.distance_out is not None:
                    outputs.append(
                        write_tree_distance_matrix(
                            args.distance_out,
                            compute_tree_distance_matrix(report.retained_tree_set_path),
                        )
                    )
                if args.topology_out is not None:
                    outputs.append(
                        write_topology_cluster_table(
                            args.topology_out,
                            cluster_trees_by_topology(report.retained_tree_set_path),
                        )
                    )
                if args.unstable_clade_out is not None:
                    outputs.append(
                        write_unstable_clade_table(
                            args.unstable_clade_out,
                            detect_unstable_clades(report.retained_tree_set_path),
                        )
                    )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=outputs,
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        metrics={
                            "total_tree_count": report.total_tree_count,
                            "kept_tree_count": report.kept_tree_count,
                            "rooted_topology_count": report.rooted_topology_count,
                            "dominant_topology_frequency": (
                                report.dominant_topology_frequency
                            ),
                            "pair_count": report.pair_count,
                            "unstable_clade_count": report.unstable_clade_count,
                            "burnin_fraction": report.burnin_fraction,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-convergence":
                report = assess_beast_convergence(
                    args.input_path,
                    burnin_fraction=args.burnin_fraction,
                    ess_threshold=args.ess_threshold,
                    mean_shift_threshold=args.mean_shift_threshold,
                )
                outputs = _finalize_outputs(
                    args, command="adapter", inputs=[args.input_path]
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=[warning["message"] for warning in report.warnings],
                        metrics={
                            "warning_count": len(report.warnings),
                            "converged": report.converged,
                            "burnin_fraction": report.burnin_fraction,
                            "sample_count": report.sample_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "beast-calibration-report":
                report = render_calibration_audit_report(
                    tree_path=args.tree_path,
                    calibration_path=args.calibration_path,
                    out_path=args.out,
                    tip_dates_path=args.tip_dates,
                    alignment_path=args.alignment,
                    date_column=args.date_column,
                )
                inputs = [
                    args.tree_path,
                    args.calibration_path,
                    *([args.tip_dates] if args.tip_dates is not None else []),
                    *([args.alignment] if args.alignment is not None else []),
                ]
                outputs = _finalize_outputs(
                    args, command="adapter", inputs=inputs, outputs=[report.output_path]
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=inputs,
                        outputs=outputs,
                        warnings=method_tier_warnings(report.method_tier),
                        metrics={
                            "invalid_calibration_count": report.invalid_calibration_count,
                            "warning_count": report.warning_count
                            + len(method_tier_warnings(report.method_tier)),
                            **method_tier_metrics(report.method_tier),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "bayesian-evidence":
                report = build_bayesian_evidence_package(
                    bundle_root=args.out_dir,
                    input_paths=args.inputs,
                    config_paths=args.configs,
                    tree_paths=args.trees,
                    log_paths=args.logs,
                    diagnostic_paths=args.diagnostics,
                    report_paths=args.reports,
                )
                inputs = [
                    *args.inputs,
                    *args.configs,
                    *args.trees,
                    *args.logs,
                    *args.diagnostics,
                    *args.reports,
                ]
                outputs = _finalize_outputs(
                    args, command="adapter", inputs=inputs, outputs=[args.out_dir]
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=inputs,
                        outputs=outputs,
                        metrics={
                            "file_count": report.file_count,
                            "valid": report.valid,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "bayesian-diagnostics-table":
                report = write_supplementary_bayesian_diagnostics_table(
                    args.out,
                    posterior_tree_path=args.posterior_trees,
                    primary_log_path=args.log,
                    additional_log_paths=args.additional_logs,
                    burnin_fractions=tuple(args.burnin_fractions),
                    ess_threshold=args.ess_threshold,
                    mean_shift_threshold=args.mean_shift_threshold,
                    cross_chain_mean_shift_threshold=args.cross_chain_mean_shift_threshold,
                )
                inputs = [args.posterior_trees, args.log, *(args.additional_logs or [])]
                outputs = _finalize_outputs(
                    args, command="adapter", inputs=inputs, outputs=[args.out]
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=inputs,
                        outputs=outputs,
                        metrics={
                            "row_count": report.row_count,
                            "warning_count": report.warning_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "bayesian-methods":
                report = write_bayesian_methods_summary_text(
                    args.out,
                    posterior_tree_path=args.posterior_trees,
                    primary_log_path=args.log,
                    additional_log_paths=args.additional_logs,
                    analysis_xml_path=args.analysis_xml,
                    tree_prior=args.tree_prior,
                    clock_model=args.clock_model,
                    calibration_path=args.calibration_path,
                    tip_dates_path=args.tip_dates_path,
                    burnin_fractions=tuple(args.burnin_fractions),
                    ess_threshold=args.ess_threshold,
                    mean_shift_threshold=args.mean_shift_threshold,
                    cross_chain_mean_shift_threshold=args.cross_chain_mean_shift_threshold,
                )
                inputs = [
                    args.posterior_trees,
                    args.log,
                    *(args.additional_logs or []),
                    *([args.analysis_xml] if args.analysis_xml is not None else []),
                ]
                outputs = _finalize_outputs(
                    args, command="adapter", inputs=inputs, outputs=[args.out]
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=inputs,
                        outputs=outputs,
                        metrics={"warning_count": report.warning_count},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "infer-large":
                report = run_large_alignment_inference(
                    args.input_path,
                    out_dir=args.out_dir,
                    prefix=args.prefix,
                    sequence_type=args.sequence_type,
                    executable=args.executable or "FastTree",
                    timeout_seconds=args.timeout_seconds,
                    resume=args.resume,
                    incomplete_run_policy=args.incomplete_run_policy,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[*report.output_paths.values()],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "sequence_count": report.input_summary.sequence_count,
                            "alignment_length": report.input_summary.alignment_length,
                            "total_site_cells": report.input_summary.total_site_cells,
                            "sequence_type": report.sequence_type,
                            "resumed": report.resumed,
                            "timeout_seconds": report.timeout_seconds,
                            "peak_memory_bytes": max(
                                (
                                    row.peak_memory_bytes
                                    for row in report.resource_rows
                                    if row.peak_memory_bytes is not None
                                ),
                                default=None,
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "compare":
                report = compare_fast_and_ml_trees(
                    args.fast_tree, args.ml_tree, out_path=args.out
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.fast_tree, args.ml_tree],
                    outputs=[report.comparison_report.output_path],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.fast_tree, args.ml_tree],
                        outputs=outputs,
                        metrics={
                            "shared_taxa": len(
                                report.comparison_report.topology.shared_taxa
                            ),
                            "robinson_foulds_distance": report.comparison_report.topology.robinson_foulds_distance,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "compare-engines":
                report = run_tree_inference_comparison(
                    args.input_path,
                    out_dir=args.out_dir,
                    prefix=args.prefix,
                    sequence_type=args.sequence_type,
                    iqtree_executable=args.iqtree_executable or "iqtree2",
                    fasttree_executable=args.fasttree_executable or "FastTree",
                    iqtree_seed=args.iqtree_seed,
                    iqtree_threads=args.iqtree_threads,
                    bootstrap_replicates=args.bootstrap_replicates,
                    resume=args.resume,
                    timeout_seconds=args.timeout_seconds,
                    incomplete_run_policy=args.incomplete_run_policy,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[*report.output_paths.values()],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "selected_model": report.selected_model,
                            "shared_taxa": len(
                                report.engine_comparison.topology.shared_taxa
                            ),
                            "robinson_foulds_distance": report.engine_comparison.topology.robinson_foulds_distance,
                            "shared_clade_count": len(report.shared_clade_rows),
                            "conflicting_clade_count": len(
                                report.conflicting_clade_rows
                            ),
                            "stable_clade_count": report.conclusion_summary.stable_clade_count,
                            "unstable_clade_count": report.conclusion_summary.unstable_clade_count,
                            "engine_specific_clade_count": report.conclusion_summary.engine_specific_clade_count,
                            "support_disagreement_count": sum(
                                1
                                for row in report.conflicting_clade_rows
                                if row.conflict_kind == "support_disagreement"
                            ),
                            "high_support_conflict_count": report.conclusion_summary.high_support_conflict_count,
                            "low_support_disagreement_count": report.conclusion_summary.low_support_disagreement_count,
                            "serious_conflict_count": report.conclusion_summary.serious_conflict_count,
                            "resumed": any(
                                workflow.resumed
                                for workflow in (
                                    report.model_selection_workflow,
                                    report.iqtree_support_workflow,
                                    report.fasttree_workflow,
                                )
                            ),
                            "timeout_seconds": args.timeout_seconds,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.adapter_command == "reproducibility":
                report = run_inference_reproducibility_check(
                    args.input_path,
                    out_dir=args.out_dir,
                    prefix=args.prefix,
                    sequence_type=args.sequence_type,
                    executable=args.iqtree_executable or "iqtree2",
                    repeats=args.repeats,
                    bootstrap_replicates=args.bootstrap_replicates,
                    seed=args.iqtree_seed,
                    threads=args.iqtree_threads,
                )
                outputs = _finalize_outputs(
                    args,
                    command="adapter",
                    inputs=[args.input_path],
                    outputs=[*report.output_paths.values()],
                )
                _print_result(
                    build_command_result(
                        command="adapter",
                        inputs=[args.input_path],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "selected_model": report.selected_model,
                            "overall_status": report.overall_status,
                            "repeat_count": report.repeat_count,
                            "unstable_comparison_count": sum(
                                1
                                for row in report.comparison_rows
                                if row.classification == "unstable"
                            ),
                            "equivalent_comparison_count": sum(
                                1
                                for row in report.comparison_rows
                                if row.classification == "equivalent"
                            ),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
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
