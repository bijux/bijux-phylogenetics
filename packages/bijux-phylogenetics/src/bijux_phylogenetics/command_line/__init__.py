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
from bijux_phylogenetics.command_line.engines import (
    add_phylo_commands,
    run_phylo_command,
)
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
from bijux_phylogenetics.command_line.comparative_trait_review import (
    add_comparative_trait_review_commands,
    run_comparative_trait_review_command,
)
from bijux_phylogenetics.command_line.comparative_support import (
    add_comparative_support_commands,
    run_comparative_support_command,
)
from bijux_phylogenetics.command_line.comparative_modeling import (
    add_comparative_modeling_commands,
    run_comparative_modeling_command,
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
from bijux_phylogenetics.command_line.tree_normalization import (
    add_tree_normalization_commands,
    run_normalize_command,
    run_normalize_taxa_command,
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
    _json_requested,
    _parse_assignment_map,
    _parse_time_bin_definition,
    _parse_transition_pairs,
    _validate_ancestral_discrete_model_arguments,
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
from bijux_phylogenetics.core.environment import inspect_environment
from bijux_phylogenetics.engines import (
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

    add_phylo_commands(subparsers)

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
    add_comparative_trait_review_commands(comparative_subparsers)
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
    add_comparative_support_commands(comparative_subparsers)
    add_comparative_modeling_commands(comparative_subparsers)

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

    add_tree_normalization_commands(subparsers)

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
            review_exit_code = run_comparative_trait_review_command(
                args,
                parser=parser,
            )
            if review_exit_code is not None:
                return review_exit_code
            support_exit_code = run_comparative_support_command(
                args,
                parser=parser,
            )
            if support_exit_code is not None:
                return support_exit_code
            modeling_exit_code = run_comparative_modeling_command(
                args,
                parser=parser,
            )
            if modeling_exit_code is not None:
                return modeling_exit_code
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
            return run_normalize_command(args)
        if args.command == "normalize-taxa":
            return run_normalize_taxa_command(args)
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
