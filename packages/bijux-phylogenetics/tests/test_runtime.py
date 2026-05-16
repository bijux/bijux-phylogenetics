from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

import bijux_phylogenetics
from bijux_phylogenetics.ancestral import (
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
    reconstruct_continuous_ancestral_states_from_dataset,
    reconstruct_discrete_ancestral_states,
    reconstruct_discrete_ancestral_states_from_dataset,
    render_ancestral_state_report,
    render_ancestral_state_tree,
    render_ancestral_state_visualization,
    validate_discrete_ancestral_reference_examples,
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
    summarize_continuous_change_branches,
    summarize_continuous_change_counts,
    summarize_discrete_ancestral_confidence,
    summarize_discrete_ancestral_report,
    summarize_discrete_ancestral_tree_set,
    summarize_discrete_ancestral_tree_set_confidence,
    summarize_discrete_ancestral_tree_set_report,
    summarize_irreversible_discrete_reconstruction,
    summarize_irreversible_discrete_report,
    summarize_ordered_discrete_reconstruction,
    summarize_ordered_discrete_report,
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
    write_continuous_change_branch_table,
    write_continuous_change_count_table,
    write_discrete_ancestral_comparison_table,
    write_discrete_ancestral_confidence_table,
    write_discrete_ancestral_exclusion_table,
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
from bijux_phylogenetics.bayesian import (
    BeastAnalysisXmlReport,
    BeastCalibration,
    BeastPosteriorConsensusReport,
    BeastPosteriorTopologyDiversityReport,
    MrBayesBurninSensitivityReport,
    MrBayesBurninSensitivitySlice,
    MrBayesParameterDiagnosticsReport,
    MrBayesParameterSummary,
    assess_beast_burnin_sensitivity,
    assess_beast_chain_mixing,
    assess_beast_convergence,
    assess_mrbayes_burnin_sensitivity,
    assess_mrbayes_convergence,
    build_bayesian_evidence_package,
    build_posterior_uncertainty_figure_package,
    compare_bayesian_tree_sets,
    compare_independent_bayesian_runs,
    compare_posterior_tree_sets_by_clock,
    compare_posterior_tree_sets_by_prior,
    compute_mrbayes_effective_sample_sizes,
    detect_impossible_calibration_constraints,
    parse_beast_log,
    parse_beast_posterior_tree_samples,
    parse_mrbayes_consensus_tree,
    parse_mrbayes_mcmc_diagnostics,
    parse_mrbayes_parameter_traces,
    parse_mrbayes_posterior_tree_samples,
    prepare_beast_time_tree_analysis,
    prepare_mrbayes_analysis,
    render_bayesian_diagnostics_report,
    render_bayesian_posterior_report,
    render_bayesian_run_comparison_report,
    render_calibration_audit_report,
    run_beast_posterior_inference,
    run_mrbayes_posterior_inference,
    subsample_beast_posterior_tree_set,
    subsample_mrbayes_posterior_tree_set,
    subsample_posterior_tree_set,
    summarize_beast_analysis_xml,
    summarize_beast_log,
    summarize_beast_posterior_topology_diversity,
    summarize_beast_posterior_trees,
    summarize_maximum_clade_credibility_tree,
    summarize_mrbayes_parameter_diagnostics,
    summarize_mrbayes_posterior_trees,
    summarize_posterior_node_ages,
    thin_posterior_tree_set,
    validate_beast_analysis_xml,
    validate_beast_posterior_log,
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
from bijux_phylogenetics.benchmark import (
    benchmark_alignment_diagnostics,
    benchmark_tree_comparison,
    benchmark_tree_validation,
)
from bijux_phylogenetics.biogeography import (
    BiogeographyRegionCountRow,
    BiogeographyReportExclusionRow,
    BiogeographyReportPackageResult,
    ConstrainedGeographicFitRow,
    ConstrainedGeographicReport,
    ConstrainedGeographicSummary,
    ConstrainedGeographicTransitionRow,
    DatedBiogeographyEventRow,
    DatedBiogeographyNodeRow,
    DatedBiogeographyReport,
    DatedBiogeographySummary,
    DatedBiogeographyTimeBinRow,
    GeographicMigrationEventReport,
    GeographicMigrationEventRow,
    GeographicMigrationEventSummary,
    GeographicMigrationTreeRow,
    GeographicMigrationTreeSetEventRow,
    GeographicMigrationTreeSetEventSummaryRow,
    GeographicMigrationTreeSetReport,
    GeographicMigrationTreeSetSummary,
    GeographicSamplingBiasNodeRow,
    GeographicSamplingBiasReport,
    GeographicSamplingBiasSummary,
    GeographicSamplingBiasTransitionRow,
    GeographicSamplingCountRow,
    TimeBinDefinition,
    TimeStratifiedBranchRow,
    TimeStratifiedTransitionMatrixRow,
    TimeStratifiedTransitionReport,
    TimeStratifiedTransitionSummary,
    UnsupportedGeographicTransitionClaimRow,
    build_biogeography_report_package,
    summarize_biogeographic_transition_chronology,
    summarize_biogeography_region_counts,
    summarize_constrained_geographic_model,
    summarize_constrained_geographic_report,
    summarize_geographic_migration_event_tree_set,
    summarize_geographic_migration_events,
    summarize_geographic_sampling_bias,
    summarize_geographic_state_model,
    summarize_time_stratified_geographic_transitions,
    write_biogeography_region_count_table,
    write_biogeography_report_exclusion_table,
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
from bijux_phylogenetics.branch_lengths import (
    BranchLengthAggregate,
    BranchLengthDistributionReport,
    BranchLengthRow,
    analyze_branch_length_distribution,
    analyze_tree_set_branch_lengths,
    write_branch_length_table,
)
from bijux_phylogenetics.clades import (
    CladeMetadataObservation,
    CladeTableReport,
    CladeTableRow,
    extract_tree_clades,
    extract_tree_set_clades,
    write_clade_table,
)
from bijux_phylogenetics.cli import main
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.comparative import (
    BranchIdentityMetadata,
    BrownianRegimeBranchRow,
    BrownianRegimeExclusion,
    BrownianRegimeFitSummaryReport,
    BrownianRegimeIdentifiabilityWarning,
    BrownianRegimeProfileRow,
    BrownianRegimeRateRow,
    CladeTraitExclusion,
    CladeTraitRow,
    CladeTraitStateCount,
    CladeTraitSummaryReport,
    ComparativeAnalysisSummaryRow,
    ComparativeAuditTableRow,
    ComparativeCladeCoefficientChangeRow,
    ComparativeCladeResidualReport,
    ComparativeCladeStabilityReport,
    ComparativeCladeStabilityRow,
    ComparativeCoefficientTableRow,
    ComparativeInterpretationRow,
    ComparativeRegressionModelExclusion,
    ComparativeRegressionModelRow,
    ComparativeRegressionModelSelectionReport,
    ComparativeRegressionPairwiseComparisonRow,
    ComparativeReportPackageResult,
    ComparativeResidualCladeRow,
    ComparativeResidualTableRow,
    ComparativeResidualTaxonRow,
    ComparativeSignalTableRow,
    CorrelatedTraitComparisonRow,
    CorrelatedTraitEvolutionReport,
    CorrelatedTraitExclusion,
    CorrelatedTraitObservationRow,
    EarlyBurstIdentifiabilityWarning,
    EarlyBurstRateChangeProfileRow,
    EarlyBurstTraitEvolutionExclusion,
    EarlyBurstTraitEvolutionSummaryReport,
    IndependentContrastRegressionReport,
    IndependentContrastRegressionRow,
    MultivariateComparativeRegressionReport,
    MultivariateResidualAssociationRow,
    MultivariateResidualCorrelationRow,
    MultivariateResidualCovarianceRow,
    MultivariateResponseCoefficientRow,
    MultivariateResponseModelRow,
    MultivariateTaxonExclusion,
    PhylogeneticLogisticCoefficient,
    PhylogeneticLogisticFittedRow,
    PhylogeneticLogisticReport,
    PhylogeneticLogisticWarning,
    PhylogeneticSignalPermutation,
    PhylogeneticSignalSummaryReport,
    PosteriorTreePGLSCoefficientRow,
    PosteriorTreePGLSCoefficientSummaryRow,
    PosteriorTreePGLSReport,
    PosteriorTreePGLSTreeFitRow,
    TraitImputationExclusion,
    TraitImputationHoldoutRow,
    TraitImputationRow,
    TraitImputationSummaryReport,
    TraitOutlierExclusion,
    TraitOutlierSummaryReport,
    TraitOutlierTaxonRow,
    TraitRateThroughTimeExclusion,
    TraitRateThroughTimeIntervalRow,
    TraitRateThroughTimeSummaryReport,
    TraitRegimeBranchRow,
    TraitRegimeExclusion,
    TraitRegimeMappingReport,
    TraitRegimeNodeRow,
    analyze_comparative_clade_stability,
    analyze_comparative_residual_clades,
    assess_comparative_method_maturity,
    audit_comparative_parameter_uncertainty,
    audit_ou_identifiability_reference_examples,
    build_branch_identity_lookup,
    build_comparative_report_package,
    build_pgls_model_matrix,
    compare_comparative_regression_models,
    compute_blombergs_k,
    fit_discrete_mk_model,
    fit_discrete_mk_model_from_dataset,
    compute_phylogenetic_independent_contrasts,
    compute_phylogenetic_independent_contrasts_from_dataset,
    compute_phylogenetic_signal_test,
    evaluate_pagels_lambda_likelihood,
    evaluate_pagels_lambda_likelihood_from_dataset,
    estimate_pagels_lambda,
    DiscreteMkFitReport,
    DiscreteMkInputAudit,
    PagelLambdaLikelihoodReport,
    PagelLambdaOptimizerDiagnostics,
    PagelLambdaProfileRow,
    inspect_pgls_inputs,
    run_multivariate_comparative_regression,
    run_pgls,
    run_posterior_tree_pgls,
    summarize_brownian_covariance,
    summarize_brownian_covariance_from_tree,
    summarize_brownian_covariance_pgls,
    summarize_brownian_regime_rates,
    summarize_brownian_trait_evolution,
    summarize_clade_traits,
    summarize_comparative_analysis,
    summarize_comparative_audit,
    summarize_comparative_coefficients,
    summarize_comparative_interpretation,
    summarize_comparative_residuals,
    summarize_comparative_signal,
    summarize_correlated_trait_evolution,
    summarize_early_burst_trait_evolution,
    summarize_independent_contrast_regression,
    summarize_numeric_trait,
    summarize_numeric_trait_readiness,
    summarize_ou_covariance_pgls,
    summarize_ou_trait_evolution,
    summarize_pgls_categorical_contrasts,
    summarize_pgls_interaction_coefficients,
    summarize_pgls_lambda_fit,
    summarize_phylogenetic_logistic,
    summarize_phylogenetic_signal,
    summarize_trait_imputation,
    summarize_trait_outliers,
    summarize_trait_rate_through_time,
    summarize_trait_regime_mapping,
    write_brownian_covariance_long_table,
    write_brownian_covariance_matrix_table,
    write_brownian_covariance_table,
    write_brownian_regime_branch_table,
    write_brownian_regime_comparison_table,
    write_brownian_regime_exclusion_table,
    write_brownian_regime_profile_table,
    write_brownian_regime_rate_table,
    write_brownian_regime_summary_table,
    write_brownian_trait_evolution_exclusion_table,
    write_brownian_trait_evolution_summary_table,
    write_clade_trait_clade_table,
    write_clade_trait_exclusion_table,
    write_clade_trait_summary_table,
    write_comparative_audit_table,
    write_comparative_clade_coefficient_change_table,
    write_comparative_clade_stability_table,
    write_comparative_coefficient_table,
    write_comparative_contrast_table,
    write_comparative_interpretation_table,
    write_comparative_model_comparison_table,
    write_comparative_regression_excluded_taxa_table,
    write_comparative_regression_model_ranking_table,
    write_comparative_regression_pairwise_table,
    write_comparative_residual_clade_table,
    write_comparative_residual_table,
    write_comparative_residual_taxon_table,
    write_comparative_signal_table,
    write_comparative_summary_table,
    write_correlated_trait_comparison_table,
    write_correlated_trait_exclusion_table,
    write_correlated_trait_observation_table,
    write_correlated_trait_summary_table,
    write_early_burst_rate_change_profile_table,
    write_early_burst_trait_evolution_comparison_table,
    write_early_burst_trait_evolution_exclusion_table,
    write_early_burst_trait_evolution_summary_table,
    write_independent_contrast_regression_table,
    write_independent_contrast_table,
    write_multivariate_excluded_taxa_table,
    write_multivariate_residual_association_table,
    write_multivariate_residual_correlation_table,
    write_multivariate_residual_covariance_table,
    write_multivariate_response_coefficient_table,
    write_multivariate_response_model_table,
    write_ou_alpha_profile_table,
    write_ou_covariance_table,
    write_ou_trait_evolution_exclusion_table,
    write_ou_trait_evolution_summary_table,
    write_pgls_categorical_contrast_table,
    write_pgls_interaction_coefficient_table,
    write_pgls_lambda_profile_table,
    write_pgls_model_matrix_table,
    write_phylogenetic_logistic_coefficient_table,
    write_phylogenetic_logistic_excluded_taxa_table,
    write_phylogenetic_logistic_fitted_table,
    write_phylogenetic_signal_permutation_table,
    write_phylogenetic_signal_summary_table,
    write_posterior_tree_pgls_coefficient_table,
    write_posterior_tree_pgls_summary_table,
    write_posterior_tree_pgls_tree_table,
    write_trait_imputation_exclusion_table,
    write_trait_imputation_holdout_table,
    write_trait_imputation_summary_table,
    write_trait_imputation_table,
    write_trait_outlier_exclusion_table,
    write_trait_outlier_summary_table,
    write_trait_outlier_taxon_table,
    write_trait_rate_through_time_exclusion_table,
    write_trait_rate_through_time_interval_table,
    write_trait_rate_through_time_summary_table,
    write_trait_regime_branch_table,
    write_trait_regime_exclusion_table,
    write_trait_regime_node_table,
    write_trait_regime_summary_table,
)
from bijux_phylogenetics.comparative.evidence_contract import (
    SUPPORTED_EVIDENCE_API_LOCATORS,
    resolve_supported_evidence_api,
)
import bijux_phylogenetics.compare as compare_api
from bijux_phylogenetics.compare.reports import build_tree_comparison_report
from bijux_phylogenetics.compare.taxon_influence import (
    TaxonInfluenceReport,
    TaxonInfluenceRow,
    analyze_taxon_influence,
    write_taxon_influence_table,
)
from bijux_phylogenetics.compare.topology_distance import (
    TopologyDistanceReport,
    TopologyDistanceSplitRow,
    compare_topology_distance,
    write_topology_distance_split_table,
)
from bijux_phylogenetics.compare.topology import (
    BranchScoreComparisonReport,
    CladeOverlapComparisonReport,
    RobinsonFouldsComparisonReport,
    SharedTaxaPruningReport,
    SupportComparisonReport,
    SupportConflictRow,
    compare_branch_lengths,
    compare_branch_score_distance,
    compare_clade_overlap,
    compare_clade_sets,
    compare_robinson_foulds,
    compare_support_values,
    compare_tree_paths,
    detect_clade_changes,
    prune_trees_to_shared_taxa,
    write_clade_overlap_table,
    write_shared_taxa_pruning_table,
    write_shared_taxa_removed_taxa_table,
    write_support_comparison_table,
)
from bijux_phylogenetics.compare.support_reference import (
    SupportReferenceObservation,
    SupportReferenceValidationReport,
    validate_support_reference_examples,
)
from bijux_phylogenetics.compare.tree_distance_reference import (
    TreeDistanceReferenceObservation,
    TreeDistanceReferenceValidationReport,
    validate_tree_distance_reference_examples,
)
from bijux_phylogenetics.core.alignment import AlignmentRecord, AlignmentSummary
from bijux_phylogenetics.core.dataset import (
    audit_dataset_inputs,
    audit_dataset_taxon_ordering,
    build_dataset_completeness_matrix,
    build_dataset_crosswalk,
    build_dataset_mismatch_report,
    summarize_dataset_readiness,
)
from bijux_phylogenetics.core.demo import run_capability_demo
from bijux_phylogenetics.core.environment import inspect_environment
from bijux_phylogenetics.core.locus_occupancy import (
    LocusOccupancyFilterIteration,
    build_locus_occupancy_report,
    filter_locus_occupancy,
)
from bijux_phylogenetics.core.manifest import build_run_manifest, write_run_manifest
from bijux_phylogenetics.core.metadata import inspect_metadata_table, join_table_to_taxa
from bijux_phylogenetics.core.partitions import (
    build_partition_summary_report,
    parse_locus_partitions,
    write_locus_partitions,
    write_partition_summary_table,
)
from bijux_phylogenetics.core.pruning import (
    drop_tree_taxa,
    prune_alignment_to_tree,
    prune_tree_to_alignment,
    prune_tree_to_requested_taxa,
    prune_tree_to_taxa,
)
from bijux_phylogenetics.core.taxonomy import (
    build_taxon_audit_report,
    build_taxon_mapping_conflict_report,
    detect_duplicate_biological_identities,
    export_tree_accepted_names,
    infer_taxon_rank,
    inspect_tree_taxa_safety,
    inspect_tree_taxon_identity,
    inspect_tree_taxon_rank_consistency,
    normalize_tree_taxa,
    write_accepted_name_mapping,
    write_taxon_mapping,
)
from bijux_phylogenetics.core.topology import (
    TreeRootingReport,
    assess_tree_monophyly,
    collapse_branches_below_length,
    extract_named_clade,
    extract_tree_clade_by_descendant_taxa,
    extract_tree_clade_by_node_id,
    find_tree_mrca,
    ladderize_tree,
    reroot_tree_by_midpoint,
    root_tree_on_outgroup,
    rotate_all_internal_nodes,
    rotate_named_node,
    sort_tree_tips_alphabetically,
    unroot_tree,
    write_tree_rooting_report,
)
from bijux_phylogenetics.core.branching_times import (
    TreeBranchingTimeReport,
    TreeBranchingTimeRow,
    compute_tree_branching_times,
    write_tree_branching_time_table,
)
from bijux_phylogenetics.core.node_depth import (
    TreeNodeDepthReport,
    TreeNodeDepthRow,
    compute_tree_node_depths,
    write_tree_node_depth_table,
)
from bijux_phylogenetics.core.ultrametric import (
    APE_ULTRAMETRIC_TOLERANCE,
    TreeUltrametricReport,
    TreeUltrametricTipRow,
    assess_tree_ultrametricity,
    write_tree_ultrametric_table,
)
from bijux_phylogenetics.core.tree_distance import (
    TipDistanceMatrixReport,
    TipDistanceMatrixRow,
    compute_tree_tip_distance_matrix,
    write_tree_tip_distance_long_table,
    write_tree_tip_distance_matrix,
)
from bijux_phylogenetics.core.traits import (
    detect_missing_trait_values,
    detect_unusable_trait_columns,
    link_tree_to_traits,
    prune_traits_to_tree,
    validate_traits_table,
)
from bijux_phylogenetics.diagnostics.assumptions import (
    assess_tree_assumptions,
    inspect_branch_length_units,
    standardize_support_labels,
)
from bijux_phylogenetics.diagnostics.root_to_tip import (
    compute_root_to_tip_distances,
    diagnose_ultrametricity,
)
from bijux_phylogenetics.diagnostics.validation import (
    diagnose_tree_path,
    forensic_tree_path,
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.discrete_evolution import (
    compare_discrete_state_models,
    count_discrete_stochastic_map_transitions,
    detect_state_imbalance_problems,
    estimate_ancestral_geographic_states,
    load_stochastic_map_collection,
    render_stochastic_map_density_artifact,
    render_discrete_state_evolution_report,
    render_tree_with_geographic_states,
    run_discrete_state_transition_model,
    simulate_discrete_stochastic_maps,
    simulate_discrete_stochastic_maps_from_fit_report,
    summarize_discrete_stochastic_map_density,
    summarize_discrete_stochastic_maps,
    validate_discrete_state_coding,
    write_discrete_model_comparison_table,
    write_node_state_probability_table,
    write_stochastic_map_aggregate_transition_matrix,
    write_stochastic_map_branch_probability_table,
    write_stochastic_map_branch_transition_count_table,
    write_stochastic_map_collection,
    write_stochastic_map_branch_occupancy_table,
    write_stochastic_map_density_branch_table,
    write_stochastic_map_density_slice_table,
    write_stochastic_map_event_table,
    write_stochastic_map_segment_table,
    write_stochastic_map_state_time_table,
    write_stochastic_map_transition_count_matrix,
    write_stochastic_map_summary_table,
    write_transition_summary_table,
)
from bijux_phylogenetics.distance import (
    assess_distance_method_maturity,
    build_distance_method_report,
    build_distance_tree,
    build_distance_tree_from_genetic_distance_matrix,
    build_tree_from_imported_distance_matrix,
    compare_distance_gap_policies,
    compare_distance_models,
    compare_distance_tree_to_reference_tree,
    compare_distance_tree_topologies,
    compute_pairwise_genetic_distance_matrix_from_dna_bin_alignment,
    compute_pairwise_genetic_distance_matrix,
    list_distance_tree_method_policies,
    load_imported_distance_matrix,
    resolve_distance_tree_method_policy,
    summarize_distance_bootstrap_support,
    validate_imported_distance_matrix,
)
from bijux_phylogenetics.diversification import (
    compare_diversification_models,
    compute_diversification_gamma_statistic,
    compute_lineage_through_time_curve,
    detect_diversification_outlier_clades,
    detect_incomplete_taxon_sampling_metadata,
    estimate_diversification_rate,
    inspect_diversification_time_tree,
    render_diversification_report,
    run_trait_dependent_diversification_analysis,
    validate_time_tree_for_diversification,
    write_clade_diversification_table,
    write_diversification_gamma_statistic_table,
    write_lineage_through_time_table,
    write_trait_dependent_diversification_table,
)
from bijux_phylogenetics.ecological_niche import (
    NicheStateNodeRow,
    NicheTransitionBranchRow,
    NicheTransitionCladeRow,
    NicheTransitionCountRow,
    NicheTransitionExclusionRow,
    NicheTransitionRateRow,
    NicheTransitionReport,
    NicheTransitionSummary,
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
    audit_alignment_inference_readiness,
    build_inference_comparison_conclusion_rows,
    build_inference_comparison_shared_clade_rows,
    build_inference_comparison_weighted_conflict_rows,
    classify_inference_workflow_failure,
    compare_fast_and_ml_trees,
    compare_inferred_tree_to_taxon_metadata,
    export_workflow_result_bundle,
    infer_unaligned_sequence_type,
    render_inference_workflow_report,
    rewrite_inference_comparison_report_html,
    run_alignment_trimming,
    run_bootstrap_consensus_tree,
    run_bootstrap_support_estimation,
    run_codon_aware_multiple_sequence_alignment,
    run_fast_tree_inference,
    run_fasta_to_tree_workflow,
    run_inference_reproducibility_check,
    run_large_alignment_inference,
    run_maximum_likelihood_tree_inference,
    run_model_selection,
    run_multiple_sequence_alignment,
    run_sh_alrt_support_estimation,
    run_tree_inference_comparison,
    validate_bootstrap_tree_set,
    validate_inference_engine_outputs,
    validate_ml_tree_contains_expected_taxa,
    validate_model_selection_against_engine_outputs,
    validate_workflow_result_bundle,
    write_inference_comparison_clade_table,
    write_inference_comparison_conclusion_table,
    write_inference_comparison_summary_table,
    write_inference_comparison_weighted_conflict_table,
)
from bijux_phylogenetics.errors import (
    AlignmentTaxonMismatchError,
    DuplicateTaxonError,
    EngineUnavailableError,
    InvalidAlignmentError,
    InvalidBranchLengthError,
    InvalidDistanceMatrixError,
    MetadataJoinError,
    NonUltrametricTreeError,
    TreeParseError,
    TreeRootingError,
    UnnamedTipError,
    UnrootedTreeError,
    UnsupportedTreeFormatError,
)
from bijux_phylogenetics.evidence.bundles import (
    bundle_directory,
    bundle_file_paths,
    validate_bundle,
)
from bijux_phylogenetics.host_association import (
    HostStateNodeRow,
    HostSwitchBranchRow,
    HostSwitchCountRow,
    HostSwitchExclusionRow,
    HostSwitchFitRow,
    HostSwitchingReport,
    HostSwitchSummary,
    UnsupportedHostSwitchClaimRow,
    summarize_host_switching,
    write_host_state_node_table,
    write_host_switch_branch_table,
    write_host_switch_count_table,
    write_host_switch_exclusion_table,
    write_host_switch_fit_table,
    write_host_switch_summary_table,
    write_unsupported_host_switch_claim_table,
)
from bijux_phylogenetics.identity import IDENTITY
from bijux_phylogenetics.io.fasta import (
    assess_alignment_low_information,
    build_alignment_forensic_report,
    build_alignment_quality_report,
    build_ambiguous_alignment_column_report,
    build_duplicate_sequence_policy_report,
    build_sequence_quality_ranking,
    classify_alignment_sequences,
    clean_alignment_with_profile,
    compare_alignment_versions,
    compute_amino_acid_composition,
    compute_nucleotide_composition,
    compute_pairwise_sequence_identity_matrix,
    detect_composition_outlier_sequences,
    detect_fasta_sequence_type,
    detect_identical_duplicate_sequences,
    detect_invalid_alignment_characters,
    detect_near_duplicate_sequences,
    detect_over_aligned_regions,
    detect_sequence_length_outliers,
    detect_sequences_with_excessive_missing_data,
    detect_sites_with_excessive_missing_data,
    detect_under_aligned_regions,
    get_alignment_filter_profile,
    infer_alignment_alphabet,
    inspect_coding_alignment,
    inspect_coding_alignment_from_dna_bin_alignment,
    link_alignment_to_tree,
    list_alignment_filter_profiles,
    load_dna_bin_alignment,
    load_fasta_alignment,
    prepare_coding_sequences_for_alignment,
    remove_all_gap_columns,
    remove_all_missing_columns,
    remove_sequences_above_missingness_threshold,
    repair_fasta_input,
    summarise_fasta,
    summarize_alignment_readiness,
    summarize_alignment_windows,
    summarize_fasta_input,
    translate_coding_alignment,
    translate_coding_alignment_from_dna_bin_alignment,
    trim_alignment,
    trim_columns_above_missingness_threshold,
    validate_fasta_input,
    write_dna_bin_alignment_fasta,
    write_fasta_alignment,
)
from bijux_phylogenetics.io.newick import dumps_newick, loads_newick, write_newick
from bijux_phylogenetics.io.nexus import load_nexus
from bijux_phylogenetics.io.phyloxml import load_phyloxml
from bijux_phylogenetics.io.roundtrip import validate_tree_roundtrip
from bijux_phylogenetics.io.trees import detect_tree_format
from bijux_phylogenetics.core.tree import PhyloTree, TaxonLabel, TreeNode
from bijux_phylogenetics.phylogeography import (
    CoordinateEstimateRow,
    CoordinateMovementBranchRow,
    CoordinateMovementExclusionRow,
    CoordinateMovementOutlierRow,
    CoordinateMovementSummary,
    CoordinateMovementVisualization,
    GeographicMapArtifact,
    GeographicMapExclusionRow,
    GeographicMapLineRow,
    GeographicMapMarkerRow,
    GeographicMapReport,
    GeographicMapSummary,
    PhylogeographicCoordinateReport,
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
from bijux_phylogenetics.reference_parity import (
    validate_reference_parity_examples,
    write_reference_parity_observation_table,
    write_reference_parity_summary_table,
)
from bijux_phylogenetics.ape_parity import (
    list_ape_parity_cases,
    run_ape_parity_cases,
    write_ape_parity_observation_table,
    write_ape_parity_summary_table,
)
from bijux_phylogenetics.phytools_parity import (
    list_phytools_parity_cases,
    run_phytools_parity_cases,
    write_phytools_parity_observation_table,
    write_phytools_parity_summary_table,
)
from bijux_phylogenetics.reference_validation import (
    build_core_workflow_failure_gallery,
    build_core_workflow_validation_report,
    build_level_one_release_gate_report,
    classify_core_workflow_maturity,
    validate_alignment_quality_reference_fixtures,
    validate_dataset_audit_reference_fixtures,
    validate_figure_reference_fixtures,
    validate_report_regression_fixtures,
    validate_taxon_naming_reference_fixtures,
    validate_tree_reference_fixtures,
)
from bijux_phylogenetics.release_truth import build_release_truth_report
from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.render.package import build_tree_figure_package
from bijux_phylogenetics.render.svg import AnnotationStrip, render_tree_svg
from bijux_phylogenetics.reports.service import (
    annotate_tree_against_table,
    distance_method_limitations,
    render_alignment_report,
    render_dataset_report,
    render_distance_report,
    render_level_one_release_gate_report,
    render_phylo_inputs_report,
    render_phylogenetics_report,
    render_release_truth_report,
    render_taxon_report,
    render_tree_report,
    render_tree_set_comparison_report,
    render_tree_uncertainty_report,
    render_workflow_validation_report,
)
from bijux_phylogenetics.simulation import (
    simulate_birth_death_trees,
    simulate_brownian_traits,
    simulate_coalescent_tree,
    simulate_coalescent_trees,
    simulate_discrete_traits,
    simulate_dna_alignment,
    simulate_early_burst_traits,
    simulate_ou_traits,
    simulate_random_tree,
    simulate_random_trees,
    simulate_protein_alignment,
    write_continuous_trait_table,
    write_discrete_trait_table,
    write_simulated_alignment,
    write_tree_simulation_envelope_table,
    write_tree_simulation_record_table,
    write_tree_set,
)
from bijux_phylogenetics.tree_set import (
    BootstrapTreeSetArtifactReport,
    BootstrapTreeSetSummaryReport,
    BootstrapUnstableBranch,
    cluster_trees_by_topology,
    compare_bootstrap_and_posterior_uncertainty,
    compare_posterior_topological_diversity,
    compare_posterior_tree_sets,
    compute_clade_frequency_table,
    compute_consensus_tree,
    compute_reference_tree_clade_support,
    compute_strict_consensus_tree,
    compute_tree_distance_matrix,
    detect_posterior_topology_multimodality,
    detect_unstable_clades,
    detect_unstable_taxa,
    load_tree_set,
    summarize_bootstrap_tree_set,
    summarize_clade_credibility_conflicts,
    summarize_uncertainty_aware_conclusions,
    write_bootstrap_tree_set_artifacts,
    write_bootstrap_tree_set_summary_table,
    write_bootstrap_unstable_branch_table,
    write_clade_credibility_conflict_table,
    write_reference_tree_clade_support_table,
    write_topology_cluster_table,
    write_uncertainty_conclusion_table,
)
from bijux_phylogenetics.tree_shape import (
    TreeShapeAggregate,
    TreeShapeReport,
    TreeShapeRow,
    summarize_tree_set_shapes,
    summarize_tree_shape,
    write_tree_shape_table,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")
REPO_ROOT = Path(__file__).resolve().parents[3]


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def _write_junit_report(
    path: Path,
    *,
    suite_name: str,
    tests: int,
    failures: int,
    skipped: int,
    errors: int = 0,
) -> Path:
    path.write_text(
        (
            "<testsuites name=\"pytest\">"
            f"<testsuite name=\"{suite_name}\" tests=\"{tests}\" failures=\"{failures}\" errors=\"{errors}\" skipped=\"{skipped}\" />"
            "</testsuites>\n"
        ),
        encoding="utf-8",
    )
    return path


def _load_robinson_foulds_reference_rows() -> list[dict[str, str]]:
    with fixture("robinson_foulds_reference.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _load_branch_score_reference_rows() -> list[dict[str, str]]:
    with fixture("branch_score_reference.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_package_identity_matches_canonical_names() -> None:
    assert bijux_phylogenetics.__name__ == "bijux_phylogenetics"
    assert IDENTITY.package_name == "bijux-phylogenetics"
    assert IDENTITY.import_name == "bijux_phylogenetics"
    assert IDENTITY.cli_name == "bijux-phylogenetics"
    assert IDENTITY.umbrella_command == "bijux phylogenetics"
    assert IDENTITY.cli_aliases == ("bijux phylo",)


def test_public_package_exports_alignment_and_topology_workflows() -> None:
    assert bijux_phylogenetics.summarise_fasta is summarise_fasta
    assert (
        bijux_phylogenetics.build_alignment_quality_report
        is build_alignment_quality_report
    )
    assert (
        bijux_phylogenetics.build_alignment_forensic_report
        is build_alignment_forensic_report
    )
    assert (
        bijux_phylogenetics.classify_alignment_sequences is classify_alignment_sequences
    )
    assert (
        bijux_phylogenetics.clean_alignment_with_profile is clean_alignment_with_profile
    )
    assert bijux_phylogenetics.compare_alignment_versions is compare_alignment_versions
    assert (
        bijux_phylogenetics.compute_pairwise_genetic_distance_matrix
        is compute_pairwise_genetic_distance_matrix
    )
    assert (
        bijux_phylogenetics.build_distance_method_report is build_distance_method_report
    )
    assert bijux_phylogenetics.build_distance_tree is build_distance_tree
    assert (
        bijux_phylogenetics.build_distance_tree_from_genetic_distance_matrix
        is build_distance_tree_from_genetic_distance_matrix
    )
    assert (
        bijux_phylogenetics.build_tree_from_imported_distance_matrix
        is build_tree_from_imported_distance_matrix
    )
    assert (
        bijux_phylogenetics.compare_distance_gap_policies
        is compare_distance_gap_policies
    )
    assert bijux_phylogenetics.compare_distance_models is compare_distance_models
    assert (
        bijux_phylogenetics.compare_distance_tree_to_reference_tree
        is compare_distance_tree_to_reference_tree
    )
    assert (
        bijux_phylogenetics.compare_distance_tree_topologies
        is compare_distance_tree_topologies
    )
    assert (
        bijux_phylogenetics.list_distance_tree_method_policies
        is list_distance_tree_method_policies
    )
    assert (
        bijux_phylogenetics.resolve_distance_tree_method_policy
        is resolve_distance_tree_method_policy
    )
    assert (
        bijux_phylogenetics.BranchScoreComparisonReport is BranchScoreComparisonReport
    )
    assert (
        bijux_phylogenetics.CladeOverlapComparisonReport is CladeOverlapComparisonReport
    )
    assert bijux_phylogenetics.CladeMetadataObservation is CladeMetadataObservation
    assert bijux_phylogenetics.BranchLengthAggregate is BranchLengthAggregate
    assert (
        bijux_phylogenetics.BranchLengthDistributionReport
        is BranchLengthDistributionReport
    )
    assert bijux_phylogenetics.BranchLengthRow is BranchLengthRow
    assert bijux_phylogenetics.CladeTableReport is CladeTableReport
    assert bijux_phylogenetics.CladeTableRow is CladeTableRow
    assert bijux_phylogenetics.TreeShapeAggregate is TreeShapeAggregate
    assert bijux_phylogenetics.TreeShapeReport is TreeShapeReport
    assert bijux_phylogenetics.TreeShapeRow is TreeShapeRow
    assert bijux_phylogenetics.SharedTaxaPruningReport is SharedTaxaPruningReport
    assert bijux_phylogenetics.SupportComparisonReport is SupportComparisonReport
    assert bijux_phylogenetics.SupportConflictRow is SupportConflictRow
    assert bijux_phylogenetics.TaxonInfluenceReport is TaxonInfluenceReport
    assert bijux_phylogenetics.TaxonInfluenceRow is TaxonInfluenceRow
    assert (
        bijux_phylogenetics.TreeDistanceReferenceObservation
        is TreeDistanceReferenceObservation
    )
    assert (
        bijux_phylogenetics.TreeDistanceReferenceValidationReport
        is TreeDistanceReferenceValidationReport
    )
    assert (
        bijux_phylogenetics.compare_branch_score_distance
        is compare_branch_score_distance
    )
    assert bijux_phylogenetics.compare_clade_overlap is compare_clade_overlap
    assert bijux_phylogenetics.analyze_taxon_influence is analyze_taxon_influence
    assert (
        bijux_phylogenetics.analyze_branch_length_distribution
        is analyze_branch_length_distribution
    )
    assert (
        bijux_phylogenetics.analyze_tree_set_branch_lengths
        is analyze_tree_set_branch_lengths
    )
    assert bijux_phylogenetics.extract_tree_clades is extract_tree_clades
    assert bijux_phylogenetics.extract_tree_set_clades is extract_tree_set_clades
    assert bijux_phylogenetics.summarize_tree_shape is summarize_tree_shape
    assert bijux_phylogenetics.summarize_tree_set_shapes is summarize_tree_set_shapes
    assert bijux_phylogenetics.extract_tree_clade_by_node_id is extract_tree_clade_by_node_id
    assert (
        bijux_phylogenetics.extract_tree_clade_by_descendant_taxa
        is extract_tree_clade_by_descendant_taxa
    )
    assert bijux_phylogenetics.assess_tree_monophyly is assess_tree_monophyly
    assert bijux_phylogenetics.find_tree_mrca is find_tree_mrca
    assert bijux_phylogenetics.TreeBranchingTimeReport is TreeBranchingTimeReport
    assert bijux_phylogenetics.TreeBranchingTimeRow is TreeBranchingTimeRow
    assert (
        bijux_phylogenetics.compute_tree_branching_times
        is compute_tree_branching_times
    )
    assert (
        bijux_phylogenetics.write_tree_branching_time_table
        is write_tree_branching_time_table
    )
    assert bijux_phylogenetics.PhyloTree is PhyloTree
    assert bijux_phylogenetics.TaxonLabel is TaxonLabel
    assert bijux_phylogenetics.TreeNode is TreeNode
    assert bijux_phylogenetics.TreeNodeDepthReport is TreeNodeDepthReport
    assert bijux_phylogenetics.TreeNodeDepthRow is TreeNodeDepthRow
    assert bijux_phylogenetics.compute_tree_node_depths is compute_tree_node_depths
    assert (
        bijux_phylogenetics.write_tree_node_depth_table
        is write_tree_node_depth_table
    )
    assert (
        bijux_phylogenetics.APE_ULTRAMETRIC_TOLERANCE
        == APE_ULTRAMETRIC_TOLERANCE
    )
    assert bijux_phylogenetics.TreeUltrametricReport is TreeUltrametricReport
    assert bijux_phylogenetics.TreeUltrametricTipRow is TreeUltrametricTipRow
    assert (
        bijux_phylogenetics.assess_tree_ultrametricity
        is assess_tree_ultrametricity
    )
    assert (
        bijux_phylogenetics.write_tree_ultrametric_table
        is write_tree_ultrametric_table
    )
    assert bijux_phylogenetics.TipDistanceMatrixReport is TipDistanceMatrixReport
    assert bijux_phylogenetics.TipDistanceMatrixRow is TipDistanceMatrixRow
    assert (
        bijux_phylogenetics.compute_tree_tip_distance_matrix
        is compute_tree_tip_distance_matrix
    )
    assert (
        bijux_phylogenetics.write_tree_tip_distance_matrix
        is write_tree_tip_distance_matrix
    )
    assert (
        bijux_phylogenetics.write_tree_tip_distance_long_table
        is write_tree_tip_distance_long_table
    )
    assert bijux_phylogenetics.prune_trees_to_shared_taxa is prune_trees_to_shared_taxa
    assert bijux_phylogenetics.compare_support_values is compare_support_values
    assert bijux_phylogenetics.write_branch_length_table is write_branch_length_table
    assert bijux_phylogenetics.write_clade_table is write_clade_table
    assert bijux_phylogenetics.write_tree_shape_table is write_tree_shape_table
    assert bijux_phylogenetics.write_clade_overlap_table is write_clade_overlap_table
    assert (
        bijux_phylogenetics.write_shared_taxa_pruning_table
        is write_shared_taxa_pruning_table
    )
    assert (
        bijux_phylogenetics.write_shared_taxa_removed_taxa_table
        is write_shared_taxa_removed_taxa_table
    )
    assert (
        bijux_phylogenetics.write_support_comparison_table
        is write_support_comparison_table
    )
    assert (
        bijux_phylogenetics.write_taxon_influence_table is write_taxon_influence_table
    )
    assert (
        bijux_phylogenetics.RobinsonFouldsComparisonReport
        is RobinsonFouldsComparisonReport
    )
    assert compare_api.CladeOverlapComparisonReport is CladeOverlapComparisonReport
    assert compare_api.SharedTaxaPruningReport is SharedTaxaPruningReport
    assert compare_api.SupportComparisonReport is SupportComparisonReport
    assert compare_api.SupportConflictRow is SupportConflictRow
    assert compare_api.TaxonInfluenceReport is TaxonInfluenceReport
    assert compare_api.TaxonInfluenceRow is TaxonInfluenceRow
    assert (
        compare_api.TreeDistanceReferenceObservation
        is TreeDistanceReferenceObservation
    )
    assert (
        compare_api.TreeDistanceReferenceValidationReport
        is TreeDistanceReferenceValidationReport
    )
    assert compare_api.compare_clade_overlap is compare_clade_overlap
    assert compare_api.prune_trees_to_shared_taxa is prune_trees_to_shared_taxa
    assert compare_api.compare_support_values is compare_support_values
    assert compare_api.analyze_taxon_influence is analyze_taxon_influence
    assert compare_api.write_clade_overlap_table is write_clade_overlap_table
    assert (
        compare_api.write_shared_taxa_pruning_table is write_shared_taxa_pruning_table
    )
    assert (
        compare_api.write_shared_taxa_removed_taxa_table
        is write_shared_taxa_removed_taxa_table
    )
    assert compare_api.write_support_comparison_table is write_support_comparison_table
    assert compare_api.write_taxon_influence_table is write_taxon_influence_table
    assert compare_api.TopologyDistanceReport is TopologyDistanceReport
    assert compare_api.TopologyDistanceSplitRow is TopologyDistanceSplitRow
    assert compare_api.compare_topology_distance is compare_topology_distance
    assert (
        compare_api.write_topology_distance_split_table
        is write_topology_distance_split_table
    )
    assert compare_api.BranchScoreComparisonReport is BranchScoreComparisonReport
    assert compare_api.compare_branch_score_distance is compare_branch_score_distance
    assert bijux_phylogenetics.TopologyDistanceReport is TopologyDistanceReport
    assert bijux_phylogenetics.TopologyDistanceSplitRow is TopologyDistanceSplitRow
    assert bijux_phylogenetics.compare_topology_distance is compare_topology_distance
    assert (
        bijux_phylogenetics.write_topology_distance_split_table
        is write_topology_distance_split_table
    )
    assert bijux_phylogenetics.compare_robinson_foulds is compare_robinson_foulds
    assert compare_api.RobinsonFouldsComparisonReport is RobinsonFouldsComparisonReport
    assert compare_api.compare_robinson_foulds is compare_robinson_foulds
    assert (
        bijux_phylogenetics.validate_support_reference_examples
        is validate_support_reference_examples
    )
    assert (
        compare_api.validate_support_reference_examples
        is validate_support_reference_examples
    )
    assert (
        bijux_phylogenetics.validate_tree_distance_reference_examples
        is validate_tree_distance_reference_examples
    )
    assert (
        compare_api.validate_tree_distance_reference_examples
        is validate_tree_distance_reference_examples
    )
    assert (
        bijux_phylogenetics.SupportReferenceObservation
        is SupportReferenceObservation
    )
    assert (
        bijux_phylogenetics.SupportReferenceValidationReport
        is SupportReferenceValidationReport
    )
    assert compare_api.SupportReferenceObservation is SupportReferenceObservation
    assert (
        compare_api.SupportReferenceValidationReport
        is SupportReferenceValidationReport
    )
    assert (
        bijux_phylogenetics.assess_distance_method_maturity
        is assess_distance_method_maturity
    )
    assert (
        bijux_phylogenetics.summarize_distance_bootstrap_support
        is summarize_distance_bootstrap_support
    )
    assert (
        bijux_phylogenetics.validate_imported_distance_matrix
        is validate_imported_distance_matrix
    )
    assert (
        bijux_phylogenetics.validate_discrete_state_coding
        is validate_discrete_state_coding
    )
    assert (
        bijux_phylogenetics.detect_state_imbalance_problems
        is detect_state_imbalance_problems
    )
    assert (
        bijux_phylogenetics.run_discrete_state_transition_model
        is run_discrete_state_transition_model
    )
    assert (
        bijux_phylogenetics.estimate_ancestral_geographic_states
        is estimate_ancestral_geographic_states
    )
    assert (
        bijux_phylogenetics.compare_discrete_state_models
        is compare_discrete_state_models
    )
    assert (
        bijux_phylogenetics.write_discrete_model_comparison_table
        is write_discrete_model_comparison_table
    )
    assert (
        bijux_phylogenetics.write_node_state_probability_table
        is write_node_state_probability_table
    )
    assert (
        bijux_phylogenetics.write_transition_summary_table
        is write_transition_summary_table
    )
    assert (
        bijux_phylogenetics.summarize_constrained_geographic_model
        is summarize_constrained_geographic_model
    )
    assert (
        bijux_phylogenetics.summarize_constrained_geographic_report
        is summarize_constrained_geographic_report
    )
    assert (
        bijux_phylogenetics.summarize_biogeographic_transition_chronology
        is summarize_biogeographic_transition_chronology
    )
    assert (
        bijux_phylogenetics.summarize_geographic_sampling_bias
        is summarize_geographic_sampling_bias
    )
    assert bijux_phylogenetics.BiogeographyRegionCountRow is BiogeographyRegionCountRow
    assert (
        bijux_phylogenetics.BiogeographyReportExclusionRow
        is BiogeographyReportExclusionRow
    )
    assert (
        bijux_phylogenetics.BiogeographyReportPackageResult
        is BiogeographyReportPackageResult
    )
    assert (
        bijux_phylogenetics.summarize_biogeography_region_counts
        is summarize_biogeography_region_counts
    )
    assert (
        bijux_phylogenetics.summarize_geographic_state_model
        is summarize_geographic_state_model
    )
    assert (
        bijux_phylogenetics.build_biogeography_report_package
        is build_biogeography_report_package
    )
    assert (
        bijux_phylogenetics.summarize_geographic_migration_events
        is summarize_geographic_migration_events
    )
    assert (
        bijux_phylogenetics.summarize_geographic_migration_event_tree_set
        is summarize_geographic_migration_event_tree_set
    )
    assert (
        bijux_phylogenetics.summarize_time_stratified_geographic_transitions
        is summarize_time_stratified_geographic_transitions
    )
    assert (
        bijux_phylogenetics.write_constrained_geographic_summary_table
        is write_constrained_geographic_summary_table
    )
    assert (
        bijux_phylogenetics.write_constrained_geographic_fit_table
        is write_constrained_geographic_fit_table
    )
    assert (
        bijux_phylogenetics.write_constrained_geographic_transition_table
        is write_constrained_geographic_transition_table
    )
    assert (
        bijux_phylogenetics.write_biogeography_region_count_table
        is write_biogeography_region_count_table
    )
    assert (
        bijux_phylogenetics.write_biogeography_report_exclusion_table
        is write_biogeography_report_exclusion_table
    )
    assert (
        bijux_phylogenetics.write_constrained_geographic_exclusion_table
        is write_constrained_geographic_exclusion_table
    )
    assert (
        bijux_phylogenetics.write_dated_biogeography_summary_table
        is write_dated_biogeography_summary_table
    )
    assert (
        bijux_phylogenetics.write_dated_biogeography_node_table
        is write_dated_biogeography_node_table
    )
    assert (
        bijux_phylogenetics.write_dated_biogeography_event_table
        is write_dated_biogeography_event_table
    )
    assert (
        bijux_phylogenetics.write_dated_biogeography_time_bin_table
        is write_dated_biogeography_time_bin_table
    )
    assert (
        bijux_phylogenetics.write_dated_biogeography_exclusion_table
        is write_dated_biogeography_exclusion_table
    )
    assert (
        bijux_phylogenetics.write_geographic_sampling_bias_summary_table
        is write_geographic_sampling_bias_summary_table
    )
    assert (
        bijux_phylogenetics.write_geographic_sampling_count_table
        is write_geographic_sampling_count_table
    )
    assert (
        bijux_phylogenetics.write_geographic_sampling_bias_node_table
        is write_geographic_sampling_bias_node_table
    )
    assert (
        bijux_phylogenetics.write_geographic_sampling_bias_transition_table
        is write_geographic_sampling_bias_transition_table
    )
    assert (
        bijux_phylogenetics.write_geographic_sampling_bias_exclusion_table
        is write_geographic_sampling_bias_exclusion_table
    )
    assert (
        bijux_phylogenetics.write_unsupported_geographic_transition_claim_table
        is write_unsupported_geographic_transition_claim_table
    )
    assert (
        bijux_phylogenetics.write_geographic_state_summary_table
        is write_geographic_state_summary_table
    )
    assert (
        bijux_phylogenetics.write_geographic_region_probability_table
        is write_geographic_region_probability_table
    )
    assert (
        bijux_phylogenetics.write_geographic_transition_rate_table
        is write_geographic_transition_rate_table
    )
    assert (
        bijux_phylogenetics.write_geographic_transition_event_table
        is write_geographic_transition_event_table
    )
    assert (
        bijux_phylogenetics.write_geographic_exclusion_table
        is write_geographic_exclusion_table
    )
    assert (
        bijux_phylogenetics.write_geographic_migration_event_summary_table
        is write_geographic_migration_event_summary_table
    )
    assert (
        bijux_phylogenetics.write_geographic_migration_event_table
        is write_geographic_migration_event_table
    )
    assert (
        bijux_phylogenetics.write_geographic_migration_exclusion_table
        is write_geographic_migration_exclusion_table
    )
    assert (
        bijux_phylogenetics.write_geographic_migration_tree_set_summary_table
        is write_geographic_migration_tree_set_summary_table
    )
    assert (
        bijux_phylogenetics.write_geographic_migration_tree_set_tree_table
        is write_geographic_migration_tree_set_tree_table
    )
    assert (
        bijux_phylogenetics.write_geographic_migration_tree_set_event_table
        is write_geographic_migration_tree_set_event_table
    )
    assert (
        bijux_phylogenetics.write_geographic_migration_tree_set_event_summary_table
        is write_geographic_migration_tree_set_event_summary_table
    )
    assert (
        bijux_phylogenetics.write_geographic_migration_tree_set_exclusion_table
        is write_geographic_migration_tree_set_exclusion_table
    )
    assert bijux_phylogenetics.TimeBinDefinition is TimeBinDefinition
    assert (
        bijux_phylogenetics.GeographicMigrationEventReport
        is GeographicMigrationEventReport
    )
    assert (
        bijux_phylogenetics.GeographicMigrationEventRow is GeographicMigrationEventRow
    )
    assert (
        bijux_phylogenetics.GeographicMigrationEventSummary
        is GeographicMigrationEventSummary
    )
    assert bijux_phylogenetics.GeographicMigrationTreeRow is GeographicMigrationTreeRow
    assert (
        bijux_phylogenetics.GeographicMigrationTreeSetEventRow
        is GeographicMigrationTreeSetEventRow
    )
    assert (
        bijux_phylogenetics.GeographicMigrationTreeSetEventSummaryRow
        is GeographicMigrationTreeSetEventSummaryRow
    )
    assert (
        bijux_phylogenetics.GeographicMigrationTreeSetReport
        is GeographicMigrationTreeSetReport
    )
    assert (
        bijux_phylogenetics.GeographicMigrationTreeSetSummary
        is GeographicMigrationTreeSetSummary
    )
    assert bijux_phylogenetics.DatedBiogeographyEventRow is DatedBiogeographyEventRow
    assert bijux_phylogenetics.DatedBiogeographyNodeRow is DatedBiogeographyNodeRow
    assert bijux_phylogenetics.DatedBiogeographyReport is DatedBiogeographyReport
    assert bijux_phylogenetics.DatedBiogeographySummary is DatedBiogeographySummary
    assert (
        bijux_phylogenetics.DatedBiogeographyTimeBinRow is DatedBiogeographyTimeBinRow
    )
    assert (
        bijux_phylogenetics.GeographicSamplingBiasNodeRow
        is GeographicSamplingBiasNodeRow
    )
    assert (
        bijux_phylogenetics.GeographicSamplingBiasReport is GeographicSamplingBiasReport
    )
    assert (
        bijux_phylogenetics.GeographicSamplingBiasSummary
        is GeographicSamplingBiasSummary
    )
    assert (
        bijux_phylogenetics.GeographicSamplingBiasTransitionRow
        is GeographicSamplingBiasTransitionRow
    )
    assert bijux_phylogenetics.GeographicSamplingCountRow is GeographicSamplingCountRow
    assert (
        bijux_phylogenetics.ConstrainedGeographicFitRow is ConstrainedGeographicFitRow
    )
    assert (
        bijux_phylogenetics.ConstrainedGeographicReport is ConstrainedGeographicReport
    )
    assert (
        bijux_phylogenetics.ConstrainedGeographicSummary is ConstrainedGeographicSummary
    )
    assert (
        bijux_phylogenetics.ConstrainedGeographicTransitionRow
        is ConstrainedGeographicTransitionRow
    )
    assert bijux_phylogenetics.TimeStratifiedBranchRow is TimeStratifiedBranchRow
    assert (
        bijux_phylogenetics.TimeStratifiedTransitionMatrixRow
        is TimeStratifiedTransitionMatrixRow
    )
    assert (
        bijux_phylogenetics.TimeStratifiedTransitionReport
        is TimeStratifiedTransitionReport
    )
    assert (
        bijux_phylogenetics.UnsupportedGeographicTransitionClaimRow
        is UnsupportedGeographicTransitionClaimRow
    )
    assert (
        bijux_phylogenetics.TimeStratifiedTransitionSummary
        is TimeStratifiedTransitionSummary
    )
    assert (
        bijux_phylogenetics.write_time_stratified_transition_summary_table
        is write_time_stratified_transition_summary_table
    )
    assert (
        bijux_phylogenetics.write_time_stratified_transition_matrix_table
        is write_time_stratified_transition_matrix_table
    )
    assert (
        bijux_phylogenetics.write_time_stratified_branch_table
        is write_time_stratified_branch_table
    )
    assert (
        bijux_phylogenetics.write_time_stratified_exclusion_table
        is write_time_stratified_exclusion_table
    )
    assert bijux_phylogenetics.summarize_host_switching is summarize_host_switching
    assert bijux_phylogenetics.HostStateNodeRow is HostStateNodeRow
    assert bijux_phylogenetics.HostSwitchBranchRow is HostSwitchBranchRow
    assert bijux_phylogenetics.HostSwitchCountRow is HostSwitchCountRow
    assert bijux_phylogenetics.HostSwitchExclusionRow is HostSwitchExclusionRow
    assert bijux_phylogenetics.HostSwitchFitRow is HostSwitchFitRow
    assert bijux_phylogenetics.HostSwitchSummary is HostSwitchSummary
    assert bijux_phylogenetics.HostSwitchingReport is HostSwitchingReport
    assert (
        bijux_phylogenetics.UnsupportedHostSwitchClaimRow
        is UnsupportedHostSwitchClaimRow
    )
    assert (
        bijux_phylogenetics.write_host_switch_summary_table
        is write_host_switch_summary_table
    )
    assert (
        bijux_phylogenetics.write_host_state_node_table is write_host_state_node_table
    )
    assert (
        bijux_phylogenetics.write_host_switch_branch_table
        is write_host_switch_branch_table
    )
    assert (
        bijux_phylogenetics.write_host_switch_count_table
        is write_host_switch_count_table
    )
    assert (
        bijux_phylogenetics.write_host_switch_fit_table is write_host_switch_fit_table
    )
    assert (
        bijux_phylogenetics.write_unsupported_host_switch_claim_table
        is write_unsupported_host_switch_claim_table
    )
    assert (
        bijux_phylogenetics.write_host_switch_exclusion_table
        is write_host_switch_exclusion_table
    )
    assert (
        bijux_phylogenetics.summarize_niche_transitions is summarize_niche_transitions
    )
    assert bijux_phylogenetics.NicheStateNodeRow is NicheStateNodeRow
    assert bijux_phylogenetics.NicheTransitionBranchRow is NicheTransitionBranchRow
    assert bijux_phylogenetics.NicheTransitionCladeRow is NicheTransitionCladeRow
    assert bijux_phylogenetics.NicheTransitionCountRow is NicheTransitionCountRow
    assert (
        bijux_phylogenetics.NicheTransitionExclusionRow is NicheTransitionExclusionRow
    )
    assert bijux_phylogenetics.NicheTransitionRateRow is NicheTransitionRateRow
    assert bijux_phylogenetics.NicheTransitionReport is NicheTransitionReport
    assert bijux_phylogenetics.NicheTransitionSummary is NicheTransitionSummary
    assert (
        bijux_phylogenetics.write_niche_transition_summary_table
        is write_niche_transition_summary_table
    )
    assert (
        bijux_phylogenetics.write_niche_state_node_table is write_niche_state_node_table
    )
    assert (
        bijux_phylogenetics.write_niche_transition_rate_table
        is write_niche_transition_rate_table
    )
    assert (
        bijux_phylogenetics.write_niche_transition_branch_table
        is write_niche_transition_branch_table
    )
    assert (
        bijux_phylogenetics.write_niche_transition_count_table
        is write_niche_transition_count_table
    )
    assert (
        bijux_phylogenetics.write_niche_transition_clade_table
        is write_niche_transition_clade_table
    )
    assert (
        bijux_phylogenetics.write_niche_transition_exclusion_table
        is write_niche_transition_exclusion_table
    )
    assert (
        bijux_phylogenetics.summarize_continuous_phylogeography
        is summarize_continuous_phylogeography
    )
    assert bijux_phylogenetics.CoordinateEstimateRow is CoordinateEstimateRow
    assert (
        bijux_phylogenetics.CoordinateMovementBranchRow is CoordinateMovementBranchRow
    )
    assert (
        bijux_phylogenetics.CoordinateMovementExclusionRow
        is CoordinateMovementExclusionRow
    )
    assert (
        bijux_phylogenetics.CoordinateMovementOutlierRow is CoordinateMovementOutlierRow
    )
    assert bijux_phylogenetics.CoordinateMovementSummary is CoordinateMovementSummary
    assert (
        bijux_phylogenetics.CoordinateMovementVisualization
        is CoordinateMovementVisualization
    )
    assert (
        bijux_phylogenetics.PhylogeographicCoordinateReport
        is PhylogeographicCoordinateReport
    )
    assert (
        bijux_phylogenetics.write_coordinate_movement_summary_table
        is write_coordinate_movement_summary_table
    )
    assert (
        bijux_phylogenetics.write_coordinate_estimate_table
        is write_coordinate_estimate_table
    )
    assert (
        bijux_phylogenetics.write_coordinate_movement_branch_table
        is write_coordinate_movement_branch_table
    )
    assert (
        bijux_phylogenetics.write_coordinate_movement_outlier_table
        is write_coordinate_movement_outlier_table
    )
    assert (
        bijux_phylogenetics.write_coordinate_movement_exclusion_table
        is write_coordinate_movement_exclusion_table
    )
    assert (
        bijux_phylogenetics.render_coordinate_movement_visualization
        is render_coordinate_movement_visualization
    )
    assert (
        bijux_phylogenetics.summarize_continuous_phylogeography_map
        is summarize_continuous_phylogeography_map
    )
    assert (
        bijux_phylogenetics.summarize_discrete_region_map
        is summarize_discrete_region_map
    )
    assert bijux_phylogenetics.GeographicMapArtifact is GeographicMapArtifact
    assert bijux_phylogenetics.GeographicMapExclusionRow is GeographicMapExclusionRow
    assert bijux_phylogenetics.GeographicMapLineRow is GeographicMapLineRow
    assert bijux_phylogenetics.GeographicMapMarkerRow is GeographicMapMarkerRow
    assert bijux_phylogenetics.GeographicMapReport is GeographicMapReport
    assert bijux_phylogenetics.GeographicMapSummary is GeographicMapSummary
    assert (
        bijux_phylogenetics.write_geographic_map_summary_table
        is write_geographic_map_summary_table
    )
    assert (
        bijux_phylogenetics.write_geographic_map_marker_table
        is write_geographic_map_marker_table
    )
    assert (
        bijux_phylogenetics.write_geographic_map_line_table
        is write_geographic_map_line_table
    )
    assert (
        bijux_phylogenetics.write_geographic_map_exclusion_table
        is write_geographic_map_exclusion_table
    )
    assert bijux_phylogenetics.render_geographic_map_html is render_geographic_map_html
    assert (
        bijux_phylogenetics.simulate_discrete_stochastic_maps
        is simulate_discrete_stochastic_maps
    )
    assert (
        bijux_phylogenetics.simulate_discrete_stochastic_maps_from_fit_report
        is simulate_discrete_stochastic_maps_from_fit_report
    )
    assert (
        bijux_phylogenetics.summarize_discrete_stochastic_maps
        is summarize_discrete_stochastic_maps
    )
    assert (
        bijux_phylogenetics.count_discrete_stochastic_map_transitions
        is count_discrete_stochastic_map_transitions
    )
    assert (
        bijux_phylogenetics.summarize_discrete_stochastic_map_density
        is summarize_discrete_stochastic_map_density
    )
    assert (
        bijux_phylogenetics.render_stochastic_map_density_artifact
        is render_stochastic_map_density_artifact
    )
    assert (
        bijux_phylogenetics.write_stochastic_map_collection
        is write_stochastic_map_collection
    )
    assert (
        bijux_phylogenetics.write_stochastic_map_branch_occupancy_table
        is write_stochastic_map_branch_occupancy_table
    )
    assert (
        bijux_phylogenetics.write_stochastic_map_branch_probability_table
        is write_stochastic_map_branch_probability_table
    )
    assert (
        bijux_phylogenetics.write_stochastic_map_transition_count_matrix
        is write_stochastic_map_transition_count_matrix
    )
    assert (
        bijux_phylogenetics.write_stochastic_map_aggregate_transition_matrix
        is write_stochastic_map_aggregate_transition_matrix
    )
    assert (
        bijux_phylogenetics.write_stochastic_map_branch_transition_count_table
        is write_stochastic_map_branch_transition_count_table
    )
    assert (
        bijux_phylogenetics.write_stochastic_map_density_branch_table
        is write_stochastic_map_density_branch_table
    )
    assert (
        bijux_phylogenetics.write_stochastic_map_density_slice_table
        is write_stochastic_map_density_slice_table
    )
    assert (
        bijux_phylogenetics.write_stochastic_map_event_table
        is write_stochastic_map_event_table
    )
    assert (
        bijux_phylogenetics.write_stochastic_map_summary_table
        is write_stochastic_map_summary_table
    )
    assert (
        bijux_phylogenetics.write_stochastic_map_state_time_table
        is write_stochastic_map_state_time_table
    )
    assert (
        bijux_phylogenetics.write_stochastic_map_segment_table
        is write_stochastic_map_segment_table
    )
    assert (
        bijux_phylogenetics.load_stochastic_map_collection
        is load_stochastic_map_collection
    )
    assert (
        bijux_phylogenetics.validate_time_tree_for_diversification
        is validate_time_tree_for_diversification
    )
    assert (
        bijux_phylogenetics.inspect_diversification_time_tree
        is inspect_diversification_time_tree
    )
    assert (
        bijux_phylogenetics.compute_lineage_through_time_curve
        is compute_lineage_through_time_curve
    )
    assert (
        bijux_phylogenetics.detect_incomplete_taxon_sampling_metadata
        is detect_incomplete_taxon_sampling_metadata
    )
    assert (
        bijux_phylogenetics.estimate_diversification_rate
        is estimate_diversification_rate
    )
    assert (
        bijux_phylogenetics.compute_diversification_gamma_statistic
        is compute_diversification_gamma_statistic
    )
    assert (
        bijux_phylogenetics.compare_diversification_models
        is compare_diversification_models
    )
    assert (
        bijux_phylogenetics.detect_diversification_outlier_clades
        is detect_diversification_outlier_clades
    )
    assert (
        bijux_phylogenetics.run_trait_dependent_diversification_analysis
        is run_trait_dependent_diversification_analysis
    )
    assert (
        bijux_phylogenetics.render_diversification_report
        is render_diversification_report
    )
    assert (
        bijux_phylogenetics.write_lineage_through_time_table
        is write_lineage_through_time_table
    )
    assert (
        bijux_phylogenetics.write_clade_diversification_table
        is write_clade_diversification_table
    )
    assert (
        bijux_phylogenetics.write_diversification_gamma_statistic_table
        is write_diversification_gamma_statistic_table
    )
    assert (
        bijux_phylogenetics.write_trait_dependent_diversification_table
        is write_trait_dependent_diversification_table
    )
    assert (
        bijux_phylogenetics.render_tree_with_geographic_states
        is render_tree_with_geographic_states
    )
    assert (
        bijux_phylogenetics.render_discrete_state_evolution_report
        is render_discrete_state_evolution_report
    )
    assert bijux_phylogenetics.assess_tree_assumptions is assess_tree_assumptions
    assert bijux_phylogenetics.inspect_coding_alignment is inspect_coding_alignment
    assert (
        bijux_phylogenetics.inspect_coding_alignment_from_dna_bin_alignment
        is inspect_coding_alignment_from_dna_bin_alignment
    )
    assert (
        bijux_phylogenetics.compute_pairwise_sequence_identity_matrix
        is compute_pairwise_sequence_identity_matrix
    )
    assert (
        bijux_phylogenetics.detect_sequence_length_outliers
        is detect_sequence_length_outliers
    )
    assert (
        bijux_phylogenetics.prepare_coding_sequences_for_alignment
        is prepare_coding_sequences_for_alignment
    )
    assert bijux_phylogenetics.detect_fasta_sequence_type is detect_fasta_sequence_type
    assert bijux_phylogenetics.validate_fasta_input is validate_fasta_input
    assert bijux_phylogenetics.repair_fasta_input is repair_fasta_input
    assert bijux_phylogenetics.summarize_fasta_input is summarize_fasta_input
    assert (
        bijux_phylogenetics.detect_over_aligned_regions is detect_over_aligned_regions
    )
    assert (
        bijux_phylogenetics.run_codon_aware_multiple_sequence_alignment
        is run_codon_aware_multiple_sequence_alignment
    )
    assert (
        bijux_phylogenetics.detect_under_aligned_regions is detect_under_aligned_regions
    )
    assert (
        bijux_phylogenetics.list_alignment_filter_profiles
        is list_alignment_filter_profiles
    )
    assert (
        bijux_phylogenetics.get_alignment_filter_profile is get_alignment_filter_profile
    )
    assert (
        bijux_phylogenetics.summarize_alignment_windows is summarize_alignment_windows
    )
    assert (
        bijux_phylogenetics.summarize_alignment_readiness
        is summarize_alignment_readiness
    )
    assert bijux_phylogenetics.audit_dataset_inputs is audit_dataset_inputs
    assert (
        bijux_phylogenetics.audit_dataset_taxon_ordering is audit_dataset_taxon_ordering
    )
    assert (
        bijux_phylogenetics.build_dataset_completeness_matrix
        is build_dataset_completeness_matrix
    )
    assert bijux_phylogenetics.build_dataset_crosswalk is build_dataset_crosswalk
    assert (
        bijux_phylogenetics.build_dataset_mismatch_report
        is build_dataset_mismatch_report
    )
    assert (
        bijux_phylogenetics.summarize_dataset_readiness is summarize_dataset_readiness
    )
    assert (
        bijux_phylogenetics.build_locus_occupancy_report is build_locus_occupancy_report
    )
    assert (
        bijux_phylogenetics.LocusOccupancyFilterIteration
        is LocusOccupancyFilterIteration
    )
    assert (
        bijux_phylogenetics.build_partition_summary_report
        is build_partition_summary_report
    )
    assert bijux_phylogenetics.filter_locus_occupancy is filter_locus_occupancy
    assert bijux_phylogenetics.parse_locus_partitions is parse_locus_partitions
    assert (
        bijux_phylogenetics.write_partition_summary_table
        is write_partition_summary_table
    )
    assert bijux_phylogenetics.write_locus_partitions is write_locus_partitions
    assert bijux_phylogenetics.build_taxon_audit_report is build_taxon_audit_report
    assert (
        bijux_phylogenetics.build_taxon_mapping_conflict_report
        is build_taxon_mapping_conflict_report
    )
    assert (
        bijux_phylogenetics.detect_duplicate_biological_identities
        is detect_duplicate_biological_identities
    )
    assert bijux_phylogenetics.export_tree_accepted_names is export_tree_accepted_names
    assert bijux_phylogenetics.infer_taxon_rank is infer_taxon_rank
    assert (
        bijux_phylogenetics.inspect_tree_taxon_rank_consistency
        is inspect_tree_taxon_rank_consistency
    )
    assert (
        bijux_phylogenetics.write_accepted_name_mapping is write_accepted_name_mapping
    )
    assert bijux_phylogenetics.load_tree_set is load_tree_set
    assert (
        bijux_phylogenetics.BootstrapTreeSetSummaryReport
        is BootstrapTreeSetSummaryReport
    )
    assert (
        bijux_phylogenetics.BootstrapTreeSetArtifactReport
        is BootstrapTreeSetArtifactReport
    )
    assert bijux_phylogenetics.BootstrapUnstableBranch is BootstrapUnstableBranch
    assert bijux_phylogenetics.compute_consensus_tree is compute_consensus_tree
    assert bijux_phylogenetics.compute_strict_consensus_tree is compute_strict_consensus_tree
    assert (
        bijux_phylogenetics.compute_clade_frequency_table
        is compute_clade_frequency_table
    )
    assert (
        bijux_phylogenetics.compute_reference_tree_clade_support
        is compute_reference_tree_clade_support
    )
    assert (
        bijux_phylogenetics.compute_tree_distance_matrix is compute_tree_distance_matrix
    )
    assert bijux_phylogenetics.cluster_trees_by_topology is cluster_trees_by_topology
    assert bijux_phylogenetics.detect_unstable_taxa is detect_unstable_taxa
    assert bijux_phylogenetics.detect_unstable_clades is detect_unstable_clades
    assert (
        bijux_phylogenetics.summarize_bootstrap_tree_set is summarize_bootstrap_tree_set
    )
    assert (
        bijux_phylogenetics.compare_posterior_topological_diversity
        is compare_posterior_topological_diversity
    )
    assert (
        bijux_phylogenetics.detect_posterior_topology_multimodality
        is detect_posterior_topology_multimodality
    )
    assert (
        bijux_phylogenetics.summarize_clade_credibility_conflicts
        is summarize_clade_credibility_conflicts
    )
    assert (
        bijux_phylogenetics.summarize_uncertainty_aware_conclusions
        is summarize_uncertainty_aware_conclusions
    )
    assert (
        bijux_phylogenetics.write_bootstrap_tree_set_summary_table
        is write_bootstrap_tree_set_summary_table
    )
    assert (
        bijux_phylogenetics.write_bootstrap_unstable_branch_table
        is write_bootstrap_unstable_branch_table
    )
    assert (
        bijux_phylogenetics.write_bootstrap_tree_set_artifacts
        is write_bootstrap_tree_set_artifacts
    )
    assert (
        bijux_phylogenetics.write_reference_tree_clade_support_table
        is write_reference_tree_clade_support_table
    )
    assert (
        bijux_phylogenetics.compare_posterior_tree_sets is compare_posterior_tree_sets
    )
    assert (
        bijux_phylogenetics.render_tree_uncertainty_report
        is render_tree_uncertainty_report
    )
    assert (
        bijux_phylogenetics.validate_tree_reference_fixtures
        is validate_tree_reference_fixtures
    )
    assert (
        bijux_phylogenetics.validate_taxon_naming_reference_fixtures
        is validate_taxon_naming_reference_fixtures
    )
    assert (
        bijux_phylogenetics.validate_alignment_quality_reference_fixtures
        is validate_alignment_quality_reference_fixtures
    )
    assert (
        bijux_phylogenetics.validate_dataset_audit_reference_fixtures
        is validate_dataset_audit_reference_fixtures
    )
    assert (
        bijux_phylogenetics.validate_figure_reference_fixtures
        is validate_figure_reference_fixtures
    )
    assert (
        bijux_phylogenetics.validate_report_regression_fixtures
        is validate_report_regression_fixtures
    )
    assert (
        bijux_phylogenetics.build_core_workflow_validation_report
        is build_core_workflow_validation_report
    )
    assert (
        bijux_phylogenetics.build_core_workflow_failure_gallery
        is build_core_workflow_failure_gallery
    )
    assert (
        bijux_phylogenetics.classify_core_workflow_maturity
        is classify_core_workflow_maturity
    )
    assert (
        bijux_phylogenetics.build_level_one_release_gate_report
        is build_level_one_release_gate_report
    )
    assert (
        bijux_phylogenetics.validate_reference_parity_examples
        is validate_reference_parity_examples
    )
    assert bijux_phylogenetics.list_ape_parity_cases is list_ape_parity_cases
    assert (
        bijux_phylogenetics.list_phytools_parity_cases is list_phytools_parity_cases
    )
    assert bijux_phylogenetics.run_ape_parity_cases is run_ape_parity_cases
    assert (
        bijux_phylogenetics.run_phytools_parity_cases
        is run_phytools_parity_cases
    )
    assert (
        bijux_phylogenetics.write_reference_parity_summary_table
        is write_reference_parity_summary_table
    )
    assert (
        bijux_phylogenetics.write_ape_parity_summary_table
        is write_ape_parity_summary_table
    )
    assert (
        bijux_phylogenetics.write_phytools_parity_summary_table
        is write_phytools_parity_summary_table
    )
    assert (
        bijux_phylogenetics.write_reference_parity_observation_table
        is write_reference_parity_observation_table
    )
    assert (
        bijux_phylogenetics.write_ape_parity_observation_table
        is write_ape_parity_observation_table
    )
    assert (
        bijux_phylogenetics.write_phytools_parity_observation_table
        is write_phytools_parity_observation_table
    )
    assert bijux_phylogenetics.render_tree_report is render_tree_report
    assert (
        bijux_phylogenetics.render_workflow_validation_report
        is render_workflow_validation_report
    )
    assert (
        bijux_phylogenetics.render_level_one_release_gate_report
        is render_level_one_release_gate_report
    )
    assert bijux_phylogenetics.simulate_birth_death_trees is simulate_birth_death_trees
    assert bijux_phylogenetics.simulate_random_tree is simulate_random_tree
    assert bijux_phylogenetics.simulate_coalescent_trees is simulate_coalescent_trees
    assert bijux_phylogenetics.simulate_coalescent_tree is simulate_coalescent_tree
    assert bijux_phylogenetics.simulate_random_trees is simulate_random_trees
    assert bijux_phylogenetics.simulate_brownian_traits is simulate_brownian_traits
    assert bijux_phylogenetics.simulate_early_burst_traits is simulate_early_burst_traits
    assert bijux_phylogenetics.simulate_ou_traits is simulate_ou_traits
    assert bijux_phylogenetics.simulate_discrete_traits is simulate_discrete_traits
    assert bijux_phylogenetics.simulate_dna_alignment is simulate_dna_alignment
    assert bijux_phylogenetics.simulate_protein_alignment is simulate_protein_alignment
    assert (
        bijux_phylogenetics.write_tree_simulation_record_table
        is write_tree_simulation_record_table
    )
    assert (
        bijux_phylogenetics.write_tree_simulation_envelope_table
        is write_tree_simulation_envelope_table
    )


def test_command_registry_exposes_discrete_evolution_surface() -> None:
    spec = get_command_spec("discrete-evolution")

    assert spec.domain == "discrete-state-evolution"
    assert spec.outputs == ("discrete-state-evolution-report",)


def test_command_registry_exposes_biogeography_surface() -> None:
    spec = get_command_spec("biogeography")

    assert spec.domain == "biogeography"
    assert spec.outputs == ("biogeography-report",)


def test_command_registry_exposes_host_association_surface() -> None:
    spec = get_command_spec("host-association")

    assert spec.domain == "host-association"
    assert spec.outputs == ("host-association-report",)


def test_command_registry_exposes_ecological_niche_surface() -> None:
    spec = get_command_spec("ecological-niche")

    assert spec.domain == "ecological-niche"
    assert spec.outputs == ("ecological-niche-report",)


def test_command_registry_exposes_phylogeography_surface() -> None:
    spec = get_command_spec("phylogeography")

    assert spec.domain == "geographic-reconstruction"
    assert spec.outputs == ("geographic-reconstruction-report",)


def test_command_registry_exposes_reference_parity_surface() -> None:
    spec = get_command_spec("parity")

    assert spec.domain == "reference-validation"
    assert spec.outputs == ("reference-parity-report",)


def test_command_registry_exposes_diversification_surface() -> None:
    spec = get_command_spec("diversification")

    assert spec.domain == "diversification-analysis"
    assert spec.outputs == ("diversification-report",)


def test_public_package_exports_comparative_and_bayesian_workflows() -> None:
    assert bijux_phylogenetics.BranchIdentityMetadata is BranchIdentityMetadata
    assert (
        bijux_phylogenetics.CorrelatedTraitComparisonRow is CorrelatedTraitComparisonRow
    )
    assert (
        bijux_phylogenetics.CorrelatedTraitEvolutionReport
        is CorrelatedTraitEvolutionReport
    )
    assert bijux_phylogenetics.CorrelatedTraitExclusion is CorrelatedTraitExclusion
    assert (
        bijux_phylogenetics.CorrelatedTraitObservationRow
        is CorrelatedTraitObservationRow
    )
    assert bijux_phylogenetics.TraitRegimeBranchRow is TraitRegimeBranchRow
    assert bijux_phylogenetics.TraitRegimeExclusion is TraitRegimeExclusion
    assert bijux_phylogenetics.TraitRegimeMappingReport is TraitRegimeMappingReport
    assert bijux_phylogenetics.TraitRegimeNodeRow is TraitRegimeNodeRow
    assert (
        bijux_phylogenetics.build_branch_identity_lookup is build_branch_identity_lookup
    )
    assert bijux_phylogenetics.benchmark_tree_validation is benchmark_tree_validation
    assert bijux_phylogenetics.benchmark_tree_comparison is benchmark_tree_comparison
    assert (
        bijux_phylogenetics.benchmark_alignment_diagnostics
        is benchmark_alignment_diagnostics
    )
    assert (
        bijux_phylogenetics.summarize_numeric_trait_readiness
        is summarize_numeric_trait_readiness
    )
    assert bijux_phylogenetics.summarize_numeric_trait is summarize_numeric_trait
    assert (
        bijux_phylogenetics.compute_phylogenetic_independent_contrasts
        is compute_phylogenetic_independent_contrasts
    )
    assert (
        bijux_phylogenetics.compute_phylogenetic_independent_contrasts_from_dataset
        is compute_phylogenetic_independent_contrasts_from_dataset
    )
    assert (
        bijux_phylogenetics.summarize_independent_contrast_regression
        is summarize_independent_contrast_regression
    )
    assert (
        bijux_phylogenetics.run_multivariate_comparative_regression
        is run_multivariate_comparative_regression
    )
    assert (
        bijux_phylogenetics.summarize_correlated_trait_evolution
        is summarize_correlated_trait_evolution
    )
    assert bijux_phylogenetics.run_posterior_tree_pgls is run_posterior_tree_pgls
    assert (
        bijux_phylogenetics.analyze_comparative_clade_stability
        is analyze_comparative_clade_stability
    )
    assert (
        bijux_phylogenetics.analyze_comparative_residual_clades
        is analyze_comparative_residual_clades
    )
    assert (
        bijux_phylogenetics.compare_comparative_regression_models
        is compare_comparative_regression_models
    )
    assert (
        bijux_phylogenetics.summarize_phylogenetic_logistic
        is summarize_phylogenetic_logistic
    )
    assert (
        bijux_phylogenetics.summarize_phylogenetic_signal
        is summarize_phylogenetic_signal
    )
    assert bijux_phylogenetics.compute_blombergs_k is compute_blombergs_k
    assert bijux_phylogenetics.fit_discrete_mk_model is fit_discrete_mk_model
    assert (
        bijux_phylogenetics.fit_discrete_mk_model_from_dataset
        is fit_discrete_mk_model_from_dataset
    )
    assert (
        bijux_phylogenetics.evaluate_pagels_lambda_likelihood
        is evaluate_pagels_lambda_likelihood
    )
    assert (
        bijux_phylogenetics.evaluate_pagels_lambda_likelihood_from_dataset
        is evaluate_pagels_lambda_likelihood_from_dataset
    )
    assert bijux_phylogenetics.estimate_pagels_lambda is estimate_pagels_lambda
    assert (
        bijux_phylogenetics.PagelLambdaLikelihoodReport
        is PagelLambdaLikelihoodReport
    )
    assert (
        bijux_phylogenetics.PagelLambdaOptimizerDiagnostics
        is PagelLambdaOptimizerDiagnostics
    )
    assert bijux_phylogenetics.PagelLambdaProfileRow is PagelLambdaProfileRow
    assert bijux_phylogenetics.DiscreteMkFitReport is DiscreteMkFitReport
    assert bijux_phylogenetics.DiscreteMkInputAudit is DiscreteMkInputAudit
    assert (
        bijux_phylogenetics.compute_phylogenetic_signal_test
        is compute_phylogenetic_signal_test
    )
    assert (
        bijux_phylogenetics.assess_comparative_method_maturity
        is assess_comparative_method_maturity
    )
    assert (
        bijux_phylogenetics.audit_comparative_parameter_uncertainty
        is audit_comparative_parameter_uncertainty
    )
    assert (
        bijux_phylogenetics.audit_ou_identifiability_reference_examples
        is audit_ou_identifiability_reference_examples
    )
    assert (
        bijux_phylogenetics.build_comparative_report_package
        is build_comparative_report_package
    )
    assert bijux_phylogenetics.build_pgls_model_matrix is build_pgls_model_matrix
    assert bijux_phylogenetics.inspect_pgls_inputs is inspect_pgls_inputs
    assert bijux_phylogenetics.run_pgls is run_pgls
    assert (
        bijux_phylogenetics.summarize_brownian_covariance
        is summarize_brownian_covariance
    )
    assert (
        bijux_phylogenetics.summarize_brownian_covariance_from_tree
        is summarize_brownian_covariance_from_tree
    )
    assert (
        bijux_phylogenetics.summarize_brownian_covariance_pgls
        is summarize_brownian_covariance_pgls
    )
    assert (
        bijux_phylogenetics.summarize_brownian_regime_rates
        is summarize_brownian_regime_rates
    )
    assert (
        bijux_phylogenetics.summarize_trait_regime_mapping
        is summarize_trait_regime_mapping
    )
    assert (
        bijux_phylogenetics.summarize_brownian_trait_evolution
        is summarize_brownian_trait_evolution
    )
    assert (
        bijux_phylogenetics.summarize_comparative_analysis
        is summarize_comparative_analysis
    )
    assert (
        bijux_phylogenetics.summarize_comparative_audit is summarize_comparative_audit
    )
    assert (
        bijux_phylogenetics.summarize_comparative_coefficients
        is summarize_comparative_coefficients
    )
    assert (
        bijux_phylogenetics.summarize_comparative_interpretation
        is summarize_comparative_interpretation
    )
    assert (
        bijux_phylogenetics.summarize_comparative_residuals
        is summarize_comparative_residuals
    )
    assert (
        bijux_phylogenetics.summarize_comparative_signal is summarize_comparative_signal
    )
    assert (
        bijux_phylogenetics.summarize_early_burst_trait_evolution
        is summarize_early_burst_trait_evolution
    )
    assert bijux_phylogenetics.summarize_clade_traits is summarize_clade_traits
    assert bijux_phylogenetics.summarize_trait_outliers is summarize_trait_outliers
    assert bijux_phylogenetics.summarize_trait_imputation is summarize_trait_imputation
    assert (
        bijux_phylogenetics.summarize_trait_rate_through_time
        is summarize_trait_rate_through_time
    )
    assert (
        bijux_phylogenetics.summarize_ou_covariance_pgls is summarize_ou_covariance_pgls
    )
    assert (
        bijux_phylogenetics.summarize_ou_trait_evolution is summarize_ou_trait_evolution
    )
    assert (
        bijux_phylogenetics.BrownianRegimeFitSummaryReport
        is BrownianRegimeFitSummaryReport
    )
    assert bijux_phylogenetics.BrownianRegimeBranchRow is BrownianRegimeBranchRow
    assert bijux_phylogenetics.BrownianRegimeRateRow is BrownianRegimeRateRow
    assert bijux_phylogenetics.BrownianRegimeProfileRow is BrownianRegimeProfileRow
    assert bijux_phylogenetics.BrownianRegimeExclusion is BrownianRegimeExclusion
    assert (
        bijux_phylogenetics.BrownianRegimeIdentifiabilityWarning
        is BrownianRegimeIdentifiabilityWarning
    )
    assert (
        bijux_phylogenetics.ComparativeAnalysisSummaryRow
        is ComparativeAnalysisSummaryRow
    )
    assert bijux_phylogenetics.ComparativeAuditTableRow is ComparativeAuditTableRow
    assert (
        bijux_phylogenetics.ComparativeCoefficientTableRow
        is ComparativeCoefficientTableRow
    )
    assert (
        bijux_phylogenetics.ComparativeInterpretationRow is ComparativeInterpretationRow
    )
    assert (
        bijux_phylogenetics.ComparativeReportPackageResult
        is ComparativeReportPackageResult
    )
    assert (
        bijux_phylogenetics.ComparativeResidualTableRow is ComparativeResidualTableRow
    )
    assert bijux_phylogenetics.ComparativeSignalTableRow is ComparativeSignalTableRow
    assert (
        bijux_phylogenetics.EarlyBurstTraitEvolutionSummaryReport
        is EarlyBurstTraitEvolutionSummaryReport
    )
    assert (
        bijux_phylogenetics.EarlyBurstTraitEvolutionExclusion
        is EarlyBurstTraitEvolutionExclusion
    )
    assert (
        bijux_phylogenetics.EarlyBurstRateChangeProfileRow
        is EarlyBurstRateChangeProfileRow
    )
    assert (
        bijux_phylogenetics.EarlyBurstIdentifiabilityWarning
        is EarlyBurstIdentifiabilityWarning
    )
    assert bijux_phylogenetics.CladeTraitSummaryReport is CladeTraitSummaryReport
    assert bijux_phylogenetics.CladeTraitRow is CladeTraitRow
    assert bijux_phylogenetics.CladeTraitStateCount is CladeTraitStateCount
    assert bijux_phylogenetics.CladeTraitExclusion is CladeTraitExclusion
    assert bijux_phylogenetics.TraitOutlierSummaryReport is TraitOutlierSummaryReport
    assert bijux_phylogenetics.TraitOutlierTaxonRow is TraitOutlierTaxonRow
    assert bijux_phylogenetics.TraitOutlierExclusion is TraitOutlierExclusion
    assert (
        bijux_phylogenetics.write_comparative_audit_table
        is write_comparative_audit_table
    )
    assert (
        bijux_phylogenetics.write_comparative_coefficient_table
        is write_comparative_coefficient_table
    )
    assert (
        bijux_phylogenetics.write_comparative_contrast_table
        is write_comparative_contrast_table
    )
    assert (
        bijux_phylogenetics.write_comparative_interpretation_table
        is write_comparative_interpretation_table
    )
    assert (
        bijux_phylogenetics.write_comparative_model_comparison_table
        is write_comparative_model_comparison_table
    )
    assert (
        bijux_phylogenetics.write_comparative_residual_table
        is write_comparative_residual_table
    )
    assert (
        bijux_phylogenetics.write_comparative_signal_table
        is write_comparative_signal_table
    )
    assert (
        bijux_phylogenetics.write_comparative_summary_table
        is write_comparative_summary_table
    )
    assert (
        bijux_phylogenetics.TraitImputationSummaryReport is TraitImputationSummaryReport
    )
    assert bijux_phylogenetics.TraitImputationRow is TraitImputationRow
    assert bijux_phylogenetics.TraitImputationHoldoutRow is TraitImputationHoldoutRow
    assert bijux_phylogenetics.TraitImputationExclusion is TraitImputationExclusion
    assert (
        bijux_phylogenetics.TraitRateThroughTimeSummaryReport
        is TraitRateThroughTimeSummaryReport
    )
    assert (
        bijux_phylogenetics.TraitRateThroughTimeIntervalRow
        is TraitRateThroughTimeIntervalRow
    )
    assert (
        bijux_phylogenetics.TraitRateThroughTimeExclusion
        is TraitRateThroughTimeExclusion
    )
    assert (
        bijux_phylogenetics.write_brownian_regime_summary_table
        is write_brownian_regime_summary_table
    )
    assert (
        bijux_phylogenetics.write_trait_regime_branch_table
        is write_trait_regime_branch_table
    )
    assert (
        bijux_phylogenetics.write_trait_regime_exclusion_table
        is write_trait_regime_exclusion_table
    )
    assert (
        bijux_phylogenetics.write_trait_regime_node_table
        is write_trait_regime_node_table
    )
    assert (
        bijux_phylogenetics.write_trait_regime_summary_table
        is write_trait_regime_summary_table
    )
    assert (
        bijux_phylogenetics.write_brownian_regime_rate_table
        is write_brownian_regime_rate_table
    )
    assert (
        bijux_phylogenetics.write_brownian_regime_profile_table
        is write_brownian_regime_profile_table
    )
    assert (
        bijux_phylogenetics.write_brownian_regime_comparison_table
        is write_brownian_regime_comparison_table
    )
    assert (
        bijux_phylogenetics.write_brownian_regime_branch_table
        is write_brownian_regime_branch_table
    )
    assert (
        bijux_phylogenetics.write_brownian_regime_exclusion_table
        is write_brownian_regime_exclusion_table
    )
    assert (
        bijux_phylogenetics.write_early_burst_trait_evolution_summary_table
        is write_early_burst_trait_evolution_summary_table
    )
    assert (
        bijux_phylogenetics.write_early_burst_trait_evolution_exclusion_table
        is write_early_burst_trait_evolution_exclusion_table
    )
    assert (
        bijux_phylogenetics.write_early_burst_trait_evolution_comparison_table
        is write_early_burst_trait_evolution_comparison_table
    )
    assert (
        bijux_phylogenetics.write_early_burst_rate_change_profile_table
        is write_early_burst_rate_change_profile_table
    )
    assert (
        bijux_phylogenetics.write_clade_trait_summary_table
        is write_clade_trait_summary_table
    )
    assert (
        bijux_phylogenetics.write_clade_trait_clade_table
        is write_clade_trait_clade_table
    )
    assert (
        bijux_phylogenetics.write_clade_trait_exclusion_table
        is write_clade_trait_exclusion_table
    )
    assert (
        bijux_phylogenetics.write_trait_outlier_summary_table
        is write_trait_outlier_summary_table
    )
    assert (
        bijux_phylogenetics.write_trait_outlier_taxon_table
        is write_trait_outlier_taxon_table
    )
    assert (
        bijux_phylogenetics.write_trait_outlier_exclusion_table
        is write_trait_outlier_exclusion_table
    )
    assert (
        bijux_phylogenetics.write_trait_imputation_summary_table
        is write_trait_imputation_summary_table
    )
    assert (
        bijux_phylogenetics.write_trait_imputation_table is write_trait_imputation_table
    )
    assert (
        bijux_phylogenetics.write_trait_imputation_holdout_table
        is write_trait_imputation_holdout_table
    )
    assert (
        bijux_phylogenetics.write_trait_imputation_exclusion_table
        is write_trait_imputation_exclusion_table
    )
    assert (
        bijux_phylogenetics.write_trait_rate_through_time_summary_table
        is write_trait_rate_through_time_summary_table
    )
    assert (
        bijux_phylogenetics.write_trait_rate_through_time_interval_table
        is write_trait_rate_through_time_interval_table
    )
    assert (
        bijux_phylogenetics.write_trait_rate_through_time_exclusion_table
        is write_trait_rate_through_time_exclusion_table
    )
    assert bijux_phylogenetics.summarize_pgls_lambda_fit is summarize_pgls_lambda_fit
    assert (
        bijux_phylogenetics.summarize_pgls_categorical_contrasts
        is summarize_pgls_categorical_contrasts
    )
    assert (
        bijux_phylogenetics.ComparativeCladeCoefficientChangeRow
        is ComparativeCladeCoefficientChangeRow
    )
    assert (
        bijux_phylogenetics.ComparativeCladeResidualReport
        is ComparativeCladeResidualReport
    )
    assert (
        bijux_phylogenetics.ComparativeCladeStabilityReport
        is ComparativeCladeStabilityReport
    )
    assert (
        bijux_phylogenetics.ComparativeCladeStabilityRow is ComparativeCladeStabilityRow
    )
    assert (
        bijux_phylogenetics.PosteriorTreePGLSTreeFitRow is PosteriorTreePGLSTreeFitRow
    )
    assert (
        bijux_phylogenetics.PosteriorTreePGLSCoefficientRow
        is PosteriorTreePGLSCoefficientRow
    )
    assert (
        bijux_phylogenetics.PosteriorTreePGLSCoefficientSummaryRow
        is PosteriorTreePGLSCoefficientSummaryRow
    )
    assert bijux_phylogenetics.PosteriorTreePGLSReport is PosteriorTreePGLSReport
    assert (
        bijux_phylogenetics.ComparativeResidualTaxonRow is ComparativeResidualTaxonRow
    )
    assert (
        bijux_phylogenetics.ComparativeResidualCladeRow is ComparativeResidualCladeRow
    )
    assert (
        bijux_phylogenetics.ComparativeRegressionModelSelectionReport
        is ComparativeRegressionModelSelectionReport
    )
    assert (
        bijux_phylogenetics.ComparativeRegressionModelRow
        is ComparativeRegressionModelRow
    )
    assert (
        bijux_phylogenetics.ComparativeRegressionPairwiseComparisonRow
        is ComparativeRegressionPairwiseComparisonRow
    )
    assert (
        bijux_phylogenetics.ComparativeRegressionModelExclusion
        is ComparativeRegressionModelExclusion
    )
    assert (
        bijux_phylogenetics.IndependentContrastRegressionReport
        is IndependentContrastRegressionReport
    )
    assert (
        bijux_phylogenetics.IndependentContrastRegressionRow
        is IndependentContrastRegressionRow
    )
    assert (
        bijux_phylogenetics.MultivariateComparativeRegressionReport
        is MultivariateComparativeRegressionReport
    )
    assert (
        bijux_phylogenetics.MultivariateResidualCorrelationRow
        is MultivariateResidualCorrelationRow
    )
    assert (
        bijux_phylogenetics.MultivariateResidualAssociationRow
        is MultivariateResidualAssociationRow
    )
    assert (
        bijux_phylogenetics.MultivariateResidualCovarianceRow
        is MultivariateResidualCovarianceRow
    )
    assert (
        bijux_phylogenetics.MultivariateResponseCoefficientRow
        is MultivariateResponseCoefficientRow
    )
    assert (
        bijux_phylogenetics.MultivariateResponseModelRow is MultivariateResponseModelRow
    )
    assert bijux_phylogenetics.MultivariateTaxonExclusion is MultivariateTaxonExclusion
    assert (
        bijux_phylogenetics.PhylogeneticLogisticCoefficient
        is PhylogeneticLogisticCoefficient
    )
    assert (
        bijux_phylogenetics.PhylogeneticLogisticFittedRow
        is PhylogeneticLogisticFittedRow
    )
    assert bijux_phylogenetics.PhylogeneticLogisticReport is PhylogeneticLogisticReport
    assert (
        bijux_phylogenetics.PhylogeneticLogisticWarning is PhylogeneticLogisticWarning
    )
    assert (
        bijux_phylogenetics.PhylogeneticSignalPermutation
        is PhylogeneticSignalPermutation
    )
    assert (
        bijux_phylogenetics.PhylogeneticSignalSummaryReport
        is PhylogeneticSignalSummaryReport
    )
    assert (
        bijux_phylogenetics.write_brownian_covariance_table
        is write_brownian_covariance_table
    )
    assert (
        bijux_phylogenetics.write_brownian_covariance_long_table
        is write_brownian_covariance_long_table
    )
    assert (
        bijux_phylogenetics.write_brownian_covariance_matrix_table
        is write_brownian_covariance_matrix_table
    )
    assert (
        bijux_phylogenetics.write_brownian_trait_evolution_summary_table
        is write_brownian_trait_evolution_summary_table
    )
    assert (
        bijux_phylogenetics.write_brownian_trait_evolution_exclusion_table
        is write_brownian_trait_evolution_exclusion_table
    )
    assert (
        bijux_phylogenetics.write_phylogenetic_logistic_coefficient_table
        is write_phylogenetic_logistic_coefficient_table
    )
    assert (
        bijux_phylogenetics.write_phylogenetic_logistic_excluded_taxa_table
        is write_phylogenetic_logistic_excluded_taxa_table
    )
    assert (
        bijux_phylogenetics.write_phylogenetic_logistic_fitted_table
        is write_phylogenetic_logistic_fitted_table
    )
    assert (
        bijux_phylogenetics.write_independent_contrast_table
        is write_independent_contrast_table
    )
    assert (
        bijux_phylogenetics.write_comparative_clade_stability_table
        is write_comparative_clade_stability_table
    )
    assert (
        bijux_phylogenetics.write_comparative_clade_coefficient_change_table
        is write_comparative_clade_coefficient_change_table
    )
    assert (
        bijux_phylogenetics.write_posterior_tree_pgls_tree_table
        is write_posterior_tree_pgls_tree_table
    )
    assert (
        bijux_phylogenetics.write_posterior_tree_pgls_coefficient_table
        is write_posterior_tree_pgls_coefficient_table
    )
    assert (
        bijux_phylogenetics.write_posterior_tree_pgls_summary_table
        is write_posterior_tree_pgls_summary_table
    )
    assert (
        bijux_phylogenetics.write_comparative_residual_taxon_table
        is write_comparative_residual_taxon_table
    )
    assert (
        bijux_phylogenetics.write_comparative_residual_clade_table
        is write_comparative_residual_clade_table
    )
    assert (
        bijux_phylogenetics.write_comparative_regression_model_ranking_table
        is write_comparative_regression_model_ranking_table
    )
    assert (
        bijux_phylogenetics.write_comparative_regression_pairwise_table
        is write_comparative_regression_pairwise_table
    )
    assert (
        bijux_phylogenetics.write_comparative_regression_excluded_taxa_table
        is write_comparative_regression_excluded_taxa_table
    )
    assert (
        bijux_phylogenetics.write_independent_contrast_regression_table
        is write_independent_contrast_regression_table
    )
    assert (
        bijux_phylogenetics.write_multivariate_residual_covariance_table
        is write_multivariate_residual_covariance_table
    )
    assert (
        bijux_phylogenetics.write_multivariate_residual_correlation_table
        is write_multivariate_residual_correlation_table
    )
    assert (
        bijux_phylogenetics.write_multivariate_residual_association_table
        is write_multivariate_residual_association_table
    )
    assert (
        bijux_phylogenetics.write_multivariate_response_model_table
        is write_multivariate_response_model_table
    )
    assert (
        bijux_phylogenetics.write_multivariate_response_coefficient_table
        is write_multivariate_response_coefficient_table
    )
    assert (
        bijux_phylogenetics.write_multivariate_excluded_taxa_table
        is write_multivariate_excluded_taxa_table
    )
    assert (
        bijux_phylogenetics.write_phylogenetic_signal_summary_table
        is write_phylogenetic_signal_summary_table
    )
    assert (
        bijux_phylogenetics.write_phylogenetic_signal_permutation_table
        is write_phylogenetic_signal_permutation_table
    )
    assert (
        bijux_phylogenetics.write_ou_alpha_profile_table is write_ou_alpha_profile_table
    )
    assert (
        bijux_phylogenetics.write_ou_trait_evolution_summary_table
        is write_ou_trait_evolution_summary_table
    )
    assert (
        bijux_phylogenetics.write_ou_trait_evolution_exclusion_table
        is write_ou_trait_evolution_exclusion_table
    )
    assert bijux_phylogenetics.write_ou_covariance_table is write_ou_covariance_table
    assert (
        bijux_phylogenetics.write_pgls_categorical_contrast_table
        is write_pgls_categorical_contrast_table
    )
    assert (
        bijux_phylogenetics.summarize_pgls_interaction_coefficients
        is summarize_pgls_interaction_coefficients
    )
    assert (
        bijux_phylogenetics.write_pgls_interaction_coefficient_table
        is write_pgls_interaction_coefficient_table
    )
    assert (
        bijux_phylogenetics.write_pgls_lambda_profile_table
        is write_pgls_lambda_profile_table
    )
    assert (
        bijux_phylogenetics.write_pgls_model_matrix_table
        is write_pgls_model_matrix_table
    )
    assert (
        bijux_phylogenetics.reconstruct_continuous_ancestral_states
        is reconstruct_continuous_ancestral_states
    )
    assert (
        bijux_phylogenetics.reconstruct_continuous_ancestral_states_from_dataset
        is reconstruct_continuous_ancestral_states_from_dataset
    )
    assert (
        bijux_phylogenetics.summarize_continuous_ancestral_report
        is summarize_continuous_ancestral_report
    )
    assert (
        bijux_phylogenetics.summarize_continuous_ancestral_confidence
        is summarize_continuous_ancestral_confidence
    )
    assert (
        bijux_phylogenetics.summarize_continuous_change_branches
        is summarize_continuous_change_branches
    )
    assert (
        bijux_phylogenetics.summarize_continuous_change_counts
        is summarize_continuous_change_counts
    )
    assert (
        bijux_phylogenetics.continuous_ancestral_exclusions
        is continuous_ancestral_exclusions
    )
    assert (
        bijux_phylogenetics.discrete_ancestral_exclusions
        is discrete_ancestral_exclusions
    )
    assert (
        bijux_phylogenetics.reconstruct_discrete_ancestral_states
        is reconstruct_discrete_ancestral_states
    )
    assert (
        bijux_phylogenetics.reconstruct_discrete_ancestral_states_from_dataset
        is reconstruct_discrete_ancestral_states_from_dataset
    )
    assert (
        bijux_phylogenetics.summarize_discrete_ancestral_report
        is summarize_discrete_ancestral_report
    )
    assert (
        bijux_phylogenetics.summarize_discrete_ancestral_confidence
        is summarize_discrete_ancestral_confidence
    )
    assert (
        bijux_phylogenetics.build_ancestral_figure_package
        is build_ancestral_figure_package
    )
    assert (
        bijux_phylogenetics.build_ancestral_report_package
        is build_ancestral_report_package
    )
    assert (
        bijux_phylogenetics.build_ancestral_sensitivity_report
        is build_ancestral_sensitivity_report
    )
    assert (
        bijux_phylogenetics.compare_continuous_ancestral_models
        is compare_continuous_ancestral_models
    )
    assert (
        bijux_phylogenetics.compare_discrete_ancestral_reconstructions
        is compare_discrete_ancestral_reconstructions
    )
    assert (
        bijux_phylogenetics.render_ancestral_state_tree is render_ancestral_state_tree
    )
    assert (
        bijux_phylogenetics.build_continuous_ancestral_confidence_rows
        is build_continuous_ancestral_confidence_rows
    )
    assert (
        bijux_phylogenetics.write_continuous_ancestral_summary_table
        is write_continuous_ancestral_summary_table
    )
    assert (
        bijux_phylogenetics.write_continuous_change_branch_table
        is write_continuous_change_branch_table
    )
    assert (
        bijux_phylogenetics.write_continuous_change_count_table
        is write_continuous_change_count_table
    )
    assert (
        bijux_phylogenetics.write_continuous_ancestral_confidence_table
        is write_continuous_ancestral_confidence_table
    )
    assert (
        bijux_phylogenetics.write_continuous_ancestral_uncertainty_table
        is write_continuous_ancestral_uncertainty_table
    )
    assert (
        bijux_phylogenetics.write_continuous_ancestral_exclusion_table
        is write_continuous_ancestral_exclusion_table
    )
    assert (
        bijux_phylogenetics.build_discrete_ancestral_confidence_rows
        is build_discrete_ancestral_confidence_rows
    )
    assert (
        bijux_phylogenetics.write_discrete_ancestral_summary_table
        is write_discrete_ancestral_summary_table
    )
    assert (
        bijux_phylogenetics.write_discrete_ancestral_confidence_table
        is write_discrete_ancestral_confidence_table
    )
    assert (
        bijux_phylogenetics.write_discrete_ancestral_probability_table
        is write_discrete_ancestral_probability_table
    )
    assert (
        bijux_phylogenetics.write_discrete_ancestral_transition_table
        is write_discrete_ancestral_transition_table
    )
    assert (
        bijux_phylogenetics.write_discrete_ancestral_comparison_table
        is write_discrete_ancestral_comparison_table
    )
    assert (
        bijux_phylogenetics.write_discrete_ancestral_exclusion_table
        is write_discrete_ancestral_exclusion_table
    )
    assert (
        bijux_phylogenetics.write_ancestral_confidence_summary_table
        is write_ancestral_confidence_summary_table
    )
    assert (
        bijux_phylogenetics.render_ancestral_state_report
        is render_ancestral_state_report
    )
    assert (
        bijux_phylogenetics.write_ancestral_state_table is write_ancestral_state_table
    )
    assert (
        bijux_phylogenetics.render_ancestral_state_visualization
        is render_ancestral_state_visualization
    )
    assert (
        bijux_phylogenetics.validate_discrete_ancestral_reference_examples
        is validate_discrete_ancestral_reference_examples
    )
    assert (
        bijux_phylogenetics.summarize_continuous_ancestral_tree_set
        is summarize_continuous_ancestral_tree_set
    )
    assert (
        bijux_phylogenetics.summarize_continuous_ancestral_tree_set_confidence
        is summarize_continuous_ancestral_tree_set_confidence
    )
    assert (
        bijux_phylogenetics.summarize_continuous_ancestral_tree_set_report
        is summarize_continuous_ancestral_tree_set_report
    )
    assert (
        bijux_phylogenetics.build_continuous_ancestral_tree_set_confidence_rows
        is build_continuous_ancestral_tree_set_confidence_rows
    )
    assert (
        bijux_phylogenetics.summarize_discrete_ancestral_tree_set
        is summarize_discrete_ancestral_tree_set
    )
    assert (
        bijux_phylogenetics.summarize_discrete_ancestral_tree_set_confidence
        is summarize_discrete_ancestral_tree_set_confidence
    )
    assert (
        bijux_phylogenetics.summarize_discrete_ancestral_tree_set_report
        is summarize_discrete_ancestral_tree_set_report
    )
    assert (
        bijux_phylogenetics.build_discrete_ancestral_tree_set_confidence_rows
        is build_discrete_ancestral_tree_set_confidence_rows
    )
    assert (
        bijux_phylogenetics.summarize_irreversible_discrete_reconstruction
        is summarize_irreversible_discrete_reconstruction
    )
    assert (
        bijux_phylogenetics.summarize_irreversible_discrete_report
        is summarize_irreversible_discrete_report
    )
    assert (
        bijux_phylogenetics.summarize_ordered_discrete_reconstruction
        is summarize_ordered_discrete_reconstruction
    )
    assert (
        bijux_phylogenetics.summarize_ordered_discrete_report
        is summarize_ordered_discrete_report
    )
    assert (
        bijux_phylogenetics.summarize_ancestral_root_sensitivity
        is summarize_ancestral_root_sensitivity
    )
    assert (
        bijux_phylogenetics.summarize_ancestral_root_sensitivity_report
        is summarize_ancestral_root_sensitivity_report
    )
    assert (
        bijux_phylogenetics.summarize_ancestral_transitions
        is summarize_ancestral_transitions
    )
    assert (
        bijux_phylogenetics.summarize_ancestral_transition_report
        is summarize_ancestral_transition_report
    )
    assert (
        bijux_phylogenetics.summarize_ancestral_transition_tree_set
        is summarize_ancestral_transition_tree_set
    )
    assert (
        bijux_phylogenetics.summarize_ancestral_transition_tree_set_report
        is summarize_ancestral_transition_tree_set_report
    )
    assert (
        bijux_phylogenetics.write_irreversible_discrete_summary_table
        is write_irreversible_discrete_summary_table
    )
    assert (
        bijux_phylogenetics.write_irreversible_discrete_fit_table
        is write_irreversible_discrete_fit_table
    )
    assert (
        bijux_phylogenetics.write_irreversible_discrete_node_table
        is write_irreversible_discrete_node_table
    )
    assert (
        bijux_phylogenetics.write_irreversible_discrete_transition_table
        is write_irreversible_discrete_transition_table
    )
    assert (
        bijux_phylogenetics.write_ordered_discrete_summary_table
        is write_ordered_discrete_summary_table
    )
    assert (
        bijux_phylogenetics.write_ordered_discrete_fit_table
        is write_ordered_discrete_fit_table
    )
    assert (
        bijux_phylogenetics.write_ordered_discrete_node_table
        is write_ordered_discrete_node_table
    )
    assert (
        bijux_phylogenetics.write_ordered_discrete_transition_table
        is write_ordered_discrete_transition_table
    )
    assert (
        bijux_phylogenetics.write_ancestral_root_sensitivity_summary_table
        is write_ancestral_root_sensitivity_summary_table
    )
    assert (
        bijux_phylogenetics.write_ancestral_root_assumption_table
        is write_ancestral_root_assumption_table
    )
    assert (
        bijux_phylogenetics.write_ancestral_root_sensitivity_node_table
        is write_ancestral_root_sensitivity_node_table
    )
    assert (
        bijux_phylogenetics.write_ancestral_tree_set_tree_table
        is write_ancestral_tree_set_tree_table
    )
    assert (
        bijux_phylogenetics.write_ancestral_transition_summary_table
        is write_ancestral_transition_summary_table
    )
    assert (
        bijux_phylogenetics.write_ancestral_transition_branch_table
        is write_ancestral_transition_branch_table
    )
    assert (
        bijux_phylogenetics.write_ancestral_transition_count_table
        is write_ancestral_transition_count_table
    )
    assert (
        bijux_phylogenetics.write_ancestral_transition_exclusion_table
        is write_ancestral_transition_exclusion_table
    )
    assert (
        bijux_phylogenetics.write_ancestral_transition_tree_set_summary_table
        is write_ancestral_transition_tree_set_summary_table
    )
    assert (
        bijux_phylogenetics.write_ancestral_transition_tree_set_tree_table
        is write_ancestral_transition_tree_set_tree_table
    )
    assert (
        bijux_phylogenetics.write_ancestral_transition_tree_set_branch_table
        is write_ancestral_transition_tree_set_branch_table
    )
    assert (
        bijux_phylogenetics.write_ancestral_transition_tree_set_count_table
        is write_ancestral_transition_tree_set_count_table
    )
    assert (
        bijux_phylogenetics.write_ancestral_tree_set_exclusion_table
        is write_ancestral_tree_set_exclusion_table
    )
    assert (
        bijux_phylogenetics.write_continuous_ancestral_tree_set_summary_table
        is write_continuous_ancestral_tree_set_summary_table
    )
    assert (
        bijux_phylogenetics.write_continuous_ancestral_tree_set_confidence_table
        is write_continuous_ancestral_tree_set_confidence_table
    )
    assert (
        bijux_phylogenetics.write_continuous_ancestral_tree_set_node_table
        is write_continuous_ancestral_tree_set_node_table
    )
    assert (
        bijux_phylogenetics.write_continuous_ancestral_tree_set_clade_table
        is write_continuous_ancestral_tree_set_clade_table
    )
    assert (
        bijux_phylogenetics.write_discrete_ancestral_tree_set_summary_table
        is write_discrete_ancestral_tree_set_summary_table
    )
    assert (
        bijux_phylogenetics.write_discrete_ancestral_tree_set_confidence_table
        is write_discrete_ancestral_tree_set_confidence_table
    )
    assert (
        bijux_phylogenetics.write_discrete_ancestral_tree_set_node_table
        is write_discrete_ancestral_tree_set_node_table
    )
    assert (
        bijux_phylogenetics.write_discrete_ancestral_tree_set_clade_table
        is write_discrete_ancestral_tree_set_clade_table
    )
    assert (
        bijux_phylogenetics.run_multiple_sequence_alignment
        is run_multiple_sequence_alignment
    )
    assert bijux_phylogenetics.run_alignment_trimming is run_alignment_trimming
    assert (
        bijux_phylogenetics.audit_alignment_inference_readiness
        is audit_alignment_inference_readiness
    )
    assert bijux_phylogenetics.run_model_selection is run_model_selection
    assert (
        bijux_phylogenetics.validate_model_selection_against_engine_outputs
        is validate_model_selection_against_engine_outputs
    )
    assert (
        bijux_phylogenetics.run_maximum_likelihood_tree_inference
        is run_maximum_likelihood_tree_inference
    )
    assert (
        bijux_phylogenetics.validate_ml_tree_contains_expected_taxa
        is validate_ml_tree_contains_expected_taxa
    )
    assert (
        bijux_phylogenetics.run_bootstrap_support_estimation
        is run_bootstrap_support_estimation
    )
    assert (
        bijux_phylogenetics.run_sh_alrt_support_estimation
        is run_sh_alrt_support_estimation
    )
    assert bijux_phylogenetics.run_fasta_to_tree_workflow is run_fasta_to_tree_workflow
    assert (
        bijux_phylogenetics.infer_unaligned_sequence_type
        is infer_unaligned_sequence_type
    )
    assert (
        bijux_phylogenetics.validate_bootstrap_tree_set is validate_bootstrap_tree_set
    )
    assert (
        bijux_phylogenetics.run_bootstrap_consensus_tree is run_bootstrap_consensus_tree
    )
    assert bijux_phylogenetics.run_fast_tree_inference is run_fast_tree_inference
    assert (
        bijux_phylogenetics.run_large_alignment_inference
        is run_large_alignment_inference
    )
    assert (
        bijux_phylogenetics.run_inference_reproducibility_check
        is run_inference_reproducibility_check
    )
    assert (
        bijux_phylogenetics.run_tree_inference_comparison
        is run_tree_inference_comparison
    )
    assert bijux_phylogenetics.compare_fast_and_ml_trees is compare_fast_and_ml_trees
    assert (
        bijux_phylogenetics.build_inference_comparison_shared_clade_rows
        is build_inference_comparison_shared_clade_rows
    )
    assert (
        bijux_phylogenetics.build_inference_comparison_weighted_conflict_rows
        is build_inference_comparison_weighted_conflict_rows
    )
    assert (
        bijux_phylogenetics.build_inference_comparison_conclusion_rows
        is build_inference_comparison_conclusion_rows
    )
    assert (
        bijux_phylogenetics.write_inference_comparison_clade_table
        is write_inference_comparison_clade_table
    )
    assert (
        bijux_phylogenetics.write_inference_comparison_weighted_conflict_table
        is write_inference_comparison_weighted_conflict_table
    )
    assert (
        bijux_phylogenetics.write_inference_comparison_conclusion_table
        is write_inference_comparison_conclusion_table
    )
    assert (
        bijux_phylogenetics.write_inference_comparison_summary_table
        is write_inference_comparison_summary_table
    )
    assert (
        bijux_phylogenetics.rewrite_inference_comparison_report_html
        is rewrite_inference_comparison_report_html
    )
    assert (
        bijux_phylogenetics.compare_inferred_tree_to_taxon_metadata
        is compare_inferred_tree_to_taxon_metadata
    )
    assert (
        bijux_phylogenetics.classify_inference_workflow_failure
        is classify_inference_workflow_failure
    )
    assert (
        bijux_phylogenetics.validate_inference_engine_outputs
        is validate_inference_engine_outputs
    )
    assert (
        bijux_phylogenetics.render_inference_workflow_report
        is render_inference_workflow_report
    )
    assert (
        bijux_phylogenetics.export_workflow_result_bundle
        is export_workflow_result_bundle
    )
    assert (
        bijux_phylogenetics.validate_workflow_result_bundle
        is validate_workflow_result_bundle
    )
    assert bijux_phylogenetics.prepare_mrbayes_analysis is prepare_mrbayes_analysis
    assert (
        bijux_phylogenetics.run_mrbayes_posterior_inference
        is run_mrbayes_posterior_inference
    )
    assert (
        bijux_phylogenetics.summarize_mrbayes_posterior_trees
        is summarize_mrbayes_posterior_trees
    )
    assert (
        bijux_phylogenetics.parse_mrbayes_parameter_traces
        is parse_mrbayes_parameter_traces
    )
    assert (
        bijux_phylogenetics.parse_mrbayes_posterior_tree_samples
        is parse_mrbayes_posterior_tree_samples
    )
    assert (
        bijux_phylogenetics.parse_mrbayes_mcmc_diagnostics
        is parse_mrbayes_mcmc_diagnostics
    )
    assert (
        bijux_phylogenetics.parse_mrbayes_consensus_tree is parse_mrbayes_consensus_tree
    )
    assert (
        bijux_phylogenetics.compute_mrbayes_effective_sample_sizes
        is compute_mrbayes_effective_sample_sizes
    )
    assert (
        bijux_phylogenetics.summarize_mrbayes_parameter_diagnostics
        is summarize_mrbayes_parameter_diagnostics
    )
    assert (
        bijux_phylogenetics.write_mrbayes_parameter_summary_table
        is write_mrbayes_parameter_summary_table
    )
    assert (
        bijux_phylogenetics.assess_mrbayes_burnin_sensitivity
        is assess_mrbayes_burnin_sensitivity
    )
    assert (
        bijux_phylogenetics.write_mrbayes_burnin_sensitivity_slice_table
        is write_mrbayes_burnin_sensitivity_slice_table
    )
    assert bijux_phylogenetics.assess_mrbayes_convergence is assess_mrbayes_convergence
    assert bijux_phylogenetics.MrBayesParameterSummary is MrBayesParameterSummary
    assert (
        bijux_phylogenetics.MrBayesParameterDiagnosticsReport
        is MrBayesParameterDiagnosticsReport
    )
    assert (
        bijux_phylogenetics.MrBayesBurninSensitivityReport
        is MrBayesBurninSensitivityReport
    )
    assert (
        bijux_phylogenetics.MrBayesBurninSensitivitySlice
        is MrBayesBurninSensitivitySlice
    )
    assert (
        bijux_phylogenetics.render_bayesian_posterior_report
        is render_bayesian_posterior_report
    )
    assert (
        bijux_phylogenetics.validate_fossil_calibration_table
        is validate_fossil_calibration_table
    )
    assert (
        bijux_phylogenetics.detect_impossible_calibration_constraints
        is detect_impossible_calibration_constraints
    )
    assert (
        bijux_phylogenetics.validate_tip_dating_metadata is validate_tip_dating_metadata
    )
    assert (
        bijux_phylogenetics.prepare_beast_time_tree_analysis
        is prepare_beast_time_tree_analysis
    )
    assert (
        bijux_phylogenetics.run_beast_posterior_inference
        is run_beast_posterior_inference
    )
    assert bijux_phylogenetics.BeastAnalysisXmlReport is BeastAnalysisXmlReport
    assert bijux_phylogenetics.BeastCalibration is BeastCalibration
    assert (
        bijux_phylogenetics.BeastPosteriorConsensusReport
        is BeastPosteriorConsensusReport
    )
    assert (
        bijux_phylogenetics.BeastPosteriorTopologyDiversityReport
        is BeastPosteriorTopologyDiversityReport
    )
    assert bijux_phylogenetics.parse_beast_log is parse_beast_log
    assert (
        bijux_phylogenetics.parse_beast_posterior_tree_samples
        is parse_beast_posterior_tree_samples
    )
    assert (
        bijux_phylogenetics.summarize_beast_analysis_xml is summarize_beast_analysis_xml
    )
    assert bijux_phylogenetics.summarize_beast_log is summarize_beast_log
    assert (
        bijux_phylogenetics.summarize_beast_posterior_trees
        is summarize_beast_posterior_trees
    )
    assert (
        bijux_phylogenetics.summarize_beast_posterior_topology_diversity
        is summarize_beast_posterior_topology_diversity
    )
    assert bijux_phylogenetics.assess_beast_convergence is assess_beast_convergence
    assert (
        bijux_phylogenetics.validate_beast_analysis_xml is validate_beast_analysis_xml
    )
    assert (
        bijux_phylogenetics.validate_beast_posterior_log is validate_beast_posterior_log
    )
    assert (
        bijux_phylogenetics.assess_beast_burnin_sensitivity
        is assess_beast_burnin_sensitivity
    )
    assert bijux_phylogenetics.assess_beast_chain_mixing is assess_beast_chain_mixing
    assert (
        bijux_phylogenetics.write_beast_burnin_sensitivity_slice_table
        is write_beast_burnin_sensitivity_slice_table
    )
    assert (
        bijux_phylogenetics.write_beast_log_summary_table
        is write_beast_log_summary_table
    )
    assert (
        bijux_phylogenetics.write_beast_posterior_tree_set
        is write_beast_posterior_tree_set
    )
    assert (
        bijux_phylogenetics.subsample_beast_posterior_tree_set
        is subsample_beast_posterior_tree_set
    )
    assert (
        bijux_phylogenetics.subsample_mrbayes_posterior_tree_set
        is subsample_mrbayes_posterior_tree_set
    )
    assert (
        bijux_phylogenetics.subsample_posterior_tree_set is subsample_posterior_tree_set
    )
    assert (
        bijux_phylogenetics.summarize_maximum_clade_credibility_tree
        is summarize_maximum_clade_credibility_tree
    )
    assert bijux_phylogenetics.thin_posterior_tree_set is thin_posterior_tree_set
    assert (
        bijux_phylogenetics.summarize_posterior_node_ages
        is summarize_posterior_node_ages
    )
    assert bijux_phylogenetics.compare_bayesian_tree_sets is compare_bayesian_tree_sets
    assert (
        bijux_phylogenetics.compare_independent_bayesian_runs
        is compare_independent_bayesian_runs
    )
    assert (
        bijux_phylogenetics.compare_posterior_tree_sets_by_prior
        is compare_posterior_tree_sets_by_prior
    )
    assert (
        bijux_phylogenetics.compare_posterior_tree_sets_by_clock
        is compare_posterior_tree_sets_by_clock
    )
    assert (
        bijux_phylogenetics.render_bayesian_run_comparison_report
        is render_bayesian_run_comparison_report
    )
    assert (
        bijux_phylogenetics.render_bayesian_diagnostics_report
        is render_bayesian_diagnostics_report
    )
    assert (
        bijux_phylogenetics.build_posterior_uncertainty_figure_package
        is build_posterior_uncertainty_figure_package
    )
    assert (
        bijux_phylogenetics.write_supplementary_bayesian_diagnostics_table
        is write_supplementary_bayesian_diagnostics_table
    )
    assert (
        bijux_phylogenetics.write_posterior_tree_subsample
        is write_posterior_tree_subsample
    )
    assert (
        bijux_phylogenetics.write_posterior_tree_subsample_table
        is write_posterior_tree_subsample_table
    )
    assert (
        bijux_phylogenetics.write_bayesian_methods_summary_text
        is write_bayesian_methods_summary_text
    )
    assert (
        bijux_phylogenetics.render_calibration_audit_report
        is render_calibration_audit_report
    )
    assert (
        bijux_phylogenetics.build_bayesian_evidence_package
        is build_bayesian_evidence_package
    )


def test_simulate_birth_death_trees_returns_requested_tree_and_tip_counts(
    tmp_path: Path,
) -> None:
    trees, report = simulate_birth_death_trees(tree_count=2, tip_count=4, seed=7)
    assert report.model == "birth-death"
    assert report.tree_count == 2
    assert report.tip_count == 4
    assert [tree.tip_count for tree in trees] == [4, 4]
    assert [row.newick for row in report.records] == [
        "(((Taxon1:0,Taxon2:0):0.204697722139132,Taxon3:0.204697722139132):0.200894048820103,Taxon4:0.405591770959235);",
        "((Taxon1:0.075806896024214,(Taxon2:0,Taxon3:0):0.075806896024214):0.054021431036002,Taxon4:0.129828327060217);",
    ]
    output_path = tmp_path / "birth-death.trees"
    write_tree_set(output_path, trees)
    assert output_path.read_text(encoding="utf-8") == (
        "(((Taxon1:0,Taxon2:0):0.204697722139132,Taxon3:0.204697722139132):0.200894048820103,Taxon4:0.405591770959235);\n"
        "((Taxon1:0.075806896024214,(Taxon2:0,Taxon3:0):0.075806896024214):0.054021431036002,Taxon4:0.129828327060217);\n"
    )


def test_simulate_coalescent_trees_returns_requested_sample_size() -> None:
    trees, report = simulate_coalescent_trees(tree_count=1, tip_count=4, seed=7)
    assert report.model == "coalescent"
    assert report.tree_count == 1
    assert report.tip_count == 4
    assert sorted(trees[0].tip_names) == ["Taxon1", "Taxon2", "Taxon3", "Taxon4"]
    assert report.records[0].newick == (
        "(((Taxon1:0.065219140705801,Taxon4:0.065219140705801):0.054506152335159,"
        "Taxon2:0.11972529304096):1.05249561807093,Taxon3:1.17222091111189);"
    )


def test_simulate_brownian_traits_generates_one_value_per_tip(tmp_path: Path) -> None:
    report = simulate_brownian_traits(
        fixture("example_tree.nwk"), seed=7, root_state=1.0, sigma=0.5
    )
    assert report.model == "brownian-motion"
    assert report.tip_count == 4
    assert [(row.taxon, row.value) for row in report.traits] == [
        ("A", 1.023647850429746),
        ("B", 0.907034485545723),
        ("C", 0.742224918944575),
        ("D", 0.90248754291519),
    ]
    output_path = tmp_path / "brownian.tsv"
    write_continuous_trait_table(output_path, report)
    assert output_path.read_text(encoding="utf-8") == (
        "taxon\tvalue\n"
        "A\t1.02364785042975\n"
        "B\t0.907034485545723\n"
        "C\t0.742224918944575\n"
        "D\t0.90248754291519\n"
    )
    assert [row.node for row in report.node_values if not row.is_tip] == [
        "A|B|C|D",
        "A|B",
        "C|D",
    ]


def test_simulate_ou_traits_uses_declared_parameters() -> None:
    report = simulate_ou_traits(
        fixture("example_tree.nwk"),
        seed=7,
        root_state=1.0,
        sigma=0.5,
        alpha=1.25,
        theta=0.25,
    )
    assert report.model == "ornstein-uhlenbeck"
    assert report.alpha == 1.25
    assert report.theta == 0.25
    assert [(row.taxon, row.value) for row in report.traits] == [
        ("A", 0.796738462243473),
        ("B", 0.687047684603513),
        ("C", 0.544493844861481),
        ("D", 0.68666212038887),
    ]


def test_simulate_early_burst_traits_uses_declared_rate_change() -> None:
    report = simulate_early_burst_traits(
        fixture("example_tree.nwk"),
        seed=7,
        root_state=1.0,
        sigma=0.5,
        rate_change=4.0,
    )
    assert report.model == "early-burst"
    assert report.rate_change == 4.0
    assert [(row.taxon, row.value) for row in report.traits] == [
        ("A", 1.062949925373677),
        ("B", 0.870045899658049),
        ("C", 0.630380931701498),
        ("D", 0.872656745957104),
    ]


def test_cli_simulate_early_burst_traits_reports_rate_change(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "early-burst.tsv"
    exit_code = main(
        [
            "simulate",
            "traits-early-burst",
            str(fixture("example_tree.nwk")),
            "--root-state",
            "1.0",
            "--sigma",
            "0.5",
            "--rate-change",
            "4.0",
            "--seed",
            "7",
            "--out",
            str(output),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["command"] == "simulate"
    assert payload["metrics"]["rate_change"] == 4.0
    assert payload["metrics"]["trait_count"] == 4
    assert payload["data"]["model"] == "early-burst"
    assert payload["data"]["rate_change"] == 4.0
    assert output.read_text(encoding="utf-8").splitlines()[0] == "taxon\tvalue"


def test_simulate_discrete_traits_assigns_a_state_to_every_tip(tmp_path: Path) -> None:
    report = simulate_discrete_traits(
        fixture("example_tree.nwk"),
        states=["wet", "dry", "mixed"],
        transition_rate=8.0,
        root_state="wet",
        seed=3,
    )
    assert report.model == "symmetric-discrete"
    assert report.tip_count == 4
    assert [(row.taxon, row.state) for row in report.traits] == [
        ("A", "wet"),
        ("B", "dry"),
        ("C", "wet"),
        ("D", "mixed"),
    ]
    output_path = tmp_path / "discrete.tsv"
    write_discrete_trait_table(output_path, report)
    assert output_path.read_text(encoding="utf-8") == (
        "taxon\tstate\nA\twet\nB\tdry\nC\twet\nD\tmixed\n"
    )
    assert [row.node for row in report.node_states if not row.is_tip] == [
        "A|B|C|D",
        "A|B",
        "C|D",
    ]
    assert [
        (
            row.parent_node,
            row.child_node,
            row.start_state,
            row.end_state,
            row.event_count,
        )
        for row in report.branch_histories
    ] == [
        ("A|B|C|D", "A|B", "wet", "mixed", 1),
        ("A|B", "A", "mixed", "wet", 1),
        ("A|B", "B", "mixed", "dry", 2),
        ("A|B|C|D", "C|D", "wet", "dry", 2),
        ("C|D", "C", "dry", "wet", 1),
        ("C|D", "D", "dry", "mixed", 4),
    ]
    assert [
        (event.source_state, event.target_state)
        for row in report.branch_histories
        for event in row.events
    ] == [
        ("wet", "mixed"),
        ("mixed", "wet"),
        ("mixed", "wet"),
        ("wet", "dry"),
        ("wet", "mixed"),
        ("mixed", "dry"),
        ("dry", "wet"),
        ("dry", "wet"),
        ("wet", "mixed"),
        ("mixed", "wet"),
        ("wet", "mixed"),
    ]


def test_simulate_dna_alignment_returns_requested_taxa_and_length(
    tmp_path: Path,
) -> None:
    report = simulate_dna_alignment(
        fixture("example_tree.nwk"), sequence_length=8, substitution_rate=1.2, seed=7
    )
    assert report.model == "jukes-cantor-like"
    assert report.tip_count == 4
    assert report.sequence_length == 8
    assert [(row.identifier, row.sequence) for row in report.records] == [
        ("A", "ACTAACGA"),
        ("B", "ACTAACGA"),
        ("C", "GCTAGAAA"),
        ("D", "GCGAAGAA"),
    ]
    output_path = tmp_path / "simulated-dna.fasta"
    write_simulated_alignment(output_path, report)
    assert output_path.read_text(encoding="utf-8") == (
        ">A\nACTAACGA\n>B\nACTAACGA\n>C\nGCTAGAAA\n>D\nGCGAAGAA\n"
    )


def test_simulate_protein_alignment_returns_requested_length_and_alphabet() -> None:
    report = simulate_protein_alignment(
        fixture("example_tree.nwk"), sequence_length=6, substitution_rate=0.8, seed=7
    )
    assert report.model == "symmetric-protein"
    assert report.inferred_alphabet == "protein"
    assert report.sequence_length == 6
    assert [(row.identifier, row.sequence) for row in report.records] == [
        ("A", "MFDCDP"),
        ("B", "MFDCDV"),
        ("C", "MFPCDV"),
        ("D", "MFICDV"),
    ]


def test_benchmark_tree_validation_reports_runtime_and_memory_by_size_class() -> None:
    report = benchmark_tree_validation(
        replicates=1,
        size_classes=[("small", 4), ("medium", 8), ("large", 16)],
    )
    assert report.replicates == 1
    assert [(row.label, row.item_count) for row in report.observations] == [
        ("small", 4),
        ("medium", 8),
        ("large", 16),
    ]
    assert all(row.runtime_seconds >= 0.0 for row in report.observations)
    assert all(row.peak_memory_bytes >= 0 for row in report.observations)


def test_benchmark_tree_comparison_reports_scaling_curve() -> None:
    report = benchmark_tree_comparison(replicates=1, taxon_counts=[4, 8, 16])
    assert report.replicates == 1
    assert [(row.label, row.item_count) for row in report.observations] == [
        ("taxa-4", 4),
        ("taxa-8", 8),
        ("taxa-16", 16),
    ]
    assert all(row.runtime_seconds >= 0.0 for row in report.observations)
    assert all(row.peak_memory_bytes >= 0 for row in report.observations)


def test_benchmark_alignment_diagnostics_reports_runtime_and_memory() -> None:
    report = benchmark_alignment_diagnostics(
        replicates=1, sequence_counts=[4, 8, 16], sequence_length=24
    )
    assert report.replicates == 1
    assert [(row.label, row.item_count) for row in report.observations] == [
        ("sequences-4", 4),
        ("sequences-8", 8),
        ("sequences-16", 16),
    ]
    assert all(row.runtime_seconds >= 0.0 for row in report.observations)
    assert all(row.peak_memory_bytes >= 0 for row in report.observations)


def test_load_tree_set_reports_tree_count_and_topology_diversity() -> None:
    report = load_tree_set(fixture("example_tree_set_left.nwk"))
    assert report.tree_count == 3
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.rooted_topology_count == 2
    assert report.unrooted_topology_count == 2
    assert [(row.index, row.tip_count) for row in report.records] == [
        (1, 4),
        (2, 4),
        (3, 4),
    ]


def test_compute_clade_frequency_table_counts_informative_clades() -> None:
    report = compute_clade_frequency_table(fixture("example_tree_set_left.nwk"))
    assert [
        (row.clade, row.tree_count, row.frequency) for row in report.clade_frequencies
    ] == [
        ("A|B", 2, 0.666666666666667),
        ("A|C", 1, 0.333333333333333),
        ("B|D", 1, 0.333333333333333),
        ("C|D", 2, 0.666666666666667),
    ]


def test_compute_consensus_tree_returns_majority_rule_consensus() -> None:
    tree, report = compute_consensus_tree(fixture("example_tree_set_left.nwk"))
    assert dumps_newick(tree) == (
        "((A:0.1,B:0.1)66.6666666666667:0.2,(C:0.1,D:0.1)66.6666666666667:0.2);"
    )
    assert report.tree_count == 3
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.consensus_method == "majority-rule"
    assert report.consensus_threshold == 0.5
    assert report.included_clade_count == 2


def test_compute_consensus_tree_stays_native_without_biopython_bridge() -> None:
    import bijux_phylogenetics.tree_set as tree_set_module

    assert not hasattr(tree_set_module, "Phylo")
    assert not hasattr(tree_set_module, "tree_from_biophylo")

    tree, report = tree_set_module.compute_consensus_tree(
        fixture("example_tree_set_left.nwk")
    )

    assert dumps_newick(tree) == (
        "((A:0.1,B:0.1)66.6666666666667:0.2,(C:0.1,D:0.1)66.6666666666667:0.2);"
    )
    assert report.tree_count == 3


def test_compute_strict_consensus_tree_returns_star_when_no_clade_is_unanimous() -> None:
    tree, report = compute_strict_consensus_tree(fixture("example_tree_set_left.nwk"))

    assert dumps_newick(tree) == "(A:0.1,B:0.1,C:0.1,D:0.1);"
    assert report.tree_count == 3
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.consensus_method == "strict"
    assert report.consensus_threshold == 1.0
    assert report.included_clade_count == 0


def test_compute_consensus_tree_requires_identical_taxon_sets() -> None:
    with pytest.raises(
        InvalidAlignmentError,
        match="share the exact same taxon set",
    ):
        compute_consensus_tree(fixture("example_tree_set_mismatched.nwk"))
    assert (
        bijux_phylogenetics.trim_columns_above_missingness_threshold
        is trim_columns_above_missingness_threshold
    )
    assert bijux_phylogenetics.trim_alignment is trim_alignment
    assert bijux_phylogenetics.translate_coding_alignment is translate_coding_alignment
    assert (
        bijux_phylogenetics.translate_coding_alignment_from_dna_bin_alignment
        is translate_coding_alignment_from_dna_bin_alignment
    )
    assert bijux_phylogenetics.load_dna_bin_alignment is load_dna_bin_alignment
    assert (
        bijux_phylogenetics.write_dna_bin_alignment_fasta
        is write_dna_bin_alignment_fasta
    )
    assert (
        bijux_phylogenetics.compute_pairwise_genetic_distance_matrix_from_dna_bin_alignment
        is compute_pairwise_genetic_distance_matrix_from_dna_bin_alignment
    )
    assert bijux_phylogenetics.TreeRootingReport is TreeRootingReport
    assert bijux_phylogenetics.root_tree_on_outgroup is root_tree_on_outgroup
    assert bijux_phylogenetics.reroot_tree_by_midpoint is reroot_tree_by_midpoint
    assert bijux_phylogenetics.rotate_named_node is rotate_named_node
    assert bijux_phylogenetics.rotate_all_internal_nodes is rotate_all_internal_nodes
    assert bijux_phylogenetics.unroot_tree is unroot_tree
    assert bijux_phylogenetics.write_tree_rooting_report is write_tree_rooting_report
    assert bijux_phylogenetics.render_phylo_inputs_report is render_phylo_inputs_report


def test_compute_tree_distance_matrix_reports_symmetric_rf_pairs() -> None:
    report = compute_tree_distance_matrix(fixture("example_tree_set_left.nwk"))
    assert report.tree_count == 3
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert [
        (
            row.left_index,
            row.right_index,
            row.robinson_foulds_distance,
            row.normalized_robinson_foulds,
        )
        for row in report.pairs
    ] == [
        (1, 1, 0, 0.0),
        (1, 2, 0, 0.0),
        (1, 3, 4, 1.0),
        (2, 2, 0, 0.0),
        (2, 3, 4, 1.0),
        (3, 3, 0, 0.0),
    ]


def test_cluster_trees_by_topology_groups_identical_rooted_signatures() -> None:
    report = cluster_trees_by_topology(fixture("example_tree_set_left.nwk"))
    assert report.tree_count == 3
    assert report.rooted_topology_count == 2
    assert [
        (
            cluster.rooted_topology_id,
            cluster.tree_indices,
            cluster.tree_count,
            cluster.frequency,
            cluster.representative_index,
            cluster.representative_newick,
        )
        for cluster in report.clusters
    ] == [
        (
            "A|B||C|D",
            [1, 2],
            2,
            0.666666666666667,
            1,
            "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);",
        ),
        (
            "A|C||B|D",
            [3],
            1,
            0.333333333333333,
            3,
            "((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);",
        ),
    ]


def test_detect_unstable_taxa_reports_inconsistent_placements() -> None:
    report = detect_unstable_taxa(fixture("example_tree_set_left.nwk"))
    assert [
        (
            row.taxon,
            row.unique_placements,
            row.dominant_frequency,
            row.instability_score,
            [
                (placement.signature, placement.tree_count, placement.frequency)
                for placement in row.placements
            ],
        )
        for row in report.taxa
    ] == [
        (
            "A",
            2,
            0.666666666666667,
            0.333333333333333,
            [("A|B", 2, 0.666666666666667), ("A|C", 1, 0.333333333333333)],
        ),
        (
            "B",
            2,
            0.666666666666667,
            0.333333333333333,
            [("A|B", 2, 0.666666666666667), ("B|D", 1, 0.333333333333333)],
        ),
        (
            "C",
            2,
            0.666666666666667,
            0.333333333333333,
            [("C|D", 2, 0.666666666666667), ("A|C", 1, 0.333333333333333)],
        ),
        (
            "D",
            2,
            0.666666666666667,
            0.333333333333333,
            [("C|D", 2, 0.666666666666667), ("B|D", 1, 0.333333333333333)],
        ),
    ]


def test_detect_unstable_clades_reports_frequencies_and_conflicts() -> None:
    report = detect_unstable_clades(fixture("example_tree_set_left.nwk"))
    assert [
        (
            row.clade,
            row.tree_count,
            row.frequency,
            row.conflict_count,
            row.instability_score,
            row.support_classification,
            row.conflicting_clades,
        )
        for row in report.clades
    ] == [
        (
            "A|B",
            2,
            0.666666666666667,
            2,
            0.333333333333333,
            "intermediate-support",
            ["A|C", "B|D"],
        ),
        (
            "A|C",
            1,
            0.333333333333333,
            2,
            0.333333333333333,
            "intermediate-support",
            ["A|B", "C|D"],
        ),
        (
            "B|D",
            1,
            0.333333333333333,
            2,
            0.333333333333333,
            "intermediate-support",
            ["A|B", "C|D"],
        ),
        (
            "C|D",
            2,
            0.666666666666667,
            2,
            0.333333333333333,
            "intermediate-support",
            ["A|C", "B|D"],
        ),
    ]


def test_compare_posterior_topological_diversity_reports_relative_dispersion() -> None:
    report = compare_posterior_topological_diversity(
        fixture("example_tree_set_left.nwk"),
        fixture("example_tree_set_right.nwk"),
    )

    assert report.left_summary.rooted_topology_count == 2
    assert report.right_summary.rooted_topology_count == 2
    assert report.left_summary.dominant_topology_frequency == 0.666666666666667
    assert report.right_summary.effective_topology_count == 1.889881574842309
    assert (
        report.left_summary.mean_within_set_normalized_robinson_foulds
        == 0.666666666666667
    )
    assert report.warnings == []


def test_detect_posterior_topology_multimodality_reports_multiple_modes() -> None:
    report = detect_posterior_topology_multimodality(
        fixture("example_tree_set_left.nwk"),
        min_mode_frequency=0.3,
    )

    assert report.multimodal is True
    assert report.mode_count == 2
    assert report.dominant_mode_frequency == 0.666666666666667
    assert [mode.tree_indices for mode in report.modes] == [[1, 2], [3]]


def test_summarize_clade_credibility_conflicts_reports_incompatible_high_frequency_clades() -> (
    None
):
    report = summarize_clade_credibility_conflicts(
        fixture("example_tree_set_left.nwk"),
        credibility_threshold=0.3,
    )

    assert report.high_credibility_clade_count == 4
    assert report.conflict_count == 2
    assert report.conflicts[0].combined_frequency == 1.0


def test_summarize_uncertainty_aware_conclusions_separates_robust_uncertain_and_conflicting_clades() -> (
    None
):
    report = summarize_uncertainty_aware_conclusions(
        fixture("example_tree_set_left.nwk"),
        robust_threshold=0.95,
        credibility_threshold=0.3,
    )

    assert report.robust_clade_count == 0
    assert report.uncertain_clade_count == 1
    assert report.conflicting_clade_count == 3
    assert report.uncertain_clades[0].clade == "C|D"
    assert {row.conclusion for row in report.conflicting_clades} == {"conflict-prone"}


def test_write_uncertainty_tables_emit_clusters_conflicts_and_conclusions(
    tmp_path: Path,
) -> None:
    cluster_path = tmp_path / "topology-clusters.tsv"
    conflict_path = tmp_path / "clade-conflicts.tsv"
    conclusion_path = tmp_path / "uncertainty-conclusions.tsv"

    write_topology_cluster_table(
        cluster_path, cluster_trees_by_topology(fixture("example_tree_set_left.nwk"))
    )
    write_clade_credibility_conflict_table(
        conflict_path,
        summarize_clade_credibility_conflicts(
            fixture("example_tree_set_left.nwk"), credibility_threshold=0.3
        ),
    )
    write_uncertainty_conclusion_table(
        conclusion_path,
        summarize_uncertainty_aware_conclusions(
            fixture("example_tree_set_left.nwk"),
            robust_threshold=0.95,
            credibility_threshold=0.3,
        ),
    )

    assert "rooted_topology_id\ttree_indices" in cluster_path.read_text(
        encoding="utf-8"
    )
    assert "left_clade\tleft_frequency\tright_clade" in conflict_path.read_text(
        encoding="utf-8"
    )
    assert "clade\tfrequency\tconclusion" in conclusion_path.read_text(encoding="utf-8")


def test_compare_posterior_tree_sets_reports_clade_deltas_and_cross_set_distance() -> (
    None
):
    report = compare_posterior_tree_sets(
        fixture("example_tree_set_left.nwk"),
        fixture("example_tree_set_right.nwk"),
    )
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.left_tree_count == 3
    assert report.right_tree_count == 3
    assert report.left_rooted_topology_count == 2
    assert report.right_rooted_topology_count == 2
    assert report.shared_rooted_topology_count == 1
    assert report.mean_between_set_robinson_foulds == 3.111111111111111
    assert report.mean_between_set_normalized_robinson_foulds == 0.777777777777778
    assert [
        (row.clade, row.left_frequency, row.right_frequency, row.delta)
        for row in report.clade_frequency_deltas
    ] == [
        ("A|B", 0.666666666666667, 0.333333333333333, -0.333333333333333),
        ("A|C", 0.333333333333333, 0.0, -0.333333333333333),
        ("A|D", 0.0, 0.666666666666667, 0.666666666666667),
        ("B|C", 0.0, 0.666666666666667, 0.666666666666667),
        ("B|D", 0.333333333333333, 0.0, -0.333333333333333),
        ("C|D", 0.666666666666667, 0.333333333333333, -0.333333333333333),
    ]


def test_compare_bootstrap_and_posterior_uncertainty_reports_conflicting_clades(
    tmp_path: Path,
) -> None:
    bootstrap_tree = tmp_path / "bootstrap-support.nwk"
    bootstrap_tree.write_text(
        "((A:0.1,B:0.1)95:0.2,(C:0.1,D:0.1)60:0.2);\n", encoding="utf-8"
    )

    report = compare_bootstrap_and_posterior_uncertainty(
        bootstrap_tree,
        fixture("example_tree_set_right.nwk"),
    )

    rows = {row.clade: row for row in report.rows}
    assert report.posterior_tree_count == 3
    assert report.topology_mismatch_detected is True
    assert report.topology_mismatch_clade_count == 2
    assert rows["A|B"].agreement == "strong_conflict"
    assert rows["C|D"].posterior_frequency == 0.333333333333333


def test_compare_bootstrap_and_posterior_uncertainty_uses_ufboot_from_composite_labels(
    tmp_path: Path,
) -> None:
    bootstrap_tree = tmp_path / "bootstrap-support.treefile"
    bootstrap_tree.write_text(
        "((A:0.1,B:0.1)82/97:0.2,(C:0.1,D:0.1)79/96:0.2);\n",
        encoding="utf-8",
    )

    report = compare_bootstrap_and_posterior_uncertainty(
        bootstrap_tree,
        fixture("example_tree_set_right.nwk"),
    )

    rows = {row.clade: row for row in report.rows}
    assert rows["A|B"].bootstrap_support == 0.97
    assert rows["C|D"].bootstrap_support == 0.96
    assert report.topology_mismatch_detected is True


def test_render_tree_set_comparison_report_embeds_tree_set_differences(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "tree-set-comparison.html"
    result = render_tree_set_comparison_report(
        left_tree_set_path=fixture("example_tree_set_left.nwk"),
        right_tree_set_path=fixture("example_tree_set_right.nwk"),
        out_path=output_path,
    )

    assert result.output_path.exists()
    assert result.report_kind == "tree-set-comparison"
    assert result.shared_rooted_topology_count == 1
    assert "left-topology-clusters" in result.machine_manifest["sections"]
    assert "topological-diversity-comparison" in result.machine_manifest["sections"]
    assert "left-topology-multimodality" in result.machine_manifest["sections"]
    assert "right-uncertainty-aware-conclusions" in result.machine_manifest["sections"]


def test_taxon_labels_preserve_raw_names_and_normalized_keys() -> None:
    tree = loads_newick("('Homo sapiens':0.1,'NCBI|123/45':0.2,'A.B-1':0.3);")
    assert [(taxon.raw, taxon.key) for taxon in tree.tip_taxa] == [
        ("Homo sapiens", "Homo_sapiens"),
        ("NCBI|123/45", "NCBI_123_45"),
        ("A.B-1", "A.B-1"),
    ]
    assert tree.branch_lengths() == [0.1, 0.2, 0.3]
    assert tree.terminal_branch_lengths() == [
        ("Homo sapiens", 0.1),
        ("NCBI|123/45", 0.2),
        ("A.B-1", 0.3),
    ]
    assert dumps_newick(tree) == "(A.B-1:0.3,'Homo sapiens':0.1,'NCBI|123/45':0.2);"


def test_normalize_tree_taxa_reports_rename_mapping() -> None:
    tree = loads_newick("('Homo sapiens':0.1,'Mus musculus':0.2,A:0.3);")
    normalized_tree, report = normalize_tree_taxa(tree, policy="spaces-to-underscores")
    assert normalized_tree.tip_names == ["Homo_sapiens", "Mus_musculus", "A"]
    assert [
        (rename.raw_label, rename.normalized_label) for rename in report.renamed_taxa
    ] == [
        ("Homo sapiens", "Homo_sapiens"),
        ("Mus musculus", "Mus_musculus"),
    ]
    assert report.original_tip_count == 3
    assert report.normalized_tip_count == 3
    assert report.unchanged_taxa == ["A"]
    assert report.topology_preserved is True
    assert report.branch_lengths_preserved is True


def test_taxon_safety_reports_unsafe_labels_and_normalization_collisions(
    tmp_path: Path,
) -> None:
    tree = loads_newick(
        "('Homo sapiens':0.1,Homo_sapiens:0.2,'NCBI/123':0.3,'Quoted''Name':0.4,A:0.5);"
    )
    report = inspect_tree_taxa_safety(tree, policy="spaces-to-underscores")
    assert [
        (entry.raw_label, entry.normalized_label, entry.reasons)
        for entry in report.unsafe_taxa
    ] == [
        (
            "Homo sapiens",
            "Homo_sapiens",
            ["contains whitespace", "collides with another label after normalization"],
        ),
        (
            "Homo_sapiens",
            "Homo_sapiens",
            ["collides with another label after normalization"],
        ),
        ("NCBI/123", "NCBI/123", ["contains slash characters"]),
        ("Quoted'Name", "Quoted'Name", ["contains quote characters"]),
    ]
    assert [
        (entry.normalized_label, entry.raw_labels) for entry in report.collisions
    ] == [("Homo_sapiens", ["Homo sapiens", "Homo_sapiens"])]

    mapping_path = tmp_path / "taxon-mapping.tsv"
    write_taxon_mapping(
        mapping_path,
        normalize_tree_taxa(tree, policy="spaces-to-underscores")[1].renamed_taxa,
    )
    assert mapping_path.read_text(encoding="utf-8") == (
        "raw_label\tnormalized_label\nHomo sapiens\tHomo_sapiens\n"
    )


def test_metadata_inspect_reports_taxon_contract() -> None:
    report = inspect_metadata_table(fixture("example_metadata.tsv"))
    assert report.format == "tsv"
    assert report.row_count == 4
    assert report.column_count == 3
    assert report.taxon_column == "taxon"
    assert report.taxa == ["A", "B", "C", "D"]
    assert [
        (row.name, row.missing_count, row.completeness_fraction)
        for row in report.column_completeness
    ] == [
        ("taxon", 0, 1.0),
        ("species", 0, 1.0),
        ("location", 0, 1.0),
    ]


def test_join_table_to_taxa_returns_tip_by_tip_metadata_rows() -> None:
    report = join_table_to_taxa(["A", "B", "Z"], fixture("example_metadata.tsv"))
    assert [
        (row.taxon, row.matched, row.values.get("species", ""))
        for row in report.joined_rows
    ] == [
        ("A", True, "Alpha species"),
        ("B", True, "Beta species"),
        ("Z", False, ""),
    ]
    assert report.missing_from_metadata == ["Z"]
    assert report.extra_metadata_taxa == ["C", "D"]


def test_inspect_environment_reports_available_and_optional_dependencies() -> None:
    report = inspect_environment()
    status_by_name = {item.name: item for item in report.dependencies}
    assert report.python_version
    assert report.host_platform
    assert status_by_name["biopython"].available is True
    assert "dendropy" in status_by_name


def test_metadata_inspect_rejects_duplicate_taxa() -> None:
    try:
        inspect_metadata_table(fixture("example_metadata_duplicate.tsv"))
    except MetadataJoinError as error:
        assert error.code == "metadata_join_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected MetadataJoinError")


def test_metadata_inspect_rejects_missing_requested_taxon_column() -> None:
    try:
        inspect_metadata_table(
            fixture("example_metadata_missing_taxon.csv"), taxon_column="taxon"
        )
    except MetadataJoinError as error:
        assert error.code == "metadata_join_error"
        assert error.message == "missing taxon column 'taxon'"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected MetadataJoinError")


def test_traits_validate_infers_numeric_and_categorical_schema() -> None:
    report = validate_traits_table(fixture("example_traits_validate.tsv"))
    assert report.taxon_column == "taxon"
    assert [
        (column.name, column.kind, column.missing_count, column.missing_fraction)
        for column in report.trait_columns
    ] == [
        ("height_cm", "numeric", 0, 0.0),
        ("habitat", "categorical", 0, 0.0),
        ("status", "categorical", 1, 0.25),
    ]


def test_traits_validate_can_distinguish_binary_and_text_columns() -> None:
    report = validate_traits_table(fixture("example_traits_schema.tsv"))
    assert [(column.name, column.kind) for column in report.trait_columns] == [
        ("height_cm", "numeric"),
        ("presence", "binary"),
        ("comment", "text"),
        ("habitat", "categorical"),
    ]


def test_traits_detect_unusable_columns_by_missingness() -> None:
    columns = detect_unusable_trait_columns(
        fixture("example_traits_validate.tsv"),
        missingness_threshold=0.2,
    )
    assert [(column.name, column.missing_fraction) for column in columns] == [
        ("status", 0.25)
    ]


def test_detect_missing_trait_values_reports_taxon_and_column() -> None:
    report = detect_missing_trait_values(fixture("example_traits_validate.tsv"))
    assert [(item.taxon, item.trait) for item in report.missing_values] == [
        ("C", "status")
    ]


def test_traits_link_reports_mismatch_and_usable_taxa() -> None:
    report = link_tree_to_traits(
        fixture("example_tree.nwk"), fixture("example_traits.tsv")
    )
    assert report.tree_taxa == 4
    assert report.trait_taxa == 4
    assert report.linked_taxa == 3
    assert report.usable_taxa == ["A", "B", "C"]
    assert report.missing_from_traits == ["D"]
    assert report.extra_trait_taxa == ["E"]


def test_traits_link_strict_mode_rejects_mismatch() -> None:
    try:
        link_tree_to_traits(
            fixture("example_tree.nwk"), fixture("example_traits.tsv"), strict=True
        )
    except MetadataJoinError as error:
        assert error.code == "metadata_join_error"
        assert error.details["failure_reason"] == "tree_trait_taxon_mismatch"
        assert error.details["missing_from_traits"] == ["D"]
        assert error.details["extra_trait_taxa"] == ["E"]
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected MetadataJoinError")


def test_prune_traits_to_tree_keeps_tree_order_for_overlapping_taxa() -> None:
    rows, report = prune_traits_to_tree(
        fixture("example_tree.nwk"),
        fixture("example_traits.tsv"),
    )
    assert [row["taxon"] for row in rows] == ["A", "B", "C"]
    assert report.original_row_count == 4
    assert report.kept_taxa == ["A", "B", "C"]
    assert report.removed_taxa == ["E"]


def test_traits_prune_cli_writes_pruned_table_in_tree_order(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "pruned-traits.tsv"
    exit_code = main(
        [
            "traits",
            "prune",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits.tsv")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "taxon\tvalue\nA\t1.2\nB\t1.4\nC\t1.8\n"
    )
    assert payload["data"]["kept_taxa"] == ["A", "B", "C"]
    assert payload["data"]["removed_taxa"] == ["E"]


def test_traits_missing_cli_reports_taxon_and_column(capsys) -> None:
    exit_code = main(
        ["traits", "missing", str(fixture("example_traits_validate.tsv")), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["missing_value_count"] == 1
    assert payload["data"]["missing_values"] == [{"taxon": "C", "trait": "status"}]


def test_dataset_readiness_reports_ready_comparative_inputs() -> None:
    report = summarize_dataset_readiness(
        fixture("example_tree.nwk"),
        fixture("example_metadata.tsv"),
        fixture("example_traits_validate.tsv"),
    )
    assert report.ready_for_comparative_analysis is True
    assert report.analysis_taxa == ["A", "B", "C", "D"]
    assert report.missing_metadata_taxa == []
    assert report.missing_trait_taxa == []
    assert report.metadata_only_taxa == []
    assert report.trait_only_taxa == []
    assert report.unusable_trait_columns == []
    assert report.blockers == []
    assert report.warnings == []


def test_dataset_readiness_reports_linkage_blockers() -> None:
    report = summarize_dataset_readiness(
        fixture("example_tree.nwk"),
        fixture("example_metadata.tsv"),
        fixture("example_traits.tsv"),
    )
    assert report.ready_for_comparative_analysis is False
    assert report.analysis_taxa == ["A", "B", "C"]
    assert report.missing_trait_taxa == ["D"]
    assert report.trait_only_taxa == ["E"]
    assert report.blockers == ["trait table is missing one or more tree taxa"]
    assert report.warnings == ["trait table contains taxa absent from the tree"]


def test_prune_tree_to_taxa_writes_expected_tip_set() -> None:
    tree, report = prune_tree_to_taxa(
        fixture("example_tree.nwk"), fixture("example_traits.tsv")
    )
    assert tree.tip_names == ["A", "B", "C"]
    assert dumps_newick(tree) == "((A:0.1,B:0.1):0.2,C:0.3);"
    assert report.kept_taxa == ["A", "B", "C"]
    assert report.removed_taxa == ["D"]
    assert [(row.taxon, row.reason) for row in report.removed_taxa_with_reasons] == [
        ("D", "absent_from_keep_table")
    ]
    assert report.summary.removed_taxa == ["D"]


def test_prune_tree_to_requested_taxa_reports_absent_requests() -> None:
    tree, report = prune_tree_to_requested_taxa(
        fixture("example_tree.nwk"),
        ["A", "C", "Z"],
    )
    assert tree.tip_names == ["A", "C"]
    assert dumps_newick(tree) == "(A:0.3,C:0.3);"
    assert tree.rooted is True
    assert report.requested_taxa == ["A", "C", "Z"]
    assert report.kept_taxa == ["A", "C"]
    assert report.removed_taxa == ["B", "D"]
    assert report.absent_requested_taxa == ["Z"]
    assert [(row.taxon, row.reason) for row in report.removed_taxa_with_reasons] == [
        ("B", "not_requested"),
        ("D", "not_requested"),
    ]
    assert report.pruning_audit.root_to_tip_complete is True
    assert report.pruning_audit.unary_internal_nodes == []


def test_prune_tree_to_requested_taxa_ignores_input_order() -> None:
    tree, report = prune_tree_to_requested_taxa(
        fixture("example_tree.nwk"),
        ["C", "A"],
    )

    assert tree.tip_names == ["A", "C"]
    assert dumps_newick(tree) == "(A:0.3,C:0.3);"
    assert tree.rooted is True
    assert report.requested_taxa == ["A", "C"]
    assert report.kept_taxa == ["A", "C"]
    assert report.removed_taxa == ["B", "D"]


def test_prune_tree_to_requested_taxa_marks_two_tip_outputs_rooted_like_ape() -> None:
    tree, report = prune_tree_to_requested_taxa(
        fixture("example_tree_unrooted.nwk"),
        ["A", "B"],
    )

    assert dumps_newick(tree) == "(A:0.1,B:0.2);"
    assert tree.rooted is True
    assert report.kept_taxa == ["A", "B"]
    assert report.removed_taxa == ["C", "D"]


def test_prune_tree_to_requested_taxa_fails_when_fewer_than_two_taxa_remain() -> None:
    with pytest.raises(ValueError, match="at least two retained taxa"):
        prune_tree_to_requested_taxa(
            fixture("example_tree.nwk"),
            ["A"],
        )


def test_drop_tree_taxa_excludes_exact_requested_tips() -> None:
    tree, report = drop_tree_taxa(
        fixture("example_tree.nwk"),
        ["B", "D", "Z"],
    )
    assert tree.tip_names == ["A", "C"]
    assert dumps_newick(tree) == "(A:0.3,C:0.3);"
    assert tree.rooted is True
    assert report.requested_taxa == ["B", "D", "Z"]
    assert report.kept_taxa == ["A", "C"]
    assert report.removed_taxa == ["B", "D"]
    assert report.absent_requested_taxa == ["Z"]
    assert [(row.taxon, row.reason) for row in report.removed_taxa_with_reasons] == [
        ("B", "excluded_explicitly"),
        ("D", "excluded_explicitly"),
    ]


def test_drop_tree_taxa_keeps_unrooted_three_tip_outputs_unrooted() -> None:
    tree, report = drop_tree_taxa(
        fixture("example_tree_unrooted.nwk"),
        ["D"],
    )

    assert dumps_newick(tree) == "(A:0.1,B:0.2,C:0.3);"
    assert tree.rooted is False
    assert report.kept_taxa == ["A", "B", "C"]
    assert report.removed_taxa == ["D"]


def test_drop_tree_taxa_marks_two_tip_outputs_rooted_like_ape() -> None:
    tree, report = drop_tree_taxa(
        fixture("example_tree_unrooted.nwk"),
        ["C", "D"],
    )

    assert dumps_newick(tree) == "(A:0.1,B:0.2);"
    assert tree.rooted is True
    assert report.kept_taxa == ["A", "B"]
    assert report.removed_taxa == ["C", "D"]


def test_drop_tree_taxa_fails_when_fewer_than_two_taxa_remain() -> None:
    with pytest.raises(ValueError, match="at least two retained taxa"):
        drop_tree_taxa(
            fixture("example_tree.nwk"),
            ["B", "C", "D"],
        )


def test_prune_cli_accepts_explicit_taxon_keep_lists(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "pruned-tree.nwk"
    pruned_taxa_path = tmp_path / "removed.tsv"
    exit_code = main(
        [
            "prune",
            str(fixture("example_tree.nwk")),
            "--taxa",
            "A",
            "C",
            "Z",
            "--out",
            str(output_path),
            "--pruned-taxa-out",
            str(pruned_taxa_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == "(A:0.3,C:0.3);\n"
    assert pruned_taxa_path.read_text(encoding="utf-8") == "taxon\nB\nD\n"
    assert payload["data"]["absent_requested_taxa"] == ["Z"]
    assert payload["data"]["kept_taxa"] == ["A", "C"]


def test_prune_cli_accepts_explicit_taxon_exclusion_lists(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "pruned-tree.nwk"
    pruned_taxa_path = tmp_path / "removed.tsv"
    exit_code = main(
        [
            "prune",
            str(fixture("example_tree.nwk")),
            "--exclude-taxa",
            "B",
            "D",
            "Z",
            "--out",
            str(output_path),
            "--pruned-taxa-out",
            str(pruned_taxa_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == "(A:0.3,C:0.3);\n"
    assert pruned_taxa_path.read_text(encoding="utf-8") == "taxon\nB\nD\n"
    assert payload["data"]["absent_requested_taxa"] == ["Z"]
    assert payload["data"]["removed_taxa"] == ["B", "D"]


def test_extract_named_clade_returns_exact_descendant_subtree() -> None:
    tree, report = extract_named_clade(
        fixture("example_tree_named_clades.nwk"),
        clade_name="Mammals",
    )
    assert tree.tip_names == ["A", "B"]
    assert dumps_newick(tree) == "(A:0.1,B:0.1)Mammals;"
    assert report.clade_name == "Mammals"
    assert report.tip_count == 2
    assert report.taxa == ["A", "B"]
    assert report.retained_all_requested_descendants is True
    assert report.missing_requested_descendants == []
    assert report.unexpected_retained_taxa == []
    assert report.summary.removed_taxa == ["C", "D"]


def test_extract_tree_clade_by_node_id_returns_internal_subtree() -> None:
    tree, report = extract_tree_clade_by_node_id(
        fixture("example_tree_named_clades.nwk"),
        node_id=6,
    )

    assert tree.tip_names == ["A", "B"]
    assert dumps_newick(tree) == "(A:0.1,B:0.1)Mammals;"
    assert report.selector_kind == "node-id"
    assert report.requested_node_id == 6
    assert report.matched_node_id == 6
    assert report.matched_node_name == "Mammals"
    assert report.taxa == ["A", "B"]
    assert report.summary.removed_taxa == ["C", "D"]


def test_extract_tree_clade_by_node_id_returns_root_clade() -> None:
    tree, report = extract_tree_clade_by_node_id(
        fixture("example_tree_named_clades.nwk"),
        node_id=5,
    )

    assert tree.tip_names == ["A", "B", "C", "D"]
    assert dumps_newick(tree) == (
        "((A:0.1,B:0.1)Mammals:0.2,(C:0.2,D:0.2)Birds:0.1)Root;"
    )
    assert report.matched_node_name == "Root"
    assert report.summary.removed_taxa == []


def test_extract_tree_clade_by_descendant_taxa_returns_matching_subtree() -> None:
    tree, report = extract_tree_clade_by_descendant_taxa(
        fixture("example_tree_named_clades.nwk"),
        descendant_taxa=["B", "A"],
    )

    assert tree.tip_names == ["A", "B"]
    assert dumps_newick(tree) == "(A:0.1,B:0.1)Mammals;"
    assert report.selector_kind == "descendant-taxa"
    assert report.requested_taxa == ["A", "B"]
    assert report.matched_node_id == 6
    assert report.matched_node_name == "Mammals"


def test_extract_tree_clade_by_descendant_taxa_returns_root_clade() -> None:
    tree, report = extract_tree_clade_by_descendant_taxa(
        fixture("example_tree_named_clades.nwk"),
        descendant_taxa=["A", "B", "C", "D"],
    )

    assert tree.tip_names == ["A", "B", "C", "D"]
    assert report.matched_node_id == 5
    assert report.matched_node_name == "Root"


def test_extract_tree_clade_by_node_id_rejects_tip_nodes() -> None:
    with pytest.raises(ValueError, match="greater than the number of tips"):
        extract_tree_clade_by_node_id(
            fixture("example_tree_named_clades.nwk"),
            node_id=1,
        )


def test_extract_tree_clade_by_node_id_rejects_out_of_bounds_nodes() -> None:
    with pytest.raises(IndexError, match="out of bounds"):
        extract_tree_clade_by_node_id(
            fixture("example_tree_named_clades.nwk"),
            node_id=8,
        )


def test_extract_tree_clade_by_descendant_taxa_rejects_non_monophyletic_sets() -> None:
    with pytest.raises(ValueError, match="do not define an internal clade"):
        extract_tree_clade_by_descendant_taxa(
            fixture("example_tree_named_clades.nwk"),
            descendant_taxa=["A", "C"],
        )


def test_find_tree_mrca_returns_internal_node_for_two_tips() -> None:
    report = find_tree_mrca(
        fixture("example_tree.nwk"),
        taxa=["A", "B"],
    )

    assert report.matched_node_id == 6
    assert report.matched_node_name is None
    assert report.matched_taxa == ["A", "B"]
    assert report.matched_extra_taxa == []
    assert report.is_root is False


def test_find_tree_mrca_returns_root_for_full_tip_set() -> None:
    report = find_tree_mrca(
        fixture("example_tree.nwk"),
        taxa=["A", "B", "C", "D"],
    )

    assert report.matched_node_id == 5
    assert report.matched_taxa == ["A", "B", "C", "D"]
    assert report.is_root is True


def test_find_tree_mrca_handles_duplicate_requested_tips_explicitly() -> None:
    report = find_tree_mrca(
        fixture("example_tree.nwk"),
        taxa=["A", "A", "B"],
    )

    assert report.requested_taxa == ["A", "A", "B"]
    assert report.unique_requested_taxa == ["A", "B"]
    assert report.duplicate_requested_taxa == ["A"]
    assert report.matched_node_id == 6


def test_find_tree_mrca_matches_many_tip_request_on_pectinate_tree() -> None:
    report = find_tree_mrca(
        fixture("example_tree_ladderized.nwk"),
        taxa=["A", "B", "C"],
    )

    assert report.matched_node_id == 6
    assert report.matched_taxa == ["A", "B", "C"]
    assert report.matched_extra_taxa == []


def test_find_tree_mrca_matches_polytomy_case() -> None:
    report = find_tree_mrca(
        fixture("example_tree_polytomy.nwk"),
        taxa=["A", "B", "C"],
    )

    assert report.matched_node_id == 6
    assert report.matched_taxa == ["A", "B", "C"]
    assert report.is_root is False


def test_find_tree_mrca_rejects_missing_tips() -> None:
    with pytest.raises(ValueError, match="requested taxa are not present in the tree: Z"):
        find_tree_mrca(
            fixture("example_tree.nwk"),
            taxa=["A", "Z"],
        )


def test_find_tree_mrca_requires_two_distinct_tips() -> None:
    with pytest.raises(ValueError, match="mrca requires at least two distinct taxa"):
        find_tree_mrca(
            fixture("example_tree.nwk"),
            taxa=["A", "A"],
        )


def test_find_tree_mrca_works_after_pruning(tmp_path: Path) -> None:
    pruned_tree, _report = prune_tree_to_requested_taxa(
        fixture("example_tree_rooted_on_d.nwk"),
        ["A", "B", "C"],
    )
    pruned_path = tmp_path / "pruned-tree.nwk"
    write_newick(pruned_path, pruned_tree)

    report = find_tree_mrca(pruned_path, taxa=["A", "B"])

    assert report.matched_node_id == 5
    assert report.matched_taxa == ["A", "B"]


def test_find_tree_mrca_works_after_rooting(tmp_path: Path) -> None:
    rooted_tree, _report = root_tree_on_outgroup(
        fixture("example_tree_rootable.nwk"),
        outgroup_taxa=["D"],
    )
    rooted_path = tmp_path / "rooted-tree.nwk"
    write_newick(rooted_path, rooted_tree)

    report = find_tree_mrca(rooted_path, taxa=["A", "B", "C"])

    assert report.matched_node_id == 6
    assert report.matched_taxa == ["A", "B", "C"]


def test_assess_tree_monophyly_matches_rooted_two_tip_clade() -> None:
    report = assess_tree_monophyly(
        fixture("example_tree.nwk"),
        taxa=["A", "B"],
    )

    assert report.monophyletic is True
    assert report.complementary_clade_used is False
    assert report.matched_node_id == 6
    assert report.matched_taxa == ["A", "B"]
    assert report.matched_extra_taxa == []
    assert report.matched_tip_count == 2
    assert report.is_root is False


def test_assess_tree_monophyly_reports_root_extra_taxa_for_non_monophyletic_group() -> None:
    report = assess_tree_monophyly(
        fixture("example_tree.nwk"),
        taxa=["A", "C"],
    )

    assert report.monophyletic is False
    assert report.complementary_clade_used is False
    assert report.matched_node_id == 5
    assert report.matched_taxa == ["A", "B", "C", "D"]
    assert report.matched_extra_taxa == ["B", "D"]
    assert report.is_root is True


def test_assess_tree_monophyly_treats_full_tip_set_as_monophyletic() -> None:
    report = assess_tree_monophyly(
        fixture("example_tree.nwk"),
        taxa=["A", "B", "C", "D"],
    )

    assert report.monophyletic is True
    assert report.matched_node_id == 5
    assert report.is_root is True


def test_assess_tree_monophyly_treats_singletons_as_monophyletic() -> None:
    report = assess_tree_monophyly(
        fixture("example_tree.nwk"),
        taxa=["A"],
    )

    assert report.monophyletic is True
    assert report.matched_node_id == 1
    assert report.matched_taxa == ["A"]
    assert report.matched_extra_taxa == []
    assert report.matched_tip_count == 1
    assert report.is_root is False


def test_assess_tree_monophyly_reports_missing_taxa_explicitly() -> None:
    report = assess_tree_monophyly(
        fixture("example_tree.nwk"),
        taxa=["A", "Z"],
    )

    assert report.monophyletic is True
    assert report.present_requested_taxa == ["A"]
    assert report.missing_requested_taxa == ["Z"]
    assert report.matched_node_id == 1


def test_assess_tree_monophyly_explicitly_controls_reroot_policy() -> None:
    default_report = assess_tree_monophyly(
        fixture("example_tree.nwk"),
        taxa=["A", "B", "C"],
    )
    rerooted_report = assess_tree_monophyly(
        fixture("example_tree.nwk"),
        taxa=["A", "B", "C"],
        reroot=True,
    )

    assert default_report.monophyletic is False
    assert default_report.complementary_clade_used is False
    assert rerooted_report.monophyletic is True
    assert rerooted_report.complementary_clade_used is True
    assert rerooted_report.matched_extra_taxa == ["D"]


def test_assess_tree_monophyly_matches_unrooted_default_reroot_false_policy() -> None:
    report = assess_tree_monophyly(
        fixture("example_tree_unrooted.nwk"),
        taxa=["A", "B"],
    )

    assert report.monophyletic is False
    assert report.rooted is False
    assert report.matched_extra_taxa == ["C", "D"]


def test_assess_tree_monophyly_matches_unrooted_reroot_true_policy() -> None:
    report = assess_tree_monophyly(
        fixture("example_tree_unrooted.nwk"),
        taxa=["A", "B", "C"],
        reroot=True,
    )

    assert report.monophyletic is True
    assert report.complementary_clade_used is True
    assert report.matched_extra_taxa == ["D"]


def test_assess_tree_monophyly_rejects_all_missing_taxa_when_rerooting() -> None:
    with pytest.raises(ValueError, match="specified outgroup not in labels of the tree"):
        assess_tree_monophyly(
            fixture("example_tree_unrooted.nwk"),
            taxa=["Z"],
            reroot=True,
        )


def test_assess_tree_monophyly_works_after_pruning(tmp_path: Path) -> None:
    pruned_tree, _report = prune_tree_to_requested_taxa(
        fixture("example_tree_rooted_on_d.nwk"),
        ["A", "B", "C"],
    )
    pruned_path = tmp_path / "pruned-tree.nwk"
    write_newick(pruned_path, pruned_tree)

    report = assess_tree_monophyly(pruned_path, taxa=["A", "B"])

    assert report.monophyletic is True
    assert report.matched_taxa == ["A", "B"]


def test_assess_tree_monophyly_works_after_rooting(tmp_path: Path) -> None:
    rooted_tree, _report = root_tree_on_outgroup(
        fixture("example_tree_rootable.nwk"),
        outgroup_taxa=["D"],
    )
    rooted_path = tmp_path / "rooted-tree.nwk"
    write_newick(rooted_path, rooted_tree)

    report = assess_tree_monophyly(
        rooted_path,
        taxa=["A", "B", "C"],
    )

    assert report.monophyletic is True
    assert report.matched_taxa == ["A", "B", "C"]


def test_rotate_named_node_reverses_child_order_without_changing_topology() -> None:
    tree, report = rotate_named_node(
        fixture("example_tree_named_clades.nwk"),
        clade_name="Mammals",
    )
    assert tree.tip_names == ["B", "A", "C", "D"]
    assert report.strategy == "rotate:Mammals"
    assert report.tip_order == ["B", "A", "C", "D"]
    assert report.rooted_topology_preserved is True
    assert report.unrooted_topology_preserved is True


def test_rotate_all_internal_nodes_reverses_every_internal_child_order() -> None:
    tree, report = rotate_all_internal_nodes(fixture("example_tree_named_clades.nwk"))
    assert tree.tip_names == ["D", "C", "B", "A"]
    assert report.strategy == "rotate-all"
    assert report.tip_order == ["D", "C", "B", "A"]
    assert report.rooted_topology_preserved is True
    assert report.unrooted_topology_preserved is True


def test_collapse_branches_below_length_turns_short_internal_edges_into_polytomies() -> (
    None
):
    tree, report = collapse_branches_below_length(
        fixture("example_tree_collapse_threshold.nwk"),
        threshold=0.05,
    )
    assert tree.tip_names == ["A", "B", "C", "D"]
    assert dumps_newick(tree) == "((A:0.1,B:0.1,C:0.2):0.3,D:0.4);"
    assert report.threshold == 0.05
    assert report.collapsed_clades == ["A|B"]
    assert report.topology_preserved is False
    assert report.summary.branch_lengths_affected == ["A|B"]


def test_ladderize_tree_orders_larger_subtrees_first() -> None:
    tree, report = ladderize_tree(fixture("example_tree_ordering.nwk"))
    assert tree.tip_names == ["X", "Y", "Z", "A", "B"]
    assert [
        len(child.children) if child.children else 1 for child in tree.root.children
    ] == [3, 2]
    assert report.strategy == "ladderize"
    assert report.tip_order == ["X", "Y", "Z", "A", "B"]
    assert report.rooted_topology_preserved is True
    assert report.unrooted_topology_preserved is True


def test_sort_tree_tips_alphabetically_preserves_topology_with_stable_tip_order() -> (
    None
):
    tree, report = sort_tree_tips_alphabetically(fixture("example_tree_ordering.nwk"))
    assert tree.tip_names == ["A", "B", "X", "Y", "Z"]
    assert [
        len(child.children) if child.children else 1 for child in tree.root.children
    ] == [2, 3]
    assert report.strategy == "alphabetical"
    assert report.tip_order == ["A", "B", "X", "Y", "Z"]
    assert report.rooted_topology_preserved is True
    assert report.unrooted_topology_preserved is True


def test_root_tree_on_outgroup_reports_absent_taxa_and_roots_tree() -> None:
    tree, report = root_tree_on_outgroup(
        fixture("example_tree_rootable.nwk"),
        outgroup_taxa=["D", "Z"],
    )
    assert tree.tip_names == ["A", "B", "C", "D"]
    assert tree.root.children
    assert report.strategy == "outgroup"
    assert report.requested_taxa == ["D", "Z"]
    assert report.matched_taxa == ["D"]
    assert report.absent_taxa == ["Z"]
    assert report.ingroup_taxa == ["A", "B", "C"]
    assert report.outgroup_monophyletic is True
    assert report.outgroup_mrca_taxa == ["D"]
    assert report.outgroup_mrca_extra_taxa == []
    assert report.rooted_outgroup_taxa == ["D"]
    assert report.rooted_ingroup_taxa == ["A", "B", "C"]
    assert report.warnings == [
        "one or more requested outgroup taxa were absent from the input tree"
    ]
    assert dumps_newick(tree) == "(((A:0.2,B:0.2):0.7,C:0.1):0,D:0.1);"
    assert report.summary.retained_taxa == ["A", "B", "C", "D"]
    assert report.summary.removed_taxa == []


def test_root_tree_on_outgroup_reports_monophyletic_outgroup_clade() -> None:
    tree, report = root_tree_on_outgroup(
        fixture("example_tree_rootable.nwk"),
        outgroup_taxa=["C", "D"],
    )

    assert dumps_newick(tree) == "((A:0.2,B:0.2):0.7,(C:0.1,D:0.1):0);"
    assert report.matched_taxa == ["C", "D"]
    assert report.absent_taxa == []
    assert report.ingroup_taxa == ["A", "B"]
    assert report.outgroup_monophyletic is True
    assert report.outgroup_mrca_taxa == ["C", "D"]
    assert report.outgroup_mrca_extra_taxa == []
    assert report.rooted_outgroup_taxa == ["C", "D"]
    assert report.rooted_ingroup_taxa == ["A", "B"]
    assert report.warnings == []


def test_root_tree_on_outgroup_rejects_non_monophyletic_outgroup_clade() -> None:
    with pytest.raises(TreeRootingError) as error:
        root_tree_on_outgroup(
            fixture("example_tree_rootable.nwk"),
            outgroup_taxa=["B", "D"],
        )

    assert error.value.code == "outgroup_not_monophyletic"
    assert str(error.value) == "requested outgroup taxa are not monophyletic in the input tree"
    assert error.value.details["matched_taxa"] == ["B", "D"]
    assert error.value.details["outgroup_mrca_taxa"] == ["A", "B", "C", "D"]
    assert error.value.details["outgroup_mrca_extra_taxa"] == ["A", "C"]


def test_root_tree_on_outgroup_preserves_already_rooted_outgroup_tree() -> None:
    tree, report = root_tree_on_outgroup(
        fixture("example_tree_rooted_on_d.nwk"),
        outgroup_taxa=["D"],
    )

    assert dumps_newick(tree) == "(((A:0.2,B:0.2):0.7,C:0.1):0,D:0.1);"
    assert report.matched_taxa == ["D"]
    assert report.rooted_outgroup_taxa == ["D"]
    assert report.rooted_ingroup_taxa == ["A", "B", "C"]
    assert report.warnings == []


def test_root_tree_on_outgroup_stays_native_without_biopython_bridge(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import bijux_phylogenetics.core.topology as topology

    def _bridge_unavailable(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("biopython bridge should not be used for outgroup rooting")

    monkeypatch.setattr(topology, "tree_to_biophylo", _bridge_unavailable)
    monkeypatch.setattr(topology, "tree_from_biophylo", _bridge_unavailable)

    tree, report = topology.root_tree_on_outgroup(
        fixture("example_tree_rootable.nwk"),
        outgroup_taxa=["D"],
    )

    assert dumps_newick(tree) == "(((A:0.2,B:0.2):0.7,C:0.1):0,D:0.1);"
    assert report.rooted_outgroup_taxa == ["D"]


def test_write_tree_rooting_report_writes_monophyly_and_root_partition_fields(
    tmp_path: Path,
) -> None:
    _tree, report = root_tree_on_outgroup(
        fixture("example_tree_rootable.nwk"),
        outgroup_taxa=["C", "D"],
    )
    output_path = tmp_path / "rooting.tsv"

    write_tree_rooting_report(output_path, report)

    text = output_path.read_text(encoding="utf-8")
    assert "outgroup_monophyletic" in text
    assert "rooted_outgroup_taxa" in text
    assert "\ttrue\tC,D\t\tC,D\tA,B\t" in text


def test_write_tree_rooting_report_writes_midpoint_position_fields(
    tmp_path: Path,
) -> None:
    _tree, report = reroot_tree_by_midpoint(fixture("example_tree_rootable.nwk"))
    output_path = tmp_path / "midpoint.tsv"

    write_tree_rooting_report(output_path, report)

    text = output_path.read_text(encoding="utf-8")
    assert "midpoint_anchor_taxa" in text
    assert "midpoint_position_kind" in text
    assert "\tA,C\t1.0\t0.5\tnode\tA,B\tC,D\ttrue\n" in text


def test_reroot_tree_by_midpoint_preserves_taxa_and_branch_lengths() -> None:
    tree, report = reroot_tree_by_midpoint(fixture("example_tree_rootable.nwk"))
    assert sorted(tree.tip_names) == ["A", "B", "C", "D"]
    assert report.strategy == "midpoint"
    assert report.tip_order == ["C", "D", "B", "A"]
    assert report.midpoint_anchor_taxa == ["A", "C"]
    assert report.midpoint_path_length == 1.0
    assert report.midpoint_distance_from_anchor == 0.5
    assert report.midpoint_position_kind == "node"
    assert report.midpoint_anchor_side_taxa == ["A", "B"]
    assert report.midpoint_opposite_side_taxa == ["C", "D"]
    assert report.midpoint_suitable is True
    assert dumps_newick(tree) == "((A:0.2,B:0.2):0.3,(C:0.1,D:0.1):0.4);"
    assert report.summary.branch_lengths_affected != []


def test_reroot_tree_by_midpoint_warns_when_tree_is_not_strictly_bifurcating() -> None:
    tree, report = reroot_tree_by_midpoint(fixture("example_tree_polytomy.nwk"))

    assert sorted(tree.tip_names) == ["A", "B", "C", "D"]
    assert report.midpoint_suitable is False
    assert report.midpoint_anchor_taxa == ["C", "D"]
    assert report.midpoint_path_length == 1.2
    assert report.midpoint_distance_from_anchor == 0.6
    assert report.midpoint_position_kind == "branch"
    assert report.midpoint_anchor_side_taxa == ["A", "B", "C"]
    assert report.midpoint_opposite_side_taxa == ["D"]
    assert report.warnings == [
        "midpoint rooting is exploratory because the input tree is not strictly bifurcating"
    ]


def test_unroot_tree_converts_rooted_binary_tree_into_trifurcation() -> None:
    tree, report = unroot_tree(fixture("example_tree_rootable.nwk"))
    assert sorted(tree.tip_names) == ["A", "B", "C", "D"]
    assert len(tree.root.children) == 3
    assert report.strategy == "unroot"
    assert dumps_newick(tree) == "(A:0.2,B:0.2,(C:0.1,D:0.1):0.7);"
    assert report.warnings == [
        "unrooting merged the removed root-edge length into the retained sibling branch to match ape::unroot"
    ]
    assert report.summary.nodes_changed != []


def test_unroot_tree_returns_already_unrooted_tree_unchanged() -> None:
    tree, report = unroot_tree(fixture("example_tree_unrooted.nwk"))

    assert dumps_newick(tree) == "(A:0.1,B:0.2,C:0.3,D:0.4);"
    assert tree.rooted is False
    assert report.warnings == [
        "input tree already behaves as an unrooted representation; returned unchanged"
    ]
    assert report.summary.nodes_changed == []


def test_unroot_tree_fails_clearly_on_invalid_tree() -> None:
    with pytest.raises(TreeParseError):
        unroot_tree(fixture("example_tree_malformed_unbalanced_parentheses.nwk"))


def test_prune_alignment_to_tree_keeps_exact_tree_taxa() -> None:
    records, report = prune_alignment_to_tree(
        fixture("example_alignment_extra_taxon.fasta"),
        fixture("example_tree.nwk"),
    )
    assert [record.identifier for record in records] == ["A", "B", "C", "D"]
    assert report.original_sequence_count == 5
    assert report.kept_ids == ["A", "B", "C", "D"]
    assert report.removed_ids == ["E"]
    assert [(row.taxon, row.reason) for row in report.removed_ids_with_reasons] == [
        ("E", "absent_from_tree")
    ]


def test_prune_tree_to_alignment_keeps_exact_alignment_taxa() -> None:
    tree, report = prune_tree_to_alignment(
        fixture("example_tree.nwk"),
        fixture("example_alignment_mismatch.fasta"),
    )
    assert tree.tip_names == ["A", "B", "C"]
    assert dumps_newick(tree) == "((A:0.1,B:0.1):0.2,C:0.3);"
    assert report.kept_taxa == ["A", "B", "C"]
    assert report.removed_taxa == ["D"]
    assert report.taxon_column == "identifier"
    assert [(row.taxon, row.reason) for row in report.removed_taxa_with_reasons] == [
        ("D", "absent_from_alignment")
    ]


def test_alignment_inspect_reports_core_diagnostics() -> None:
    report = summarise_fasta(fixture("example_alignment.fasta"))
    assert isinstance(report, AlignmentSummary)
    assert report.sequence_count == 4
    assert report.alignment_length == 8
    assert report.ids == ["A", "B", "C", "D"]
    assert report.missing_data_fraction == 0.0
    assert report.gap_fraction == 0.0
    assert report.constant_site_count == 6
    assert report.variable_site_count == 2
    assert report.parsimony_informative_site_count == 2
    assert report.inferred_alphabet == "dna"
    assert report.nucleotide_composition == {
        "A": 0.3125,
        "C": 0.25,
        "G": 0.25,
        "T": 0.1875,
    }
    assert report.whole_alignment_gc_content == 0.5
    assert [
        (row.identifier, row.missing_fraction)
        for row in report.per_sequence_missingness
    ] == [
        ("A", 0.0),
        ("B", 0.0),
        ("C", 0.0),
        ("D", 0.0),
    ]
    assert [
        (row.identifier, row.gc_fraction) for row in report.per_sequence_gc_content
    ] == [
        ("A", 0.5),
        ("B", 0.375),
        ("C", 0.625),
        ("D", 0.5),
    ]


def test_alignment_detects_sequence_alphabet_types() -> None:
    dna_records = load_fasta_alignment(fixture("example_alignment.fasta"))
    protein_records = load_fasta_alignment(fixture("example_alignment_protein.fasta"))
    assert infer_alignment_alphabet(dna_records) == "dna"
    assert infer_alignment_alphabet(protein_records) == "protein"


def test_alignment_detects_invalid_characters_for_declared_alphabet() -> None:
    invalid = detect_invalid_alignment_characters(
        fixture("example_alignment_invalid_dna.fasta"),
        alphabet="dna",
    )
    assert [(row.identifier, row.position, row.character) for row in invalid] == [
        ("A", 5, "Z")
    ]


def test_alignment_reports_nucleotide_and_amino_acid_composition() -> None:
    dna_records = load_fasta_alignment(fixture("example_alignment.fasta"))
    protein_records = load_fasta_alignment(fixture("example_alignment_protein.fasta"))
    assert compute_nucleotide_composition(dna_records, alphabet="dna") == {
        "A": 0.3125,
        "C": 0.25,
        "G": 0.25,
        "T": 0.1875,
    }
    assert compute_amino_acid_composition(protein_records, alphabet="protein") == {
        "F": 0.083333333333333,
        "I": 0.083333333333333,
        "K": 0.083333333333333,
        "L": 0.125,
        "M": 0.25,
        "R": 0.041666666666667,
        "T": 0.125,
        "V": 0.041666666666667,
        "W": 0.125,
        "Y": 0.041666666666667,
    }


def test_alignment_detects_gc_composition_outliers() -> None:
    outliers = detect_composition_outlier_sequences(
        fixture("example_alignment_gc_outlier.fasta"),
        deviation_threshold=0.2,
    )
    assert [(row.identifier, row.deviation) for row in outliers] == [("C", 1.0)]
    assert outliers[0].robust_z_score is None or outliers[0].robust_z_score > 0


def test_alignment_detects_identical_and_near_duplicate_sequences() -> None:
    duplicates = detect_identical_duplicate_sequences(
        fixture("example_alignment_duplicates.fasta")
    )
    near_duplicates = detect_near_duplicate_sequences(
        fixture("example_alignment_duplicates.fasta"),
        identity_threshold=0.875,
    )
    assert [(group.identifiers, group.sequence) for group in duplicates] == [
        (["A", "B"], "ACTGACTG")
    ]
    assert [
        (
            pair.left_identifier,
            pair.right_identifier,
            pair.identity,
            pair.comparable_sites,
        )
        for pair in near_duplicates
    ] == [
        ("A", "C", 0.875, 8),
        ("A", "D", 0.875, 8),
        ("B", "C", 0.875, 8),
        ("B", "D", 0.875, 8),
        ("C", "D", 0.875, 8),
    ]


def test_alignment_quality_report_collects_composition_duplicates_and_warnings() -> (
    None
):
    report = build_alignment_quality_report(
        fixture("example_alignment_duplicates.fasta")
    )
    assert report.sequence_count == 4
    assert report.alignment_length == 8
    assert report.variable_site_count == 1
    assert report.inferred_alphabet == "dna"
    assert report.invalid_characters == []
    assert report.duplicate_sequence_groups
    assert report.sequence_length_outliers == []
    assert report.near_duplicate_pairs == []
    assert report.warnings == ["alignment contains identical duplicate sequences"]


def test_alignment_classifies_raw_and_ambiguous_fastas() -> None:
    raw = classify_alignment_sequences(fixture("example_sequences_raw.fasta"))
    ambiguous = classify_alignment_sequences(fixture("example_alignment.fasta"))
    assert raw.state == "raw_sequence_fasta"
    assert ambiguous.state == "ambiguous_equal_length_fasta"


def test_alignment_detects_sequence_length_outliers_before_alignment() -> None:
    outliers = detect_sequence_length_outliers(fixture("example_sequences_raw.fasta"))
    assert [(row.identifier, row.raw_length, row.note) for row in outliers] == [
        ("B", 8, "shorter than baseline"),
        ("C", 16, "longer than baseline"),
    ]


def test_alignment_inspect_reports_per_sequence_missingness() -> None:
    report = summarise_fasta(fixture("example_alignment_missingness.fasta"))
    assert report.sequence_count == 3
    assert report.alignment_length == 6
    assert report.constant_site_count == 6
    assert report.variable_site_count == 0
    assert report.parsimony_informative_site_count == 0
    assert report.missing_data_fraction == 4 / 18
    assert [
        (row.identifier, row.missing_fraction)
        for row in report.per_sequence_missingness
    ] == [
        ("A", 2 / 6),
        ("B", 2 / 6),
        ("C", 0.0),
    ]
    assert [
        (row.position, row.missing_fraction) for row in report.per_site_missingness
    ] == [
        (1, 0.0),
        (2, 0.0),
        (3, 0.0),
        (4, 0.0),
        (5, 2 / 3),
        (6, 2 / 3),
    ]
    assert report.all_gap_columns == []
    assert report.all_missing_columns == []


def test_alignment_inspect_reports_site_missingness_and_empty_columns() -> None:
    report = summarise_fasta(fixture("example_alignment_site_missingness.fasta"))
    assert report.alignment_length == 5
    assert [
        (row.position, row.missing_fraction) for row in report.per_site_missingness
    ] == [
        (1, 0.0),
        (2, 0.0),
        (3, 1.0),
        (4, 1.0),
        (5, 0.5),
    ]
    assert report.all_gap_columns == [2]
    assert report.all_missing_columns == [3, 4]


def test_alignment_inspect_separates_gaps_missingness_and_ambiguity() -> None:
    report = summarise_fasta(fixture("example_alignment_ambiguity.fasta"))
    assert report.missing_data_fraction == 5 / 18
    assert report.gap_fraction == 1 / 18
    assert report.ambiguity_fraction == 2 / 18
    assert [
        (row.identifier, row.gap_fraction, row.missing_fraction, row.ambiguity_fraction)
        for row in report.per_sequence_uncertainty
    ] == [
        ("A", 0.0, 2 / 6, 1 / 6),
        ("B", 0.0, 2 / 6, 1 / 6),
        ("C", 1 / 6, 1 / 6, 0.0),
    ]


def test_alignment_windows_report_over_and_under_aligned_regions() -> None:
    windows = summarize_alignment_windows(
        fixture("example_alignment_regions.fasta"),
        window_size=4,
        step_size=4,
    )
    over_aligned = detect_over_aligned_regions(
        fixture("example_alignment_regions.fasta"),
        window_size=4,
        step_size=4,
    )
    under_aligned = detect_under_aligned_regions(
        fixture("example_alignment_regions.fasta"),
        window_size=4,
        step_size=4,
    )
    assert [
        (
            window.start,
            window.end,
            window.gap_fraction,
            window.missing_fraction,
            window.ambiguity_fraction,
        )
        for window in windows
    ] == [
        (1, 4, 0.25, 0.5, 0.25),
        (5, 8, 0.0, 0.0, 0.0),
    ]
    assert [(row.start, row.end, row.kind) for row in over_aligned] == [
        (1, 4, "over_aligned")
    ]
    assert [(row.start, row.end, row.kind) for row in under_aligned] == [
        (5, 8, "under_aligned")
    ]


def test_alignment_readiness_reports_method_specific_decisions() -> None:
    raw = summarize_alignment_readiness(fixture("example_sequences_raw.fasta"))
    coding = summarize_alignment_readiness(fixture("example_alignment_coding.fasta"))
    raw_methods = {row.analysis: row for row in raw.methods}
    coding_methods = {row.analysis: row for row in coding.methods}

    assert raw_methods["distance"].ready is False
    assert raw_methods["maximum_likelihood"].blockers == [
        "input sequences are not yet aligned"
    ]
    assert coding_methods["distance"].ready is True
    assert coding_methods["coding"].ready is False
    assert (
        "one or more sequences contain premature stop codons"
        in coding_methods["coding"].blockers
    )
    assert (
        "one or more sequences contain partial codons after gaps and missing data are removed"
        in coding_methods["coding"].blockers
    )


def test_alignment_filter_profiles_include_coding_safe_profile() -> None:
    profiles = {profile.name: profile for profile in list_alignment_filter_profiles()}
    assert sorted(profiles) == [
        "aggressive",
        "coding-safe",
        "conservative",
        "moderate",
        "phylogenomics-scale",
    ]
    assert get_alignment_filter_profile("coding-safe").preserve_codon_structure is True


def test_coding_alignment_reports_mixed_coding_behaviors() -> None:
    diagnostics = inspect_coding_alignment(fixture("example_alignment_coding.fasta"))
    assert diagnostics.mixed_coding_signals is True
    assert [
        (row.identifier, row.coding_like, row.premature_stop_count)
        for row in diagnostics.coding_behaviors
    ] == [
        ("A", True, 0),
        ("B", True, 0),
        ("C", False, 0),
        ("D", False, 1),
    ]


def test_alignment_cleaning_profile_reports_comparison_and_bias_warnings() -> None:
    records, report = clean_alignment_with_profile(
        fixture("example_alignment_filtering.fasta"),
        profile_name="moderate",
        group_table_path=fixture("example_alignment_groups.tsv"),
        group_columns=["region"],
    )
    assert [record.identifier for record in records] == ["A", "B", "C"]
    assert [record.sequence for record in records] == [
        "ACTGACTG",
        "ACTGACTG",
        "ACTGACTA",
    ]
    assert report.profile.name == "moderate"
    assert report.trim.trimmed_alignment_length == 8
    assert [(row.position, row.reason) for row in report.trim.removed_columns] == [
        (5, "missingness-threshold"),
        (6, "missingness-threshold"),
        (7, "missingness-threshold"),
        (8, "missingness-threshold"),
    ]
    assert [(row.identifier, row.reason) for row in report.trim.removed_sequences] == [
        ("D", "missingness-threshold")
    ]
    assert report.comparison.left_alignment_length == 12
    assert report.comparison.right_alignment_length == 8
    assert any(
        group.column == "region"
        and group.value == "island"
        and group.removed_fraction == 1.0
        for group in report.group_retention
    )
    assert (
        "cleaning removed most taxa from one or more metadata or trait groups"
        in report.warnings
    )


def test_alignment_version_comparison_reports_taxa_length_and_signal_changes() -> None:
    report = compare_alignment_versions(
        fixture("example_alignment_filtering.fasta"),
        fixture("example_alignment_filtering_cleaned_moderate.fasta"),
    )
    assert report.shared_taxa == ["A", "B", "C"]
    assert report.left_only_taxa == ["D"]
    assert report.right_only_taxa == []
    assert report.left_alignment_length == 12
    assert report.right_alignment_length == 8
    assert (
        report.left_parsimony_informative_site_count
        >= report.right_parsimony_informative_site_count
    )


def test_alignment_quality_report_includes_transparent_score_components() -> None:
    report = build_alignment_quality_report(
        fixture("example_alignment_filtering.fasta")
    )
    assert set(report.quality_components) == {
        "composition_outliers",
        "duplicates",
        "gap_burden",
        "informative_density",
        "missingness",
    }
    assert 0.0 <= report.quality_score <= 100.0


def test_alignment_quality_report_exposes_per_sequence_and_per_column_gap_profiles() -> (
    None
):
    report = build_alignment_quality_report(
        fixture("example_alignment_ambiguity.fasta")
    )

    assert report.invariant_site_count == 4
    assert [
        (row.identifier, row.gap_fraction, row.missing_fraction, row.ambiguity_fraction)
        for row in report.per_sequence_uncertainty
    ] == [
        ("A", 0.0, 2 / 6, 1 / 6),
        ("B", 0.0, 2 / 6, 1 / 6),
        ("C", 1 / 6, 1 / 6, 0.0),
    ]
    assert [
        (row.position, row.gap_fraction, row.missing_fraction, row.ambiguity_fraction)
        for row in report.per_site_uncertainty
    ] == [
        (1, 0.0, 0.0, 0.0),
        (2, 0.0, 0.0, 0.0),
        (3, 0.0, 0.0, 0.0),
        (4, 0.0, 0.0, 0.0),
        (5, 1 / 3, 2 / 3, 2 / 3),
        (6, 0.0, 1.0, 0.0),
    ]


def test_alignment_quality_report_flags_missingness_concentration_and_suspicious_windows() -> (
    None
):
    missingness_report = build_alignment_quality_report(
        fixture("example_alignment_missingness.fasta")
    )
    suspicious_window_report = build_alignment_quality_report(
        fixture("example_alignment_filtering.fasta")
    )

    assert missingness_report.missing_data_concentration.concentrated_column_count == 2
    assert missingness_report.missing_data_concentration.longest_concentrated_run == 2
    assert (
        missingness_report.missing_data_concentration.longest_concentrated_run_start
        == 5
    )
    assert (
        missingness_report.missing_data_concentration.longest_concentrated_run_end == 6
    )
    assert missingness_report.suspicious_alignment is True
    assert (
        "alignment concentrates missing data into adjacent columns"
        in missingness_report.suspicious_reasons
    )
    assert (
        "alignment has low information content for defensible inference"
        in missingness_report.suspicious_reasons
    )

    assert suspicious_window_report.suspicious_alignment is True
    assert (
        "alignment contains suspiciously over-aligned windows"
        in suspicious_window_report.suspicious_reasons
    )


def test_alignment_low_information_detection_blocks_sparse_inference_inputs() -> None:
    report = assess_alignment_low_information(
        fixture("example_alignment_missingness.fasta")
    )
    readiness = summarize_alignment_readiness(
        fixture("example_alignment_missingness.fasta")
    )
    methods = {row.analysis: row for row in readiness.methods}

    assert report.low_information is True
    assert report.parsimony_informative_site_count == 0
    assert any("parsimony-informative sites" in reason for reason in report.reasons)
    assert (
        "alignment has too few parsimony-informative sites for defensible inference"
        in methods["distance"].blockers
    )
    assert (
        "alignment has too few parsimony-informative sites for defensible inference"
        in methods["maximum_likelihood"].blockers
    )
    assert (
        "alignment has too few parsimony-informative sites for defensible inference"
        in methods["bayesian"].blockers
    )


def test_duplicate_sequence_policy_report_recommends_collapse_and_review() -> None:
    report = build_duplicate_sequence_policy_report(
        fixture("example_alignment_duplicates.fasta"),
        near_duplicate_threshold=0.875,
    )

    assert [
        (group.identifiers, group.sequence) for group in report.exact_duplicate_groups
    ] == [(["A", "B"], "ACTGACTG")]
    assert any(
        action.action == "collapse_exact_duplicates"
        and action.affected_identifiers == ["A", "B"]
        for action in report.policy_actions
    )
    assert any(
        action.action == "review_near_duplicates"
        and action.affected_identifiers == ["A", "C"]
        for action in report.policy_actions
    )
    assert any("deduplicated" in warning for warning in report.warnings)


def test_ambiguous_alignment_column_report_lists_uncertainty_heavy_sites() -> None:
    report = build_ambiguous_alignment_column_report(
        fixture("example_alignment_site_missingness.fasta"),
        threshold=0.5,
    )

    assert [
        (row.position, row.gap_fraction, row.missing_fraction) for row in report.rows
    ] == [
        (2, 1.0, 0.0),
        (3, 0.0, 1.0),
        (4, 0.0, 1.0),
        (5, 0.0, 0.5),
    ]
    assert report.warnings == [
        "alignment contains ambiguity-heavy columns that may be unsuitable for inference without masking"
    ]


def test_sequence_quality_ranking_orders_sequences_by_burden() -> None:
    report = build_sequence_quality_ranking(
        fixture("example_alignment_ambiguity.fasta")
    )

    assert [(row.identifier, row.rank, row.score) for row in report.rows] == [
        ("A", 1, 78.333),
        ("B", 2, 78.333),
        ("C", 3, 84.167),
    ]
    assert report.warnings == [
        "lower-ranked sequences should be reviewed before publication or inference"
    ]


def test_alignment_forensic_report_integrates_alignment_risks() -> None:
    report = build_alignment_forensic_report(fixture("example_alignment_coding.fasta"))
    assert report.safe_for_distance_analysis is True
    assert report.safe_for_coding_analysis is False
    assert (
        "alignment mixes coding-like and noncoding-like sequence behavior"
        in report.warnings
    )
    assert report.low_information.low_information is False
    assert any(
        "near-duplicate sequences" in warning
        for warning in report.duplicate_policy.warnings
    )
    assert report.ambiguous_columns.rows == []
    assert report.sequence_ranking.rows


def test_dataset_audit_integrates_alignment_and_time_tree_surfaces() -> None:
    report = audit_dataset_inputs(
        fixture("example_tree_named_clades.nwk"),
        fixture("example_metadata.tsv"),
        fixture("example_traits_validate.tsv"),
        alignment_path=fixture("example_alignment.fasta"),
        tip_dates_path=fixture("example_tip_dates.tsv"),
        calibration_path=fixture("example_calibrations.tsv"),
    )
    assert report.readiness_decision == "ready_with_warnings"
    assert "comparative" in report.allowed_analyses
    assert "distance" in report.allowed_analyses
    assert "time_tree" in report.allowed_analyses
    assert report.alignment_forensic is not None
    assert report.tip_dates is not None
    assert report.calibrations is not None
    assert "alignment" in report.warning_categories
    assert any(
        row.analysis == "distance" and row.decision == "risky"
        for row in report.analysis_decisions
    )
    assert any(
        level.level == "publication_ready" and level.decision == "risky"
        for level in report.readiness_levels
    )


def test_dataset_audit_blocks_invalid_time_tree_inputs() -> None:
    report = audit_dataset_inputs(
        fixture("example_tree_named_clades.nwk"),
        fixture("example_metadata.tsv"),
        fixture("example_traits_validate.tsv"),
        alignment_path=fixture("example_alignment.fasta"),
        tip_dates_path=fixture("example_tip_dates_invalid.tsv"),
        calibration_path=fixture("example_calibrations_invalid.tsv"),
    )
    assert report.readiness_decision == "blocked"
    assert "time_tree" in report.blocked_analyses
    assert set(report.blocker_categories) >= {"calibration", "tip_dates"}


def test_dataset_crosswalk_and_completeness_matrix_report_surface_presence() -> None:
    crosswalk = build_dataset_crosswalk(
        fixture("example_tree_named_clades.nwk"),
        fixture("example_metadata.tsv"),
        fixture("example_traits_validate.tsv"),
        alignment_path=fixture("example_alignment.fasta"),
        tip_dates_path=fixture("example_tip_dates.tsv"),
        calibration_path=fixture("example_calibrations.tsv"),
    )
    matrix = build_dataset_completeness_matrix(
        fixture("example_tree_named_clades.nwk"),
        fixture("example_metadata.tsv"),
        fixture("example_traits_validate.tsv"),
        alignment_path=fixture("example_alignment.fasta"),
        tip_dates_path=fixture("example_tip_dates.tsv"),
        calibration_path=fixture("example_calibrations.tsv"),
    )
    crosswalk_by_taxon = {row.taxon: row for row in crosswalk.rows}
    matrix_by_taxon = {row.taxon: row for row in matrix.rows}

    assert crosswalk_by_taxon["A"].tree_tip == "A"
    assert crosswalk_by_taxon["A"].alignment_id == "A"
    assert crosswalk_by_taxon["A"].metadata_id == "A"
    assert crosswalk_by_taxon["A"].trait_id == "A"
    assert crosswalk_by_taxon["A"].tip_date_id == "A"
    assert crosswalk_by_taxon["A"].calibration_targets == ["cal-mammals"]
    assert crosswalk_by_taxon["C"].calibration_targets == ["cal-birds"]
    assert matrix_by_taxon["A"].in_geography is True
    assert matrix.surface_counts["calibrations"] == 4


def test_dataset_mismatch_report_lists_taxa_missing_requested_surfaces() -> None:
    report = build_dataset_mismatch_report(
        fixture("example_tree.nwk"),
        fixture("example_alignment_groups.tsv"),
        fixture("example_traits.tsv"),
        alignment_path=fixture("example_alignment_filtering_cleaned_moderate.fasta"),
    )

    rows = {row.taxon: row for row in report.rows}
    assert "D" in rows
    assert "alignment" in rows["D"].missing_surfaces
    assert "tree" in rows["D"].present_surfaces
    assert report.mismatch_counts["alignment"] >= 1


def test_dataset_ordering_audit_detects_reordered_metadata_rows() -> None:
    report = audit_dataset_taxon_ordering(
        fixture("example_tree.nwk"),
        fixture("example_metadata_reordered.tsv"),
        fixture("example_traits_validate.tsv"),
    )
    assert report.consistent is False
    assert report.drifted_surfaces == ["metadata"]
    assert any(
        row.surface == "metadata"
        and row.taxon == "C"
        and row.expected_index == 3
        and row.observed_index == 1
        for row in report.conflicts
    )


def test_dataset_audit_reports_group_imbalance_and_exclusion_rows() -> None:
    report = audit_dataset_inputs(
        fixture("example_tree.nwk"),
        fixture("example_alignment_groups.tsv"),
        fixture("example_traits.tsv"),
        alignment_path=fixture("example_alignment_filtering_cleaned_moderate.fasta"),
    )
    assert any(
        row.surface == "metadata"
        and row.group_column == "region"
        and row.group == "island"
        and row.removed_fraction == 1.0
        for row in report.group_imbalance_warnings
    )
    exclusion_by_taxon = {row.taxon: row for row in report.exclusion_table.rows}
    assert exclusion_by_taxon["D"].causes == [
        "absent_from_alignment",
        "absent_from_traits",
    ]
    assert exclusion_by_taxon["D"].first_failed_surface == "alignment"
    assert report.mismatch_report.rows
    assert report.risk_score.total_score > 0.0
    assert any(
        component.component == "alignment" for component in report.risk_score.components
    )
    assert report.minimal_fix_plan.recommendations
    assert any(
        item.section == "dataset_risk" for item in report.reviewer_checklist.items
    )


def test_alignment_inspect_rejects_unequal_lengths() -> None:
    try:
        summarise_fasta(fixture("example_alignment_invalid_lengths.fasta"))
    except InvalidAlignmentError as error:
        assert error.code == "invalid_alignment_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected InvalidAlignmentError")


def test_alignment_link_reports_exact_mismatch() -> None:
    report = link_alignment_to_tree(
        fixture("example_tree.nwk"), fixture("example_alignment.fasta")
    )
    assert report.tree_taxa == 4
    assert report.alignment_ids == 4
    assert report.linked_taxa == 4
    assert report.missing_from_alignment == []
    assert report.extra_alignment_ids == []


def test_alignment_link_strict_mode_rejects_mismatch() -> None:
    try:
        link_alignment_to_tree(
            fixture("example_tree.nwk"),
            fixture("example_alignment_mismatch.fasta"),
            strict=True,
        )
    except AlignmentTaxonMismatchError as error:
        assert error.code == "alignment_taxon_mismatch_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected AlignmentTaxonMismatchError")


def test_write_fasta_alignment_preserves_record_order_and_sequences(
    tmp_path: Path,
) -> None:
    records = load_fasta_alignment(fixture("example_alignment.fasta"))
    output = tmp_path / "alignment.fasta"
    write_fasta_alignment(output, records)
    assert output.read_text(encoding="utf-8") == (
        ">A\nACTGACTG\n>B\nACTGACTA\n>C\nACTGACGG\n>D\nACTGACGA\n"
    )
    assert load_fasta_alignment(output) == records


def test_write_dna_bin_alignment_fasta_preserves_normalized_nucleotide_states(
    tmp_path: Path,
) -> None:
    alignment = load_dna_bin_alignment(fixture("example_alignment_ambiguity.fasta"))
    output = tmp_path / "dnabin-alignment.fasta"

    write_dna_bin_alignment_fasta(output, alignment)

    assert output.read_text(encoding="utf-8") == (
        ">A\nacgtn?\n>B\nacgtr?\n>C\nacgt-?\n"
    )
    assert load_dna_bin_alignment(output).records == alignment.records


def test_alignment_detects_sequences_with_excessive_missing_data() -> None:
    rows = detect_sequences_with_excessive_missing_data(
        fixture("example_alignment_missingness.fasta"),
        threshold=0.3,
    )
    assert [(row.identifier, row.missing_fraction) for row in rows] == [
        ("A", 2 / 6),
        ("B", 2 / 6),
    ]


def test_alignment_detects_sites_with_excessive_missing_data() -> None:
    rows = detect_sites_with_excessive_missing_data(
        fixture("example_alignment_site_missingness.fasta"),
        threshold=0.4,
    )
    assert [(row.position, row.missing_fraction) for row in rows] == [
        (3, 1.0),
        (4, 1.0),
        (5, 0.5),
    ]


def test_alignment_removes_all_gap_columns() -> None:
    records, report = remove_all_gap_columns(
        fixture("example_alignment_site_missingness.fasta")
    )
    assert [record.sequence for record in records] == ["AN?T", "CN?N", "GN?A", "TN?N"]
    assert report.original_alignment_length == 5
    assert report.trimmed_alignment_length == 4
    assert [(column.position, column.reason) for column in report.removed_columns] == [
        (2, "all-gap")
    ]


def test_alignment_removes_all_missing_columns() -> None:
    records, report = remove_all_missing_columns(
        fixture("example_alignment_site_missingness.fasta")
    )
    assert [record.sequence for record in records] == ["A-T", "C-N", "G-A", "T-N"]
    assert report.original_alignment_length == 5
    assert report.trimmed_alignment_length == 3
    assert [(column.position, column.reason) for column in report.removed_columns] == [
        (3, "all-missing"),
        (4, "all-missing"),
    ]


def test_alignment_removes_sequences_above_missingness_threshold() -> None:
    records, report = remove_sequences_above_missingness_threshold(
        fixture("example_alignment_missingness.fasta"),
        threshold=0.3,
    )
    assert [record.identifier for record in records] == ["C"]
    assert report.original_sequence_count == 3
    assert report.trimmed_sequence_count == 1
    assert [
        (row.identifier, row.missing_fraction, row.reason)
        for row in report.removed_sequences
    ] == [
        ("A", 2 / 6, "missingness-threshold"),
        ("B", 2 / 6, "missingness-threshold"),
    ]


def test_alignment_trims_columns_above_missingness_threshold() -> None:
    records, report = trim_columns_above_missingness_threshold(
        fixture("example_alignment_site_missingness.fasta"),
        threshold=0.4,
    )
    assert [record.sequence for record in records] == ["A-", "C-", "G-", "T-"]
    assert report.original_alignment_length == 5
    assert report.trimmed_alignment_length == 2
    assert [(column.position, column.reason) for column in report.removed_columns] == [
        (3, "missingness-threshold"),
        (4, "missingness-threshold"),
        (5, "missingness-threshold"),
    ]


def test_alignment_trimming_report_combines_sequence_and_column_transforms() -> None:
    records, report = trim_alignment(
        fixture("example_alignment_trim.fasta"),
        site_missingness_threshold=0.4,
        sequence_missingness_threshold=0.3,
    )
    assert [(record.identifier, record.sequence) for record in records] == [
        ("B", "ACG")
    ]
    assert report.original_sequence_count == 3
    assert report.trimmed_sequence_count == 1
    assert report.original_alignment_length == 6
    assert report.trimmed_alignment_length == 3
    assert [(column.position, column.reason) for column in report.removed_columns] == [
        (3, "all-gap"),
        (4, "all-missing"),
        (5, "missingness-threshold"),
    ]
    assert [(row.identifier, row.reason) for row in report.removed_sequences] == [
        ("A", "missingness-threshold"),
        ("C", "missingness-threshold"),
    ]


def test_alignment_identity_matrix_reports_pairs_and_comparable_sites() -> None:
    report = compute_pairwise_sequence_identity_matrix(
        fixture("example_alignment_duplicates.fasta")
    )
    assert report.identifiers == ["A", "B", "C", "D"]
    assert [
        (
            pair.left_identifier,
            pair.right_identifier,
            pair.identity,
            pair.comparable_sites,
        )
        for pair in report.pairs
    ] == [
        ("A", "A", 1.0, 8),
        ("A", "B", 1.0, 8),
        ("A", "C", 0.875, 8),
        ("A", "D", 0.875, 8),
        ("B", "B", 1.0, 8),
        ("B", "C", 0.875, 8),
        ("B", "D", 0.875, 8),
        ("C", "C", 1.0, 8),
        ("C", "D", 0.875, 8),
        ("D", "D", 1.0, 8),
    ]


def test_coding_alignment_reports_frameshift_like_sequences_and_stop_codons() -> None:
    diagnostics = inspect_coding_alignment(fixture("example_alignment_coding.fasta"))
    assert diagnostics.genetic_code_id == 1
    assert diagnostics.genetic_code_name == "Standard"
    assert diagnostics.sequence_count == 4
    assert diagnostics.alignment_length_multiple_of_three is True
    assert [
        (row.identifier, row.comparable_length, row.remainder)
        for row in diagnostics.frameshift_like_sequences
    ] == [("C", 8, 2)]
    assert [
        (row.identifier, row.comparable_length, row.trailing_bases)
        for row in diagnostics.partial_codon_sequences
    ] == [("C", 8, 2)]
    assert [
        (row.identifier, row.codon_index, row.nucleotide_start, row.codon, row.terminal)
        for row in diagnostics.stop_codons
    ] == [
        ("A", 3, 7, "TAA", True),
        ("D", 2, 4, "TAG", False),
    ]
    assert diagnostics.invalid_codons == []


def test_translate_coding_alignment_emits_amino_acid_records() -> None:
    records, report = translate_coding_alignment(
        fixture("example_alignment_coding.fasta")
    )
    assert [(record.identifier, record.sequence) for record in records] == [
        ("A", "ME*"),
        ("B", "MEW"),
        ("C", "MXW"),
        ("D", "M*W"),
    ]
    assert report.genetic_code_id == 1
    assert report.genetic_code_name == "Standard"
    assert report.translated_sequence_count == 4
    assert report.source_alignment_length == 9
    assert report.translated_alignment_length == 3
    assert report.invalid_codon_count == 1
    assert report.stop_codon_count == 2
    assert report.internal_stop_sequence_count == 1
    assert report.terminal_stop_sequence_count == 1
    assert report.trailing_partial_codon_sequence_count == 0
    assert report.dropped_trailing_nucleotide_count == 0
    assert report.warnings == []
    assert [
        (row.identifier, row.codon_index, row.codon, row.amino_acid, row.translation_status)
        for row in report.codon_observations[:5]
    ] == [
        ("A", 1, "ATG", "M", "translated"),
        ("A", 2, "GAA", "E", "translated"),
        ("A", 3, "TAA", "*", "terminal-stop-codon"),
        ("B", 1, "ATG", "M", "translated"),
        ("B", 2, "GAA", "E", "translated"),
    ]


def test_prepare_coding_sequences_for_alignment_excludes_frame_and_stop_failures(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "coding-raw.fasta"
    write_fasta_alignment(
        input_path,
        [
            AlignmentRecord(identifier="good", sequence="ATGGAATGG"),
            AlignmentRecord(identifier="terminal_stop", sequence="ATGGAATAA"),
            AlignmentRecord(identifier="frameshift", sequence="ATGGAATG"),
            AlignmentRecord(identifier="internal_stop", sequence="ATGTAGTGG"),
        ],
    )

    records, report = prepare_coding_sequences_for_alignment(input_path)

    assert [(record.identifier, record.sequence) for record in records] == [
        ("good", "ATGGAATGG"),
        ("terminal_stop", "ATGGAATAA"),
    ]
    assert report.sequence_type == "dna"
    assert report.genetic_code_id == 1
    assert report.genetic_code_name == "Standard"
    assert report.input_sequence_count == 4
    assert report.accepted_sequence_count == 2
    assert report.accepted_identifiers == ["good", "terminal_stop"]
    assert report.invalid_codon_sequence_count == 0
    assert report.terminal_stop_sequence_count == 1
    assert [
        (row.identifier, row.reason, row.invalid_codon_count, row.trailing_bases)
        for row in report.excluded_sequences
    ] == [
        ("frameshift", "frame-error", 0, 2),
        ("internal_stop", "internal-stop-codon", 0, 0),
    ]
    assert (
        "one or more coding sequences were excluded before codon-aware alignment"
        in report.warnings
    )
    assert (
        "terminal stop codons were retained in accepted coding sequences"
        in report.warnings
    )


def test_prepare_coding_sequences_for_alignment_preserves_rna_residues(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "coding-rna.fasta"
    write_fasta_alignment(
        input_path,
        [
            AlignmentRecord(identifier="rna_a", sequence="AUGGAAUGG"),
            AlignmentRecord(identifier="rna_b", sequence="AUGGAAUAA"),
        ],
    )

    records, report = prepare_coding_sequences_for_alignment(input_path)

    assert [(record.identifier, record.sequence) for record in records] == [
        ("rna_a", "AUGGAAUGG"),
        ("rna_b", "AUGGAAUAA"),
    ]
    assert report.sequence_type == "rna"


def test_prepare_coding_sequences_for_alignment_excludes_invalid_codons(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "coding-invalid-codon.fasta"
    write_fasta_alignment(
        input_path,
        [
            AlignmentRecord(identifier="good", sequence="ATGGAATGG"),
            AlignmentRecord(identifier="ambiguous_one", sequence="ATGNNNTGG"),
            AlignmentRecord(identifier="ambiguous_two", sequence="ATGRYATGG"),
        ],
    )

    records, report = prepare_coding_sequences_for_alignment(input_path)

    assert [(record.identifier, record.sequence) for record in records] == [
        ("good", "ATGGAATGG"),
    ]
    assert report.invalid_codon_sequence_count == 2
    assert [
        (row.identifier, row.reason, row.invalid_codon_count)
        for row in report.excluded_sequences
    ] == [
        ("ambiguous_one", "invalid-codon", 1),
        ("ambiguous_two", "invalid-codon", 1),
    ]
    assert any("ambiguous or invalid codons" in warning for warning in report.warnings)


def test_prepare_coding_sequences_for_alignment_honors_configurable_genetic_code(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "coding-mito.fasta"
    write_fasta_alignment(
        input_path,
        [
            AlignmentRecord(identifier="standard_only", sequence="ATGTGAGGG"),
            AlignmentRecord(identifier="shared_good", sequence="ATGGAATGG"),
        ],
    )

    standard_records, standard_report = prepare_coding_sequences_for_alignment(
        input_path,
        genetic_code="1",
    )
    mitochondrial_records, mitochondrial_report = (
        prepare_coding_sequences_for_alignment(
            input_path,
            genetic_code="2",
        )
    )

    assert [record.identifier for record in standard_records] == ["shared_good"]
    assert [
        (row.identifier, row.reason, row.premature_stop_count)
        for row in standard_report.excluded_sequences
    ] == [("standard_only", "internal-stop-codon", 1)]
    assert [record.identifier for record in mitochondrial_records] == [
        "shared_good",
        "standard_only",
    ]
    assert mitochondrial_report.excluded_sequences == []
    assert mitochondrial_report.genetic_code_id == 2
    assert mitochondrial_report.genetic_code_name == "Vertebrate Mitochondrial"


def test_translate_coding_alignment_honors_configurable_genetic_code_and_invalid_codons(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "coding-translate.fasta"
    write_fasta_alignment(
        input_path,
        [
            AlignmentRecord(identifier="mito_triplet", sequence="ATGTGAGGG"),
            AlignmentRecord(identifier="ambiguous_triplet", sequence="ATGNNNGGG"),
        ],
    )

    standard_records, standard_report = translate_coding_alignment(
        input_path,
        genetic_code="1",
    )
    mitochondrial_records, mitochondrial_report = translate_coding_alignment(
        input_path,
        genetic_code="2",
    )

    assert [(record.identifier, record.sequence) for record in standard_records] == [
        ("mito_triplet", "M*G"),
        ("ambiguous_triplet", "MXG"),
    ]
    assert [
        (record.identifier, record.sequence) for record in mitochondrial_records
    ] == [
        ("mito_triplet", "MWG"),
        ("ambiguous_triplet", "MXG"),
    ]
    assert standard_report.stop_codon_count == 1
    assert mitochondrial_report.stop_codon_count == 0
    assert mitochondrial_report.invalid_codon_count == 1


def test_translate_coding_alignment_truncates_trailing_partial_codon_like_ape() -> None:
    records, report = translate_coding_alignment(
        fixture("example_alignment_coding_frame_error.fasta")
    )

    assert [(record.identifier, record.sequence) for record in records] == [
        ("frame_error", "ME"),
    ]
    assert report.translated_alignment_length == 2
    assert report.dropped_trailing_nucleotide_count == 2
    assert report.trailing_partial_codon_sequence_count == 1
    assert report.warnings == [
        "sequence length not a multiple of 3: 2 nucleotides dropped"
    ]


def test_cli_alignment_trim_writes_trimmed_fasta_and_report(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "trimmed.fasta"
    exit_code = main(
        [
            "alignment",
            "trim",
            str(fixture("example_alignment_trim.fasta")),
            "--out",
            str(output_path),
            "--site-missingness-threshold",
            "0.4",
            "--sequence-missingness-threshold",
            "0.3",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == ">B\nACG\n"
    assert payload["metrics"]["removed_column_count"] == 3
    assert payload["metrics"]["removed_sequence_count"] == 2


def test_cli_alignment_identity_matrix_writes_tsv(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "identity.tsv"
    exit_code = main(
        [
            "alignment",
            "identity-matrix",
            str(fixture("example_alignment_duplicates.fasta")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8").splitlines()[:4] == [
        "left_identifier\tright_identifier\tidentity\tcomparable_sites",
        "A\tA\t1\t8",
        "A\tB\t1\t8",
        "A\tC\t0.875\t8",
    ]
    assert payload["metrics"]["pair_count"] == 10


def test_cli_alignment_distance_matrix_writes_tsv(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "distance.tsv"
    exit_code = main(
        [
            "alignment",
            "distance-matrix",
            str(fixture("example_alignment_distance_gaps.fasta")),
            "--model",
            "jukes-cantor",
            "--gap-handling",
            "complete-deletion",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8").splitlines()[:4] == [
        "left_identifier\tright_identifier\tdistance\tcomparable_sites",
        "A\tA\t0\t4",
        "A\tB\t0\t4",
        "A\tC\t0.304098831081123\t4",
    ]
    assert payload["metrics"]["model"] == "jukes-cantor"
    assert payload["metrics"]["gap_handling"] == "complete-deletion"


def test_compute_pairwise_genetic_distance_matrix_reports_p_distance() -> None:
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance.fasta")
    )
    assert report.model == "p-distance"
    assert report.gap_handling == "pairwise-deletion"
    assert report.identifiers == ["A", "B", "C", "D"]
    assert [
        (
            pair.left_identifier,
            pair.right_identifier,
            pair.distance,
            pair.comparable_sites,
        )
        for pair in report.pairs
    ] == [
        ("A", "A", 0.0, 8),
        ("A", "B", 0.125, 8),
        ("A", "C", 0.5, 8),
        ("A", "D", 0.625, 8),
        ("B", "B", 0.0, 8),
        ("B", "C", 0.625, 8),
        ("B", "D", 0.5, 8),
        ("C", "C", 0.0, 8),
        ("C", "D", 0.125, 8),
        ("D", "D", 0.0, 8),
    ]


def test_compute_pairwise_genetic_distance_matrix_uses_pairwise_deletion() -> None:
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance_gaps.fasta")
    )
    assert [
        (
            pair.left_identifier,
            pair.right_identifier,
            pair.distance,
            pair.comparable_sites,
        )
        for pair in report.pairs
    ] == [
        ("A", "A", 0.0, 6),
        ("A", "B", 0.166666666666667, 6),
        ("A", "C", 0.5, 6),
        ("A", "D", 0.0, 4),
        ("B", "B", 0.0, 6),
        ("B", "C", 0.5, 6),
        ("B", "D", 0.0, 4),
        ("C", "C", 0.0, 6),
        ("C", "D", 0.25, 4),
        ("D", "D", 0.0, 4),
    ]


def test_compute_pairwise_genetic_distance_matrix_supports_complete_deletion() -> None:
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance_gaps.fasta"),
        gap_handling="complete-deletion",
    )
    assert report.gap_handling == "complete-deletion"
    assert [
        (
            pair.left_identifier,
            pair.right_identifier,
            pair.distance,
            pair.comparable_sites,
        )
        for pair in report.pairs
    ] == [
        ("A", "A", 0.0, 4),
        ("A", "B", 0.0, 4),
        ("A", "C", 0.25, 4),
        ("A", "D", 0.0, 4),
        ("B", "B", 0.0, 4),
        ("B", "C", 0.25, 4),
        ("B", "D", 0.0, 4),
        ("C", "C", 0.0, 4),
        ("C", "D", 0.25, 4),
        ("D", "D", 0.0, 4),
    ]


def test_compute_pairwise_genetic_distance_matrix_supports_jukes_cantor() -> None:
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance.fasta"),
        model="jukes-cantor",
    )
    assert report.model == "jukes-cantor"
    assert [
        (
            pair.left_identifier,
            pair.right_identifier,
            pair.distance,
            pair.comparable_sites,
        )
        for pair in report.pairs
    ] == [
        ("A", "A", 0.0, 8),
        ("A", "B", 0.136741167595466, 8),
        ("A", "C", 0.823959216501082, 8),
        ("A", "D", 1.343819601921041, 8),
        ("B", "B", 0.0, 8),
        ("B", "C", 1.343819601921041, 8),
        ("B", "D", 0.823959216501082, 8),
        ("C", "C", 0.0, 8),
        ("C", "D", 0.136741167595466, 8),
        ("D", "D", 0.0, 8),
    ]


def test_compute_pairwise_genetic_distance_matrix_marks_saturated_jukes_cantor_pairs_undefined() -> (
    None
):
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance_saturated.fasta"),
        model="jukes-cantor",
    )
    assert [
        (
            pair.left_identifier,
            pair.right_identifier,
            pair.distance,
            pair.comparable_sites,
        )
        for pair in report.pairs
    ] == [
        ("A", "A", 0.0, 4),
        ("A", "B", None, 4),
        ("A", "C", 0.304098831081123, 4),
        ("B", "B", 0.0, 4),
        ("B", "C", None, 4),
        ("C", "C", 0.0, 4),
    ]


def test_load_imported_distance_matrix_reads_exported_long_form() -> None:
    entries = load_imported_distance_matrix(fixture("example_distance_matrix.tsv"))
    assert [
        (
            entry.left_identifier,
            entry.right_identifier,
            entry.distance,
            entry.comparable_sites,
        )
        for entry in entries[:4]
    ] == [
        ("A", "A", 0.0, 8),
        ("A", "B", 0.125, 8),
        ("A", "C", 0.5, 8),
        ("B", "A", 0.125, 8),
    ]


def test_validate_imported_distance_matrix_reports_clean_matrix() -> None:
    report = validate_imported_distance_matrix(fixture("example_distance_matrix.tsv"))
    assert report.complete is True
    assert report.zero_diagonal is True
    assert report.symmetric is True
    assert report.nonnegative is True
    assert report.missing_pairs == []
    assert report.nonmetric_observations == []
    assert report.warnings == []


def test_validate_imported_distance_matrix_detects_nonmetric_violations() -> None:
    report = validate_imported_distance_matrix(
        fixture("example_distance_matrix_nonmetric.tsv")
    )
    assert [
        (
            row.left_identifier,
            row.middle_identifier,
            row.right_identifier,
            row.direct_distance,
            row.indirect_distance,
        )
        for row in report.nonmetric_observations
    ] == [
        ("A", "B", "C", 5.0, 2.0),
    ]
    assert report.warnings == [
        "distance matrix violates triangle inequality for one or more taxon triples"
    ]


def test_validate_imported_distance_matrix_detects_asymmetry() -> None:
    report = validate_imported_distance_matrix(
        fixture("example_distance_matrix_asymmetric.tsv")
    )
    assert report.complete is True
    assert report.symmetric is False
    assert [
        (
            row.left_identifier,
            row.right_identifier,
            row.left_to_right_distance,
            row.right_to_left_distance,
        )
        for row in report.asymmetric_pairs
    ] == [
        ("A", "B", 1.0, 2.0),
    ]


def test_build_tree_from_imported_distance_matrix_constructs_neighbor_joining_tree() -> (
    None
):
    tree, report = build_tree_from_imported_distance_matrix(
        fixture("example_distance_matrix.tsv"),
        method="neighbor-joining",
    )
    assert dumps_newick(tree) == "(A:0,B:0.125,C:0.5)Inner1;"
    assert report.method == "neighbor-joining"
    assert report.taxon_count == 3


def test_build_tree_from_imported_distance_matrix_rejects_asymmetric_input() -> None:
    try:
        build_tree_from_imported_distance_matrix(
            fixture("example_distance_matrix_asymmetric.tsv"),
            method="upgma",
        )
    except InvalidDistanceMatrixError as error:
        assert error.code == "invalid_distance_matrix_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected InvalidDistanceMatrixError")


def test_distance_method_limitations_explain_approximate_methods() -> None:
    limitations = distance_method_limitations()
    assert len(limitations) == 5
    assert limitations[0].startswith("distance methods collapse")
    assert "BIONJ is explicitly excluded" in limitations[-1]


def test_render_distance_report_embeds_limitations_and_validation(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "distance-report.html"
    result = render_distance_report(
        out_path=output_path,
        matrix_path=fixture("example_distance_matrix.tsv"),
    )
    html = output_path.read_text(encoding="utf-8")
    assert result.source_kind == "imported-distance-matrix"
    assert "distance-method-limitations" in html
    assert "imported-distance-matrix-quality" in html
    assert "neighbor-joining-tree" in html


def test_render_tree_uncertainty_report_embeds_consensus_and_instability_sections(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "tree-uncertainty-report.html"
    result = render_tree_uncertainty_report(
        tree_set_path=fixture("example_tree_set_left.nwk"),
        out_path=output_path,
    )
    html = output_path.read_text(encoding="utf-8")
    assert result.report_kind == "tree-uncertainty"
    assert result.tree_count == 3
    assert result.rooted_topology_count == 2
    assert result.machine_manifest["report_kind"] == "tree-uncertainty"
    assert "consensus-tree" in html
    assert "rf-distance-distribution" in html
    assert "topology-multimodality" in html
    assert "clade-credibility-conflicts" in html
    assert "unstable-clades" in html


def test_build_distance_tree_constructs_neighbor_joining_tree() -> None:
    tree, report = build_distance_tree(
        fixture("example_alignment_distance.fasta"),
        method="neighbor-joining",
    )
    assert (
        dumps_newick(tree)
        == "((A:0.0625,B:0.0625)Inner1:0.4375,C:0.0625,D:0.0625)Inner2;"
    )
    assert report.method == "neighbor-joining"
    assert report.taxon_count == 4


def test_build_distance_tree_from_genetic_distance_matrix_matches_path_surface() -> None:
    matrix = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance.fasta")
    )
    direct_tree, direct_report = build_distance_tree_from_genetic_distance_matrix(
        matrix,
        method="neighbor-joining",
    )
    path_tree, path_report = build_distance_tree(
        fixture("example_alignment_distance.fasta"),
        method="neighbor-joining",
    )

    assert dumps_newick(direct_tree) == dumps_newick(path_tree)
    assert direct_report.alignment_path == path_report.alignment_path
    assert direct_report.method == path_report.method
    assert direct_report.method_policy == path_report.method_policy
    assert direct_report.assumptions == path_report.assumptions


def test_build_distance_tree_constructs_upgma_tree() -> None:
    tree, report = build_distance_tree(
        fixture("example_alignment_distance.fasta"),
        method="upgma",
    )
    assert (
        dumps_newick(tree)
        == "((A:0.0625,B:0.0625)Inner2:0.21875,(C:0.0625,D:0.0625)Inner1:0.21875)Inner3;"
    )
    assert report.method == "upgma"
    assert report.taxon_count == 4


def test_compare_distance_tree_topologies_reports_rooting_difference() -> None:
    report = compare_distance_tree_topologies(
        fixture("example_alignment_distance.fasta")
    )
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.topology_equal is False
    assert report.same_unrooted_topology is True
    assert report.same_taxa_different_rooting is True
    assert report.robinson_foulds_distance == 1


def test_build_distance_tree_rejects_undefined_corrected_distances() -> None:
    try:
        build_distance_tree(
            fixture("example_alignment_distance_saturated.fasta"),
            method="neighbor-joining",
            model="jukes-cantor",
        )
    except InvalidAlignmentError as error:
        assert "undefined entries" in error.message
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected InvalidAlignmentError")


def test_cli_alignment_build_tree_writes_newick(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "distance-tree.nwk"
    exit_code = main(
        [
            "alignment",
            "build-tree",
            str(fixture("example_alignment_distance.fasta")),
            "--method",
            "upgma",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "((A:0.0625,B:0.0625)Inner2:0.21875,(C:0.0625,D:0.0625)Inner1:0.21875)Inner3;\n"
    )
    assert payload["metrics"]["method"] == "upgma"


def test_cli_alignment_compare_distance_trees_reports_rooting_difference(
    capsys,
) -> None:
    exit_code = main(
        [
            "alignment",
            "compare-distance-trees",
            str(fixture("example_alignment_distance.fasta")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["robinson_foulds_distance"] == 1
    assert payload["metrics"]["same_unrooted_topology"] is True


def test_cli_distance_validate_reports_imported_matrix_status(capsys) -> None:
    exit_code = main(
        [
            "distance",
            "validate",
            str(fixture("example_distance_matrix_nonmetric.tsv")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["nonmetric_observation_count"] == 1
    assert payload["data"]["warnings"] == [
        "distance matrix violates triangle inequality for one or more taxon triples"
    ]


def test_cli_distance_build_tree_writes_newick(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "imported-tree.nwk"
    exit_code = main(
        [
            "distance",
            "build-tree",
            str(fixture("example_distance_matrix.tsv")),
            "--method",
            "upgma",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert (
        output_path.read_text(encoding="utf-8")
        == "((A:0.0625,B:0.0625)Inner1:0.21875,C:0.28125)Inner2;\n"
    )
    assert payload["metrics"]["method"] == "upgma"


def test_cli_distance_report_writes_html(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "distance-report.html"
    exit_code = main(
        [
            "distance",
            "report",
            str(fixture("example_distance_matrix.tsv")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert "distance-method-limitations" in output_path.read_text(encoding="utf-8")
    assert payload["metrics"]["section_count"] >= 3


def test_cli_distance_explain_reports_limitations(capsys) -> None:
    exit_code = main(
        ["distance", "explain", str(fixture("example_distance_matrix.tsv")), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["limitation_count"] == 4


def test_cli_tree_set_consensus_writes_newick(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "consensus.nwk"
    exit_code = main(
        [
            "tree-set",
            "consensus",
            str(fixture("example_tree_set_left.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "((A:0.1,B:0.1)66.6666666666667:0.2,(C:0.1,D:0.1)66.6666666666667:0.2);\n"
    )
    assert payload["metrics"]["tree_count"] == 3
    assert payload["metrics"]["consensus_method"] == "majority-rule"
    assert payload["metrics"]["consensus_threshold"] == 0.5
    assert payload["metrics"]["included_clade_count"] == 2


def test_cli_tree_set_consensus_supports_strict_mode_and_frequency_ledger(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "strict-consensus.nwk"
    frequency_path = tmp_path / "clade-frequencies.tsv"

    exit_code = main(
        [
            "tree-set",
            "consensus",
            str(fixture("example_tree_set_left.nwk")),
            "--out",
            str(output_path),
            "--method",
            "strict",
            "--clade-frequencies-out",
            str(frequency_path),
            "--json",
        ]
    )

    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == "(A:0.1,B:0.1,C:0.1,D:0.1);\n"
    assert frequency_path.read_text(encoding="utf-8") == (
        "clade\ttree_count\tfrequency\n"
        "A|B\t2\t0.666666666666667\n"
        "A|C\t1\t0.333333333333333\n"
        "B|D\t1\t0.333333333333333\n"
        "C|D\t2\t0.666666666666667\n"
    )
    assert payload["metrics"]["consensus_method"] == "strict"
    assert payload["metrics"]["consensus_threshold"] == 1.0
    assert payload["metrics"]["included_clade_count"] == 0


def test_cli_tree_set_support_map_writes_reference_support_table(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "reference-support.tsv"

    exit_code = main(
        [
            "tree-set",
            "support-map",
            str(fixture("example_tree.nwk")),
            str(fixture("example_tree_set_left.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )

    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert output_path.exists()
    assert "node_id\tnode_kind\tnode_label\tdescendant_taxa" in output_path.read_text(
        encoding="utf-8"
    )
    assert payload["metrics"]["tree_count"] == 3
    assert payload["metrics"]["supported_clade_count"] == 2
    assert payload["metrics"]["absent_clade_count"] == 0
    assert payload["metrics"]["unscored_clade_count"] == 0


def test_cli_tree_set_compare_reports_shared_topologies(capsys) -> None:
    exit_code = main(
        [
            "tree-set",
            "compare",
            str(fixture("example_tree_set_left.nwk")),
            str(fixture("example_tree_set_right.nwk")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["shared_rooted_topology_count"] == 1
    assert (
        payload["data"]["mean_between_set_normalized_robinson_foulds"]
        == 0.777777777777778
    )


def test_cli_tree_set_report_writes_html(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "tree-set-report.html"
    exit_code = main(
        [
            "tree-set",
            "report",
            str(fixture("example_tree_set_left.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert "unstable-taxa" in output_path.read_text(encoding="utf-8")
    assert payload["metrics"]["section_count"] == 11


def test_cli_simulate_birth_death_writes_tree_set(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "simulated.trees"
    exit_code = main(
        [
            "simulate",
            "tree-birth-death",
            "--tree-count",
            "2",
            "--tip-count",
            "4",
            "--seed",
            "7",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["tree_count"] == 2
    assert output_path.read_text(encoding="utf-8").count(";\n") == 2


def test_cli_simulate_random_tree_writes_envelope_ledgers(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "simulated-random.trees"
    record_path = tmp_path / "simulation-records.tsv"
    envelope_path = tmp_path / "simulation-envelope.tsv"
    exit_code = main(
        [
            "simulate",
            "tree-random",
            "--tree-count",
            "2",
            "--tip-count",
            "4",
            "--seed",
            "7",
            "--out",
            str(output_path),
            "--record-table-out",
            str(record_path),
            "--envelope-table-out",
            str(envelope_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["tree_count"] == 2
    assert payload["metrics"]["branch_length_model"] == "uniform"
    assert payload["metrics"]["envelope_metric_count"] == 6
    assert output_path.read_text(encoding="utf-8").count(";\n") == 2
    assert "normalized_colless_imbalance" in record_path.read_text(encoding="utf-8")
    assert "branch_length\tedge\t12\t" in envelope_path.read_text(encoding="utf-8")


def test_cli_simulate_coalescent_writes_envelope_ledgers(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "simulated-coalescent.trees"
    record_path = tmp_path / "simulation-records.tsv"
    envelope_path = tmp_path / "simulation-envelope.tsv"
    exit_code = main(
        [
            "simulate",
            "tree-coalescent",
            "--tree-count",
            "2",
            "--tip-count",
            "4",
            "--population-size",
            "2.5",
            "--seed",
            "7",
            "--out",
            str(output_path),
            "--record-table-out",
            str(record_path),
            "--envelope-table-out",
            str(envelope_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["tree_count"] == 2
    assert payload["metrics"]["pooled_branch_count"] == 12
    assert payload["metrics"]["envelope_metric_count"] == 6
    assert output_path.read_text(encoding="utf-8").count(";\n") == 2
    assert "tree_height_branch_length" in record_path.read_text(encoding="utf-8")
    assert "total_branch_length\ttree\t2\t" in envelope_path.read_text(
        encoding="utf-8"
    )


def test_cli_simulate_dna_alignment_writes_fasta(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "simulated.fasta"
    exit_code = main(
        [
            "simulate",
            "alignment-dna",
            str(fixture("example_tree.nwk")),
            "--sequence-length",
            "8",
            "--substitution-rate",
            "1.2",
            "--seed",
            "7",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["sequence_length"] == 8
    assert ">A\nACTAACGA\n" in output_path.read_text(encoding="utf-8")


@pytest.mark.slow
def test_cli_benchmark_tree_validation_reports_observations(capsys) -> None:
    exit_code = main(["benchmark", "tree-validation", "--replicates", "1", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["observation_count"] == 3
    assert payload["data"]["replicates"] == 1


def test_cli_alignment_coding_reports_frameshifts_and_stops(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "coding",
            str(fixture("example_alignment_coding.fasta")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["frameshift_like_sequence_count"] == 1
    assert payload["metrics"]["stop_codon_count"] == 2


def test_cli_alignment_translate_writes_amino_acid_alignment(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "translated.fasta"
    codon_table_path = tmp_path / "codon-validation.tsv"
    exclusion_table_path = tmp_path / "excluded-sequences.tsv"
    exit_code = main(
        [
            "alignment",
            "translate",
            str(fixture("example_alignment_coding.fasta")),
            "--out",
            str(output_path),
            "--codon-validation-out",
            str(codon_table_path),
            "--excluded-sequences-out",
            str(exclusion_table_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        ">A\nME*\n>B\nMEW\n>C\nMXW\n>D\nM*W\n"
    )
    assert codon_table_path.read_text(encoding="utf-8").splitlines()[0] == (
        "identifier\tcodon_index\tnucleotide_start\tcodon\tamino_acid\ttranslation_status"
    )
    assert exclusion_table_path.read_text(encoding="utf-8") == "identifier\treason\tnote\n"
    assert payload["metrics"]["translated_sequence_count"] == 4
    assert payload["metrics"]["invalid_codon_count"] == 1
    assert payload["metrics"]["stop_codon_count"] == 2
    assert payload["metrics"]["internal_stop_sequence_count"] == 1
    assert payload["metrics"]["terminal_stop_sequence_count"] == 1


def test_cli_alignment_translate_truncates_trailing_partial_codon(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "translated.fasta"
    exit_code = main(
        [
            "alignment",
            "translate",
            str(fixture("example_alignment_coding_frame_error.fasta")),
            "--out",
            str(output_path),
            "--json",
        ]
    )

    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == ">frame_error\nME\n"
    assert payload["metrics"]["translated_alignment_length"] == 2
    assert payload["metrics"]["dropped_trailing_nucleotide_count"] == 2
    assert payload["warnings"] == [
        "sequence length not a multiple of 3: 2 nucleotides dropped"
    ]


def test_cli_topology_root_outgroup_writes_rooted_tree_and_report(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "rooted.nwk"
    report_path = tmp_path / "rooting.tsv"
    exit_code = main(
        [
            "topology",
            "root-outgroup",
            str(fixture("example_tree_rootable.nwk")),
            "--taxa",
            "D",
            "Z",
            "--out",
            str(output_path),
            "--report-out",
            str(report_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert (
        output_path.read_text(encoding="utf-8")
        == "(((A:0.2,B:0.2):0.7,C:0.1):0,D:0.1);\n"
    )
    report_text = report_path.read_text(encoding="utf-8")
    assert "outgroup_monophyletic" in report_text
    assert "\ttrue\tD\t\tD\tA,B,C\t" in report_text
    assert payload["metrics"]["matched_taxa"] == 1
    assert payload["metrics"]["absent_taxa"] == 1
    assert payload["metrics"]["ingroup_taxa"] == 3
    assert payload["metrics"]["outgroup_monophyletic"] is True
    assert payload["metrics"]["outgroup_mrca_extra_taxa"] == 0
    assert payload["metrics"]["rooted_outgroup_taxa"] == 1
    assert payload["metrics"]["rooted_ingroup_taxa"] == 3
    assert payload["metrics"]["warning_count"] == 1
    assert payload["data"]["warnings"] == [
        "one or more requested outgroup taxa were absent from the input tree"
    ]


def test_cli_topology_root_outgroup_reports_non_monophyletic_error(
    tmp_path: Path, capsys
) -> None:
    exit_code = main(
        [
            "topology",
            "root-outgroup",
            str(fixture("example_tree_rootable.nwk")),
            "--taxa",
            "B",
            "D",
            "--out",
            str(tmp_path / "rooted.nwk"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"] == [
        {
            "code": "outgroup_not_monophyletic",
            "details": {
                "matched_taxa": ["B", "D"],
                "outgroup_mrca_extra_taxa": ["A", "C"],
                "outgroup_mrca_taxa": ["A", "B", "C", "D"],
            },
            "message": "requested outgroup taxa are not monophyletic in the input tree",
        }
    ]


def test_cli_topology_reroot_midpoint_writes_tree_and_report(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "midpoint.nwk"
    report_path = tmp_path / "midpoint.tsv"
    exit_code = main(
        [
            "topology",
            "reroot-midpoint",
            str(fixture("example_tree_rootable.nwk")),
            "--out",
            str(output_path),
            "--report-out",
            str(report_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert (
        output_path.read_text(encoding="utf-8")
        == "((A:0.2,B:0.2):0.3,(C:0.1,D:0.1):0.4);\n"
    )
    report_text = report_path.read_text(encoding="utf-8")
    assert "midpoint_anchor_taxa" in report_text
    assert "\tA,C\t1.0\t0.5\tnode\tA,B\tC,D\ttrue\n" in report_text
    assert payload["data"]["strategy"] == "midpoint"
    assert payload["metrics"]["tip_count"] == 4
    assert payload["metrics"]["midpoint_anchor_taxa"] == 2
    assert payload["metrics"]["midpoint_path_length"] == 1.0
    assert payload["metrics"]["midpoint_position_kind"] == "node"
    assert payload["metrics"]["midpoint_anchor_side_taxa"] == 2
    assert payload["metrics"]["midpoint_opposite_side_taxa"] == 2
    assert payload["metrics"]["midpoint_suitable"] is True


def test_cli_topology_reroot_midpoint_reports_exploratory_warning(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "midpoint-polytomy.nwk"
    report_path = tmp_path / "midpoint-polytomy.tsv"
    exit_code = main(
        [
            "topology",
            "reroot-midpoint",
            str(fixture("example_tree_polytomy.nwk")),
            "--out",
            str(output_path),
            "--report-out",
            str(report_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.exists()
    assert "midpoint_suitable" in report_path.read_text(encoding="utf-8")
    assert payload["metrics"]["midpoint_suitable"] is False
    assert payload["metrics"]["warning_count"] == 1
    assert payload["data"]["warnings"] == [
        "midpoint rooting is exploratory because the input tree is not strictly bifurcating"
    ]


def test_cli_topology_unroot_writes_trifurcating_tree(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "unrooted.nwk"
    exit_code = main(
        [
            "topology",
            "unroot",
            str(fixture("example_tree_rootable.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert (
        output_path.read_text(encoding="utf-8") == "(A:0.2,B:0.2,(C:0.1,D:0.1):0.7);\n"
    )
    assert payload["data"]["strategy"] == "unroot"
    assert payload["metrics"]["tip_count"] == 4
    assert payload["data"]["warnings"] == [
        "unrooting merged the removed root-edge length into the retained sibling branch to match ape::unroot"
    ]


def test_cli_topology_unroot_returns_already_unrooted_tree_unchanged(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "already-unrooted.nwk"
    exit_code = main(
        [
            "topology",
            "unroot",
            str(fixture("example_tree_unrooted.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == "(A:0.1,B:0.2,C:0.3,D:0.4);\n"
    assert payload["data"]["warnings"] == [
        "input tree already behaves as an unrooted representation; returned unchanged"
    ]


def test_build_run_manifest_captures_checksums_and_environment(tmp_path: Path) -> None:
    input_path = tmp_path / "input.txt"
    output_path = tmp_path / "output.txt"
    input_path.write_text("input\n", encoding="utf-8")
    output_path.write_text("output\n", encoding="utf-8")
    manifest = build_run_manifest(
        command="inspect",
        arguments=["inspect", str(input_path), "--json"],
        input_paths=[input_path],
        output_paths=[output_path],
    )
    manifest_path = write_run_manifest(tmp_path / "run.manifest.json", manifest)
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert payload["command"] == "inspect"
    assert payload["arguments"] == ["inspect", str(input_path), "--json"]
    assert payload["input_paths"] == [str(input_path)]
    assert payload["output_paths"] == [str(output_path)]
    assert payload["input_checksums"][str(input_path)]
    assert payload["output_checksums"][str(output_path)]
    assert payload["python_version"]
    assert payload["host_platform"]


def test_validate_tree_path_reports_expected_counts() -> None:
    report = validate_tree_path(fixture("example_tree.nwk"))
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert report.rooted is True
    assert report.ultrametric is True
    assert report.source_format == "newick"


def test_validate_tree_path_rejects_duplicate_tip_labels_by_default() -> None:
    try:
        validate_tree_path(fixture("example_tree_duplicate.nwk"))
    except DuplicateTaxonError as error:
        assert error.code == "duplicate_taxon_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected DuplicateTaxonError")


def test_validate_tree_path_warns_for_unnamed_tips_in_non_strict_mode() -> None:
    report = validate_tree_path(fixture("example_tree_unnamed_tip.nwk"))
    assert report.missing_taxa == 1
    assert "tree contains unnamed tips" in report.warnings


def test_validate_tree_path_rejects_unnamed_tips_in_strict_mode() -> None:
    try:
        validate_tree_path(fixture("example_tree_unnamed_tip.nwk"), strict=True)
    except UnnamedTipError as error:
        assert error.code == "unnamed_tip_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected UnnamedTipError")


def test_validate_tree_path_rejects_negative_branch_lengths_by_default() -> None:
    try:
        validate_tree_path(fixture("example_tree_negative_length.nwk"))
    except InvalidBranchLengthError as error:
        assert error.code == "invalid_branch_length_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected InvalidBranchLengthError")


def test_validate_tree_path_can_require_rooted_tree() -> None:
    try:
        validate_tree_path(fixture("example_tree_unrooted.nwk"), require_rooted=True)
    except UnrootedTreeError as error:
        assert error.code == "unrooted_tree_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected UnrootedTreeError")


def test_validate_tree_path_can_require_ultrametric_tree() -> None:
    try:
        validate_tree_path(
            fixture("example_tree_ladderized.nwk"), require_ultrametric=True
        )
    except NonUltrametricTreeError as error:
        assert error.code == "non_ultrametric_tree_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected NonUltrametricTreeError")


def test_validate_tree_path_accepts_near_ultrametric_tree_when_required() -> None:
    report = validate_tree_path(
        fixture("example_tree_near_ultrametric.nwk"), require_ultrametric=True
    )

    assert report.ultrametric is True


def test_validate_tree_path_warns_for_zero_length_branches() -> None:
    report = validate_tree_path(fixture("example_tree_zero_lengths.nwk"))
    assert report.zero_length_branches == 3
    assert "tree contains zero-length branches" in report.warnings


def test_validate_tree_path_localizes_missing_internal_and_terminal_branch_lengths() -> (
    None
):
    internal = validate_tree_path(fixture("example_tree_missing_internal_length.nwk"))
    terminal = validate_tree_path(fixture("example_tree_partial_lengths.nwk"))
    assert internal.missing_internal_branch_nodes == ["A|B"]
    assert internal.missing_terminal_branch_taxa == []
    assert "tree contains internal branches without lengths" in internal.warnings
    assert terminal.missing_internal_branch_nodes == []
    assert terminal.missing_terminal_branch_taxa == ["B"]
    assert "tree contains terminal branches without lengths" in terminal.warnings


def test_validate_tree_path_detects_singleton_internal_nodes() -> None:
    report = validate_tree_path(fixture("example_tree_singleton.nwk"))
    assert report.singleton_internal_nodes == ["A"]
    assert "tree contains singleton internal nodes" in report.warnings


def test_inspect_tree_path_returns_normalized_json_summary_contract() -> None:
    report = inspect_tree_path(fixture("example_tree.nwk"))
    assert report.tip_count == 4
    assert report.node_count == 7
    assert report.internal_node_count == 3
    assert report.edge_count == 6
    assert report.clade_count == 3
    assert report.has_branch_lengths is True
    assert report.is_binary is True
    assert [(row.node, row.child_count) for row in report.internal_child_counts] == [
        ("A|B|C|D", 2),
        ("A|B", 2),
        ("C|D", 2),
    ]
    assert report.singleton_internal_nodes == []
    assert report.missing_internal_branch_nodes == []
    assert report.missing_terminal_branch_taxa == []
    assert report.is_ultrametric is True
    assert report.branch_length_summary is not None
    assert (
        report.branch_length_summary.count,
        report.branch_length_summary.minimum,
        report.branch_length_summary.maximum,
        report.branch_length_summary.mean,
        report.branch_length_summary.median,
        report.branch_length_summary.first_quartile,
        report.branch_length_summary.third_quartile,
    ) == (6, 0.1, 0.2, 0.15, 0.15, 0.1, 0.2)
    assert report.tree_diameter == 0.6
    assert report.zero_length_branch_count == 0
    assert report.max_depth == 2
    assert report.mean_depth == 2.0
    assert report.colless_imbalance_index == 0.0
    assert report.normalized_colless_imbalance == 0.0
    assert report.sackin_imbalance_index == 8
    assert report.unusually_imbalanced is False
    assert report.long_branch_taxa == []
    assert report.star_like is False
    assert report.comb_like is False
    assert report.tree_quality_score == 100.0
    assert report.tree_quality_warnings == []
    assert report.imbalance_summary == "balanced"
    assert report.cherry_count == 2
    assert report.warnings == []
    assert report.taxa == ["A", "B", "C", "D"]


def test_inspect_tree_path_reports_structural_and_missing_branch_diagnostics() -> None:
    singleton = inspect_tree_path(fixture("example_tree_singleton.nwk"))
    missing_internal = inspect_tree_path(
        fixture("example_tree_missing_internal_length.nwk")
    )
    missing_terminal = inspect_tree_path(fixture("example_tree_partial_lengths.nwk"))
    assert [(row.node, row.child_count) for row in singleton.internal_child_counts] == [
        ("A|B|C", 2),
        ("A|B", 2),
        ("A", 1),
    ]
    assert singleton.singleton_internal_nodes == ["A"]
    assert missing_internal.missing_internal_branch_nodes == ["A|B"]
    assert missing_internal.missing_terminal_branch_taxa == []
    assert missing_terminal.missing_internal_branch_nodes == []
    assert missing_terminal.missing_terminal_branch_taxa == ["B"]


def test_inspect_tree_path_distinguishes_ladderized_shape() -> None:
    report = inspect_tree_path(fixture("example_tree_ladderized.nwk"))
    assert report.tree_diameter == 0.4
    assert report.max_depth == 3
    assert report.mean_depth == 2.25
    assert report.colless_imbalance_index == 3.0
    assert report.normalized_colless_imbalance == 1.0
    assert report.sackin_imbalance_index == 9
    assert report.unusually_imbalanced is True
    assert report.long_branch_taxa == []
    assert report.star_like is False
    assert report.comb_like is True
    assert report.tree_quality_score == 75.0
    assert [warning.code for warning in report.tree_quality_warnings] == [
        "unusually_imbalanced",
        "comb_like",
    ]
    assert report.imbalance_summary == "ladderized"
    assert report.cherry_count == 1


def test_inspect_tree_path_distinguishes_rooted_and_unrooted_fixtures() -> None:
    rooted = inspect_tree_path(fixture("example_tree.nwk"))
    unrooted = inspect_tree_path(fixture("example_tree_unrooted.nwk"))
    assert rooted.rooted is True
    assert unrooted.rooted is False


def test_inspect_tree_path_reports_exact_polytomy_nodes() -> None:
    report = inspect_tree_path(fixture("example_tree_polytomy.nwk"))
    assert report.is_binary is False
    assert report.colless_imbalance_index is None
    assert report.normalized_colless_imbalance is None
    assert report.unusually_imbalanced is None
    assert report.long_branch_taxa == []
    assert report.star_like is False
    assert report.comb_like is False
    assert report.tree_quality_score == 90.0
    assert [warning.code for warning in report.tree_quality_warnings] == ["polytomies"]
    assert report.polytomy_count == 1
    assert report.polytomy_nodes == ["A|B|C"]


def test_inspect_tree_path_detects_long_branch_taxa() -> None:
    report = inspect_tree_path(fixture("example_tree_long_branch.nwk"))
    assert report.long_branch_taxa == ["A"]
    assert [
        (row.node, row.branch_length, row.branch_type)
        for row in report.long_branch_outliers
    ] == [("A", 1.0, "terminal")]
    assert report.short_branch_outliers == []
    assert report.star_like is False
    assert report.tree_quality_score == 85.0
    assert [warning.code for warning in report.tree_quality_warnings] == [
        "long_branches"
    ]


def test_inspect_tree_path_detects_internal_long_and_short_branch_outliers() -> None:
    long_report = inspect_tree_path(fixture("example_tree_internal_long_branch.nwk"))
    short_report = inspect_tree_path(fixture("example_tree_short_branch.nwk"))
    assert [
        (row.node, row.branch_length, row.branch_type)
        for row in long_report.long_branch_outliers
    ] == [("A|B", 1.0, "internal")]
    assert long_report.long_branch_taxa == []
    assert long_report.short_branch_outliers == []
    assert [
        (row.node, row.branch_length, row.branch_type)
        for row in short_report.short_branch_outliers
    ] == [("B", 0.001, "terminal")]
    assert [warning.code for warning in short_report.tree_quality_warnings] == [
        "short_branches"
    ]


def test_inspect_tree_path_classifies_internal_support_and_name_labels() -> None:
    support = inspect_tree_path(fixture("example_tree_support_mixed.nwk"))
    names = inspect_tree_path(fixture("example_tree_named_clades.nwk"))
    assert [
        (row.node, row.label, row.interpretation, row.numeric_value)
        for row in support.likely_support_labels
    ] == [
        ("A|B", "0.95", "fractional_support", 0.95),
        ("A|B|C|D", "99", "percentage_support", 99.0),
        ("C|D", "88", "percentage_support", 88.0),
    ]
    assert support.likely_named_internal_labels == []
    assert [
        (row.node, row.label, row.interpretation)
        for row in names.likely_named_internal_labels
    ] == [
        ("A|B", "Mammals", "named_internal_label"),
        ("A|B|C|D", "Root", "named_internal_label"),
        ("C|D", "Birds", "named_internal_label"),
    ]
    assert names.likely_support_labels == []

    invalid = inspect_tree_path(fixture("example_tree_support_invalid.nwk"))
    assert [
        (row.label, row.interpretation) for row in invalid.likely_support_labels
    ] == [
        ("120", "out_of_range_support"),
        ("101", "out_of_range_support"),
        ("-5", "out_of_range_support"),
    ]


def test_inspect_tree_path_detects_suspicious_and_mixed_support_scales() -> None:
    invalid = inspect_tree_path(fixture("example_tree_support_invalid.nwk"))
    mixed = inspect_tree_path(fixture("example_tree_support_mixed.nwk"))
    assert invalid.suspicious_support_value_ranges == [
        "support value 101 at A|B|C|D exceeds 100",
        "support value 120 at A|B exceeds 100",
        "support value -5 at C|D is negative",
    ]
    assert invalid.mixed_support_scales is False
    assert mixed.suspicious_support_value_ranges == []
    assert mixed.mixed_support_scales is True
    assert [warning.code for warning in mixed.tree_quality_warnings] == [
        "mixed_support_scales"
    ]


def test_standardize_support_labels_normalizes_fraction_and_percentage_scales() -> None:
    rows = standardize_support_labels(fixture("example_tree_support_mixed.nwk"))
    assert [
        (row.node, row.raw_value, row.scale, row.support_fraction, row.support_percent)
        for row in rows
    ] == [
        ("A|B", 0.95, "fraction", 0.95, 95.0),
        ("A|B|C|D", 99.0, "percentage", 0.99, 99.0),
        ("C|D", 88.0, "percentage", 0.88, 88.0),
    ]
    assert [
        (row.node, row.normalized_probability, row.confidence_of_inference)
        for row in rows
    ] == [
        ("A|B", 0.95, "medium"),
        ("A|B|C|D", 0.99, "medium"),
        ("C|D", 0.88, "medium"),
    ]


def test_validate_tree_roundtrip_preserves_tree_structure_across_formats() -> None:
    nexus_report = validate_tree_roundtrip(
        fixture("example_tree.nwk"), target_format="nexus"
    )
    phyloxml_report = validate_tree_roundtrip(
        fixture("example_tree.nwk"), target_format="phyloxml"
    )
    assert nexus_report.preserved_taxa is True
    assert nexus_report.preserved_topology is True
    assert nexus_report.preserved_branch_lengths is True
    assert nexus_report.preserved_support_labels is True
    assert phyloxml_report.preserved_taxa is True
    assert phyloxml_report.preserved_topology is True


def test_validate_tree_roundtrip_reports_semantic_loss_for_nexus_and_phyloxml() -> None:
    nexus_report = validate_tree_roundtrip(
        fixture("example_tree.nex"), target_format="newick"
    )
    phyloxml_report = validate_tree_roundtrip(
        fixture("example_tree_annotated.phyloxml"), target_format="newick"
    )
    assert [item.feature for item in nexus_report.semantic_loss] == [
        "nexus-metadata-blocks",
        "nexus-translate",
        "target-newick-structure",
    ]
    assert sorted(item.feature for item in phyloxml_report.semantic_loss) == [
        "property",
        "target-newick-structure",
        "taxonomy",
    ]


def test_taxon_identity_audit_reports_collisions_and_near_duplicates() -> None:
    report = inspect_tree_taxon_identity(load_nexus(fixture("example_tree.nex")))
    assert report.spelling_variants == []
    assert report.whitespace_variants == []
    identity_report = inspect_tree_taxon_identity(
        load_phyloxml(fixture("example_tree_annotated.phyloxml"))
    )
    assert identity_report.case_collisions == []

    ambiguous = inspect_tree_taxon_identity(
        loads_newick(fixture("example_tree_identity.nwk").read_text(encoding="utf-8"))
    )
    assert {
        (row.left_label, row.right_label)
        for row in ambiguous.underscore_space_collisions
    } == {
        ("Homo sapiens", "Homo_sapiens"),
        ("Homo sapiens", "homo sapiens"),
        ("Homo_sapiens", "homo sapiens"),
    }
    assert [(row.left_label, row.right_label) for row in ambiguous.case_collisions] == [
        ("Homo sapiens", "homo sapiens"),
    ]
    assert ("Homo sapiens", "Hoomo sapiens") in [
        (row.left_label, row.right_label)
        for row in ambiguous.suspicious_near_duplicates
    ]


def test_inspect_branch_length_units_reports_time_and_substitution_metadata() -> None:
    time_report = inspect_branch_length_units(
        fixture("example_metadata_branch_units_time.tsv")
    )
    substitution_report = inspect_branch_length_units(
        fixture("example_metadata_branch_units_substitution.tsv")
    )
    conflict_report = inspect_branch_length_units(
        fixture("example_metadata_branch_units_conflict.tsv")
    )
    assert (
        time_report.declared_unit,
        time_report.compatible_with_time_tree,
        time_report.compatible_with_substitution_tree,
    ) == (
        "years",
        True,
        False,
    )
    assert (
        substitution_report.declared_unit,
        substitution_report.compatible_with_time_tree,
        substitution_report.compatible_with_substitution_tree,
    ) == (
        "substitutions_per_site",
        False,
        True,
    )
    assert conflict_report.conflicting_values == ["substitutions_per_site", "years"]


def test_assess_tree_assumptions_reports_time_and_substitution_compatibility() -> None:
    time_report = assess_tree_assumptions(
        fixture("example_tree.nwk"),
        metadata_path=fixture("example_metadata_branch_units_time.tsv"),
    )
    substitution_report = assess_tree_assumptions(
        fixture("example_tree.nwk"),
        metadata_path=fixture("example_metadata_branch_units_substitution.tsv"),
    )
    partial_report = assess_tree_assumptions(
        fixture("example_tree_partial_lengths.nwk"),
        metadata_path=fixture("example_metadata_branch_units_substitution.tsv"),
    )
    assert time_report.time_tree_compatible is True
    assert time_report.substitution_tree_compatible is False
    assert substitution_report.time_tree_compatible is False
    assert substitution_report.substitution_tree_compatible is True
    assert partial_report.time_tree_compatible is False
    assert "tree requires complete branch lengths" in partial_report.blockers


def test_inspect_tree_path_detects_star_like_tree() -> None:
    report = inspect_tree_path(fixture("example_tree_star.nwk"))
    assert report.star_like is True
    assert report.long_branch_taxa == []
    assert report.tree_quality_score == 80.0
    assert [warning.code for warning in report.tree_quality_warnings] == [
        "polytomies",
        "star_like",
    ]


def test_inspect_tree_path_classifies_branch_length_completeness() -> None:
    complete = inspect_tree_path(fixture("example_tree.nwk"))
    partial = inspect_tree_path(fixture("example_tree_partial_lengths.nwk"))
    absent = inspect_tree_path(fixture("example_tree_no_lengths.nwk"))
    assert complete.branch_length_status == "complete"
    assert partial.branch_length_status == "partial"
    assert absent.branch_length_status == "absent"
    assert complete.branch_length_summary is not None
    assert partial.branch_length_summary is not None
    assert absent.branch_length_summary is None
    assert complete.tree_diameter == 0.6
    assert partial.tree_diameter is None
    assert absent.tree_diameter is None
    assert partial.has_branch_lengths is True
    assert partial.warnings == [
        "tree contains terminal branches without lengths",
        "tree contains partial branch lengths",
    ]
    assert absent.has_branch_lengths is False
    assert absent.warnings == [
        "tree contains internal branches without lengths",
        "tree contains terminal branches without lengths",
        "tree contains no branch lengths",
    ]
    assert complete.root_state_confidence.classification == "apparently_rooted"


def test_validate_tree_path_reports_validity_decision_contexts_and_repair_guidance() -> (
    None
):
    report = validate_tree_path(
        fixture("example_tree_partial_lengths.nwk"), allow_duplicates=True
    )
    assert report.syntax_valid is True
    assert report.biologically_safe is False
    assert report.validity_decision == "invalid"
    assert report.root_state_confidence.classification == "apparently_rooted"
    assert {
        context.context: context.allowed for context in report.branch_length_contexts
    } == {
        "topology_only": True,
        "substitution_tree": False,
        "time_tree": False,
        "comparative_methods": False,
    }
    assert [item.issue_code for item in report.branch_length_repair_suggestions] == [
        "partial_branch_lengths"
    ]


def test_forensic_tree_path_reports_unsafe_external_labels_and_downstream_safety() -> (
    None
):
    report = forensic_tree_path(fixture("example_tree_labels.nwk"))
    assert report.validity_decision == "valid_with_warnings"
    assert report.safe_for_topology_comparison is True
    assert report.safe_for_time_tree_analysis is False
    assert report.safe_for_comparative_methods is False
    assert report.unsafe_external_labels[0].engines == [
        "iqtree",
        "raxml",
        "mrbayes",
        "beast",
        "r",
        "shell",
    ]
    assert any(
        item.raw_label == "Homo sapiens" for item in report.unsafe_external_labels
    )


def test_newick_loader_raises_invalid_branch_length_error() -> None:
    try:
        loads_newick("((A:abc,B:0.2):0.3,C:0.4);")
    except InvalidBranchLengthError as error:
        assert error.code == "invalid_branch_length_error"
        assert error.details == {"position": 4, "line": 1, "column": 5}
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected InvalidBranchLengthError")


def test_newick_loader_preserves_numeric_support_metadata() -> None:
    tree = loads_newick("((A:0.1,B:0.1)95:0.2,C:0.3)Root;")

    support_node = next(
        node
        for node in tree.iter_internal_nodes()
        if node is not tree.root and set(node.descendant_taxa) == {"A", "B"}
    )

    assert support_node.name == "95"
    assert support_node.metadata["confidence"] == 95.0


def test_newick_loader_reports_location_for_malformed_structure() -> None:
    try:
        loads_newick("((A:0.1,B:0.2):0.3,C:0.4")
    except TreeParseError as error:
        assert error.details == {"position": 24, "line": 1, "column": 25}
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected TreeParseError")


def test_nexus_loader_reads_translation_block_fixture() -> None:
    tree = load_nexus(fixture("example_tree.nex"))
    assert tree.source_format == "nexus"
    assert tree.tip_names == ["A", "B", "C", "D"]
    assert tree.tip_count == 4


def test_phyloxml_loader_reads_annotated_tree_fixture() -> None:
    tree = load_phyloxml(fixture("example_tree.phyloxml"))
    assert tree.source_format == "phyloxml"
    assert tree.tip_names == ["A", "B", "C"]
    assert tree.tip_count == 3


def test_detect_tree_format_uses_filename_suffixes() -> None:
    assert detect_tree_format(Path("x.nwk")) == "newick"
    assert detect_tree_format(Path("x.nex")) == "nexus"
    assert detect_tree_format(Path("x.phyloxml")) == "phyloxml"


def test_detect_tree_format_prefers_file_content_over_misleading_suffix(
    tmp_path: Path,
) -> None:
    path = tmp_path / "tree.nex"
    path.write_text("(A:1,B:1);\n", encoding="utf-8")
    assert detect_tree_format(path) == "newick"


def test_validate_cli_reports_unsupported_format_error(tmp_path: Path, capsys) -> None:
    path = tmp_path / "tree.txt"
    path.write_text("not a tree\n", encoding="utf-8")
    exit_code = main(["validate", str(path), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"] == [
        {
            "code": UnsupportedTreeFormatError.code,
            "message": f"unsupported tree format for {path}",
        }
    ]


def test_validate_cli_can_allow_duplicate_tip_labels(capsys) -> None:
    exit_code = main(
        [
            "validate",
            str(fixture("example_tree_duplicate.nwk")),
            "--allow-duplicates",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["data"]["duplicate_taxa"] == ["A"]


def test_validate_cli_strict_mode_rejects_unnamed_tips(capsys) -> None:
    exit_code = main(
        ["validate", str(fixture("example_tree_unnamed_tip.nwk")), "--strict", "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"] == [
        {
            "code": UnnamedTipError.code,
            "message": "tree contains 1 unnamed tip labels",
        }
    ]


def test_validate_cli_can_allow_negative_branch_lengths(capsys) -> None:
    exit_code = main(
        [
            "validate",
            str(fixture("example_tree_negative_length.nwk")),
            "--allow-negative-branches",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["data"]["negative_branch_lengths"] == 1


def test_validate_cli_can_require_rooted_and_ultrametric_typed_errors(capsys) -> None:
    exit_code = main(
        [
            "validate",
            str(fixture("example_tree_unrooted.nwk")),
            "--require-rooted",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["errors"][0]["code"] == UnrootedTreeError.code

    exit_code = main(
        [
            "validate",
            str(fixture("example_tree_ladderized.nwk")),
            "--require-ultrametric",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["errors"][0]["code"] == NonUltrametricTreeError.code


def test_cli_inspect_accepts_explicit_tree_format(capsys) -> None:
    exit_code = main(
        ["inspect", str(fixture("example_tree.nex")), "--format", "nexus", "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["data"]["source_format"] == "nexus"
    assert payload["data"]["node_count"] == 7
    assert payload["data"]["edge_count"] == 6
    assert payload["data"]["clade_count"] == 3
    assert payload["data"]["tree_diameter"] == 1.6
    assert payload["data"]["mean_depth"] == 2.0
    assert payload["data"]["colless_imbalance_index"] == 0.0
    assert payload["data"]["sackin_imbalance_index"] == 8
    assert payload["data"]["imbalance_summary"] == "balanced"
    assert payload["data"]["tree_quality_score"] == 100.0
    assert payload["metrics"]["cherry_count"] == 2
    assert payload["metrics"]["tree_diameter"] == 1.6
    assert payload["metrics"]["tree_quality_score"] == 100.0
    assert payload["metrics"]["likely_support_label_count"] == 0
    assert payload["data"]["taxa"] == ["A", "B", "C", "D"]
    assert payload["metrics"]["tip_count"] == 4


def test_cli_validate_and_inspect_surface_structural_and_support_diagnostics(
    capsys,
) -> None:
    validate_exit = main(
        ["validate", str(fixture("example_tree_singleton.nwk")), "--json"]
    )
    validate_payload = json.loads(capsys.readouterr().out)
    assert validate_exit == 0
    assert validate_payload["metrics"]["singleton_internal_node_count"] == 1

    inspect_exit = main(
        ["inspect", str(fixture("example_tree_support_mixed.nwk")), "--json"]
    )
    inspect_payload = json.loads(capsys.readouterr().out)
    assert inspect_exit == 0
    assert inspect_payload["metrics"]["likely_support_label_count"] == 3
    assert inspect_payload["metrics"]["suspicious_support_range_count"] == 0
    assert inspect_payload["data"]["mixed_support_scales"] is True


def test_cli_normalize_writes_canonical_newick(tmp_path: Path, capsys) -> None:
    output = tmp_path / "normalized.nwk"
    exit_code = main(
        [
            "normalize",
            str(fixture("example_tree.nex")),
            "--format",
            "nexus",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["outputs"] == [str(output)]
    assert (
        output.read_text(encoding="utf-8").strip()
        == "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);"
    )


def test_cli_normalize_taxa_writes_mapping_file(tmp_path: Path, capsys) -> None:
    output = tmp_path / "normalized-taxa.nwk"
    mapping = tmp_path / "normalized-taxa.tsv"
    exit_code = main(
        [
            "normalize-taxa",
            str(fixture("example_tree_labels.nwk")),
            "--policy",
            "spaces-to-underscores",
            "--out",
            str(output),
            "--mapping-out",
            str(mapping),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["renamed_taxa"] == 2
    assert (
        output.read_text(encoding="utf-8").strip()
        == "(A.B-1:0.3,Homo_sapiens:0.1,Mus_musculus:0.2);"
    )
    assert mapping.read_text(encoding="utf-8") == (
        "raw_label\tnormalized_label\n"
        "Homo sapiens\tHomo_sapiens\n"
        "Mus musculus\tMus_musculus\n"
    )


def test_compare_tree_paths_reports_nonzero_distance() -> None:
    report = compare_tree_paths(
        fixture("example_tree.nwk"), fixture("example_tree_alt.nwk")
    )
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.robinson_foulds_distance > 0


def test_compare_robinson_foulds_matches_reference_fixture_cases() -> None:
    for row in _load_robinson_foulds_reference_rows():
        report = compare_robinson_foulds(
            fixture(row["left_tree"]),
            fixture(row["right_tree"]),
            rf_mode=row["rf_mode"],
            taxon_overlap_policy=row["taxon_overlap_policy"],
        )
        assert report.left_split_count == int(row["left_split_count"])
        assert report.right_split_count == int(row["right_split_count"])
        assert report.robinson_foulds_distance == int(row["robinson_foulds_distance"])
        assert report.normalized_robinson_foulds == pytest.approx(
            row["normalized_robinson_foulds"]
            if isinstance(row["normalized_robinson_foulds"], float)
            else float(row["normalized_robinson_foulds"]),
            abs=1e-12,
        )


def test_prune_trees_to_shared_taxa_keeps_identical_tip_sets() -> None:
    left, right, report = prune_trees_to_shared_taxa(
        fixture("example_tree.nwk"),
        fixture("example_tree_overlap.nwk"),
    )
    assert left.tip_names == ["A", "B", "C"]
    assert right.tip_names == ["A", "B", "C"]
    assert report.shared_taxa == ["A", "B", "C"]
    assert report.left_only_taxa == ["D"]
    assert report.right_only_taxa == ["E"]


def test_compare_tree_paths_reports_identical_topology_boolean() -> None:
    report = compare_tree_paths(
        fixture("example_tree.nwk"), fixture("example_tree.nwk")
    )
    assert report.topology_equal is True
    assert report.same_unrooted_topology is True
    assert report.same_taxa_different_rooting is False
    assert report.same_topology_different_branch_lengths is False
    assert report.robinson_foulds_distance == 0


def test_compare_tree_paths_supports_explicit_unrooted_rf_mode() -> None:
    report = compare_tree_paths(
        fixture("example_tree.nwk"),
        fixture("example_tree_rooting_diff.nwk"),
        rf_mode="unrooted",
    )
    assert report.rf_mode == "unrooted"
    assert report.robinson_foulds_distance == 0
    assert report.unrooted_robinson_foulds_distance == 0
    assert report.rooted_robinson_foulds_distance == 2
    assert report.same_unrooted_topology is True
    assert report.same_taxa_different_rooting is True


def test_compare_tree_paths_reports_different_topology_boolean() -> None:
    report = compare_tree_paths(
        fixture("example_tree.nwk"), fixture("example_tree_topology_diff.nwk")
    )
    assert report.topology_equal is False
    assert report.same_unrooted_topology is False
    assert report.same_taxa_different_rooting is False


def test_compare_tree_paths_enforces_identical_taxa_when_requested() -> None:
    try:
        compare_tree_paths(
            fixture("example_tree.nwk"),
            fixture("example_tree_overlap.nwk"),
            taxon_overlap_policy="require-identical",
        )
    except ValueError as error:
        assert "requires identical taxon sets" in str(error)
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected identical-taxon RF comparison to fail")


def test_compare_clade_sets_reports_shared_and_unique_clades() -> None:
    report = compare_clade_sets(
        fixture("example_tree.nwk"), fixture("example_tree_alt.nwk")
    )
    assert report.shared_clades == ["A|B"]
    assert report.left_only_clades == ["C|D"]
    assert report.right_only_clades == ["A|B|C"]


def test_compare_clade_overlap_reports_multi_tree_presence_and_support() -> None:
    report = compare_clade_overlap(
        [
            fixture("example_tree.nwk"),
            fixture("example_tree_support_left.nwk"),
            fixture("example_tree_alt.nwk"),
        ]
    )
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.shared_clades == ["A|B"]
    assert report.conflicting_clades == ["C|D", "A|B|C"]
    assert [
        (
            summary.tree_path.name,
            summary.clade_count,
            summary.support_clade_count,
            summary.unique_clades,
        )
        for summary in report.tree_summaries
    ] == [
        ("example_tree.nwk", 2, 0, []),
        ("example_tree_support_left.nwk", 2, 2, []),
        ("example_tree_alt.nwk", 2, 0, ["A|B|C"]),
    ]
    by_clade = {row.clade_id: row for row in report.clade_rows}
    assert [
        (item.tree_path.name, item.present, item.support)
        for item in by_clade["A|B"].observations
    ] == [
        ("example_tree.nwk", True, None),
        ("example_tree_support_left.nwk", True, 95.0),
        ("example_tree_alt.nwk", True, None),
    ]
    assert [
        (item.tree_path.name, item.present, item.support)
        for item in by_clade["C|D"].observations
    ] == [
        ("example_tree.nwk", True, None),
        ("example_tree_support_left.nwk", True, 88.0),
        ("example_tree_alt.nwk", False, None),
    ]


def test_write_clade_overlap_table_writes_one_row_per_clade_per_tree(
    tmp_path: Path,
) -> None:
    output = tmp_path / "clade-overlap.tsv"
    write_clade_overlap_table(
        output,
        [
            fixture("example_tree.nwk"),
            fixture("example_tree_support_left.nwk"),
            fixture("example_tree_alt.nwk"),
        ],
    )
    assert output.read_text(encoding="utf-8") == (
        "clade_id\ttree_path\tpresent\tsupport\tpresent_in_all_trees\tpresent_tree_count\tabsent_tree_count\n"
        f"A|B\t{fixture('example_tree.nwk')}\ttrue\t\ttrue\t3\t0\n"
        f"A|B\t{fixture('example_tree_support_left.nwk')}\ttrue\t95.0\ttrue\t3\t0\n"
        f"A|B\t{fixture('example_tree_alt.nwk')}\ttrue\t\ttrue\t3\t0\n"
        f"C|D\t{fixture('example_tree.nwk')}\ttrue\t\tfalse\t2\t1\n"
        f"C|D\t{fixture('example_tree_support_left.nwk')}\ttrue\t88.0\tfalse\t2\t1\n"
        f"C|D\t{fixture('example_tree_alt.nwk')}\tfalse\t\tfalse\t2\t1\n"
        f"A|B|C\t{fixture('example_tree.nwk')}\tfalse\t\tfalse\t1\t2\n"
        f"A|B|C\t{fixture('example_tree_support_left.nwk')}\tfalse\t\tfalse\t1\t2\n"
        f"A|B|C\t{fixture('example_tree_alt.nwk')}\ttrue\t\tfalse\t1\t2\n"
    )


def test_compare_tree_paths_detects_same_topology_with_different_branch_lengths() -> (
    None
):
    report = compare_tree_paths(
        fixture("example_tree.nwk"), fixture("example_tree_branch_lengths_right.nwk")
    )
    assert report.topology_equal is True
    assert report.same_topology_different_branch_lengths is True


def test_compare_tree_paths_detects_same_taxa_with_different_rooting() -> None:
    report = compare_tree_paths(
        fixture("example_tree.nwk"), fixture("example_tree_rooting_diff.nwk")
    )
    assert report.topology_equal is False
    assert report.same_unrooted_topology is True
    assert report.same_taxa_different_rooting is True


def test_detect_clade_changes_reports_lost_and_gained_sets() -> None:
    report = detect_clade_changes(
        fixture("example_tree.nwk"), fixture("example_tree_alt.nwk")
    )
    assert report.lost_clades == ["C|D"]
    assert report.gained_clades == ["A|B|C"]


def test_compare_support_values_pairs_shared_clades() -> None:
    report = compare_support_values(
        fixture("example_tree_support_left.nwk"),
        fixture("example_tree_support_right.nwk"),
    )
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert [
        (row.split_id, row.left_support, row.right_support)
        for row in report.shared_clades
    ] == [
        ("A|B", 95.0, 90.0),
        ("C|D", 88.0, 85.0),
    ]


def test_compare_support_values_prefers_ufboot_when_iqtree_dual_support_is_present(
    tmp_path: Path,
) -> None:
    left_path = tmp_path / "left.nwk"
    right_path = tmp_path / "right.nwk"
    left_path.write_text(
        "((A:0.1,B:0.1)82/97:0.2,(C:0.1,D:0.1)79/96:0.2);\n",
        encoding="utf-8",
    )
    right_path.write_text(
        "((A:0.1,B:0.1)81/93:0.2,(C:0.1,D:0.1)78/94:0.2);\n",
        encoding="utf-8",
    )

    report = compare_support_values(left_path, right_path)

    assert [
        (row.split_id, row.left_support, row.right_support)
        for row in report.shared_clades
    ] == [
        ("A|B", 97.0, 93.0),
        ("C|D", 96.0, 94.0),
    ]


def test_compare_branch_lengths_reports_delta_ratio_and_missing_lengths() -> None:
    scaled = compare_branch_lengths(
        fixture("example_tree.nwk"), fixture("example_tree_branch_lengths_right.nwk")
    )
    missing = compare_branch_lengths(
        fixture("example_tree.nwk"), fixture("example_tree_branch_lengths_missing.nwk")
    )
    assert [(row.split_id, row.delta, row.ratio) for row in scaled.shared_splits] == [
        ("A|B", 0.2, 2.0),
        ("C|D", 0.1, 2.0),
    ]
    assert scaled.same_taxon_set is True
    assert scaled.branch_score.branch_score_distance == 0.33166247903554
    assert [
        (row.split_id, row.comparison_status, row.left_length, row.right_length)
        for row in scaled.branch_score.splits
    ] == [
        ("A", "shared", 0.1, 0.2),
        ("B", "shared", 0.1, 0.2),
        ("C", "shared", 0.2, 0.2),
        ("D", "shared", 0.2, 0.2),
        ("A|B", "shared", 0.30000000000000004, 0.6000000000000001),
    ]
    assert [
        (row.split_id, row.left_length, row.right_length, row.delta, row.ratio)
        for row in missing.shared_splits
    ] == [
        ("A|B", 0.2, None, None, None),
        ("C|D", 0.1, 0.2, 0.1, 2.0),
    ]
    assert missing.branch_score.branch_score_distance is None
    assert missing.branch_score.missing_length_split_count == 1


def test_compare_branch_score_distance_matches_reference_fixture_cases() -> None:
    for row in _load_branch_score_reference_rows():
        report = compare_branch_score_distance(
            fixture(row["left_tree"]),
            fixture(row["right_tree"]),
            taxon_overlap_policy=row["taxon_overlap_policy"],
        )
        assert report.same_taxon_set is (row["same_taxon_set"] == "true")
        assert report.branch_score_distance == pytest.approx(
            float(row["branch_score_distance"]),
            abs=1e-12,
        )


def test_compare_branch_score_distance_enforces_identical_taxa_when_requested() -> None:
    try:
        compare_branch_score_distance(
            fixture("example_tree.nwk"),
            fixture("example_tree_overlap.nwk"),
            taxon_overlap_policy="require-identical",
        )
    except ValueError as error:
        assert "identical taxon sets" in str(error)
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected identical-taxon branch-score comparison to fail")


def test_build_tree_comparison_report_writes_html_with_checksums(
    tmp_path: Path,
) -> None:
    output = tmp_path / "compare.html"
    result = build_tree_comparison_report(
        fixture("example_tree_support_left.nwk"),
        fixture("example_tree_support_right.nwk"),
        out_path=output,
    )
    html = output.read_text(encoding="utf-8")
    assert result.output_path == output
    assert "Bijux Tree Comparison Report" in html
    assert "input-checksums" in html
    assert "clade-comparison" in html
    assert "clade-changes" in html
    assert "support-comparison" in html
    assert "conflicting_clades" in html
    assert "unique_clades" in html


def test_render_tree_svg_writes_static_tree_image(tmp_path: Path) -> None:
    output = tmp_path / "tree.svg"
    result = render_tree_svg(fixture("example_tree.nwk"), out_path=output)
    svg = output.read_text(encoding="utf-8")
    assert result.output_path == output
    assert result.format == "svg"
    assert result.layout == "cladogram"
    assert result.has_scale_bar is False
    assert "<svg" in svg
    assert "A" in svg and "D" in svg


def test_render_tree_svg_can_render_phylogram_with_scale_bar(tmp_path: Path) -> None:
    output = tmp_path / "phylogram.svg"
    result = render_tree_svg(
        fixture("example_tree.nwk"), out_path=output, layout="phylogram"
    )
    svg = output.read_text(encoding="utf-8")
    assert result.layout == "phylogram"
    assert result.has_scale_bar is True
    assert result.scale_bar_length == 0.1
    assert result.max_branch_distance == 0.30000000000000004
    assert 'class="scale-bar"' in svg
    assert 'class="scale-label"' in svg


def test_render_tree_svg_can_render_circular_layout(tmp_path: Path) -> None:
    output = tmp_path / "circular.svg"
    result = render_tree_svg(
        fixture("example_tree.nwk"), out_path=output, layout="circular"
    )
    svg = output.read_text(encoding="utf-8")
    assert result.layout == "circular"
    assert result.has_scale_bar is False
    assert '<path d="M ' in svg
    assert "text-anchor=" in svg


def test_render_tree_svg_can_render_internal_support_values(tmp_path: Path) -> None:
    output = tmp_path / "support.svg"
    result = render_tree_svg(
        fixture("example_tree_support_left.nwk"),
        out_path=output,
        layout="phylogram",
        show_support_values=True,
    )
    svg = output.read_text(encoding="utf-8")
    assert result.rendered_support_count == 2
    assert 'class="support-label"' in svg
    assert ">95<" in svg
    assert ">88<" in svg


def test_render_tree_svg_can_render_branch_color_overlays(tmp_path: Path) -> None:
    output = tmp_path / "branch-colors.svg"
    result = render_tree_svg(
        fixture("example_tree.nwk"),
        out_path=output,
        branch_colors={
            "A": "#dc2626",
            "C|D": "#0f766e",
        },
    )
    svg = output.read_text(encoding="utf-8")
    assert result.rendered_branch_color_count == 2
    assert 'class="branch branch-colored"' in svg
    assert "#dc2626" in svg
    assert "#0f766e" in svg


def test_render_tree_svg_can_render_internal_pie_markers(tmp_path: Path) -> None:
    output = tmp_path / "internal-pies.svg"
    result = render_tree_svg(
        fixture("example_tree.nwk"),
        out_path=output,
        internal_pies={
            "A|B": {"forest": 1.0},
            "C|D": {"forest": 0.25, "desert": 0.75},
        },
        internal_pie_colors={
            "forest": "#0f766e",
            "desert": "#c2410c",
        },
    )
    svg = output.read_text(encoding="utf-8")
    assert result.rendered_internal_pie_count == 2
    assert 'class="internal-pie-slice"' in svg
    assert "#0f766e" in svg
    assert "#c2410c" in svg


def test_render_tree_svg_can_render_categorical_tip_traits(tmp_path: Path) -> None:
    output = tmp_path / "categorical.svg"
    result = render_tree_svg(
        fixture("example_tree.nwk"),
        out_path=output,
        categorical_traits={
            "A": "Sweden",
            "B": "Norway",
            "C": "Denmark",
            "D": "Finland",
        },
    )
    svg = output.read_text(encoding="utf-8")
    assert result.rendered_categorical_trait_count == 4
    assert "categorical trait" in svg
    assert "<circle" in svg
    assert "Sweden" in svg


def test_render_tree_svg_can_render_continuous_tip_traits(tmp_path: Path) -> None:
    output = tmp_path / "continuous.svg"
    result = render_tree_svg(
        fixture("example_tree.nwk"),
        out_path=output,
        continuous_traits={"A": 1.2, "B": 1.4, "C": 1.8, "D": 2.0},
    )
    svg = output.read_text(encoding="utf-8")
    assert result.rendered_continuous_trait_count == 4
    assert "continuous trait" in svg
    assert "continuous-trait-gradient" in svg
    assert 'class="trait-bar-fill"' in svg


def test_render_tree_svg_can_render_collapsed_clades(tmp_path: Path) -> None:
    output = tmp_path / "collapsed.svg"
    result = render_tree_svg(
        fixture("example_tree_named_clades.nwk"),
        out_path=output,
        collapsed_clades=["Mammals"],
    )
    svg = output.read_text(encoding="utf-8")
    assert result.collapsed_clade_count == 1
    assert "Mammals (2 tips)" in svg
    assert 'class="collapsed-clade"' in svg
    assert ">A<" not in svg
    assert ">B<" not in svg


def test_render_tree_svg_can_render_metadata_strips(tmp_path: Path) -> None:
    output = tmp_path / "metadata-strips.svg"
    result = render_tree_svg(
        fixture("example_tree.nwk"),
        out_path=output,
        metadata_strips=[
            AnnotationStrip(
                "species",
                {
                    "A": "Alpha species",
                    "B": "Beta species",
                    "C": "Gamma species",
                    "D": "Delta species",
                },
            ),
            AnnotationStrip(
                "location",
                {"A": "Sweden", "B": "Norway", "C": "Denmark", "D": "Finland"},
            ),
        ],
    )
    svg = output.read_text(encoding="utf-8")
    assert result.rendered_metadata_strip_count == 2
    assert 'class="metadata-strip-cell"' in svg
    assert "species" in svg
    assert "location" in svg


def test_render_tree_svg_can_render_trait_heatmap_columns(tmp_path: Path) -> None:
    output = tmp_path / "heatmap.svg"
    result = render_tree_svg(
        fixture("example_tree.nwk"),
        out_path=output,
        heatmap_columns=[
            AnnotationStrip(
                "temperature", {"A": "1.2", "B": "1.4", "C": "1.8", "D": "2.0"}
            ),
            AnnotationStrip(
                "status", {"A": "high", "B": "medium", "C": "medium", "D": "low"}
            ),
        ],
    )
    svg = output.read_text(encoding="utf-8")
    assert result.rendered_heatmap_column_count == 2
    assert 'class="heatmap-cell"' in svg
    assert "temperature" in svg
    assert "status" in svg


def test_build_tree_figure_package_writes_publication_bundle(tmp_path: Path) -> None:
    output_dir = tmp_path / "tree-figure"
    result = build_tree_figure_package(
        fixture("example_tree_support_left.nwk"),
        out_dir=output_dir,
        show_support_values=True,
        categorical_traits={
            "A": "Sweden",
            "B": "Norway",
            "C": "Denmark",
            "D": "Finland",
        },
        metadata_strips=[
            AnnotationStrip(
                "location",
                {"A": "Sweden", "B": "Norway", "C": "Denmark", "D": "Finland"},
            )
        ],
        heatmap_columns=[
            AnnotationStrip(
                "temperature", {"A": "1.2", "B": "1.4", "C": "1.8", "D": "2.0"}
            )
        ],
    )
    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert result.figure_path.exists()
    assert result.caption_path.exists()
    assert result.annotations_path.exists()
    assert manifest["layout"] == "phylogram"
    assert manifest["render"]["rendered_support_count"] == 2
    assert manifest["audit"]["scale_bar_valid"] is True
    assert manifest["audit"]["support_audit"]["validated"] is True
    assert (
        "taxon\tlabel\tcategorical_trait\tcontinuous_trait\tlocation\ttemperature"
        in result.annotations_path.read_text(encoding="utf-8")
    )
    caption = result.caption_path.read_text(encoding="utf-8")
    assert "## Reviewer Summary" in caption
    assert "## Legend" in caption
    assert "## Limitations" in caption


def test_build_tree_figure_package_records_collapsed_clade_and_annotation_audits(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "tree-figure"
    result = build_tree_figure_package(
        fixture("example_tree_named_clades.nwk"),
        out_dir=output_dir,
        labels={"A": "Alpha", "B": "Beta", "C": "Gamma", "D": "Delta"},
        layout="phylogram",
        metadata_strips=[
            AnnotationStrip(
                "location",
                {"A": "Sweden", "B": "Sweden", "C": "Denmark", "D": "Finland"},
            )
        ],
        collapsed_clades=["Mammals"],
    )
    assert result.audit.collapsed_clades[0].clade_name == "Mammals"
    assert result.audit.collapsed_clades[0].descendant_count == 2
    assert result.audit.collapsed_clades[0].metadata_summaries == ["location: Sweden=2"]
    assert result.audit.annotation_coverage[0].aligned is False
    assert result.audit.annotation_coverage[0].missing_taxa == ["Mammals"]
    assert result.audit.table_consistency.consistent is True


def test_build_tree_figure_package_withholds_unvalidated_support_labels(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "tree-figure"
    result = build_tree_figure_package(
        fixture("example_tree_support_invalid.nwk"),
        out_dir=output_dir,
        layout="phylogram",
        show_support_values=True,
    )
    assert result.render.rendered_support_count == 0
    assert result.audit.support_audit.validated is False
    assert "support labels were withheld" in result.audit.reviewer_summary[1]


def test_render_tree_svg_can_use_metadata_labels(tmp_path: Path) -> None:
    output = tmp_path / "annotated.svg"
    result = render_tree_svg(
        fixture("example_tree.nwk"),
        out_path=output,
        labels={"A": "Alpha species", "B": "Beta species", "C": "Gamma species"},
    )
    svg = output.read_text(encoding="utf-8")
    assert "Alpha species" in svg
    assert "Beta species" in svg
    assert result.missing_metadata_labels == ["D"]


def test_diagnose_tree_path_combines_inspection_and_validation() -> None:
    report = diagnose_tree_path(fixture("example_tree.nwk"))
    assert report.inspection.tip_count == 4
    assert report.validation.tip_count == 4


def test_compute_root_to_tip_distances_reports_one_row_per_tip() -> None:
    report = compute_root_to_tip_distances(fixture("example_tree.nwk"))
    assert [(row.tip, row.distance) for row in report.distances] == [
        ("A", 0.30000000000000004),
        ("B", 0.30000000000000004),
        ("C", 0.30000000000000004),
        ("D", 0.30000000000000004),
    ]


def test_diagnose_ultrametricity_reports_max_deviation() -> None:
    ultrametric = diagnose_ultrametricity(fixture("example_tree.nwk"), tolerance=1e-6)
    non_ultrametric = diagnose_ultrametricity(
        fixture("example_tree_ladderized.nwk"), tolerance=1e-6
    )
    assert ultrametric.ultrametric is True
    assert ultrametric.max_deviation == 0.0
    assert non_ultrametric.ultrametric is False
    assert non_ultrametric.max_deviation == 0.2


def test_diagnose_ultrametricity_uses_shared_default_tolerance() -> None:
    near_default = diagnose_ultrametricity(fixture("example_tree_near_ultrametric.nwk"))
    near_tight = diagnose_ultrametricity(
        fixture("example_tree_near_ultrametric.nwk"),
        tolerance=1e-12,
    )

    assert near_default.ultrametric is True
    assert near_default.max_deviation == pytest.approx(1e-9, abs=1e-15)
    assert near_tight.ultrametric is False


def test_annotate_tree_against_table_finds_missing_and_extra_taxa() -> None:
    report = annotate_tree_against_table(
        fixture("example_tree.nwk"), fixture("example_traits.tsv")
    )
    assert report.linked_taxa == 3
    assert report.annotated_taxa == ["A", "B", "C"]
    assert report.missing_from_table == ["D"]
    assert report.extra_table_entries == ["E"]
    assert [(row.taxon, row.matched) for row in report.joined_rows] == [
        ("A", True),
        ("B", True),
        ("C", True),
        ("D", False),
    ]


def test_cli_metadata_inspect_json_output(capsys) -> None:
    exit_code = main(
        ["metadata", "inspect", str(fixture("example_metadata.tsv")), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "metadata"
    assert payload["data"]["taxon_column"] == "taxon"
    assert payload["metrics"]["taxon_count"] == 4


def test_cli_annotate_can_write_joined_tip_rows(tmp_path: Path, capsys) -> None:
    report_path = tmp_path / "annotation.json"
    joined_path = tmp_path / "annotation.tsv"
    exit_code = main(
        [
            "annotate",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_metadata.tsv")),
            "--out",
            str(report_path),
            "--joined-out",
            str(joined_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["data"]["joined_rows"][3]["taxon"] == "D"
    assert joined_path.read_text(encoding="utf-8") == (
        "taxon\tmatched\tspecies\tlocation\n"
        "A\ttrue\tAlpha species\tSweden\n"
        "B\ttrue\tBeta species\tNorway\n"
        "C\ttrue\tGamma species\tDenmark\n"
        "D\ttrue\tDelta species\tFinland\n"
    )


def test_cli_env_inspect_json_output(capsys) -> None:
    exit_code = main(["env", "inspect", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "env"
    assert payload["metrics"]["dependency_count"] >= 1


def test_cli_annotate_writes_annotation_json(tmp_path: Path, capsys) -> None:
    output = tmp_path / "tree.annotated.json"
    exit_code = main(
        [
            "annotate",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_traits.tsv")),
            "--taxon-column",
            "taxon",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    artifact = json.loads(output.read_text(encoding="utf-8"))
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["outputs"] == [str(output)]
    assert artifact["annotated_taxa"] == ["A", "B", "C"]


def test_cli_traits_validate_json_output(capsys) -> None:
    exit_code = main(
        ["traits", "validate", str(fixture("example_traits_validate.tsv")), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "traits"
    assert payload["metrics"]["trait_column_count"] == 3
    assert payload["data"]["trait_columns"][0]["kind"] == "numeric"


def test_cli_traits_link_json_output(capsys) -> None:
    exit_code = main(
        [
            "traits",
            "link",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits.tsv")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["linked_taxa"] == 3
    assert payload["data"]["missing_from_traits"] == ["D"]


def test_cli_traits_link_strict_mode_returns_typed_error(capsys) -> None:
    exit_code = main(
        [
            "traits",
            "link",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits.tsv")),
            "--strict",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == MetadataJoinError.code
    assert payload["errors"][0]["details"]["failure_reason"] == "tree_trait_taxon_mismatch"
    assert payload["errors"][0]["details"]["evidence"]["missing_from_traits"] == ["D"]
    assert payload["errors"][0]["details"]["evidence"]["extra_trait_taxa"] == ["E"]


def test_cli_prune_writes_tree_and_pruned_taxa_report(tmp_path: Path, capsys) -> None:
    output = tmp_path / "tree.pruned.nwk"
    exit_code = main(
        [
            "prune",
            str(fixture("example_tree.nwk")),
            "--keep-from",
            str(fixture("example_traits.tsv")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert output.read_text(encoding="utf-8").strip() == "((A:0.1,B:0.1):0.2,C:0.3);"
    assert payload["data"]["removed_taxa"] == ["D"]
    assert (tmp_path / "pruned_taxa.tsv").read_text(encoding="utf-8") == "taxon\nD\n"


def test_cli_alignment_inspect_json_output(capsys) -> None:
    exit_code = main(
        ["alignment", "inspect", str(fixture("example_alignment.fasta")), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "alignment"
    assert payload["metrics"]["alignment_length"] == 8
    assert payload["metrics"]["alphabet"] == "dna"
    assert payload["metrics"]["duplicate_group_count"] == 0
    assert payload["data"]["variable_site_count"] == 2


def test_cli_alignment_link_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "link",
            str(fixture("example_tree.nwk")),
            str(fixture("example_alignment.fasta")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["linked_taxa"] == 4


def test_cli_alignment_link_strict_mode_returns_typed_error(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "link",
            str(fixture("example_tree.nwk")),
            str(fixture("example_alignment_mismatch.fasta")),
            "--strict",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == AlignmentTaxonMismatchError.code


def test_cli_alignment_quality_json_output_surfaces_suspicion_metrics(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "quality",
            str(fixture("example_alignment_duplicates.fasta")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["duplicate_group_count"] == 1
    assert payload["metrics"]["sequence_length_outlier_count"] == 0
    assert payload["metrics"]["near_duplicate_count"] == 0
    assert payload["warnings"] == ["alignment contains identical duplicate sequences"]


def test_cli_alignment_classify_json_output(capsys) -> None:
    exit_code = main(
        ["alignment", "classify", str(fixture("example_sequences_raw.fasta")), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["state"] == "raw_sequence_fasta"
    assert payload["data"]["sequence_count"] == 4


def test_cli_alignment_sequence_type_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "sequence-type",
            str(fixture("example_sequences_raw.fasta")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["detected_type"] == "dna"
    assert payload["metrics"]["selected_type"] == "dna"
    assert payload["metrics"]["confidence"] == "medium"
    assert payload["data"]["compatible_types"] == ["dna", "protein"]


def test_cli_alignment_validate_input_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "validate-input",
            str(fixture("example_sequences_invalid_input.fasta")),
            "--sequence-type",
            "dna",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["duplicate_identifier_count"] == 1
    assert payload["metrics"]["illegal_character_count"] == 1
    assert payload["metrics"]["empty_sequence_count"] == 1
    assert payload["metrics"]["sequence_length_outlier_count"] == 2
    assert payload["metrics"]["detected_type"] == "invalid"
    assert payload["metrics"]["selected_type"] is None
    assert payload["metrics"]["sequence_type_confidence"] == "blocked"


def test_cli_alignment_repair_input_writes_repaired_fasta(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "repaired.fasta"
    exit_code = main(
        [
            "alignment",
            "repair-input",
            str(fixture("example_sequences_invalid_input.fasta")),
            "--out",
            str(output),
            "--sequence-type",
            "dna",
            "--normalize-identifiers",
            "--remove-invalid-records",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output.read_text(encoding="utf-8") == (
        ">Alpha_sample\nACTGACTG\n>rare_taxon\nACTGACTGACTGACTGACTGACTG\n"
    )
    assert payload["metrics"]["normalized_identifier_count"] == 2
    assert payload["metrics"]["removed_record_count"] == 2
    assert payload["metrics"]["remaining_duplicate_identifier_count"] == 0
    assert payload["metrics"]["remaining_illegal_character_count"] == 0
    assert payload["metrics"]["remaining_empty_sequence_count"] == 0


def test_cli_alignment_windows_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "windows",
            str(fixture("example_alignment_regions.fasta")),
            "--window-size",
            "4",
            "--step-size",
            "4",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["window_count"] == 2
    assert payload["metrics"]["over_aligned_region_count"] == 1
    assert payload["metrics"]["under_aligned_region_count"] == 1


def test_cli_alignment_readiness_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "readiness",
            str(fixture("example_alignment_coding.fasta")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["ready_method_count"] == 3
    assert payload["metrics"]["blocked_method_count"] == 2
    assert any(
        method["analysis"] == "coding" and method["ready"] is False
        for method in payload["data"]["methods"]
    )


def test_cli_alignment_profiles_json_output(capsys) -> None:
    exit_code = main(["alignment", "profiles", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["profile_count"] == 5
    assert any(profile["name"] == "coding-safe" for profile in payload["data"])


def test_cli_alignment_occupancy_writes_tables_and_filtered_outputs(
    tmp_path: Path, capsys
) -> None:
    taxa_out = tmp_path / "taxa.tsv"
    loci_out = tmp_path / "loci.tsv"
    matrix_out = tmp_path / "matrix.tsv"
    filtered_alignment_out = tmp_path / "filtered_alignment.fasta"
    filtered_partitions_out = tmp_path / "filtered_partitions.txt"

    exit_code = main(
        [
            "alignment",
            "occupancy",
            str(fixture("example_multilocus_alignment.fasta")),
            str(fixture("example_multilocus_partitions.txt")),
            "--taxon-coverage-threshold",
            "0.6",
            "--locus-coverage-threshold",
            "0.6",
            "--taxa-out",
            str(taxa_out),
            "--loci-out",
            str(loci_out),
            "--matrix-out",
            str(matrix_out),
            "--filtered-alignment-out",
            str(filtered_alignment_out),
            "--filtered-partitions-out",
            str(filtered_partitions_out),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["taxon_count"] == 5
    assert payload["metrics"]["locus_count"] == 3
    assert payload["metrics"]["minimum_locus_occupancy"] == 0.0
    assert payload["metrics"]["filtered_taxon_count"] == 2
    assert payload["metrics"]["filtered_locus_count"] == 2
    assert payload["data"]["filter_report"]["removed_taxa"] == [
        "TaxonC",
        "TaxonD",
        "TaxonE",
    ]
    assert payload["data"]["filter_report"]["removed_loci"] == ["gene_beta"]

    assert taxa_out.read_text(encoding="utf-8").startswith(
        "taxon\tcovered_locus_count\ttotal_locus_count\t"
        "locus_coverage_fraction\tobserved_site_count\t"
        "total_site_count\tsite_coverage_fraction\tlow_coverage"
    )
    assert "gene_gamma\t2\t5\t0.4\t6\t15\t0.4\ttrue\n" in loci_out.read_text(
        encoding="utf-8"
    )
    assert matrix_out.read_text(encoding="utf-8").splitlines()[0] == (
        "taxon\tgene_alpha\tgene_beta\tgene_gamma\tcovered_locus_count\t"
        "total_locus_count\tlocus_coverage_fraction\tobserved_site_count\t"
        "total_site_count\tsite_coverage_fraction\tlow_coverage"
    )
    assert filtered_alignment_out.read_text(encoding="utf-8") == (
        ">TaxonA\nAAAAGGG\n>TaxonB\nAAAAGGG\n"
    )
    assert filtered_partitions_out.read_text(encoding="utf-8") == (
        "DNA,gene_alpha = 1-4\nDNA,gene_gamma = 5-7\n"
    )


def test_cli_alignment_occupancy_supports_minimum_locus_occupancy(
    tmp_path: Path, capsys
) -> None:
    taxa_out = tmp_path / "partial-taxa.tsv"
    loci_out = tmp_path / "partial-loci.tsv"
    matrix_out = tmp_path / "partial-matrix.tsv"
    filtered_alignment_out = tmp_path / "partial-filtered.fasta"
    filtered_partitions_out = tmp_path / "partial-filtered.partitions.txt"

    exit_code = main(
        [
            "alignment",
            "occupancy",
            str(fixture("example_multilocus_partial_occupancy.fasta")),
            str(fixture("example_multilocus_partial_occupancy_partitions.txt")),
            "--taxon-coverage-threshold",
            "0.5",
            "--locus-coverage-threshold",
            "0.75",
            "--minimum-locus-occupancy",
            "0.75",
            "--taxa-out",
            str(taxa_out),
            "--loci-out",
            str(loci_out),
            "--matrix-out",
            str(matrix_out),
            "--filtered-alignment-out",
            str(filtered_alignment_out),
            "--filtered-partitions-out",
            str(filtered_partitions_out),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["minimum_locus_occupancy"] == 0.75
    assert payload["metrics"]["filtered_taxon_count"] == 2
    assert payload["metrics"]["filtered_locus_count"] == 2
    assert payload["data"]["report"]["minimum_locus_occupancy"] == 0.75
    assert payload["data"]["report"]["low_coverage_loci"] == [
        "gene_alpha",
        "gene_gamma",
    ]
    assert payload["data"]["filter_report"]["retained_taxa"] == ["TaxonA", "TaxonB"]
    assert payload["data"]["filter_report"]["retained_loci"] == [
        "gene_alpha",
        "gene_beta",
    ]
    assert payload["data"]["filter_report"]["filter_iterations"][0] == {
        "input_loci": ["gene_alpha", "gene_beta", "gene_gamma"],
        "input_taxa": ["TaxonA", "TaxonB", "TaxonC", "TaxonD", "TaxonE"],
        "iteration": 1,
        "low_coverage_loci": ["gene_gamma"],
        "low_coverage_taxa": ["TaxonC", "TaxonD", "TaxonE"],
        "removed_loci": ["gene_gamma"],
        "removed_taxa": ["TaxonC", "TaxonD", "TaxonE"],
        "retained_loci": ["gene_alpha", "gene_beta"],
        "retained_taxa": ["TaxonA", "TaxonB"],
    }
    assert filtered_alignment_out.read_text(encoding="utf-8") == (
        ">TaxonA\nAAAACCCC\n>TaxonB\nAAAACCCC\n"
    )
    assert filtered_partitions_out.read_text(encoding="utf-8") == (
        "DNA,gene_alpha = 1-4\nDNA,gene_beta = 5-8\n"
    )
    assert "TaxonB\t2\t3\t0.666666666667\t9\t12\t0.75\tfalse\n" in taxa_out.read_text(
        encoding="utf-8"
    )
    assert "gene_gamma\t2\t5\t0.4\t9\t20\t0.45\ttrue\n" in loci_out.read_text(
        encoding="utf-8"
    )
    assert matrix_out.read_text(encoding="utf-8").splitlines()[0] == (
        "taxon\tgene_alpha\tgene_beta\tgene_gamma\tcovered_locus_count\t"
        "total_locus_count\tlocus_coverage_fraction\tobserved_site_count\t"
        "total_site_count\tsite_coverage_fraction\tlow_coverage"
    )


def test_cli_alignment_forensic_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "forensic",
            str(fixture("example_alignment_coding.fasta")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["safe_for_coding_analysis"] is False
    assert payload["data"]["coding"]["mixed_coding_signals"] is True


def test_cli_alignment_quality_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "quality",
            str(fixture("example_alignment_missingness.fasta")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["quality_score"] == payload["data"]["quality_score"]
    assert payload["metrics"]["invariant_site_count"] == 6
    assert payload["metrics"]["parsimony_informative_site_count"] == 0
    assert payload["metrics"]["suspicious_alignment"] is True
    assert payload["metrics"]["suspicious_reason_count"] >= 2
    assert payload["metrics"]["concentrated_column_count"] == 2
    assert (
        payload["data"]["missing_data_concentration"]["longest_concentrated_run"] == 2
    )


def test_cli_alignment_low_information_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "low-information",
            str(fixture("example_alignment_missingness.fasta")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["low_information"] is True
    assert payload["metrics"]["parsimony_informative_site_count"] == 0


def test_cli_alignment_duplicate_policy_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "duplicate-policy",
            str(fixture("example_alignment_duplicates.fasta")),
            "--identity-threshold",
            "0.875",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["exact_duplicate_group_count"] == 1
    assert payload["metrics"]["policy_action_count"] >= 2


def test_cli_alignment_ambiguous_columns_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "ambiguous-columns",
            str(fixture("example_alignment_site_missingness.fasta")),
            "--threshold",
            "0.5",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["ambiguous_column_count"] == 4


def test_cli_alignment_sequence_ranking_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "sequence-ranking",
            str(fixture("example_alignment_ambiguity.fasta")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["sequence_count"] == 3
    assert payload["data"]["rows"][0]["identifier"] == "A"


def test_cli_alignment_filter_json_output(tmp_path: Path, capsys) -> None:
    output = tmp_path / "cleaned.fasta"
    exit_code = main(
        [
            "alignment",
            "filter",
            str(fixture("example_alignment_filtering.fasta")),
            "--profile",
            "moderate",
            "--group-table",
            str(fixture("example_alignment_groups.tsv")),
            "--group-columns",
            "region",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output.read_text(encoding="utf-8") == fixture(
        "example_alignment_filtering_cleaned_moderate.fasta"
    ).read_text(encoding="utf-8")
    assert payload["metrics"]["profile"] == "moderate"
    assert payload["metrics"]["trimmed_alignment_length"] == 8


def test_cli_alignment_compare_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "compare",
            str(fixture("example_alignment_filtering.fasta")),
            str(fixture("example_alignment_filtering_cleaned_moderate.fasta")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["shared_taxa"] == 3
    assert payload["metrics"]["left_alignment_length"] == 12
    assert payload["metrics"]["right_alignment_length"] == 8


def test_cli_report_dataset_json_output_includes_audit(capsys, tmp_path: Path) -> None:
    output = tmp_path / "dataset-report.html"
    exit_code = main(
        [
            "report",
            "dataset",
            "--tree",
            str(fixture("example_tree_named_clades.nwk")),
            "--metadata",
            str(fixture("example_metadata.tsv")),
            "--traits",
            str(fixture("example_traits_validate.tsv")),
            "--alignment",
            str(fixture("example_alignment.fasta")),
            "--tip-dates",
            str(fixture("example_tip_dates.tsv")),
            "--calibrations",
            str(fixture("example_calibrations.tsv")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["readiness_decision"] == "ready_with_warnings"
    assert payload["metrics"]["excluded_taxa"] == 0
    assert payload["metrics"]["risky_analysis_count"] >= 1
    assert payload["data"]["dataset_audit"]["alignment_forensic"] is not None


def test_cli_alignment_length_outliers_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "length-outliers",
            str(fixture("example_sequences_raw.fasta")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["sequence_length_outlier_count"] == 2
    assert [row["identifier"] for row in payload["data"]] == ["B", "C"]


def test_cli_alignment_composition_json_output(capsys) -> None:
    exit_code = main(
        ["alignment", "composition", str(fixture("example_alignment.fasta")), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["alphabet"] == "dna"
    assert payload["data"]["whole_alignment_gc_content"] == 0.5
    assert payload["data"]["nucleotide_composition"] == {
        "A": 0.3125,
        "C": 0.25,
        "G": 0.25,
        "T": 0.1875,
    }


def test_cli_alignment_alphabet_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "alphabet",
            str(fixture("example_alignment_protein.fasta")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["alphabet"] == "protein"
    assert payload["data"]["inferred_alphabet"] == "protein"


def test_cli_alignment_gc_json_output(capsys) -> None:
    exit_code = main(
        ["alignment", "gc", str(fixture("example_alignment.fasta")), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["data"]["whole_alignment_gc_content"] == 0.5
    assert payload["data"]["per_sequence_gc_content"][1] == {
        "gc_fraction": 0.375,
        "identifier": "B",
    }


def test_cli_alignment_invalid_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "invalid",
            str(fixture("example_alignment_invalid_dna.fasta")),
            "--alphabet",
            "dna",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["invalid_character_count"] == 1
    assert payload["data"] == [{"character": "Z", "identifier": "A", "position": 5}]


def test_cli_alignment_duplicates_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "duplicates",
            str(fixture("example_alignment_duplicates.fasta")),
            "--identity-threshold",
            "0.875",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["duplicate_group_count"] == 1
    assert payload["metrics"]["near_duplicate_count"] == 5
    assert payload["data"]["duplicate_sequence_groups"][0]["identifiers"] == ["A", "B"]


def test_cli_alignment_outliers_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "outliers",
            str(fixture("example_alignment_gc_outlier.fasta")),
            "--deviation-threshold",
            "0.2",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["composition_outlier_count"] == 1
    assert payload["data"] == [
        {"deviation": 1.0, "identifier": "C", "robust_z_score": None}
    ]


def test_cli_compare_support_json_output(capsys) -> None:
    exit_code = main(
        [
            "compare",
            "support",
            str(fixture("example_tree_support_left.nwk")),
            str(fixture("example_tree_support_right.nwk")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["shared_clades"] == 2
    assert payload["data"]["shared_clades"][0]["split_id"] == "A|B"


def test_cli_compare_branch_lengths_json_output(capsys) -> None:
    exit_code = main(
        [
            "compare",
            "branch-lengths",
            str(fixture("example_tree.nwk")),
            str(fixture("example_tree_branch_lengths_right.nwk")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["shared_splits"] == 2
    assert payload["metrics"]["branch_score_distance"] == 0.33166247903554
    assert payload["metrics"]["same_taxon_set"] is True
    assert payload["data"]["branch_score"]["split_count"] == 5


def test_cli_compare_branch_lengths_reports_pruned_taxon_set_status(capsys) -> None:
    exit_code = main(
        [
            "compare",
            "branch-lengths",
            str(fixture("example_tree.nwk")),
            str(fixture("example_tree_overlap.nwk")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["same_taxon_set"] is False
    assert payload["data"]["left_only_taxa"] == ["D"]
    assert payload["data"]["right_only_taxa"] == ["E"]
    assert payload["data"]["branch_score"]["branch_score_distance"] == 0.0


def test_cli_compare_branch_lengths_requires_identical_taxa_when_requested() -> None:
    try:
        main(
            [
                "compare",
                "branch-lengths",
                str(fixture("example_tree.nwk")),
                str(fixture("example_tree_overlap.nwk")),
                "--taxon-overlap-policy",
                "require-identical",
            ]
        )
    except SystemExit as error:
        assert error.code == 2
    else:  # pragma: no cover - defensive assertion
        raise AssertionError(
            "expected CLI branch-length comparison to reject mismatched taxa"
        )


def test_cli_compare_clades_json_output(capsys) -> None:
    exit_code = main(
        [
            "compare",
            "clades",
            str(fixture("example_tree.nwk")),
            str(fixture("example_tree_alt.nwk")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["tree_count"] == 2
    assert payload["data"]["shared_clades"] == ["A|B"]
    assert payload["data"]["conflicting_clades"] == ["C|D", "A|B|C"]
    assert [
        (row["tree_path"].split("/")[-1], row["unique_clades"])
        for row in payload["data"]["tree_summaries"]
    ] == [
        ("example_tree.nwk", ["C|D"]),
        ("example_tree_alt.nwk", ["A|B|C"]),
    ]


def test_cli_compare_clades_supports_multiple_trees_and_table_output(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "clade-overlap.tsv"
    exit_code = main(
        [
            "compare",
            "clades",
            str(fixture("example_tree.nwk")),
            str(fixture("example_tree_support_left.nwk")),
            "--tree",
            str(fixture("example_tree_alt.nwk")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["tree_count"] == 3
    assert payload["metrics"]["shared_clades"] == 1
    assert payload["metrics"]["conflicting_clades"] == 2
    assert output.exists()
    assert "clade_id\ttree_path\tpresent\tsupport" in output.read_text(encoding="utf-8")


def test_cli_compare_json_output_supports_unrooted_rf_mode(capsys) -> None:
    exit_code = main(
        [
            "compare",
            str(fixture("example_tree.nwk")),
            str(fixture("example_tree_rooting_diff.nwk")),
            "--rf-mode",
            "unrooted",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["robinson_foulds_distance"] == 0
    assert payload["metrics"]["rf_mode"] == "unrooted"
    assert payload["data"]["rooted_robinson_foulds_distance"] == 2
    assert payload["data"]["same_unrooted_topology"] is True


def test_cli_compare_writes_topology_distance_split_table(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "topology-distance.tsv"

    exit_code = main(
        [
            "compare",
            str(fixture("example_tree.nwk")),
            str(fixture("example_tree_alt.nwk")),
            "--split-table-out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert output_path.exists()
    assert payload["outputs"] == [str(output_path)]
    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "split_id\tsplit_kind\tcomparison_status\ttaxon_count\tdescendant_taxa\tleft_present\tright_present",
        "A|B\tclade\tshared\t2\tA|B\ttrue\ttrue",
        "C|D\tclade\tleft_only\t2\tC|D\ttrue\tfalse",
        "A|B|C\tclade\tright_only\t3\tA|B|C\tfalse\ttrue",
    ]


def test_cli_compare_requires_identical_taxa_when_policy_requests_it() -> None:
    try:
        main(
            [
                "compare",
                str(fixture("example_tree.nwk")),
                str(fixture("example_tree_overlap.nwk")),
                "--taxon-overlap-policy",
                "require-identical",
            ]
        )
    except SystemExit as error:
        assert error.code == 2
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected CLI compare to reject mismatched taxa")


def test_cli_compare_changes_json_output(capsys) -> None:
    exit_code = main(
        [
            "compare",
            "changes",
            str(fixture("example_tree.nwk")),
            str(fixture("example_tree_alt.nwk")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["data"]["lost_clades"] == ["C|D"]
    assert payload["data"]["gained_clades"] == ["A|B|C"]


def test_cli_compare_prune_writes_shared_taxon_trees(tmp_path: Path, capsys) -> None:
    output_dir = tmp_path / "shared"
    exit_code = main(
        [
            "compare",
            "prune",
            str(fixture("example_tree.nwk")),
            str(fixture("example_tree_overlap.nwk")),
            "--out",
            str(output_dir),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["data"]["shared_taxa"] == ["A", "B", "C"]
    assert (output_dir / "left-shared.nwk").read_text(
        encoding="utf-8"
    ) == "((A:0.1,B:0.1):0.2,C:0.3);\n"
    assert (output_dir / "right-shared.nwk").read_text(
        encoding="utf-8"
    ) == "((A:0.1,B:0.1):0.2,C:0.3);\n"


def test_cli_compare_table_writes_tsv_output(tmp_path: Path, capsys) -> None:
    output = tmp_path / "comparison.tsv"
    exit_code = main(
        [
            "compare",
            "table",
            str(fixture("example_tree.nwk")),
            str(fixture("example_tree_alt.nwk")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["table_rows"] == 7
    assert output.read_text(encoding="utf-8").startswith(
        "split_id\tcomparison_status\tshared_clade\t"
    )
    assert payload["data"]["table_path"] == str(output)


def test_cli_compare_report_json_output(tmp_path: Path, capsys) -> None:
    output = tmp_path / "compare.html"
    exit_code = main(
        [
            "compare",
            "report",
            str(fixture("example_tree_support_left.nwk")),
            str(fixture("example_tree_support_right.nwk")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["outputs"] == [str(output)]
    assert "Bijux Tree Comparison Report" in output.read_text(encoding="utf-8")


def test_cli_render_writes_svg_output(tmp_path: Path, capsys) -> None:
    output = tmp_path / "tree.svg"
    exit_code = main(
        ["render", str(fixture("example_tree.nwk")), "--out", str(output), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["outputs"] == [str(output)]
    assert "<svg" in output.read_text(encoding="utf-8")


def test_cli_render_with_metadata_labels_reports_missing_taxa(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "annotated.svg"
    exit_code = main(
        [
            "render",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_metadata.tsv")),
            "--label-column",
            "species",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["data"]["render"]["missing_metadata_labels"] == []
    assert "Alpha species" in output.read_text(encoding="utf-8")


def test_cli_render_with_partial_metadata_warns_for_missing_labels(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "partial.svg"
    exit_code = main(
        [
            "render",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_traits.tsv")),
            "--label-column",
            "value",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["warnings"] == ["D"]
    assert payload["data"]["render"]["missing_metadata_labels"] == ["D"]


def test_cli_render_can_build_annotated_figure_package(tmp_path: Path, capsys) -> None:
    output = tmp_path / "annotated.svg"
    package_dir = tmp_path / "figure-package"
    exit_code = main(
        [
            "render",
            str(fixture("example_tree_support_left.nwk")),
            "--metadata",
            str(fixture("example_metadata.tsv")),
            "--label-column",
            "species",
            "--metadata-strip-columns",
            "location",
            "--traits",
            str(fixture("example_traits_validate.tsv")),
            "--layout",
            "phylogram",
            "--support-labels",
            "--categorical-column",
            "habitat",
            "--continuous-column",
            "height_cm",
            "--heatmap-columns",
            "height_cm,status",
            "--package-dir",
            str(package_dir),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["rendered_support_count"] == 2
    assert payload["metrics"]["rendered_metadata_strip_count"] == 1
    assert payload["metrics"]["rendered_heatmap_column_count"] == 2
    assert payload["data"]["render"]["layout"] == "phylogram"
    assert payload["data"]["figure_package_dir"] == str(package_dir)
    assert payload["data"]["figure_package_audit"]["scale_bar_valid"] is True
    assert (package_dir / "figure.svg").exists()
    assert "Alpha species" in output.read_text(encoding="utf-8")


def test_render_phylogenetics_report_writes_html(tmp_path: Path) -> None:
    output = tmp_path / "report.html"
    result = render_phylogenetics_report(
        tree_path=fixture("example_tree.nwk"),
        alignment_path=fixture("example_alignment.fasta"),
        traits_path=fixture("example_traits.tsv"),
        metadata_path=fixture("example_traits.tsv"),
        out_path=output,
    )
    assert result.output_path == output
    assert output.exists()
    assert "Bijux Phylogenetics Report" in output.read_text(encoding="utf-8")


def test_render_tree_report_writes_embedded_manifest(tmp_path: Path) -> None:
    output = tmp_path / "tree-report.html"
    result = render_tree_report(tree_path=fixture("example_tree.nwk"), out_path=output)
    text = output.read_text(encoding="utf-8")
    assert result.report_kind == "tree"
    assert result.machine_manifest["report_kind"] == "tree"
    assert result.machine_manifest["input_paths"] == [str(fixture("example_tree.nwk"))]
    assert result.machine_manifest_path.exists()
    assert 'id="bijux-report-manifest"' in text
    assert "Bijux Tree Report" in text


def test_write_html_report_renders_summary_metrics_and_artifact_links(
    tmp_path: Path,
) -> None:
    output = tmp_path / "summary-report.html"
    write_html_report(
        title="Summary Report",
        sections=[("overview", "compact reviewer summary")],
        out_path=output,
        embedded_json={"report_kind": "summary"},
        summary_metrics=[("tree count", 12), ("report mode", "scaled-summary")],
        artifact_links=[
            ("clade-frequencies", "summary-report.artifacts/clade-frequencies.tsv", "128 bytes"),
            ("report-manifest", "summary-report.artifacts/summary.manifest.json", None),
        ],
    )

    text = output.read_text(encoding="utf-8")
    assert "<h2>summary</h2>" in text
    assert "tree count" in text
    assert "scaled-summary" in text
    assert 'href="summary-report.artifacts/clade-frequencies.tsv"' in text
    assert "128 bytes" in text
    assert 'href="summary-report.artifacts/summary.manifest.json"' in text


def test_render_tree_report_preserves_support_and_branch_diagnostics(
    tmp_path: Path,
) -> None:
    output = tmp_path / "tree-support-report.html"
    result = render_tree_report(
        tree_path=fixture("example_tree_support_invalid.nwk"), out_path=output
    )
    assert result.validation.missing_internal_branch_nodes == []
    assert result.inspection.suspicious_support_value_ranges == [
        "support value 101 at A|B|C|D exceeds 100",
        "support value 120 at A|B exceeds 100",
        "support value -5 at C|D is negative",
    ]
    assert [warning.code for warning in result.inspection.tree_quality_warnings] == [
        "suspicious_support_ranges"
    ]
    assert (
        output.read_text(encoding="utf-8").count("suspicious_support_value_ranges") >= 1
    )


def test_render_dataset_report_writes_metadata_sections(tmp_path: Path) -> None:
    output = tmp_path / "dataset-report.html"
    result = render_dataset_report(
        tree_path=fixture("example_tree.nwk"),
        metadata_path=fixture("example_metadata.tsv"),
        traits_path=fixture("example_traits.tsv"),
        out_path=output,
    )
    text = output.read_text(encoding="utf-8")
    assert result.report_kind == "dataset"
    assert result.metadata_linkage is not None
    assert result.dataset_readiness is not None
    assert result.dataset_audit is not None
    assert result.machine_manifest["sections"] == [
        "reviewer-summary",
        "tree-validation",
        "tree-inspection",
        "tree-forensic",
        "metadata-linkage",
        "traits-linkage",
        "trait-missing-values",
        "dataset-readiness",
        "dataset-audit",
        "dataset-findings",
        "dataset-analysis-decisions",
        "dataset-readiness-levels",
        "dataset-crosswalk",
        "dataset-completeness",
        "dataset-exclusions",
        "dataset-mismatches",
        "dataset-risk-score",
        "dataset-minimal-fix-plan",
        "dataset-reviewer-checklist",
        "dataset-ordering",
        "dataset-pruning",
        "dataset-group-imbalance",
        "dataset-input-ledger",
        "limitations",
    ]
    assert result.trait_missing_values is not None
    assert result.trait_missing_values.missing_values == []
    assert result.machine_manifest_path.exists()
    assert "Bijux Dataset Report" in text
    assert len(result.input_ledger) == 3
    assert result.input_ledger[0].role == "tree"


def test_render_phylo_inputs_report_writes_alignment_sections(tmp_path: Path) -> None:
    output = tmp_path / "phylo-inputs-report.html"
    result = render_phylo_inputs_report(
        tree_path=fixture("example_tree.nwk"),
        alignment_path=fixture("example_alignment_missingness.fasta"),
        out_path=output,
    )
    text = output.read_text(encoding="utf-8")
    assert result.report_kind == "phylo-inputs"
    assert result.alignment is not None
    assert result.alignment_quality is not None
    assert result.alignment_forensic is not None
    assert result.alignment_coding is not None
    assert result.alignment_identity_matrix is not None
    assert result.alignment_linkage is not None
    assert result.machine_manifest["sections"] == [
        "reviewer-summary",
        "tree-validation",
        "tree-inspection",
        "tree-forensic",
        "alignment-summary",
        "alignment-quality",
        "alignment-low-information",
        "alignment-duplicate-policy",
        "alignment-ambiguous-columns",
        "alignment-sequence-ranking",
        "alignment-forensic",
        "alignment-coding",
        "alignment-identity-matrix",
        "alignment-linkage",
        "limitations",
    ]
    assert result.machine_manifest_path.exists()
    assert "Bijux Phylo Inputs Report" in text
    assert "alignment suspicious diagnostics: flagged" in text
    assert result.alignment_low_information is not None
    assert result.alignment_duplicate_policy is not None
    assert result.alignment_ambiguous_columns is not None
    assert result.alignment_sequence_ranking is not None


def test_render_alignment_report_writes_alignment_sections(tmp_path: Path) -> None:
    output = tmp_path / "alignment-report.html"
    result = render_alignment_report(
        alignment_path=fixture("example_alignment_filtering.fasta"),
        out_path=output,
    )

    text = output.read_text(encoding="utf-8")
    assert result.report_kind == "alignment"
    assert result.machine_manifest["sections"] == [
        "reviewer-summary",
        "alignment-summary",
        "alignment-quality",
        "alignment-readiness",
        "alignment-low-information",
        "alignment-duplicate-policy",
        "alignment-ambiguous-columns",
        "alignment-sequence-ranking",
        "alignment-filter-profiles",
        "alignment-suspicious-windows",
        "alignment-forensic",
        "alignment-coding",
        "alignment-identity-matrix",
        "limitations",
    ]
    assert result.machine_manifest_path.exists()
    assert "Bijux Alignment Report" in text
    assert "alignment suspicious diagnostics: flagged" in text
    assert "longest concentrated missing-data run: 4" in text


def test_render_taxon_report_writes_taxon_sections(tmp_path: Path) -> None:
    output = tmp_path / "taxonomy-report.html"
    result = render_taxon_report(
        tree_path=fixture("example_taxonomy_tree.nwk"),
        synonym_table_path=fixture("example_taxon_synonyms.tsv"),
        out_path=output,
    )

    text = output.read_text(encoding="utf-8")
    assert result.report_kind == "taxonomy"
    assert result.machine_manifest["sections"] == [
        "reviewer-summary",
        "taxon-audit",
        "taxon-identity",
        "taxon-safety",
        "taxon-namespaces",
        "taxon-rank-consistency",
        "taxon-synonyms",
        "taxon-duplicate-identities",
        "taxon-mapping-conflicts",
        "taxon-accepted-names",
        "limitations",
    ]
    assert result.machine_manifest_path.exists()
    assert "Bijux Taxon Audit Report" in text


def test_render_taxon_report_includes_crosswalk_exclusions_loss_stability_and_input_ledger(
    tmp_path: Path,
) -> None:
    output = tmp_path / "taxonomy-report-workflow.html"
    result = render_taxon_report(
        tree_path=fixture("example_taxon_workflow_tree.nwk"),
        metadata_path=fixture("example_taxon_workflow_metadata.csv"),
        traits_path=fixture("example_taxon_workflow_traits.csv"),
        alignment_path=fixture("example_taxon_workflow_alignment.fasta"),
        filtered_alignment_path=fixture(
            "example_taxon_workflow_filtered_alignment.fasta"
        ),
        inference_tree_path=fixture("example_taxon_workflow_inference.nwk"),
        reported_taxa_path=fixture("example_taxon_workflow_reported.csv"),
        out_path=output,
    )

    sidecar = json.loads(result.machine_manifest_path.read_text(encoding="utf-8"))
    assert result.taxon_crosswalk is not None
    assert result.taxon_exclusions is not None
    assert result.taxon_workflow_loss is not None
    assert result.taxon_stability is not None
    assert result.machine_manifest["sections"] == [
        "reviewer-summary",
        "taxon-audit",
        "taxon-identity",
        "taxon-safety",
        "taxon-namespaces",
        "taxon-rank-consistency",
        "taxon-duplicate-identities",
        "taxon-mapping-conflicts",
        "taxon-crosswalk",
        "taxon-exclusions",
        "taxon-loss",
        "taxon-loss-events",
        "taxon-stability",
        "limitations",
    ]
    assert sidecar["metrics"]["crosswalk_rows"] == 4
    assert sidecar["metrics"]["excluded_taxa"] == 2
    assert sidecar["metrics"]["loss_stage_count"] == 3
    assert sidecar["metrics"]["unstable_taxa"] == 3
    assert [row["role"] for row in sidecar["input_ledger"]] == [
        "tree",
        "metadata",
        "traits",
        "alignment",
        "filtered_alignment",
        "inference_tree",
        "reported_taxa",
    ]
    assert any(
        line == "excluded taxa with explicit causes: 2"
        for line in sidecar["reviewer_summary"]
    )
    assert any(
        line == "unstable taxa across linked sources: 3"
        for line in sidecar["reviewer_summary"]
    )


def test_render_reports_write_machine_readable_sidecars_and_reviewer_sections(
    tmp_path: Path,
) -> None:
    tree_output = tmp_path / "tree-report.html"
    dataset_output = tmp_path / "dataset-report.html"
    phylo_output = tmp_path / "phylo-report.html"
    tree_result = render_tree_report(
        tree_path=fixture("example_tree_support_invalid.nwk"), out_path=tree_output
    )
    dataset_result = render_dataset_report(
        tree_path=fixture("example_tree.nwk"),
        metadata_path=fixture("example_metadata.tsv"),
        traits_path=fixture("example_traits.tsv"),
        alignment_path=fixture("example_alignment.fasta"),
        out_path=dataset_output,
    )
    phylo_result = render_phylo_inputs_report(
        tree_path=fixture("example_tree.nwk"),
        alignment_path=fixture("example_alignment_coding.fasta"),
        out_path=phylo_output,
    )
    tree_sidecar = json.loads(
        tree_result.machine_manifest_path.read_text(encoding="utf-8")
    )
    dataset_sidecar = json.loads(
        dataset_result.machine_manifest_path.read_text(encoding="utf-8")
    )
    phylo_sidecar = json.loads(
        phylo_result.machine_manifest_path.read_text(encoding="utf-8")
    )
    assert (
        tree_sidecar["reviewer_summary"][0]
        == "tree validity decision: valid_with_warnings"
    )
    assert "limitations" in dataset_sidecar
    assert dataset_sidecar["input_ledger"][0]["role"] == "tree"
    assert "reviewer-summary" in dataset_result.machine_manifest["sections"]
    assert len(phylo_sidecar["limitations"]) >= 1


def test_reference_validation_suites_cover_goals_91_through_96() -> None:
    suites = [
        validate_tree_reference_fixtures(),
        validate_taxon_naming_reference_fixtures(),
        validate_alignment_quality_reference_fixtures(),
        validate_dataset_audit_reference_fixtures(),
        validate_figure_reference_fixtures(),
        validate_report_regression_fixtures(),
    ]

    assert [suite.goal_id for suite in suites] == [91, 92, 93, 94, 95, 96]
    assert all(suite.passed for suite in suites)
    assert sum(suite.fixture_count for suite in suites) >= 18


def test_alignment_quality_reference_fixtures_lock_suspicion_and_concentration() -> (
    None
):
    suite = validate_alignment_quality_reference_fixtures()
    fixtures = {fixture.name: fixture for fixture in suite.fixtures}

    assert fixtures["clean_alignment_quality"].observed["suspicious_alignment"] is False
    assert (
        fixtures["missingness_heavy_alignment"].observed["concentrated_missing_run"]
        == 2
    )
    assert fixtures["ambiguity_heavy_alignment"].observed["warning_count"] == 5
    assert suite.passed is True


def test_build_core_workflow_validation_report_aggregates_suites_failures_and_maturity() -> (
    None
):
    report = build_core_workflow_validation_report()

    assert report.total_fixture_count == 18
    assert report.passed_fixture_count == 18
    assert report.failed_fixture_count == 0
    assert [suite.goal_id for suite in report.suites] == [91, 92, 93, 94, 95, 96]
    assert len(report.workflows) == 6
    assert len(report.failure_gallery) == 4
    assert all(case.passed for case in report.failure_gallery)
    assert {row.workflow for row in report.maturity_classifications} == {
        "tree-review",
        "taxon-audit",
        "alignment-review",
        "dataset-audit",
        "figure-package",
        "phylo-inputs-review",
    }


def test_render_workflow_validation_report_writes_fixture_sections(
    tmp_path: Path,
) -> None:
    output = tmp_path / "workflow-validation.html"
    result = render_workflow_validation_report(out_path=output)

    text = output.read_text(encoding="utf-8")
    assert result.report_kind == "workflow-validation"
    assert result.machine_manifest["sections"] == [
        "reviewer-summary",
        "validation-overview",
        "validation-suites",
        "workflow-coverage",
        "failure-gallery",
        "maturity-classification",
        "limitations",
    ]
    assert result.machine_manifest["metrics"]["passed_fixture_count"] == 18
    assert result.machine_manifest_path.exists()
    assert "Bijux Core Workflow Validation Report" in text


def test_build_level_one_release_gate_report_tracks_taxon_loss_and_blocked_analysis() -> (
    None
):
    report = build_level_one_release_gate_report()

    assert report.gate.decision == "blocked"
    assert report.dataset_readiness_decision == "blocked"
    assert report.gate.retained_taxa == ["A"]
    assert report.gate.excluded_taxa == ["B", "C", "D"]
    assert report.taxon_first_loss_stage == {
        "A": None,
        "B": "alignment_filtering",
        "C": "trait_missingness",
        "D": "alignment",
    }
    assert report.exclusion_causes["B"][0] == "alignment_filtering"
    assert "maximum_likelihood" in report.gate.blocked_analyses


def test_render_level_one_release_gate_report_writes_traceability_sections(
    tmp_path: Path,
) -> None:
    output = tmp_path / "release-gate.html"
    result = render_level_one_release_gate_report(out_path=output)

    text = output.read_text(encoding="utf-8")
    assert result.report_kind == "release-gate"
    assert result.machine_manifest["sections"] == [
        "reviewer-summary",
        "gate-decision",
        "dataset-readiness",
        "taxon-loss-traceability",
        "workflow-validation",
        "limitations",
    ]
    assert result.machine_manifest["metrics"]["excluded_taxa"] == 3
    assert result.machine_manifest_path.exists()
    assert "Bijux Level 1 Release Gate Report" in text


def test_build_release_truth_report_aggregates_test_counts_and_release_surfaces(
    tmp_path: Path,
) -> None:
    total_report = _write_junit_report(
        tmp_path / "full.xml",
        suite_name="full",
        tests=14,
        failures=1,
        skipped=2,
    )
    real_engine_report = _write_junit_report(
        tmp_path / "real-engine.xml",
        suite_name="real-engine",
        tests=6,
        failures=0,
        skipped=1,
    )

    report = build_release_truth_report(
        test_report_paths=[total_report],
        real_engine_test_report_paths=[real_engine_report],
        include_extended_parity=False,
        stress_tier="small",
    )

    assert report.total_tests.total_tests == 14
    assert report.total_tests.passed_tests == 11
    assert report.real_engine_tests.total_tests == 6
    assert report.real_engine_tests.passed_tests == 5
    assert any(
        workflow.surface == "fasta-to-tree"
        for workflow in report.supported_workflows
    )
    assert any(
        dataset.demo_command == "rabies-cross-host-geography-panel"
        for dataset in report.flagship_datasets
    )
    assert report.reference_parity.case_count > 0
    assert len(report.stress_suite.observations) == 5


def test_render_release_truth_report_writes_release_sections(
    tmp_path: Path,
) -> None:
    total_report = _write_junit_report(
        tmp_path / "full.xml",
        suite_name="full",
        tests=12,
        failures=1,
        skipped=1,
    )
    real_engine_report = _write_junit_report(
        tmp_path / "real-engine.xml",
        suite_name="real-engine",
        tests=5,
        failures=0,
        skipped=2,
    )
    output = tmp_path / "release-truth.html"

    result = render_release_truth_report(
        out_path=output,
        test_report_paths=[total_report],
        real_engine_test_report_paths=[real_engine_report],
        include_extended_parity=False,
        stress_tier="small",
    )

    text = output.read_text(encoding="utf-8")
    assert result.report_kind == "release-truth"
    assert result.machine_manifest["sections"] == [
        "reviewer-summary",
        "total-tests",
        "real-engine-tests",
        "supported-workflows",
        "experimental-workflows",
        "advisory-workflows",
        "parser-only-workflows",
        "flagship-datasets",
        "workflow-validation",
        "release-gate",
        "reference-parity",
        "stress-suite",
        "known-limitations",
    ]
    assert result.machine_manifest["metrics"]["total_tests"] == 12
    assert result.machine_manifest["metrics"]["real_engine_tests"] == 5
    assert result.machine_manifest["metrics"]["flagship_dataset_count"] >= 1
    assert result.machine_manifest_path.exists()
    assert "Bijux Release Truth Report" in text


def test_bundle_directory_copies_files_and_manifest(tmp_path: Path) -> None:
    inputs_root = tmp_path / "inputs"
    outputs_root = tmp_path / "outputs"
    inputs_root.mkdir()
    outputs_root.mkdir()
    (inputs_root / "tree.nwk").write_text("(A:0.1,B:0.1);\n", encoding="utf-8")
    (outputs_root / "report.html").write_text("<html></html>\n", encoding="utf-8")
    report = bundle_directory([inputs_root], [outputs_root], tmp_path / "bundle")
    assert report.file_count == 2
    assert report.input_file_count == 1
    assert report.output_file_count == 1
    manifest = (tmp_path / "bundle" / "manifest.json").read_text(encoding="utf-8")
    assert "tree.nwk" in manifest
    assert (tmp_path / "bundle" / "checksums.tsv").exists()
    assert (tmp_path / "bundle" / "environment.json").exists()
    assert (tmp_path / "bundle" / "README.md").exists()


def test_validate_bundle_accepts_matching_checksums(tmp_path: Path) -> None:
    inputs_root = tmp_path / "inputs"
    outputs_root = tmp_path / "outputs"
    inputs_root.mkdir()
    outputs_root.mkdir()
    (inputs_root / "artifact.txt").write_text("evidence\n", encoding="utf-8")
    (outputs_root / "summary.txt").write_text("result\n", encoding="utf-8")
    bundle_root = tmp_path / "bundle"
    bundle_directory([inputs_root], [outputs_root], bundle_root)
    report = validate_bundle(bundle_root)
    assert report.valid is True
    assert report.mismatches == []


def test_validate_bundle_detects_checksum_drift(tmp_path: Path) -> None:
    inputs_root = tmp_path / "inputs"
    outputs_root = tmp_path / "outputs"
    inputs_root.mkdir()
    outputs_root.mkdir()
    (inputs_root / "artifact.txt").write_text("evidence\n", encoding="utf-8")
    (outputs_root / "summary.txt").write_text("result\n", encoding="utf-8")
    bundle_root = tmp_path / "bundle"
    bundle_directory([inputs_root], [outputs_root], bundle_root)
    (bundle_root / "outputs" / "outputs" / "summary.txt").write_text(
        "drift\n", encoding="utf-8"
    )
    report = validate_bundle(bundle_root)
    assert report.valid is False
    assert report.mismatches[0].reason in {"checksum_mismatch", "size_mismatch"}


def test_bundle_file_paths_copies_explicit_input_and_output_files(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "inputs" / "alignment.fasta"
    output_path = tmp_path / "outputs" / "analysis.xml"
    input_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    input_path.write_text(">A\nACTG\n", encoding="utf-8")
    output_path.write_text("<beast />\n", encoding="utf-8")

    report = bundle_file_paths([input_path], [output_path], tmp_path / "bundle")

    assert report.file_count == 2
    assert report.input_file_count == 1
    assert report.output_file_count == 1
    assert (report.output_root / "manifest.json").exists()


def test_cli_validate_json_output(capsys) -> None:
    exit_code = main(["validate", str(fixture("example_tree.nwk")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "validate"
    assert payload["metrics"]["tip_count"] == 4
    assert payload["data"]["tip_count"] == 4


def test_cli_inspect_reports_zero_length_branch_count(capsys) -> None:
    exit_code = main(
        ["inspect", str(fixture("example_tree_zero_lengths.nwk")), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["zero_length_branch_count"] == 3
    assert payload["data"]["zero_length_branch_count"] == 3
    assert payload["warnings"] == ["tree contains zero-length branches"]


def test_cli_diagnose_json_output(capsys) -> None:
    exit_code = main(["diagnose", str(fixture("example_tree.nwk")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "diagnose"
    assert payload["metrics"]["cherry_count"] == 2
    assert payload["metrics"]["tree_diameter"] == 0.6
    assert payload["metrics"]["tree_quality_score"] == 100.0
    assert payload["data"]["inspection"]["imbalance_summary"] == "balanced"


def test_cli_diagnose_distances_writes_tsv(tmp_path: Path, capsys) -> None:
    output = tmp_path / "distances.tsv"
    exit_code = main(
        [
            "diagnose",
            "distances",
            str(fixture("example_tree.nwk")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["outputs"] == [str(output)]
    assert output.read_text(encoding="utf-8") == (
        "tip\tdistance\nA\t0.3\nB\t0.3\nC\t0.3\nD\t0.3\n"
    )


def test_cli_diagnose_assumptions_reports_unit_aware_compatibility(capsys) -> None:
    exit_code = main(
        [
            "diagnose",
            "assumptions",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_metadata_branch_units_time.tsv")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["standardized_support_count"] == 0
    assert payload["metrics"]["time_tree_compatible"] is True
    assert payload["metrics"]["substitution_tree_compatible"] is False


def test_cli_diagnose_ultrametric_reports_tolerance_and_deviation(capsys) -> None:
    exit_code = main(
        [
            "diagnose",
            "ultrametric",
            str(fixture("example_tree_ladderized.nwk")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["data"]["ultrametric"] is False
    assert payload["metrics"]["max_deviation"] == 0.2


def test_cli_validate_writes_run_manifest(tmp_path: Path, capsys) -> None:
    manifest = tmp_path / "validate.manifest.json"
    exit_code = main(
        [
            "validate",
            str(fixture("example_tree.nwk")),
            "--json",
            "--manifest",
            str(manifest),
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert exit_code == 0
    assert payload["outputs"] == [str(manifest)]
    assert manifest_payload["command"] == "validate"
    assert manifest_payload["arguments"] == [
        "validate",
        str(fixture("example_tree.nwk")),
        "--json",
        "--manifest",
        str(manifest),
    ]
    assert manifest_payload["input_checksums"][str(fixture("example_tree.nwk"))]


def test_cli_normalize_includes_manifest_in_output_list(tmp_path: Path, capsys) -> None:
    output = tmp_path / "normalized.nwk"
    manifest = tmp_path / "normalize.manifest.json"
    exit_code = main(
        [
            "normalize",
            str(fixture("example_tree.nex")),
            "--format",
            "nexus",
            "--out",
            str(output),
            "--json",
            "--manifest",
            str(manifest),
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["outputs"] == [str(output), str(manifest)]


def test_cli_commands_json_lists_registered_taxonomy(capsys) -> None:
    exit_code = main(["commands", "--format", "json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    command_names = [item["name"] for item in payload["data"]["commands"]]
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert command_names == [
        "env",
        "metadata",
        "traits",
        "prune",
        "alignment",
        "comparative",
        "ancestral",
        "biogeography",
        "host-association",
        "ecological-niche",
        "phylogeography",
        "discrete-evolution",
        "diversification",
        "distance",
        "tree-set",
        "simulate",
        "benchmark",
        "parity",
        "inspect",
        "validate",
        "normalize",
        "normalize-taxa",
        "taxonomy",
        "topology",
        "compare",
        "annotate",
        "diagnose",
        "render",
        "report",
        "demo",
        "evidence",
        "adapter",
    ]


def test_supported_evidence_api_contract_resolves_public_comparative_entrypoints() -> (
    None
):
    resolved = {
        locator: resolve_supported_evidence_api(locator)
        for locator in SUPPORTED_EVIDENCE_API_LOCATORS
    }

    assert (
        resolved["bijux_phylogenetics.comparative:inspect_pgls_inputs"]
        is inspect_pgls_inputs
    )
    assert (
        resolved["bijux_phylogenetics.comparative:build_pgls_model_matrix"]
        is build_pgls_model_matrix
    )
    assert (
        resolved["bijux_phylogenetics.comparative:summarize_pgls_lambda_fit"]
        is summarize_pgls_lambda_fit
    )
    assert (
        resolved["bijux_phylogenetics.comparative:summarize_brownian_covariance"]
        is summarize_brownian_covariance
    )
    assert (
        resolved["bijux_phylogenetics.comparative:summarize_brownian_covariance_pgls"]
        is summarize_brownian_covariance_pgls
    )
    assert (
        resolved["bijux_phylogenetics.comparative:summarize_brownian_regime_rates"]
        is summarize_brownian_regime_rates
    )
    assert (
        resolved["bijux_phylogenetics.comparative:summarize_trait_regime_mapping"]
        is summarize_trait_regime_mapping
    )
    assert (
        resolved["bijux_phylogenetics.comparative:summarize_brownian_trait_evolution"]
        is summarize_brownian_trait_evolution
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:summarize_early_burst_trait_evolution"
        ]
        is summarize_early_burst_trait_evolution
    )
    assert (
        resolved["bijux_phylogenetics.comparative:summarize_clade_traits"]
        is summarize_clade_traits
    )
    assert (
        resolved["bijux_phylogenetics.comparative:summarize_trait_outliers"]
        is summarize_trait_outliers
    )
    assert (
        resolved["bijux_phylogenetics.comparative:summarize_trait_imputation"]
        is summarize_trait_imputation
    )
    assert (
        resolved["bijux_phylogenetics.comparative:summarize_trait_rate_through_time"]
        is summarize_trait_rate_through_time
    )
    assert (
        resolved["bijux_phylogenetics.comparative:summarize_ou_covariance_pgls"]
        is summarize_ou_covariance_pgls
    )
    assert (
        resolved["bijux_phylogenetics.comparative:summarize_ou_trait_evolution"]
        is summarize_ou_trait_evolution
    )
    assert (
        resolved["bijux_phylogenetics.comparative:summarize_correlated_trait_evolution"]
        is summarize_correlated_trait_evolution
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_early_burst_trait_evolution_summary_table"
        ]
        is write_early_burst_trait_evolution_summary_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_correlated_trait_summary_table"]
        is write_correlated_trait_summary_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_correlated_trait_comparison_table"
        ]
        is write_correlated_trait_comparison_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_correlated_trait_observation_table"
        ]
        is write_correlated_trait_observation_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_correlated_trait_exclusion_table"
        ]
        is write_correlated_trait_exclusion_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_early_burst_trait_evolution_exclusion_table"
        ]
        is write_early_burst_trait_evolution_exclusion_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_early_burst_trait_evolution_comparison_table"
        ]
        is write_early_burst_trait_evolution_comparison_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_early_burst_rate_change_profile_table"
        ]
        is write_early_burst_rate_change_profile_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_clade_trait_summary_table"]
        is write_clade_trait_summary_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_clade_trait_clade_table"]
        is write_clade_trait_clade_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_clade_trait_exclusion_table"]
        is write_clade_trait_exclusion_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_trait_outlier_summary_table"]
        is write_trait_outlier_summary_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_trait_outlier_taxon_table"]
        is write_trait_outlier_taxon_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_trait_outlier_exclusion_table"]
        is write_trait_outlier_exclusion_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_trait_imputation_summary_table"]
        is write_trait_imputation_summary_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_trait_imputation_table"]
        is write_trait_imputation_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_trait_imputation_holdout_table"]
        is write_trait_imputation_holdout_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_trait_imputation_exclusion_table"
        ]
        is write_trait_imputation_exclusion_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_trait_rate_through_time_summary_table"
        ]
        is write_trait_rate_through_time_summary_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_trait_rate_through_time_interval_table"
        ]
        is write_trait_rate_through_time_interval_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_trait_rate_through_time_exclusion_table"
        ]
        is write_trait_rate_through_time_exclusion_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:summarize_pgls_categorical_contrasts"]
        is summarize_pgls_categorical_contrasts
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:summarize_independent_contrast_regression"
        ]
        is summarize_independent_contrast_regression
    )
    assert (
        resolved["bijux_phylogenetics.comparative:summarize_phylogenetic_signal"]
        is summarize_phylogenetic_signal
    )
    assert (
        resolved["bijux_phylogenetics.comparative:summarize_phylogenetic_logistic"]
        is summarize_phylogenetic_logistic
    )
    assert (
        resolved["bijux_phylogenetics.comparative:run_posterior_tree_pgls"]
        is run_posterior_tree_pgls
    )
    assert (
        resolved["bijux_phylogenetics.comparative:analyze_comparative_clade_stability"]
        is analyze_comparative_clade_stability
    )
    assert (
        resolved["bijux_phylogenetics.comparative:analyze_comparative_residual_clades"]
        is analyze_comparative_residual_clades
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:compare_comparative_regression_models"
        ]
        is compare_comparative_regression_models
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:run_multivariate_comparative_regression"
        ]
        is run_multivariate_comparative_regression
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:summarize_pgls_interaction_coefficients"
        ]
        is summarize_pgls_interaction_coefficients
    )
    assert resolved["bijux_phylogenetics.comparative:run_pgls"] is run_pgls
    assert (
        resolved["bijux_phylogenetics.comparative:write_brownian_covariance_table"]
        is write_brownian_covariance_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_brownian_covariance_long_table"
        ]
        is write_brownian_covariance_long_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_brownian_covariance_matrix_table"
        ]
        is write_brownian_covariance_matrix_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_brownian_regime_summary_table"]
        is write_brownian_regime_summary_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_trait_regime_branch_table"]
        is write_trait_regime_branch_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_trait_regime_exclusion_table"]
        is write_trait_regime_exclusion_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_trait_regime_node_table"]
        is write_trait_regime_node_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_trait_regime_summary_table"]
        is write_trait_regime_summary_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_brownian_regime_rate_table"]
        is write_brownian_regime_rate_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_brownian_regime_profile_table"]
        is write_brownian_regime_profile_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_brownian_regime_comparison_table"
        ]
        is write_brownian_regime_comparison_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_brownian_regime_branch_table"]
        is write_brownian_regime_branch_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_brownian_regime_exclusion_table"
        ]
        is write_brownian_regime_exclusion_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_brownian_trait_evolution_summary_table"
        ]
        is write_brownian_trait_evolution_summary_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_brownian_trait_evolution_exclusion_table"
        ]
        is write_brownian_trait_evolution_exclusion_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_comparative_clade_stability_table"
        ]
        is write_comparative_clade_stability_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_comparative_clade_coefficient_change_table"
        ]
        is write_comparative_clade_coefficient_change_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_posterior_tree_pgls_tree_table"]
        is write_posterior_tree_pgls_tree_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_posterior_tree_pgls_coefficient_table"
        ]
        is write_posterior_tree_pgls_coefficient_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_posterior_tree_pgls_summary_table"
        ]
        is write_posterior_tree_pgls_summary_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_comparative_residual_taxon_table"
        ]
        is write_comparative_residual_taxon_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_comparative_residual_clade_table"
        ]
        is write_comparative_residual_clade_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_comparative_regression_model_ranking_table"
        ]
        is write_comparative_regression_model_ranking_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_comparative_regression_pairwise_table"
        ]
        is write_comparative_regression_pairwise_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_comparative_regression_excluded_taxa_table"
        ]
        is write_comparative_regression_excluded_taxa_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_phylogenetic_logistic_coefficient_table"
        ]
        is write_phylogenetic_logistic_coefficient_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_phylogenetic_logistic_excluded_taxa_table"
        ]
        is write_phylogenetic_logistic_excluded_taxa_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_phylogenetic_logistic_fitted_table"
        ]
        is write_phylogenetic_logistic_fitted_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_independent_contrast_table"]
        is write_independent_contrast_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_multivariate_residual_covariance_table"
        ]
        is write_multivariate_residual_covariance_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_multivariate_residual_correlation_table"
        ]
        is write_multivariate_residual_correlation_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_multivariate_residual_association_table"
        ]
        is write_multivariate_residual_association_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_multivariate_response_model_table"
        ]
        is write_multivariate_response_model_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_multivariate_response_coefficient_table"
        ]
        is write_multivariate_response_coefficient_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_multivariate_excluded_taxa_table"
        ]
        is write_multivariate_excluded_taxa_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_independent_contrast_regression_table"
        ]
        is write_independent_contrast_regression_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_phylogenetic_signal_summary_table"
        ]
        is write_phylogenetic_signal_summary_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_phylogenetic_signal_permutation_table"
        ]
        is write_phylogenetic_signal_permutation_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_ou_alpha_profile_table"]
        is write_ou_alpha_profile_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_ou_trait_evolution_summary_table"
        ]
        is write_ou_trait_evolution_summary_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_ou_trait_evolution_exclusion_table"
        ]
        is write_ou_trait_evolution_exclusion_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_ou_covariance_table"]
        is write_ou_covariance_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_pgls_categorical_contrast_table"
        ]
        is write_pgls_categorical_contrast_table
    )
    assert (
        resolved[
            "bijux_phylogenetics.comparative:write_pgls_interaction_coefficient_table"
        ]
        is write_pgls_interaction_coefficient_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_pgls_lambda_profile_table"]
        is write_pgls_lambda_profile_table
    )
    assert (
        resolved["bijux_phylogenetics.comparative:write_pgls_model_matrix_table"]
        is write_pgls_model_matrix_table
    )


def test_run_capability_demo_creates_expected_artifacts(tmp_path: Path) -> None:
    result = run_capability_demo(tmp_path / "demo")
    assert result.tree_report.exists()
    assert result.dataset_report.exists()
    assert result.phylo_inputs_report.exists()
    assert result.comparison_report.exists()
    assert result.capability_summary.exists()


def test_cli_evidence_bundle_and_validate_json_output(tmp_path: Path, capsys) -> None:
    inputs_root = tmp_path / "inputs"
    outputs_root = tmp_path / "outputs"
    bundle_root = tmp_path / "bundle"
    inputs_root.mkdir()
    outputs_root.mkdir()
    (inputs_root / "artifact.txt").write_text("evidence\n", encoding="utf-8")
    (outputs_root / "summary.txt").write_text("result\n", encoding="utf-8")

    exit_code = main(
        [
            "evidence",
            "bundle",
            "--inputs",
            str(inputs_root),
            "--outputs",
            str(outputs_root),
            "--out",
            str(bundle_root),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    bundle_payload = json.loads(captured.out)
    assert exit_code == 0
    assert bundle_payload["status"] == "ok"
    assert bundle_payload["metrics"]["input_file_count"] == 1
    assert bundle_payload["metrics"]["output_file_count"] == 1

    exit_code = main(["evidence", "validate", str(bundle_root), "--json"])
    captured = capsys.readouterr()
    validate_payload = json.loads(captured.out)
    assert exit_code == 0
    assert validate_payload["status"] == "ok"
    assert validate_payload["data"]["valid"] is True


def test_cli_evidence_book_studies_json_output(capsys, monkeypatch) -> None:
    monkeypatch.chdir(REPO_ROOT)

    exit_code = main(["evidence", "book", "studies", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["study_count"] == 2
    assert payload["metrics"]["partial_rerun_capable_count"] == 2
    assert payload["data"]["studies"][0]["study_id"] == "primate-longevity-signal"


@pytest.mark.slow
def test_cli_evidence_book_validate_json_output(capsys, monkeypatch) -> None:
    monkeypatch.chdir(REPO_ROOT)

    exit_code = main(["evidence", "book", "validate", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["data"]["valid"] is True
    assert payload["metrics"]["bundle_count"] == 19
    assert payload["metrics"]["coverage_gap_count"] >= 1
    assert payload["metrics"]["downgraded_claim_count"] >= 1
    assert payload["metrics"]["foundational_numerical_trust_status"] == "bounded"
    assert payload["metrics"]["maturity_tier"] == "reviewable_but_incomplete"
    assert payload["metrics"]["completion_not_ready_count"] == 2


@pytest.mark.slow
def test_cli_evidence_book_build_json_output(capsys, monkeypatch) -> None:
    monkeypatch.chdir(REPO_ROOT)

    exit_code = main(["evidence", "book", "build", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["reviewer_summary_count"] == 19
    assert payload["metrics"]["updated_path_count"] >= 1
    assert payload["metrics"]["bundle_count"] == 19
    assert payload["metrics"]["reviewer_readiness_status"] == "bounded"


@pytest.mark.slow
def test_cli_evidence_book_build_selected_evidence_json_output(
    capsys, monkeypatch
) -> None:
    monkeypatch.chdir(REPO_ROOT)

    exit_code = main(
        [
            "evidence",
            "book",
            "build",
            "primate-pgls-and-signal",
            "--evidence-id",
            "evidence-002",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["selected_study_count"] == 1
    assert payload["metrics"]["selected_evidence_count"] == 1
    assert payload["data"]["study_id"] == "primate-pgls-and-signal"
    assert payload["data"]["selected_evidence_ids"] == ["evidence-002"]


@pytest.mark.slow
def test_cli_evidence_book_rerun_json_output(capsys, monkeypatch) -> None:
    monkeypatch.chdir(REPO_ROOT)

    exit_code = main(
        [
            "evidence",
            "book",
            "rerun",
            "primate-pgls-and-signal",
            "evidence-002",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["selected_evidence_count"] == 1
    assert payload["metrics"]["downgraded_claim_count"] >= 1
    assert payload["data"]["rerun_report"]["study_id"] == "primate-pgls-and-signal"
    assert payload["data"]["rerun_report"]["selected_evidence_ids"] == ["evidence-002"]


def test_cli_demo_run_json_output_reports_generated_artifacts(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "demo"
    exit_code = main(["demo", "run", "--out", str(output), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["command"] == "demo"
    assert payload["metrics"]["artifact_count"] == 5
    assert payload["data"]["tree_report"] == str(
        output / "reports" / "tree-report.html"
    )
    assert payload["data"]["capability_summary"] == str(
        output / "capability-summary.md"
    )


def test_cli_report_json_output_uses_result_envelope(tmp_path: Path, capsys) -> None:
    output = tmp_path / "report.html"
    exit_code = main(
        [
            "report",
            "tree",
            str(fixture("example_tree.nwk")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "report"
    assert payload["outputs"] == [str(output), str(output.with_suffix(".json"))]
    assert payload["data"]["report_kind"] == "tree"


def test_cli_report_dataset_json_output_uses_dataset_contract(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "dataset-report.html"
    exit_code = main(
        [
            "report",
            "dataset",
            "--tree",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_metadata.tsv")),
            "--traits",
            str(fixture("example_traits.tsv")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["command"] == "report"
    assert payload["outputs"] == [str(output), str(output.with_suffix(".json"))]
    assert payload["data"]["report_kind"] == "dataset"
    assert payload["metrics"]["linked_taxa"] == 4
    assert payload["data"]["input_ledger"][0]["role"] == "tree"


def test_cli_report_phylo_inputs_json_output_uses_alignment_contract(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "phylo-inputs-report.html"
    exit_code = main(
        [
            "report",
            "phylo-inputs",
            "--tree",
            str(fixture("example_tree.nwk")),
            "--alignment",
            str(fixture("example_alignment.fasta")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["command"] == "report"
    assert payload["outputs"] == [str(output), str(output.with_suffix(".json"))]
    assert payload["data"]["report_kind"] == "phylo-inputs"
    assert payload["metrics"]["alignment_length"] == 8
    assert payload["metrics"]["linked_taxa"] == 4


def test_cli_report_alignment_json_output_uses_alignment_contract(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "alignment-report.html"
    exit_code = main(
        [
            "report",
            "alignment",
            "--alignment",
            str(fixture("example_alignment_coding.fasta")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["command"] == "report"
    assert payload["outputs"] == [str(output), str(output.with_suffix(".json"))]
    assert payload["data"]["report_kind"] == "alignment"
    assert payload["metrics"]["sequence_count"] == 4


def test_cli_report_taxonomy_json_output_uses_taxon_contract(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "taxonomy-report.html"
    exit_code = main(
        [
            "report",
            "taxonomy",
            "--tree",
            str(fixture("example_taxonomy_tree.nwk")),
            "--synonym-table",
            str(fixture("example_taxon_synonyms.tsv")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["command"] == "report"
    assert payload["outputs"] == [str(output), str(output.with_suffix(".json"))]
    assert payload["data"]["report_kind"] == "taxonomy"
    assert payload["metrics"]["tree_tip_count"] == 6


def test_cli_report_workflow_validation_json_output_uses_reference_contract(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "workflow-validation.html"
    exit_code = main(
        [
            "report",
            "workflow-validation",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["command"] == "report"
    assert payload["outputs"] == [str(output), str(output.with_suffix(".json"))]
    assert payload["data"]["report_kind"] == "workflow-validation"
    assert payload["metrics"]["passed_fixture_count"] == 18


def test_cli_report_release_gate_json_output_uses_gate_contract(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "release-gate.html"
    exit_code = main(
        [
            "report",
            "release-gate",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["command"] == "report"
    assert payload["outputs"] == [str(output), str(output.with_suffix(".json"))]
    assert payload["data"]["report_kind"] == "release-gate"
    assert payload["metrics"]["decision"] == "blocked"
    assert payload["metrics"]["excluded_taxa"] == 3


def test_cli_report_release_truth_json_output_uses_release_contract(
    tmp_path: Path, capsys
) -> None:
    total_report = _write_junit_report(
        tmp_path / "full.xml",
        suite_name="full",
        tests=10,
        failures=1,
        skipped=2,
    )
    real_engine_report = _write_junit_report(
        tmp_path / "real-engine.xml",
        suite_name="real-engine",
        tests=4,
        failures=0,
        skipped=1,
    )
    output = tmp_path / "release-truth.html"
    exit_code = main(
        [
            "report",
            "release-truth",
            "--test-report",
            str(total_report),
            "--real-engine-test-report",
            str(real_engine_report),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["command"] == "report"
    assert payload["outputs"] == [str(output), str(output.with_suffix(".json"))]
    assert payload["data"]["report_kind"] == "release-truth"
    assert payload["metrics"]["total_tests"] == 10
    assert payload["metrics"]["real_engine_tests"] == 4
    assert payload["metrics"]["supported_workflow_count"] >= 1


def test_cli_adapter_returns_typed_engine_error(capsys) -> None:
    exit_code = main(
        [
            "adapter",
            "inspect",
            "iqtree",
            "--executable",
            "definitely-not-installed-engine",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == EngineUnavailableError.code


def test_public_package_exports_python_workflow_api_surface() -> None:
    from bijux_phylogenetics.api import (
        AlignmentWorkflowResult,
        AncestralReconstructionWorkflowResult,
        ComparativeModelWorkflowResult,
        ConfiguredPhyloWorkflowResult,
        FastaValidationResult,
        InferenceWorkflowResult,
        ReportWorkflowResult,
        SequenceToTreeWorkflowResult,
        SupportWorkflowResult,
        TreeComparisonWorkflowResult,
        TrimmingWorkflowResult,
        render_report_workflow,
        run_alignment_workflow,
        run_ancestral_reconstruction_workflow,
        run_comparative_model_workflow,
        run_configured_phylo_workflow,
        run_fasta_validation_workflow,
        run_sequence_to_tree_workflow,
        run_support_workflow,
        run_tree_comparison_workflow,
        run_tree_inference_workflow,
        run_trimming_workflow,
    )

    assert bijux_phylogenetics.AlignmentWorkflowResult is AlignmentWorkflowResult
    assert (
        bijux_phylogenetics.AncestralReconstructionWorkflowResult
        is AncestralReconstructionWorkflowResult
    )
    assert (
        bijux_phylogenetics.ComparativeModelWorkflowResult
        is ComparativeModelWorkflowResult
    )
    assert (
        bijux_phylogenetics.ConfiguredPhyloWorkflowResult
        is ConfiguredPhyloWorkflowResult
    )
    assert bijux_phylogenetics.FastaValidationResult is FastaValidationResult
    assert bijux_phylogenetics.InferenceWorkflowResult is InferenceWorkflowResult
    assert bijux_phylogenetics.ReportWorkflowResult is ReportWorkflowResult
    assert (
        bijux_phylogenetics.SequenceToTreeWorkflowResult
        is SequenceToTreeWorkflowResult
    )
    assert bijux_phylogenetics.SupportWorkflowResult is SupportWorkflowResult
    assert (
        bijux_phylogenetics.TreeComparisonWorkflowResult
        is TreeComparisonWorkflowResult
    )
    assert bijux_phylogenetics.TrimmingWorkflowResult is TrimmingWorkflowResult
    assert bijux_phylogenetics.run_alignment_workflow is run_alignment_workflow
    assert (
        bijux_phylogenetics.run_ancestral_reconstruction_workflow
        is run_ancestral_reconstruction_workflow
    )
    assert (
        bijux_phylogenetics.run_comparative_model_workflow
        is run_comparative_model_workflow
    )
    assert (
        bijux_phylogenetics.run_configured_phylo_workflow
        is run_configured_phylo_workflow
    )
    assert (
        bijux_phylogenetics.run_fasta_validation_workflow
        is run_fasta_validation_workflow
    )
    assert (
        bijux_phylogenetics.run_sequence_to_tree_workflow
        is run_sequence_to_tree_workflow
    )
    assert bijux_phylogenetics.run_support_workflow is run_support_workflow
    assert (
        bijux_phylogenetics.run_tree_comparison_workflow
        is run_tree_comparison_workflow
    )
    assert (
        bijux_phylogenetics.run_tree_inference_workflow
        is run_tree_inference_workflow
    )
    assert bijux_phylogenetics.run_trimming_workflow is run_trimming_workflow
    assert bijux_phylogenetics.render_report_workflow is render_report_workflow
