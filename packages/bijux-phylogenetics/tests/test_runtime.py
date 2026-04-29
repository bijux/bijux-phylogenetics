from __future__ import annotations

import json
from pathlib import Path

import bijux_phylogenetics
from bijux_phylogenetics.ancestral import (
    build_ancestral_figure_package,
    build_ancestral_sensitivity_report,
    compare_continuous_ancestral_models,
    reconstruct_continuous_ancestral_states,
    reconstruct_discrete_ancestral_states,
    render_ancestral_state_report,
    render_ancestral_state_tree,
    write_ancestral_state_table,
)
from bijux_phylogenetics.bayesian import (
    assess_beast_burnin_sensitivity,
    assess_beast_chain_mixing,
    assess_beast_convergence,
    assess_mrbayes_convergence,
    build_posterior_uncertainty_figure_package,
    build_bayesian_evidence_package,
    compare_bayesian_tree_sets,
    compare_independent_bayesian_runs,
    compare_posterior_tree_sets_by_clock,
    compare_posterior_tree_sets_by_prior,
    compute_mrbayes_effective_sample_sizes,
    detect_impossible_calibration_constraints,
    parse_beast_log,
    parse_mrbayes_parameter_traces,
    prepare_beast_time_tree_analysis,
    prepare_mrbayes_analysis,
    render_bayesian_diagnostics_report,
    render_bayesian_posterior_report,
    render_bayesian_run_comparison_report,
    render_calibration_audit_report,
    run_mrbayes_posterior_inference,
    summarize_maximum_clade_credibility_tree,
    summarize_mrbayes_posterior_trees,
    summarize_posterior_node_ages,
    thin_posterior_tree_set,
    validate_beast_posterior_log,
    validate_fossil_calibration_table,
    validate_tip_dating_metadata,
    write_bayesian_methods_summary_text,
    write_supplementary_bayesian_diagnostics_table,
)
from bijux_phylogenetics.benchmark import (
    benchmark_alignment_diagnostics,
    benchmark_tree_comparison,
    benchmark_tree_validation,
)
from bijux_phylogenetics.cli import main
from bijux_phylogenetics.compare.topology import (
    compare_branch_lengths,
    compare_clade_sets,
    compare_support_values,
    compare_tree_paths,
    detect_clade_changes,
    prune_trees_to_shared_taxa,
    write_tree_comparison_table,
)
from bijux_phylogenetics.comparative import (
    compute_blombergs_k,
    compute_phylogenetic_independent_contrasts,
    compute_phylogenetic_signal_test,
    estimate_pagels_lambda,
    inspect_pgls_inputs,
    run_pgls,
    summarize_numeric_trait,
    summarize_numeric_trait_readiness,
)
from bijux_phylogenetics.compare.reports import build_tree_comparison_report
from bijux_phylogenetics.core.alignment import AlignmentSummary
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
from bijux_phylogenetics.core.manifest import build_run_manifest, write_run_manifest
from bijux_phylogenetics.core.metadata import inspect_metadata_table, join_table_to_taxa
from bijux_phylogenetics.distance import (
    build_distance_tree,
    build_tree_from_imported_distance_matrix,
    compare_distance_tree_topologies,
    compute_pairwise_genetic_distance_matrix,
    load_imported_distance_matrix,
    validate_imported_distance_matrix,
)
from bijux_phylogenetics.discrete_evolution import (
    compare_discrete_state_models,
    detect_state_imbalance_problems,
    estimate_ancestral_geographic_states,
    load_stochastic_map_collection,
    render_discrete_state_evolution_report,
    render_tree_with_geographic_states,
    run_discrete_state_transition_model,
    simulate_discrete_stochastic_maps,
    summarize_discrete_stochastic_maps,
    validate_discrete_state_coding,
    write_discrete_model_comparison_table,
    write_node_state_probability_table,
    write_stochastic_map_collection,
    write_stochastic_map_summary_table,
    write_transition_summary_table,
)
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.diversification import (
    compare_diversification_models,
    compute_lineage_through_time_curve,
    detect_diversification_outlier_clades,
    detect_incomplete_taxon_sampling_metadata,
    estimate_diversification_rate,
    inspect_diversification_time_tree,
    render_diversification_report,
    run_trait_dependent_diversification_analysis,
    validate_time_tree_for_diversification,
    write_clade_diversification_table,
    write_lineage_through_time_table,
    write_trait_dependent_diversification_table,
)
from bijux_phylogenetics.engines import (
    audit_alignment_inference_readiness,
    classify_inference_workflow_failure,
    compare_inferred_tree_to_taxon_metadata,
    compare_fast_and_ml_trees,
    render_inference_workflow_report,
    run_alignment_trimming,
    run_bootstrap_consensus_tree,
    run_bootstrap_support_estimation,
    run_fast_tree_inference,
    run_maximum_likelihood_tree_inference,
    run_model_selection,
    run_multiple_sequence_alignment,
    validate_bootstrap_tree_set,
    validate_inference_engine_outputs,
    validate_ml_tree_contains_expected_taxa,
    validate_model_selection_against_engine_outputs,
)
from bijux_phylogenetics.core.pruning import (
    drop_tree_taxa,
    prune_alignment_to_tree,
    prune_tree_to_alignment,
    prune_tree_to_requested_taxa,
    prune_tree_to_taxa,
)
from bijux_phylogenetics.core.topology import (
    collapse_branches_below_length,
    extract_named_clade,
    ladderize_tree,
    reroot_tree_by_midpoint,
    root_tree_on_outgroup,
    sort_tree_tips_alphabetically,
    unroot_tree,
)
from bijux_phylogenetics.core.taxonomy import (
    inspect_tree_taxa_safety,
    inspect_tree_taxon_identity,
    normalize_tree_taxa,
    write_taxon_mapping,
)
from bijux_phylogenetics.core.traits import (
    detect_missing_trait_values,
    detect_unusable_trait_columns,
    link_tree_to_traits,
    prune_traits_to_tree,
    validate_traits_table,
)
from bijux_phylogenetics.diagnostics.root_to_tip import compute_root_to_tip_distances
from bijux_phylogenetics.diagnostics.root_to_tip import diagnose_ultrametricity
from bijux_phylogenetics.diagnostics.assumptions import (
    assess_tree_assumptions,
    inspect_branch_length_units,
    standardize_support_labels,
)
from bijux_phylogenetics.diagnostics.validation import diagnose_tree_path, forensic_tree_path, inspect_tree_path, validate_tree_path
from bijux_phylogenetics.evidence.bundles import bundle_directory, bundle_file_paths, validate_bundle
from bijux_phylogenetics.errors import (
    AlignmentTaxonMismatchError,
    DuplicateTaxonError,
    EngineUnavailableError,
    InvalidBranchLengthError,
    InvalidAlignmentError,
    InvalidDistanceMatrixError,
    MetadataJoinError,
    NonUltrametricTreeError,
    UnnamedTipError,
    UnsupportedTreeFormatError,
    UnrootedTreeError,
)
from bijux_phylogenetics.identity import IDENTITY
from bijux_phylogenetics.io.newick import dumps_newick, loads_newick
from bijux_phylogenetics.io.nexus import load_nexus
from bijux_phylogenetics.io.phyloxml import load_phyloxml
from bijux_phylogenetics.io.roundtrip import validate_tree_roundtrip
from bijux_phylogenetics.io.trees import detect_tree_format
from bijux_phylogenetics.io.fasta import link_alignment_to_tree, load_fasta_alignment, summarise_fasta
from bijux_phylogenetics.io.fasta import (
    build_alignment_forensic_report,
    build_alignment_quality_report,
    clean_alignment_with_profile,
    classify_alignment_sequences,
    compare_alignment_versions,
    compute_pairwise_sequence_identity_matrix,
    compute_amino_acid_composition,
    compute_nucleotide_composition,
    detect_composition_outlier_sequences,
    detect_identical_duplicate_sequences,
    detect_invalid_alignment_characters,
    detect_near_duplicate_sequences,
    detect_over_aligned_regions,
    detect_sequence_length_outliers,
    detect_under_aligned_regions,
    get_alignment_filter_profile,
    infer_alignment_alphabet,
    list_alignment_filter_profiles,
    detect_sequences_with_excessive_missing_data,
    detect_sites_with_excessive_missing_data,
    inspect_coding_alignment,
    summarize_alignment_readiness,
    summarize_alignment_windows,
    remove_all_gap_columns,
    remove_all_missing_columns,
    remove_sequences_above_missingness_threshold,
    trim_columns_above_missingness_threshold,
    translate_coding_alignment,
    trim_alignment,
    write_fasta_alignment,
)
from bijux_phylogenetics.render.package import build_tree_figure_package
from bijux_phylogenetics.render.svg import AnnotationStrip, render_tree_svg
from bijux_phylogenetics.reports.service import (
    annotate_tree_against_table,
    distance_method_limitations,
    render_dataset_report,
    render_distance_report,
    render_phylo_inputs_report,
    render_phylogenetics_report,
    render_tree_set_comparison_report,
    render_tree_uncertainty_report,
    render_tree_report,
)
from bijux_phylogenetics.simulation import (
    simulate_birth_death_trees,
    simulate_brownian_traits,
    simulate_coalescent_trees,
    simulate_discrete_traits,
    simulate_dna_alignment,
    simulate_ou_traits,
    simulate_protein_alignment,
    write_continuous_trait_table,
    write_discrete_trait_table,
    write_simulated_alignment,
    write_tree_set,
)
from bijux_phylogenetics.tree_set import (
    cluster_trees_by_topology,
    compare_bootstrap_and_posterior_uncertainty,
    compare_posterior_topological_diversity,
    compare_posterior_tree_sets,
    compute_clade_frequency_table,
    compute_consensus_tree,
    compute_tree_distance_matrix,
    detect_posterior_topology_multimodality,
    detect_unstable_clades,
    detect_unstable_taxa,
    load_tree_set,
    summarize_clade_credibility_conflicts,
    summarize_uncertainty_aware_conclusions,
    write_clade_credibility_conflict_table,
    write_topology_cluster_table,
    write_uncertainty_conclusion_table,
)


FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_package_identity_matches_canonical_names() -> None:
    assert bijux_phylogenetics.__name__ == "bijux_phylogenetics"
    assert IDENTITY.package_name == "bijux-phylogenetics"
    assert IDENTITY.import_name == "bijux_phylogenetics"
    assert IDENTITY.cli_name == "bijux-phylogenetics"
    assert "bijux phylogenetics" == IDENTITY.umbrella_command


def test_public_package_exports_alignment_and_topology_workflows() -> None:
    assert bijux_phylogenetics.summarise_fasta is summarise_fasta
    assert bijux_phylogenetics.build_alignment_quality_report is build_alignment_quality_report
    assert bijux_phylogenetics.build_alignment_forensic_report is build_alignment_forensic_report
    assert bijux_phylogenetics.classify_alignment_sequences is classify_alignment_sequences
    assert bijux_phylogenetics.clean_alignment_with_profile is clean_alignment_with_profile
    assert bijux_phylogenetics.compare_alignment_versions is compare_alignment_versions
    assert bijux_phylogenetics.compute_pairwise_genetic_distance_matrix is compute_pairwise_genetic_distance_matrix
    assert bijux_phylogenetics.build_distance_tree is build_distance_tree
    assert bijux_phylogenetics.build_tree_from_imported_distance_matrix is build_tree_from_imported_distance_matrix
    assert bijux_phylogenetics.compare_distance_tree_topologies is compare_distance_tree_topologies
    assert bijux_phylogenetics.validate_imported_distance_matrix is validate_imported_distance_matrix
    assert bijux_phylogenetics.validate_discrete_state_coding is validate_discrete_state_coding
    assert bijux_phylogenetics.detect_state_imbalance_problems is detect_state_imbalance_problems
    assert bijux_phylogenetics.run_discrete_state_transition_model is run_discrete_state_transition_model
    assert bijux_phylogenetics.estimate_ancestral_geographic_states is estimate_ancestral_geographic_states
    assert bijux_phylogenetics.compare_discrete_state_models is compare_discrete_state_models
    assert bijux_phylogenetics.write_discrete_model_comparison_table is write_discrete_model_comparison_table
    assert bijux_phylogenetics.write_node_state_probability_table is write_node_state_probability_table
    assert bijux_phylogenetics.write_transition_summary_table is write_transition_summary_table
    assert bijux_phylogenetics.simulate_discrete_stochastic_maps is simulate_discrete_stochastic_maps
    assert bijux_phylogenetics.summarize_discrete_stochastic_maps is summarize_discrete_stochastic_maps
    assert bijux_phylogenetics.write_stochastic_map_collection is write_stochastic_map_collection
    assert bijux_phylogenetics.write_stochastic_map_summary_table is write_stochastic_map_summary_table
    assert bijux_phylogenetics.load_stochastic_map_collection is load_stochastic_map_collection
    assert bijux_phylogenetics.validate_time_tree_for_diversification is validate_time_tree_for_diversification
    assert bijux_phylogenetics.inspect_diversification_time_tree is inspect_diversification_time_tree
    assert bijux_phylogenetics.compute_lineage_through_time_curve is compute_lineage_through_time_curve
    assert bijux_phylogenetics.detect_incomplete_taxon_sampling_metadata is detect_incomplete_taxon_sampling_metadata
    assert bijux_phylogenetics.estimate_diversification_rate is estimate_diversification_rate
    assert bijux_phylogenetics.compare_diversification_models is compare_diversification_models
    assert bijux_phylogenetics.detect_diversification_outlier_clades is detect_diversification_outlier_clades
    assert bijux_phylogenetics.run_trait_dependent_diversification_analysis is run_trait_dependent_diversification_analysis
    assert bijux_phylogenetics.render_diversification_report is render_diversification_report
    assert bijux_phylogenetics.write_lineage_through_time_table is write_lineage_through_time_table
    assert bijux_phylogenetics.write_clade_diversification_table is write_clade_diversification_table
    assert bijux_phylogenetics.write_trait_dependent_diversification_table is write_trait_dependent_diversification_table
    assert bijux_phylogenetics.render_tree_with_geographic_states is render_tree_with_geographic_states
    assert bijux_phylogenetics.render_discrete_state_evolution_report is render_discrete_state_evolution_report
    assert bijux_phylogenetics.assess_tree_assumptions is assess_tree_assumptions
    assert bijux_phylogenetics.inspect_coding_alignment is inspect_coding_alignment
    assert bijux_phylogenetics.compute_pairwise_sequence_identity_matrix is compute_pairwise_sequence_identity_matrix
    assert bijux_phylogenetics.detect_sequence_length_outliers is detect_sequence_length_outliers
    assert bijux_phylogenetics.detect_over_aligned_regions is detect_over_aligned_regions
    assert bijux_phylogenetics.detect_under_aligned_regions is detect_under_aligned_regions
    assert bijux_phylogenetics.list_alignment_filter_profiles is list_alignment_filter_profiles
    assert bijux_phylogenetics.get_alignment_filter_profile is get_alignment_filter_profile
    assert bijux_phylogenetics.summarize_alignment_windows is summarize_alignment_windows
    assert bijux_phylogenetics.summarize_alignment_readiness is summarize_alignment_readiness
    assert bijux_phylogenetics.audit_dataset_inputs is audit_dataset_inputs
    assert bijux_phylogenetics.audit_dataset_taxon_ordering is audit_dataset_taxon_ordering
    assert bijux_phylogenetics.build_dataset_completeness_matrix is build_dataset_completeness_matrix
    assert bijux_phylogenetics.build_dataset_crosswalk is build_dataset_crosswalk
    assert bijux_phylogenetics.summarize_dataset_readiness is summarize_dataset_readiness
    assert bijux_phylogenetics.load_tree_set is load_tree_set
    assert bijux_phylogenetics.compute_consensus_tree is compute_consensus_tree
    assert bijux_phylogenetics.compute_clade_frequency_table is compute_clade_frequency_table
    assert bijux_phylogenetics.compute_tree_distance_matrix is compute_tree_distance_matrix
    assert bijux_phylogenetics.cluster_trees_by_topology is cluster_trees_by_topology
    assert bijux_phylogenetics.detect_unstable_taxa is detect_unstable_taxa
    assert bijux_phylogenetics.detect_unstable_clades is detect_unstable_clades
    assert bijux_phylogenetics.compare_posterior_topological_diversity is compare_posterior_topological_diversity
    assert bijux_phylogenetics.detect_posterior_topology_multimodality is detect_posterior_topology_multimodality
    assert bijux_phylogenetics.summarize_clade_credibility_conflicts is summarize_clade_credibility_conflicts
    assert bijux_phylogenetics.summarize_uncertainty_aware_conclusions is summarize_uncertainty_aware_conclusions
    assert bijux_phylogenetics.compare_posterior_tree_sets is compare_posterior_tree_sets
    assert bijux_phylogenetics.render_tree_uncertainty_report is render_tree_uncertainty_report
    assert bijux_phylogenetics.simulate_birth_death_trees is simulate_birth_death_trees
    assert bijux_phylogenetics.simulate_coalescent_trees is simulate_coalescent_trees
    assert bijux_phylogenetics.simulate_brownian_traits is simulate_brownian_traits
    assert bijux_phylogenetics.simulate_ou_traits is simulate_ou_traits
    assert bijux_phylogenetics.simulate_discrete_traits is simulate_discrete_traits
    assert bijux_phylogenetics.simulate_dna_alignment is simulate_dna_alignment
    assert bijux_phylogenetics.simulate_protein_alignment is simulate_protein_alignment


def test_command_registry_exposes_discrete_evolution_surface() -> None:
    spec = get_command_spec("discrete-evolution")

    assert spec.domain == "discrete-state-evolution"
    assert spec.outputs == ("discrete-state-evolution-report",)


def test_command_registry_exposes_diversification_surface() -> None:
    spec = get_command_spec("diversification")

    assert spec.domain == "diversification-analysis"
    assert spec.outputs == ("diversification-report",)


def test_public_package_exports_comparative_and_bayesian_workflows() -> None:
    assert bijux_phylogenetics.benchmark_tree_validation is benchmark_tree_validation
    assert bijux_phylogenetics.benchmark_tree_comparison is benchmark_tree_comparison
    assert bijux_phylogenetics.benchmark_alignment_diagnostics is benchmark_alignment_diagnostics
    assert bijux_phylogenetics.summarize_numeric_trait_readiness is summarize_numeric_trait_readiness
    assert bijux_phylogenetics.summarize_numeric_trait is summarize_numeric_trait
    assert bijux_phylogenetics.compute_phylogenetic_independent_contrasts is compute_phylogenetic_independent_contrasts
    assert bijux_phylogenetics.compute_blombergs_k is compute_blombergs_k
    assert bijux_phylogenetics.estimate_pagels_lambda is estimate_pagels_lambda
    assert bijux_phylogenetics.compute_phylogenetic_signal_test is compute_phylogenetic_signal_test
    assert bijux_phylogenetics.inspect_pgls_inputs is inspect_pgls_inputs
    assert bijux_phylogenetics.run_pgls is run_pgls
    assert bijux_phylogenetics.reconstruct_continuous_ancestral_states is reconstruct_continuous_ancestral_states
    assert bijux_phylogenetics.reconstruct_discrete_ancestral_states is reconstruct_discrete_ancestral_states
    assert bijux_phylogenetics.build_ancestral_figure_package is build_ancestral_figure_package
    assert bijux_phylogenetics.build_ancestral_sensitivity_report is build_ancestral_sensitivity_report
    assert bijux_phylogenetics.compare_continuous_ancestral_models is compare_continuous_ancestral_models
    assert bijux_phylogenetics.render_ancestral_state_tree is render_ancestral_state_tree
    assert bijux_phylogenetics.render_ancestral_state_report is render_ancestral_state_report
    assert bijux_phylogenetics.write_ancestral_state_table is write_ancestral_state_table
    assert bijux_phylogenetics.run_multiple_sequence_alignment is run_multiple_sequence_alignment
    assert bijux_phylogenetics.run_alignment_trimming is run_alignment_trimming
    assert bijux_phylogenetics.audit_alignment_inference_readiness is audit_alignment_inference_readiness
    assert bijux_phylogenetics.run_model_selection is run_model_selection
    assert bijux_phylogenetics.validate_model_selection_against_engine_outputs is validate_model_selection_against_engine_outputs
    assert bijux_phylogenetics.run_maximum_likelihood_tree_inference is run_maximum_likelihood_tree_inference
    assert bijux_phylogenetics.validate_ml_tree_contains_expected_taxa is validate_ml_tree_contains_expected_taxa
    assert bijux_phylogenetics.run_bootstrap_support_estimation is run_bootstrap_support_estimation
    assert bijux_phylogenetics.validate_bootstrap_tree_set is validate_bootstrap_tree_set
    assert bijux_phylogenetics.run_bootstrap_consensus_tree is run_bootstrap_consensus_tree
    assert bijux_phylogenetics.run_fast_tree_inference is run_fast_tree_inference
    assert bijux_phylogenetics.compare_fast_and_ml_trees is compare_fast_and_ml_trees
    assert bijux_phylogenetics.compare_inferred_tree_to_taxon_metadata is compare_inferred_tree_to_taxon_metadata
    assert bijux_phylogenetics.classify_inference_workflow_failure is classify_inference_workflow_failure
    assert bijux_phylogenetics.validate_inference_engine_outputs is validate_inference_engine_outputs
    assert bijux_phylogenetics.render_inference_workflow_report is render_inference_workflow_report
    assert bijux_phylogenetics.prepare_mrbayes_analysis is prepare_mrbayes_analysis
    assert bijux_phylogenetics.run_mrbayes_posterior_inference is run_mrbayes_posterior_inference
    assert bijux_phylogenetics.summarize_mrbayes_posterior_trees is summarize_mrbayes_posterior_trees
    assert bijux_phylogenetics.parse_mrbayes_parameter_traces is parse_mrbayes_parameter_traces
    assert bijux_phylogenetics.compute_mrbayes_effective_sample_sizes is compute_mrbayes_effective_sample_sizes
    assert bijux_phylogenetics.assess_mrbayes_convergence is assess_mrbayes_convergence
    assert bijux_phylogenetics.render_bayesian_posterior_report is render_bayesian_posterior_report
    assert bijux_phylogenetics.validate_fossil_calibration_table is validate_fossil_calibration_table
    assert bijux_phylogenetics.detect_impossible_calibration_constraints is detect_impossible_calibration_constraints
    assert bijux_phylogenetics.validate_tip_dating_metadata is validate_tip_dating_metadata
    assert bijux_phylogenetics.prepare_beast_time_tree_analysis is prepare_beast_time_tree_analysis
    assert bijux_phylogenetics.parse_beast_log is parse_beast_log
    assert bijux_phylogenetics.assess_beast_convergence is assess_beast_convergence
    assert bijux_phylogenetics.validate_beast_posterior_log is validate_beast_posterior_log
    assert bijux_phylogenetics.assess_beast_burnin_sensitivity is assess_beast_burnin_sensitivity
    assert bijux_phylogenetics.assess_beast_chain_mixing is assess_beast_chain_mixing
    assert bijux_phylogenetics.summarize_maximum_clade_credibility_tree is summarize_maximum_clade_credibility_tree
    assert bijux_phylogenetics.thin_posterior_tree_set is thin_posterior_tree_set
    assert bijux_phylogenetics.summarize_posterior_node_ages is summarize_posterior_node_ages
    assert bijux_phylogenetics.compare_bayesian_tree_sets is compare_bayesian_tree_sets
    assert bijux_phylogenetics.compare_independent_bayesian_runs is compare_independent_bayesian_runs
    assert bijux_phylogenetics.compare_posterior_tree_sets_by_prior is compare_posterior_tree_sets_by_prior
    assert bijux_phylogenetics.compare_posterior_tree_sets_by_clock is compare_posterior_tree_sets_by_clock
    assert bijux_phylogenetics.render_bayesian_run_comparison_report is render_bayesian_run_comparison_report
    assert bijux_phylogenetics.render_bayesian_diagnostics_report is render_bayesian_diagnostics_report
    assert bijux_phylogenetics.build_posterior_uncertainty_figure_package is build_posterior_uncertainty_figure_package
    assert bijux_phylogenetics.write_supplementary_bayesian_diagnostics_table is write_supplementary_bayesian_diagnostics_table
    assert bijux_phylogenetics.write_bayesian_methods_summary_text is write_bayesian_methods_summary_text
    assert bijux_phylogenetics.render_calibration_audit_report is render_calibration_audit_report
    assert bijux_phylogenetics.build_bayesian_evidence_package is build_bayesian_evidence_package


def test_simulate_birth_death_trees_returns_requested_tree_and_tip_counts(tmp_path: Path) -> None:
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
        "((Taxon1:0.41605101339611,Taxon4:0.41605101339611):0.455215777552457,"
        "(Taxon2:0.065219140705801,Taxon3:0.065219140705801):0.806047650242766);"
    )


def test_simulate_brownian_traits_generates_one_value_per_tip(tmp_path: Path) -> None:
    report = simulate_brownian_traits(fixture("example_tree.nwk"), seed=7, root_state=1.0, sigma=0.5)
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
    assert [row.node for row in report.node_values if not row.is_tip] == ["A|B|C|D", "A|B", "C|D"]


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
        "taxon\tstate\n"
        "A\twet\n"
        "B\tdry\n"
        "C\twet\n"
        "D\tmixed\n"
    )
    assert [row.node for row in report.node_states if not row.is_tip] == ["A|B|C|D", "A|B", "C|D"]


def test_simulate_dna_alignment_returns_requested_taxa_and_length(tmp_path: Path) -> None:
    report = simulate_dna_alignment(fixture("example_tree.nwk"), sequence_length=8, substitution_rate=1.2, seed=7)
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
        ">A\nACTAACGA\n"
        ">B\nACTAACGA\n"
        ">C\nGCTAGAAA\n"
        ">D\nGCGAAGAA\n"
    )


def test_simulate_protein_alignment_returns_requested_length_and_alphabet() -> None:
    report = simulate_protein_alignment(fixture("example_tree.nwk"), sequence_length=6, substitution_rate=0.8, seed=7)
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
    report = benchmark_alignment_diagnostics(replicates=1, sequence_counts=[4, 8, 16], sequence_length=24)
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
    assert [(row.index, row.tip_count) for row in report.records] == [(1, 4), (2, 4), (3, 4)]


def test_compute_clade_frequency_table_counts_informative_clades() -> None:
    report = compute_clade_frequency_table(fixture("example_tree_set_left.nwk"))
    assert [(row.clade, row.tree_count, row.frequency) for row in report.clade_frequencies] == [
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
    assert bijux_phylogenetics.trim_columns_above_missingness_threshold is trim_columns_above_missingness_threshold
    assert bijux_phylogenetics.trim_alignment is trim_alignment
    assert bijux_phylogenetics.translate_coding_alignment is translate_coding_alignment
    assert bijux_phylogenetics.root_tree_on_outgroup is root_tree_on_outgroup
    assert bijux_phylogenetics.reroot_tree_by_midpoint is reroot_tree_by_midpoint
    assert bijux_phylogenetics.unroot_tree is unroot_tree
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
        ("A|B||C|D", [1, 2], 2, 0.666666666666667, 1, "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);"),
        ("A|C||B|D", [3], 1, 0.333333333333333, 3, "((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);"),
    ]


def test_detect_unstable_taxa_reports_inconsistent_placements() -> None:
    report = detect_unstable_taxa(fixture("example_tree_set_left.nwk"))
    assert [
        (
            row.taxon,
            row.unique_placements,
            row.dominant_frequency,
            row.instability_score,
            [(placement.signature, placement.tree_count, placement.frequency) for placement in row.placements],
        )
        for row in report.taxa
    ] == [
        ("A", 2, 0.666666666666667, 0.333333333333333, [("A|B", 2, 0.666666666666667), ("A|C", 1, 0.333333333333333)]),
        ("B", 2, 0.666666666666667, 0.333333333333333, [("A|B", 2, 0.666666666666667), ("B|D", 1, 0.333333333333333)]),
        ("C", 2, 0.666666666666667, 0.333333333333333, [("C|D", 2, 0.666666666666667), ("A|C", 1, 0.333333333333333)]),
        ("D", 2, 0.666666666666667, 0.333333333333333, [("C|D", 2, 0.666666666666667), ("B|D", 1, 0.333333333333333)]),
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
        ("A|B", 2, 0.666666666666667, 2, 0.333333333333333, "intermediate-support", ["A|C", "B|D"]),
        ("A|C", 1, 0.333333333333333, 2, 0.333333333333333, "intermediate-support", ["A|B", "C|D"]),
        ("B|D", 1, 0.333333333333333, 2, 0.333333333333333, "intermediate-support", ["A|B", "C|D"]),
        ("C|D", 2, 0.666666666666667, 2, 0.333333333333333, "intermediate-support", ["A|C", "B|D"]),
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
    assert report.left_summary.mean_within_set_normalized_robinson_foulds == 0.666666666666667
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


def test_summarize_clade_credibility_conflicts_reports_incompatible_high_frequency_clades() -> None:
    report = summarize_clade_credibility_conflicts(
        fixture("example_tree_set_left.nwk"),
        credibility_threshold=0.3,
    )

    assert report.high_credibility_clade_count == 4
    assert report.conflict_count == 2
    assert report.conflicts[0].combined_frequency == 1.0


def test_summarize_uncertainty_aware_conclusions_separates_robust_uncertain_and_conflicting_clades() -> None:
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


def test_write_uncertainty_tables_emit_clusters_conflicts_and_conclusions(tmp_path: Path) -> None:
    cluster_path = tmp_path / "topology-clusters.tsv"
    conflict_path = tmp_path / "clade-conflicts.tsv"
    conclusion_path = tmp_path / "uncertainty-conclusions.tsv"

    write_topology_cluster_table(cluster_path, cluster_trees_by_topology(fixture("example_tree_set_left.nwk")))
    write_clade_credibility_conflict_table(
        conflict_path,
        summarize_clade_credibility_conflicts(fixture("example_tree_set_left.nwk"), credibility_threshold=0.3),
    )
    write_uncertainty_conclusion_table(
        conclusion_path,
        summarize_uncertainty_aware_conclusions(fixture("example_tree_set_left.nwk"), robust_threshold=0.95, credibility_threshold=0.3),
    )

    assert "rooted_topology_id\ttree_indices" in cluster_path.read_text(encoding="utf-8")
    assert "left_clade\tleft_frequency\tright_clade" in conflict_path.read_text(encoding="utf-8")
    assert "clade\tfrequency\tconclusion" in conclusion_path.read_text(encoding="utf-8")


def test_compare_posterior_tree_sets_reports_clade_deltas_and_cross_set_distance() -> None:
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


def test_compare_bootstrap_and_posterior_uncertainty_reports_conflicting_clades(tmp_path: Path) -> None:
    bootstrap_tree = tmp_path / "bootstrap-support.nwk"
    bootstrap_tree.write_text("((A:0.1,B:0.1)95:0.2,(C:0.1,D:0.1)60:0.2);\n", encoding="utf-8")

    report = compare_bootstrap_and_posterior_uncertainty(
        bootstrap_tree,
        fixture("example_tree_set_right.nwk"),
    )

    rows = {row.clade: row for row in report.rows}
    assert report.posterior_tree_count == 3
    assert rows["A|B"].agreement == "strong_conflict"
    assert rows["C|D"].posterior_frequency == 0.333333333333333


def test_render_tree_set_comparison_report_embeds_tree_set_differences(tmp_path: Path) -> None:
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
    assert [(rename.raw_label, rename.normalized_label) for rename in report.renamed_taxa] == [
        ("Homo sapiens", "Homo_sapiens"),
        ("Mus musculus", "Mus_musculus"),
    ]


def test_taxon_safety_reports_unsafe_labels_and_normalization_collisions(tmp_path: Path) -> None:
    tree = loads_newick(
        "('Homo sapiens':0.1,Homo_sapiens:0.2,'NCBI/123':0.3,'Quoted''Name':0.4,A:0.5);"
    )
    report = inspect_tree_taxa_safety(tree, policy="spaces-to-underscores")
    assert [(entry.raw_label, entry.normalized_label, entry.reasons) for entry in report.unsafe_taxa] == [
        ("Homo sapiens", "Homo_sapiens", ["contains whitespace", "collides with another label after normalization"]),
        ("Homo_sapiens", "Homo_sapiens", ["collides with another label after normalization"]),
        ("NCBI/123", "NCBI/123", ["contains slash characters"]),
        ("Quoted'Name", "Quoted'Name", ["contains quote characters"]),
    ]
    assert [(entry.normalized_label, entry.raw_labels) for entry in report.collisions] == [
        ("Homo_sapiens", ["Homo sapiens", "Homo_sapiens"])
    ]

    mapping_path = tmp_path / "taxon-mapping.tsv"
    write_taxon_mapping(mapping_path, normalize_tree_taxa(tree, policy="spaces-to-underscores")[1].renamed_taxa)
    assert mapping_path.read_text(encoding="utf-8") == (
        "raw_label\tnormalized_label\n"
        "Homo sapiens\tHomo_sapiens\n"
    )


def test_metadata_inspect_reports_taxon_contract() -> None:
    report = inspect_metadata_table(fixture("example_metadata.tsv"))
    assert report.format == "tsv"
    assert report.row_count == 4
    assert report.column_count == 3
    assert report.taxon_column == "taxon"
    assert report.taxa == ["A", "B", "C", "D"]
    assert [(row.name, row.missing_count, row.completeness_fraction) for row in report.column_completeness] == [
        ("taxon", 0, 1.0),
        ("species", 0, 1.0),
        ("location", 0, 1.0),
    ]


def test_join_table_to_taxa_returns_tip_by_tip_metadata_rows() -> None:
    report = join_table_to_taxa(["A", "B", "Z"], fixture("example_metadata.tsv"))
    assert [(row.taxon, row.matched, row.values.get("species", "")) for row in report.joined_rows] == [
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
        inspect_metadata_table(fixture("example_metadata_missing_taxon.csv"), taxon_column="taxon")
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
    assert [(column.name, column.missing_fraction) for column in columns] == [("status", 0.25)]


def test_detect_missing_trait_values_reports_taxon_and_column() -> None:
    report = detect_missing_trait_values(fixture("example_traits_validate.tsv"))
    assert [(item.taxon, item.trait) for item in report.missing_values] == [("C", "status")]


def test_traits_link_reports_mismatch_and_usable_taxa() -> None:
    report = link_tree_to_traits(fixture("example_tree.nwk"), fixture("example_traits.tsv"))
    assert report.tree_taxa == 4
    assert report.trait_taxa == 4
    assert report.linked_taxa == 3
    assert report.usable_taxa == ["A", "B", "C"]
    assert report.missing_from_traits == ["D"]
    assert report.extra_trait_taxa == ["E"]


def test_traits_link_strict_mode_rejects_mismatch() -> None:
    try:
        link_tree_to_traits(fixture("example_tree.nwk"), fixture("example_traits.tsv"), strict=True)
    except MetadataJoinError as error:
        assert error.code == "metadata_join_error"
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


def test_traits_prune_cli_writes_pruned_table_in_tree_order(tmp_path: Path, capsys) -> None:
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
        "taxon\tvalue\n"
        "A\t1.2\n"
        "B\t1.4\n"
        "C\t1.8\n"
    )
    assert payload["data"]["kept_taxa"] == ["A", "B", "C"]
    assert payload["data"]["removed_taxa"] == ["E"]


def test_traits_missing_cli_reports_taxon_and_column(capsys) -> None:
    exit_code = main(["traits", "missing", str(fixture("example_traits_validate.tsv")), "--json"])
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
    tree, report = prune_tree_to_taxa(fixture("example_tree.nwk"), fixture("example_traits.tsv"))
    assert tree.tip_names == ["A", "B", "C"]
    assert dumps_newick(tree) == "((A:0.1,B:0.1):0.2,C:0.3);"
    assert report.kept_taxa == ["A", "B", "C"]
    assert report.removed_taxa == ["D"]
    assert [(row.taxon, row.reason) for row in report.removed_taxa_with_reasons] == [("D", "absent_from_keep_table")]
    assert report.summary.removed_taxa == ["D"]


def test_prune_tree_to_requested_taxa_reports_absent_requests() -> None:
    tree, report = prune_tree_to_requested_taxa(
        fixture("example_tree.nwk"),
        ["A", "C", "Z"],
    )
    assert tree.tip_names == ["A", "C"]
    assert dumps_newick(tree) == "(A:0.3,C:0.3);"
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


def test_drop_tree_taxa_excludes_exact_requested_tips() -> None:
    tree, report = drop_tree_taxa(
        fixture("example_tree.nwk"),
        ["B", "D", "Z"],
    )
    assert tree.tip_names == ["A", "C"]
    assert dumps_newick(tree) == "(A:0.3,C:0.3);"
    assert report.requested_taxa == ["B", "D", "Z"]
    assert report.kept_taxa == ["A", "C"]
    assert report.removed_taxa == ["B", "D"]
    assert report.absent_requested_taxa == ["Z"]
    assert [(row.taxon, row.reason) for row in report.removed_taxa_with_reasons] == [
        ("B", "excluded_explicitly"),
        ("D", "excluded_explicitly"),
    ]


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


def test_prune_cli_accepts_explicit_taxon_exclusion_lists(tmp_path: Path, capsys) -> None:
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


def test_collapse_branches_below_length_turns_short_internal_edges_into_polytomies() -> None:
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
    assert [len(child.children) if child.children else 1 for child in tree.root.children] == [3, 2]
    assert report.strategy == "ladderize"
    assert report.tip_order == ["X", "Y", "Z", "A", "B"]
    assert report.rooted_topology_preserved is True
    assert report.unrooted_topology_preserved is True


def test_sort_tree_tips_alphabetically_preserves_topology_with_stable_tip_order() -> None:
    tree, report = sort_tree_tips_alphabetically(fixture("example_tree_ordering.nwk"))
    assert tree.tip_names == ["A", "B", "X", "Y", "Z"]
    assert [len(child.children) if child.children else 1 for child in tree.root.children] == [2, 3]
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
    assert dumps_newick(tree) == "(((A:0.2,B:0.2):0.7,C:0.1):0.1,D:0);"
    assert report.summary.retained_taxa == ["A", "B", "C", "D"]
    assert report.summary.removed_taxa == []


def test_reroot_tree_by_midpoint_preserves_taxa_and_branch_lengths() -> None:
    tree, report = reroot_tree_by_midpoint(fixture("example_tree_rootable.nwk"))
    assert sorted(tree.tip_names) == ["A", "B", "C", "D"]
    assert report.strategy == "midpoint"
    assert report.tip_order == ["C", "D", "B", "A"]
    assert dumps_newick(tree) == "((A:0.2,B:0.2):0.3,(C:0.1,D:0.1):0.4);"
    assert report.summary.branch_lengths_affected != []


def test_unroot_tree_converts_rooted_binary_tree_into_trifurcation() -> None:
    tree, report = unroot_tree(fixture("example_tree_rootable.nwk"))
    assert sorted(tree.tip_names) == ["A", "B", "C", "D"]
    assert len(tree.root.children) == 3
    assert report.strategy == "unroot"
    assert dumps_newick(tree) == "(A:0.5,B:0.5,(C:0.1,D:0.1):0.4);"
    assert report.summary.nodes_changed != []


def test_prune_alignment_to_tree_keeps_exact_tree_taxa() -> None:
    records, report = prune_alignment_to_tree(
        fixture("example_alignment_extra_taxon.fasta"),
        fixture("example_tree.nwk"),
    )
    assert [record.identifier for record in records] == ["A", "B", "C", "D"]
    assert report.original_sequence_count == 5
    assert report.kept_ids == ["A", "B", "C", "D"]
    assert report.removed_ids == ["E"]
    assert [(row.taxon, row.reason) for row in report.removed_ids_with_reasons] == [("E", "absent_from_tree")]


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
    assert [(row.taxon, row.reason) for row in report.removed_taxa_with_reasons] == [("D", "absent_from_alignment")]


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
    assert report.nucleotide_composition == {"A": 0.3125, "C": 0.25, "G": 0.25, "T": 0.1875}
    assert report.whole_alignment_gc_content == 0.5
    assert [(row.identifier, row.missing_fraction) for row in report.per_sequence_missingness] == [
        ("A", 0.0),
        ("B", 0.0),
        ("C", 0.0),
        ("D", 0.0),
    ]
    assert [(row.identifier, row.gc_fraction) for row in report.per_sequence_gc_content] == [
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
    assert [(row.identifier, row.position, row.character) for row in invalid] == [("A", 5, "Z")]


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
    duplicates = detect_identical_duplicate_sequences(fixture("example_alignment_duplicates.fasta"))
    near_duplicates = detect_near_duplicate_sequences(
        fixture("example_alignment_duplicates.fasta"),
        identity_threshold=0.875,
    )
    assert [(group.identifiers, group.sequence) for group in duplicates] == [
        (["A", "B"], "ACTGACTG")
    ]
    assert [(pair.left_identifier, pair.right_identifier, pair.identity, pair.comparable_sites) for pair in near_duplicates] == [
        ("A", "C", 0.875, 8),
        ("A", "D", 0.875, 8),
        ("B", "C", 0.875, 8),
        ("B", "D", 0.875, 8),
        ("C", "D", 0.875, 8),
    ]


def test_alignment_quality_report_collects_composition_duplicates_and_warnings() -> None:
    report = build_alignment_quality_report(fixture("example_alignment_duplicates.fasta"))
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
    assert [(row.identifier, row.missing_fraction) for row in report.per_sequence_missingness] == [
        ("A", 2 / 6),
        ("B", 2 / 6),
        ("C", 0.0),
    ]
    assert [(row.position, row.missing_fraction) for row in report.per_site_missingness] == [
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
    assert [(row.position, row.missing_fraction) for row in report.per_site_missingness] == [
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
    assert [(row.identifier, row.gap_fraction, row.missing_fraction, row.ambiguity_fraction) for row in report.per_sequence_uncertainty] == [
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
    assert [(window.start, window.end, window.gap_fraction, window.missing_fraction, window.ambiguity_fraction) for window in windows] == [
        (1, 4, 0.25, 0.5, 0.25),
        (5, 8, 0.0, 0.0, 0.0),
    ]
    assert [(row.start, row.end, row.kind) for row in over_aligned] == [(1, 4, "over_aligned")]
    assert [(row.start, row.end, row.kind) for row in under_aligned] == [(5, 8, "under_aligned")]


def test_alignment_readiness_reports_method_specific_decisions() -> None:
    raw = summarize_alignment_readiness(fixture("example_sequences_raw.fasta"))
    coding = summarize_alignment_readiness(fixture("example_alignment_coding.fasta"))
    raw_methods = {row.analysis: row for row in raw.methods}
    coding_methods = {row.analysis: row for row in coding.methods}

    assert raw_methods["distance"].ready is False
    assert raw_methods["maximum_likelihood"].blockers == ["input sequences are not yet aligned"]
    assert coding_methods["distance"].ready is True
    assert coding_methods["coding"].ready is False
    assert "one or more sequences contain premature stop codons" in coding_methods["coding"].blockers
    assert "one or more sequences contain partial codons after gaps and missing data are removed" in coding_methods["coding"].blockers


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
    assert [(row.identifier, row.coding_like, row.premature_stop_count) for row in diagnostics.coding_behaviors] == [
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
    assert [record.sequence for record in records] == ["ACTGACTG", "ACTGACTG", "ACTGACTA"]
    assert report.profile.name == "moderate"
    assert report.trim.trimmed_alignment_length == 8
    assert [(row.position, row.reason) for row in report.trim.removed_columns] == [
        (5, "missingness-threshold"),
        (6, "missingness-threshold"),
        (7, "missingness-threshold"),
        (8, "missingness-threshold"),
    ]
    assert [(row.identifier, row.reason) for row in report.trim.removed_sequences] == [("D", "missingness-threshold")]
    assert report.comparison.left_alignment_length == 12
    assert report.comparison.right_alignment_length == 8
    assert any(group.column == "region" and group.value == "island" and group.removed_fraction == 1.0 for group in report.group_retention)
    assert "cleaning removed most taxa from one or more metadata or trait groups" in report.warnings


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
    assert report.left_parsimony_informative_site_count >= report.right_parsimony_informative_site_count


def test_alignment_quality_report_includes_transparent_score_components() -> None:
    report = build_alignment_quality_report(fixture("example_alignment_filtering.fasta"))
    assert set(report.quality_components) == {
        "composition_outliers",
        "duplicates",
        "gap_burden",
        "informative_density",
        "missingness",
    }
    assert 0.0 <= report.quality_score <= 100.0


def test_alignment_forensic_report_integrates_alignment_risks() -> None:
    report = build_alignment_forensic_report(fixture("example_alignment_coding.fasta"))
    assert report.safe_for_distance_analysis is True
    assert report.safe_for_coding_analysis is False
    assert "alignment mixes coding-like and noncoding-like sequence behavior" in report.warnings


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
    assert any(row.analysis == "distance" and row.decision == "risky" for row in report.analysis_decisions)
    assert any(level.level == "publication_ready" and level.decision == "risky" for level in report.readiness_levels)


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
        row.surface == "metadata" and row.taxon == "C" and row.expected_index == 3 and row.observed_index == 1
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
        row.surface == "metadata" and row.group_column == "region" and row.group == "island" and row.removed_fraction == 1.0
        for row in report.group_imbalance_warnings
    )
    exclusion_by_taxon = {row.taxon: row for row in report.exclusion_table.rows}
    assert exclusion_by_taxon["D"].causes == ["absent_from_alignment", "absent_from_traits"]
    assert exclusion_by_taxon["D"].first_failed_surface == "alignment"
    assert report.mismatch_report.rows
    assert report.risk_score.total_score > 0.0
    assert any(component.component == "alignment" for component in report.risk_score.components)
    assert report.minimal_fix_plan.recommendations
    assert any(item.section == "dataset_risk" for item in report.reviewer_checklist.items)


def test_alignment_inspect_rejects_unequal_lengths() -> None:
    try:
        summarise_fasta(fixture("example_alignment_invalid_lengths.fasta"))
    except InvalidAlignmentError as error:
        assert error.code == "invalid_alignment_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected InvalidAlignmentError")


def test_alignment_link_reports_exact_mismatch() -> None:
    report = link_alignment_to_tree(fixture("example_tree.nwk"), fixture("example_alignment.fasta"))
    assert report.tree_taxa == 4
    assert report.alignment_ids == 4
    assert report.linked_taxa == 4
    assert report.missing_from_alignment == []
    assert report.extra_alignment_ids == []


def test_alignment_link_strict_mode_rejects_mismatch() -> None:
    try:
        link_alignment_to_tree(fixture("example_tree.nwk"), fixture("example_alignment_mismatch.fasta"), strict=True)
    except AlignmentTaxonMismatchError as error:
        assert error.code == "alignment_taxon_mismatch_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected AlignmentTaxonMismatchError")


def test_write_fasta_alignment_preserves_record_order_and_sequences(tmp_path: Path) -> None:
    records = load_fasta_alignment(fixture("example_alignment.fasta"))
    output = tmp_path / "alignment.fasta"
    write_fasta_alignment(output, records)
    assert output.read_text(encoding="utf-8") == (
        ">A\nACTGACTG\n"
        ">B\nACTGACTA\n"
        ">C\nACTGACGG\n"
        ">D\nACTGACGA\n"
    )
    assert load_fasta_alignment(output) == records


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
    records, report = remove_all_gap_columns(fixture("example_alignment_site_missingness.fasta"))
    assert [record.sequence for record in records] == ["AN?T", "CN?N", "GN?A", "TN?N"]
    assert report.original_alignment_length == 5
    assert report.trimmed_alignment_length == 4
    assert [(column.position, column.reason) for column in report.removed_columns] == [(2, "all-gap")]


def test_alignment_removes_all_missing_columns() -> None:
    records, report = remove_all_missing_columns(fixture("example_alignment_site_missingness.fasta"))
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
    assert [(row.identifier, row.missing_fraction, row.reason) for row in report.removed_sequences] == [
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
    assert [(record.identifier, record.sequence) for record in records] == [("B", "ACG")]
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
    report = compute_pairwise_sequence_identity_matrix(fixture("example_alignment_duplicates.fasta"))
    assert report.identifiers == ["A", "B", "C", "D"]
    assert [(pair.left_identifier, pair.right_identifier, pair.identity, pair.comparable_sites) for pair in report.pairs] == [
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
    assert diagnostics.sequence_count == 4
    assert diagnostics.alignment_length_multiple_of_three is True
    assert [(row.identifier, row.comparable_length, row.remainder) for row in diagnostics.frameshift_like_sequences] == [
        ("C", 8, 2)
    ]
    assert [(row.identifier, row.comparable_length, row.trailing_bases) for row in diagnostics.partial_codon_sequences] == [
        ("C", 8, 2)
    ]
    assert [(row.identifier, row.codon_index, row.nucleotide_start, row.codon, row.terminal) for row in diagnostics.stop_codons] == [
        ("A", 3, 7, "TAA", True),
        ("D", 2, 4, "TAG", False),
    ]


def test_translate_coding_alignment_emits_amino_acid_records() -> None:
    records, report = translate_coding_alignment(fixture("example_alignment_coding.fasta"))
    assert [(record.identifier, record.sequence) for record in records] == [
        ("A", "ME*"),
        ("B", "MEW"),
        ("C", "MXW"),
        ("D", "M*W"),
    ]
    assert report.translated_sequence_count == 4
    assert report.source_alignment_length == 9
    assert report.translated_alignment_length == 3
    assert report.stop_codon_count == 2
    assert report.frameshift_like_sequence_count == 1


def test_cli_alignment_trim_writes_trimmed_fasta_and_report(tmp_path: Path, capsys) -> None:
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
    report = compute_pairwise_genetic_distance_matrix(fixture("example_alignment_distance.fasta"))
    assert report.model == "p-distance"
    assert report.gap_handling == "pairwise-deletion"
    assert report.identifiers == ["A", "B", "C", "D"]
    assert [(pair.left_identifier, pair.right_identifier, pair.distance, pair.comparable_sites) for pair in report.pairs] == [
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
    report = compute_pairwise_genetic_distance_matrix(fixture("example_alignment_distance_gaps.fasta"))
    assert [(pair.left_identifier, pair.right_identifier, pair.distance, pair.comparable_sites) for pair in report.pairs] == [
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
    assert [(pair.left_identifier, pair.right_identifier, pair.distance, pair.comparable_sites) for pair in report.pairs] == [
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
    assert [(pair.left_identifier, pair.right_identifier, pair.distance, pair.comparable_sites) for pair in report.pairs] == [
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


def test_compute_pairwise_genetic_distance_matrix_marks_saturated_jukes_cantor_pairs_undefined() -> None:
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance_saturated.fasta"),
        model="jukes-cantor",
    )
    assert [(pair.left_identifier, pair.right_identifier, pair.distance, pair.comparable_sites) for pair in report.pairs] == [
        ("A", "A", 0.0, 4),
        ("A", "B", None, 4),
        ("A", "C", 0.304098831081123, 4),
        ("B", "B", 0.0, 4),
        ("B", "C", None, 4),
        ("C", "C", 0.0, 4),
    ]


def test_load_imported_distance_matrix_reads_exported_long_form() -> None:
    entries = load_imported_distance_matrix(fixture("example_distance_matrix.tsv"))
    assert [(entry.left_identifier, entry.right_identifier, entry.distance, entry.comparable_sites) for entry in entries[:4]] == [
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
    report = validate_imported_distance_matrix(fixture("example_distance_matrix_nonmetric.tsv"))
    assert [(row.left_identifier, row.middle_identifier, row.right_identifier, row.direct_distance, row.indirect_distance) for row in report.nonmetric_observations] == [
        ("A", "B", "C", 5.0, 2.0),
    ]
    assert report.warnings == ["distance matrix violates triangle inequality for one or more taxon triples"]


def test_validate_imported_distance_matrix_detects_asymmetry() -> None:
    report = validate_imported_distance_matrix(fixture("example_distance_matrix_asymmetric.tsv"))
    assert report.complete is True
    assert report.symmetric is False
    assert [(row.left_identifier, row.right_identifier, row.left_to_right_distance, row.right_to_left_distance) for row in report.asymmetric_pairs] == [
        ("A", "B", 1.0, 2.0),
    ]


def test_build_tree_from_imported_distance_matrix_constructs_neighbor_joining_tree() -> None:
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
    assert len(limitations) == 4
    assert limitations[0].startswith("distance methods collapse")


def test_render_distance_report_embeds_limitations_and_validation(tmp_path: Path) -> None:
    output_path = tmp_path / "distance-report.html"
    result = render_distance_report(
        out_path=output_path,
        matrix_path=fixture("example_distance_matrix.tsv"),
    )
    html = output_path.read_text(encoding="utf-8")
    assert result.source_kind == "imported-distance-matrix"
    assert "distance-method-limitations" in html
    assert "imported-distance-matrix-validation" in html
    assert "neighbor-joining-tree" in html


def test_render_tree_uncertainty_report_embeds_consensus_and_instability_sections(tmp_path: Path) -> None:
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
    assert "pairwise-tree-distances" in html
    assert "topology-multimodality" in html
    assert "clade-credibility-conflicts" in html
    assert "unstable-clades" in html


def test_build_distance_tree_constructs_neighbor_joining_tree() -> None:
    tree, report = build_distance_tree(
        fixture("example_alignment_distance.fasta"),
        method="neighbor-joining",
    )
    assert dumps_newick(tree) == "((A:0.0625,B:0.0625)Inner1:0.4375,C:0.0625,D:0.0625)Inner2;"
    assert report.method == "neighbor-joining"
    assert report.taxon_count == 4


def test_build_distance_tree_constructs_upgma_tree() -> None:
    tree, report = build_distance_tree(
        fixture("example_alignment_distance.fasta"),
        method="upgma",
    )
    assert dumps_newick(tree) == "((A:0.0625,B:0.0625)Inner2:0.21875,(C:0.0625,D:0.0625)Inner1:0.21875)Inner3;"
    assert report.method == "upgma"
    assert report.taxon_count == 4


def test_compare_distance_tree_topologies_reports_rooting_difference() -> None:
    report = compare_distance_tree_topologies(fixture("example_alignment_distance.fasta"))
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


def test_cli_alignment_compare_distance_trees_reports_rooting_difference(capsys) -> None:
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
    exit_code = main(["distance", "validate", str(fixture("example_distance_matrix_nonmetric.tsv")), "--json"])
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
    assert output_path.read_text(encoding="utf-8") == "((A:0.0625,B:0.0625)Inner1:0.21875,C:0.28125)Inner2;\n"
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
    exit_code = main(["distance", "explain", str(fixture("example_distance_matrix.tsv")), "--json"])
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
    assert payload["data"]["mean_between_set_normalized_robinson_foulds"] == 0.777777777777778


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


def test_cli_benchmark_tree_validation_reports_observations(capsys) -> None:
    exit_code = main(["benchmark", "tree-validation", "--replicates", "1", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["observation_count"] == 3
    assert payload["data"]["replicates"] == 1


def test_cli_alignment_coding_reports_frameshifts_and_stops(capsys) -> None:
    exit_code = main(["alignment", "coding", str(fixture("example_alignment_coding.fasta")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["frameshift_like_sequence_count"] == 1
    assert payload["metrics"]["stop_codon_count"] == 2


def test_cli_alignment_translate_writes_amino_acid_alignment(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "translated.fasta"
    exit_code = main(
        [
            "alignment",
            "translate",
            str(fixture("example_alignment_coding.fasta")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        ">A\nME*\n"
        ">B\nMEW\n"
        ">C\nMXW\n"
        ">D\nM*W\n"
    )
    assert payload["metrics"]["translated_sequence_count"] == 4
    assert payload["metrics"]["stop_codon_count"] == 2


def test_cli_topology_root_outgroup_writes_rooted_tree(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "rooted.nwk"
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
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == "(((A:0.2,B:0.2):0.7,C:0.1):0.1,D:0);\n"
    assert payload["metrics"]["matched_taxa"] == 1
    assert payload["metrics"]["absent_taxa"] == 1


def test_cli_topology_reroot_midpoint_writes_tree(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "midpoint.nwk"
    exit_code = main(
        [
            "topology",
            "reroot-midpoint",
            str(fixture("example_tree_rootable.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == "((A:0.2,B:0.2):0.3,(C:0.1,D:0.1):0.4);\n"
    assert payload["data"]["strategy"] == "midpoint"
    assert payload["metrics"]["tip_count"] == 4


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
    assert output_path.read_text(encoding="utf-8") == "(A:0.5,B:0.5,(C:0.1,D:0.1):0.4);\n"
    assert payload["data"]["strategy"] == "unroot"
    assert payload["metrics"]["tip_count"] == 4


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
        validate_tree_path(fixture("example_tree_ladderized.nwk"), require_ultrametric=True)
    except NonUltrametricTreeError as error:
        assert error.code == "non_ultrametric_tree_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected NonUltrametricTreeError")


def test_validate_tree_path_warns_for_zero_length_branches() -> None:
    report = validate_tree_path(fixture("example_tree_zero_lengths.nwk"))
    assert report.zero_length_branches == 3
    assert "tree contains zero-length branches" in report.warnings


def test_validate_tree_path_localizes_missing_internal_and_terminal_branch_lengths() -> None:
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
    missing_internal = inspect_tree_path(fixture("example_tree_missing_internal_length.nwk"))
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
    assert [(row.node, row.branch_length, row.branch_type) for row in report.long_branch_outliers] == [
        ("A", 1.0, "terminal")
    ]
    assert report.short_branch_outliers == []
    assert report.star_like is False
    assert report.tree_quality_score == 85.0
    assert [warning.code for warning in report.tree_quality_warnings] == ["long_branches"]


def test_inspect_tree_path_detects_internal_long_and_short_branch_outliers() -> None:
    long_report = inspect_tree_path(fixture("example_tree_internal_long_branch.nwk"))
    short_report = inspect_tree_path(fixture("example_tree_short_branch.nwk"))
    assert [(row.node, row.branch_length, row.branch_type) for row in long_report.long_branch_outliers] == [
        ("A|B", 1.0, "internal")
    ]
    assert long_report.long_branch_taxa == []
    assert long_report.short_branch_outliers == []
    assert [(row.node, row.branch_length, row.branch_type) for row in short_report.short_branch_outliers] == [
        ("B", 0.001, "terminal")
    ]
    assert [warning.code for warning in short_report.tree_quality_warnings] == ["short_branches"]


def test_inspect_tree_path_classifies_internal_support_and_name_labels() -> None:
    support = inspect_tree_path(fixture("example_tree_support_mixed.nwk"))
    names = inspect_tree_path(fixture("example_tree_named_clades.nwk"))
    assert [(row.node, row.label, row.numeric_value) for row in support.likely_support_labels] == [
        ("A|B", "0.95", 0.95),
        ("A|B|C|D", "99", 99.0),
        ("C|D", "88", 88.0),
    ]
    assert support.likely_named_internal_labels == []
    assert [(row.node, row.label, row.interpretation) for row in names.likely_named_internal_labels] == [
        ("A|B", "Mammals", "name"),
        ("A|B|C|D", "Root", "name"),
        ("C|D", "Birds", "name"),
    ]
    assert names.likely_support_labels == []


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
    assert [warning.code for warning in mixed.tree_quality_warnings] == ["mixed_support_scales"]


def test_standardize_support_labels_normalizes_fraction_and_percentage_scales() -> None:
    rows = standardize_support_labels(fixture("example_tree_support_mixed.nwk"))
    assert [(row.node, row.raw_value, row.scale, row.support_fraction, row.support_percent) for row in rows] == [
        ("A|B", 0.95, "fraction", 0.95, 95.0),
        ("A|B|C|D", 99.0, "percentage", 0.99, 99.0),
        ("C|D", 88.0, "percentage", 0.88, 88.0),
    ]
    assert [(row.node, row.normalized_probability, row.confidence_of_inference) for row in rows] == [
        ("A|B", 0.95, "medium"),
        ("A|B|C|D", 0.99, "medium"),
        ("C|D", 0.88, "medium"),
    ]


def test_validate_tree_roundtrip_preserves_tree_structure_across_formats() -> None:
    nexus_report = validate_tree_roundtrip(fixture("example_tree.nwk"), target_format="nexus")
    phyloxml_report = validate_tree_roundtrip(fixture("example_tree.nwk"), target_format="phyloxml")
    assert nexus_report.preserved_taxa is True
    assert nexus_report.preserved_topology is True
    assert nexus_report.preserved_branch_lengths is True
    assert nexus_report.preserved_support_labels is True
    assert phyloxml_report.preserved_taxa is True
    assert phyloxml_report.preserved_topology is True


def test_validate_tree_roundtrip_reports_semantic_loss_for_nexus_and_phyloxml() -> None:
    nexus_report = validate_tree_roundtrip(fixture("example_tree.nex"), target_format="newick")
    phyloxml_report = validate_tree_roundtrip(fixture("example_tree_annotated.phyloxml"), target_format="newick")
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
    identity_report = inspect_tree_taxon_identity(load_phyloxml(fixture("example_tree_annotated.phyloxml")))
    assert identity_report.case_collisions == []

    ambiguous = inspect_tree_taxon_identity(loads_newick(fixture("example_tree_identity.nwk").read_text(encoding="utf-8")))
    assert {(row.left_label, row.right_label) for row in ambiguous.underscore_space_collisions} == {
        ("Homo sapiens", "Homo_sapiens"),
        ("Homo sapiens", "homo sapiens"),
        ("Homo_sapiens", "homo sapiens"),
    }
    assert [(row.left_label, row.right_label) for row in ambiguous.case_collisions] == [
        ("Homo sapiens", "homo sapiens"),
    ]
    assert ("Homo sapiens", "Hoomo sapiens") in [
        (row.left_label, row.right_label) for row in ambiguous.suspicious_near_duplicates
    ]


def test_inspect_branch_length_units_reports_time_and_substitution_metadata() -> None:
    time_report = inspect_branch_length_units(fixture("example_metadata_branch_units_time.tsv"))
    substitution_report = inspect_branch_length_units(fixture("example_metadata_branch_units_substitution.tsv"))
    conflict_report = inspect_branch_length_units(fixture("example_metadata_branch_units_conflict.tsv"))
    assert (time_report.declared_unit, time_report.compatible_with_time_tree, time_report.compatible_with_substitution_tree) == (
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
    assert [warning.code for warning in report.tree_quality_warnings] == ["polytomies", "star_like"]


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


def test_validate_tree_path_reports_validity_decision_contexts_and_repair_guidance() -> None:
    report = validate_tree_path(fixture("example_tree_partial_lengths.nwk"), allow_duplicates=True)
    assert report.syntax_valid is True
    assert report.biologically_safe is False
    assert report.validity_decision == "invalid"
    assert report.root_state_confidence.classification == "apparently_rooted"
    assert {context.context: context.allowed for context in report.branch_length_contexts} == {
        "topology_only": True,
        "substitution_tree": False,
        "time_tree": False,
        "comparative_methods": False,
    }
    assert [item.issue_code for item in report.branch_length_repair_suggestions] == ["partial_branch_lengths"]


def test_forensic_tree_path_reports_unsafe_external_labels_and_downstream_safety() -> None:
    report = forensic_tree_path(fixture("example_tree_labels.nwk"))
    assert report.validity_decision == "valid_with_warnings"
    assert report.safe_for_topology_comparison is True
    assert report.safe_for_time_tree_analysis is False
    assert report.safe_for_comparative_methods is False
    assert report.unsafe_external_labels[0].engines == ["iqtree", "raxml", "mrbayes", "beast", "r", "shell"]
    assert any(item.raw_label == "Homo sapiens" for item in report.unsafe_external_labels)


def test_newick_loader_raises_invalid_branch_length_error() -> None:
    try:
        loads_newick("((A:abc,B:0.2):0.3,C:0.4);")
    except InvalidBranchLengthError as error:
        assert error.code == "invalid_branch_length_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected InvalidBranchLengthError")


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
    exit_code = main(["validate", str(fixture("example_tree_duplicate.nwk")), "--allow-duplicates", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["data"]["duplicate_taxa"] == ["A"]


def test_validate_cli_strict_mode_rejects_unnamed_tips(capsys) -> None:
    exit_code = main(["validate", str(fixture("example_tree_unnamed_tip.nwk")), "--strict", "--json"])
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
        ["validate", str(fixture("example_tree_negative_length.nwk")), "--allow-negative-branches", "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["data"]["negative_branch_lengths"] == 1


def test_validate_cli_can_require_rooted_and_ultrametric_typed_errors(capsys) -> None:
    exit_code = main(["validate", str(fixture("example_tree_unrooted.nwk")), "--require-rooted", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["errors"][0]["code"] == UnrootedTreeError.code

    exit_code = main(["validate", str(fixture("example_tree_ladderized.nwk")), "--require-ultrametric", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["errors"][0]["code"] == NonUltrametricTreeError.code


def test_cli_inspect_accepts_explicit_tree_format(capsys) -> None:
    exit_code = main(["inspect", str(fixture("example_tree.nex")), "--format", "nexus", "--json"])
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


def test_cli_validate_and_inspect_surface_structural_and_support_diagnostics(capsys) -> None:
    validate_exit = main(["validate", str(fixture("example_tree_singleton.nwk")), "--json"])
    validate_payload = json.loads(capsys.readouterr().out)
    assert validate_exit == 0
    assert validate_payload["metrics"]["singleton_internal_node_count"] == 1

    inspect_exit = main(["inspect", str(fixture("example_tree_support_mixed.nwk")), "--json"])
    inspect_payload = json.loads(capsys.readouterr().out)
    assert inspect_exit == 0
    assert inspect_payload["metrics"]["likely_support_label_count"] == 3
    assert inspect_payload["metrics"]["suspicious_support_range_count"] == 0
    assert inspect_payload["data"]["mixed_support_scales"] is True


def test_cli_normalize_writes_canonical_newick(tmp_path: Path, capsys) -> None:
    output = tmp_path / "normalized.nwk"
    exit_code = main(
        ["normalize", str(fixture("example_tree.nex")), "--format", "nexus", "--out", str(output), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["outputs"] == [str(output)]
    assert output.read_text(encoding="utf-8").strip() == "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);"


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
    assert output.read_text(encoding="utf-8").strip() == "(A.B-1:0.3,Homo_sapiens:0.1,Mus_musculus:0.2);"
    assert mapping.read_text(encoding="utf-8") == (
        "raw_label\tnormalized_label\n"
        "Homo sapiens\tHomo_sapiens\n"
        "Mus musculus\tMus_musculus\n"
    )


def test_compare_tree_paths_reports_nonzero_distance() -> None:
    report = compare_tree_paths(fixture("example_tree.nwk"), fixture("example_tree_alt.nwk"))
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.robinson_foulds_distance > 0


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
    report = compare_tree_paths(fixture("example_tree.nwk"), fixture("example_tree.nwk"))
    assert report.topology_equal is True
    assert report.same_unrooted_topology is True
    assert report.same_taxa_different_rooting is False
    assert report.same_topology_different_branch_lengths is False
    assert report.robinson_foulds_distance == 0


def test_compare_tree_paths_reports_different_topology_boolean() -> None:
    report = compare_tree_paths(fixture("example_tree.nwk"), fixture("example_tree_topology_diff.nwk"))
    assert report.topology_equal is False
    assert report.same_unrooted_topology is False
    assert report.same_taxa_different_rooting is False


def test_compare_clade_sets_reports_shared_and_unique_clades() -> None:
    report = compare_clade_sets(fixture("example_tree.nwk"), fixture("example_tree_alt.nwk"))
    assert report.shared_clades == ["A|B"]
    assert report.left_only_clades == ["C|D"]
    assert report.right_only_clades == ["A|B|C"]


def test_compare_tree_paths_detects_same_topology_with_different_branch_lengths() -> None:
    report = compare_tree_paths(fixture("example_tree.nwk"), fixture("example_tree_branch_lengths_right.nwk"))
    assert report.topology_equal is True
    assert report.same_topology_different_branch_lengths is True


def test_compare_tree_paths_detects_same_taxa_with_different_rooting() -> None:
    report = compare_tree_paths(fixture("example_tree.nwk"), fixture("example_tree_rooting_diff.nwk"))
    assert report.topology_equal is False
    assert report.same_unrooted_topology is True
    assert report.same_taxa_different_rooting is True


def test_detect_clade_changes_reports_lost_and_gained_sets() -> None:
    report = detect_clade_changes(fixture("example_tree.nwk"), fixture("example_tree_alt.nwk"))
    assert report.lost_clades == ["C|D"]
    assert report.gained_clades == ["A|B|C"]


def test_compare_support_values_pairs_shared_clades() -> None:
    report = compare_support_values(fixture("example_tree_support_left.nwk"), fixture("example_tree_support_right.nwk"))
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert [(row.split_id, row.left_support, row.right_support) for row in report.shared_clades] == [
        ("A|B", 95.0, 90.0),
        ("C|D", 88.0, 85.0),
    ]


def test_compare_branch_lengths_reports_delta_ratio_and_missing_lengths() -> None:
    scaled = compare_branch_lengths(fixture("example_tree.nwk"), fixture("example_tree_branch_lengths_right.nwk"))
    missing = compare_branch_lengths(fixture("example_tree.nwk"), fixture("example_tree_branch_lengths_missing.nwk"))
    assert [(row.split_id, row.delta, row.ratio) for row in scaled.shared_splits] == [
        ("A|B", 0.2, 2.0),
        ("C|D", 0.1, 2.0),
    ]
    assert [(row.split_id, row.left_length, row.right_length, row.delta, row.ratio) for row in missing.shared_splits] == [
        ("A|B", 0.2, None, None, None),
        ("C|D", 0.1, 0.2, 0.1, 2.0),
    ]


def test_write_tree_comparison_table_writes_one_row_per_compared_split(tmp_path: Path) -> None:
    output = tmp_path / "comparison.tsv"
    write_tree_comparison_table(output, fixture("example_tree.nwk"), fixture("example_tree_alt.nwk"))
    assert output.read_text(encoding="utf-8") == (
        "split_id\tcomparison_status\tshared_clade\tleft_support\tright_support\tleft_length\tright_length\tlength_delta\tlength_ratio\n"
        "A|B\tshared\ttrue\t\t\t0.2\t0.1\t-0.1\t0.5\n"
        "A|B|C\tright_only\tfalse\t\t\t\t\t\t\n"
        "C|D\tleft_only\tfalse\t\t\t\t\t\t\n"
    )


def test_build_tree_comparison_report_writes_html_with_checksums(tmp_path: Path) -> None:
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
    result = render_tree_svg(fixture("example_tree.nwk"), out_path=output, layout="phylogram")
    svg = output.read_text(encoding="utf-8")
    assert result.layout == "phylogram"
    assert result.has_scale_bar is True
    assert 'class="scale-bar"' in svg
    assert 'class="scale-label"' in svg


def test_render_tree_svg_can_render_circular_layout(tmp_path: Path) -> None:
    output = tmp_path / "circular.svg"
    result = render_tree_svg(fixture("example_tree.nwk"), out_path=output, layout="circular")
    svg = output.read_text(encoding="utf-8")
    assert result.layout == "circular"
    assert result.has_scale_bar is False
    assert "<path d=\"M " in svg
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


def test_render_tree_svg_can_render_categorical_tip_traits(tmp_path: Path) -> None:
    output = tmp_path / "categorical.svg"
    result = render_tree_svg(
        fixture("example_tree.nwk"),
        out_path=output,
        categorical_traits={"A": "Sweden", "B": "Norway", "C": "Denmark", "D": "Finland"},
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
            AnnotationStrip("species", {"A": "Alpha species", "B": "Beta species", "C": "Gamma species", "D": "Delta species"}),
            AnnotationStrip("location", {"A": "Sweden", "B": "Norway", "C": "Denmark", "D": "Finland"}),
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
            AnnotationStrip("temperature", {"A": "1.2", "B": "1.4", "C": "1.8", "D": "2.0"}),
            AnnotationStrip("status", {"A": "high", "B": "medium", "C": "medium", "D": "low"}),
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
        categorical_traits={"A": "Sweden", "B": "Norway", "C": "Denmark", "D": "Finland"},
        metadata_strips=[AnnotationStrip("location", {"A": "Sweden", "B": "Norway", "C": "Denmark", "D": "Finland"})],
        heatmap_columns=[AnnotationStrip("temperature", {"A": "1.2", "B": "1.4", "C": "1.8", "D": "2.0"})],
    )
    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert result.figure_path.exists()
    assert result.caption_path.exists()
    assert result.annotations_path.exists()
    assert manifest["layout"] == "phylogram"
    assert manifest["render"]["rendered_support_count"] == 2
    assert "taxon\tlabel\tcategorical_trait\tcontinuous_trait\tlocation\ttemperature" in result.annotations_path.read_text(encoding="utf-8")


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
    non_ultrametric = diagnose_ultrametricity(fixture("example_tree_ladderized.nwk"), tolerance=1e-6)
    assert ultrametric.ultrametric is True
    assert ultrametric.max_deviation == 0.0
    assert non_ultrametric.ultrametric is False
    assert non_ultrametric.max_deviation == 0.2


def test_annotate_tree_against_table_finds_missing_and_extra_taxa() -> None:
    report = annotate_tree_against_table(fixture("example_tree.nwk"), fixture("example_traits.tsv"))
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
    exit_code = main(["metadata", "inspect", str(fixture("example_metadata.tsv")), "--json"])
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
    exit_code = main(["traits", "validate", str(fixture("example_traits_validate.tsv")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "traits"
    assert payload["metrics"]["trait_column_count"] == 3
    assert payload["data"]["trait_columns"][0]["kind"] == "numeric"


def test_cli_traits_link_json_output(capsys) -> None:
    exit_code = main(
        ["traits", "link", str(fixture("example_tree.nwk")), str(fixture("example_traits.tsv")), "--json"]
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
    exit_code = main(["alignment", "inspect", str(fixture("example_alignment.fasta")), "--json"])
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
        ["alignment", "link", str(fixture("example_tree.nwk")), str(fixture("example_alignment.fasta")), "--json"]
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


def test_cli_alignment_quality_json_output(capsys) -> None:
    exit_code = main(["alignment", "quality", str(fixture("example_alignment_duplicates.fasta")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["duplicate_group_count"] == 1
    assert payload["metrics"]["sequence_length_outlier_count"] == 0
    assert payload["metrics"]["near_duplicate_count"] == 0
    assert payload["warnings"] == ["alignment contains identical duplicate sequences"]


def test_cli_alignment_classify_json_output(capsys) -> None:
    exit_code = main(["alignment", "classify", str(fixture("example_sequences_raw.fasta")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["state"] == "raw_sequence_fasta"
    assert payload["data"]["sequence_count"] == 4


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
    exit_code = main(["alignment", "readiness", str(fixture("example_alignment_coding.fasta")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["ready_method_count"] == 3
    assert payload["metrics"]["blocked_method_count"] == 2
    assert any(method["analysis"] == "coding" and method["ready"] is False for method in payload["data"]["methods"])


def test_cli_alignment_profiles_json_output(capsys) -> None:
    exit_code = main(["alignment", "profiles", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["profile_count"] == 5
    assert any(profile["name"] == "coding-safe" for profile in payload["data"])


def test_cli_alignment_forensic_json_output(capsys) -> None:
    exit_code = main(["alignment", "forensic", str(fixture("example_alignment_coding.fasta")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["safe_for_coding_analysis"] is False
    assert payload["data"]["coding"]["mixed_coding_signals"] is True


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
    assert output.read_text(encoding="utf-8") == fixture("example_alignment_filtering_cleaned_moderate.fasta").read_text(encoding="utf-8")
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
    exit_code = main(["alignment", "length-outliers", str(fixture("example_sequences_raw.fasta")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["sequence_length_outlier_count"] == 2
    assert [row["identifier"] for row in payload["data"]] == ["B", "C"]


def test_cli_alignment_composition_json_output(capsys) -> None:
    exit_code = main(["alignment", "composition", str(fixture("example_alignment.fasta")), "--json"])
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
    exit_code = main(["alignment", "alphabet", str(fixture("example_alignment_protein.fasta")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["alphabet"] == "protein"
    assert payload["data"]["inferred_alphabet"] == "protein"


def test_cli_alignment_gc_json_output(capsys) -> None:
    exit_code = main(["alignment", "gc", str(fixture("example_alignment.fasta")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["data"]["whole_alignment_gc_content"] == 0.5
    assert payload["data"]["per_sequence_gc_content"][1] == {"gc_fraction": 0.375, "identifier": "B"}


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
    assert payload["data"] == [{"deviation": 1.0, "identifier": "C", "robust_z_score": None}]


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


def test_cli_compare_clades_json_output(capsys) -> None:
    exit_code = main(
        ["compare", "clades", str(fixture("example_tree.nwk")), str(fixture("example_tree_alt.nwk")), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["data"]["shared_clades"] == ["A|B"]
    assert payload["data"]["left_only_clades"] == ["C|D"]
    assert payload["data"]["right_only_clades"] == ["A|B|C"]


def test_cli_compare_changes_json_output(capsys) -> None:
    exit_code = main(
        ["compare", "changes", str(fixture("example_tree.nwk")), str(fixture("example_tree_alt.nwk")), "--json"]
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
    assert (output_dir / "left-shared.nwk").read_text(encoding="utf-8") == "((A:0.1,B:0.1):0.2,C:0.3);\n"
    assert (output_dir / "right-shared.nwk").read_text(encoding="utf-8") == "((A:0.1,B:0.1):0.2,C:0.3);\n"


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
    assert payload["metrics"]["table_rows"] == 3
    assert output.read_text(encoding="utf-8").startswith("split_id\tcomparison_status\tshared_clade\t")
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
    exit_code = main(["render", str(fixture("example_tree.nwk")), "--out", str(output), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["outputs"] == [str(output)]
    assert "<svg" in output.read_text(encoding="utf-8")


def test_cli_render_with_metadata_labels_reports_missing_taxa(tmp_path: Path, capsys) -> None:
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


def test_cli_render_with_partial_metadata_warns_for_missing_labels(tmp_path: Path, capsys) -> None:
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
    assert 'id="bijux-report-manifest"' in text
    assert "Bijux Tree Report" in text


def test_render_tree_report_preserves_support_and_branch_diagnostics(tmp_path: Path) -> None:
    output = tmp_path / "tree-support-report.html"
    result = render_tree_report(tree_path=fixture("example_tree_support_invalid.nwk"), out_path=output)
    assert result.validation.missing_internal_branch_nodes == []
    assert result.inspection.suspicious_support_value_ranges == [
        "support value 101 at A|B|C|D exceeds 100",
        "support value 120 at A|B exceeds 100",
        "support value -5 at C|D is negative",
    ]
    assert [warning.code for warning in result.inspection.tree_quality_warnings] == ["suspicious_support_ranges"]
    assert output.read_text(encoding="utf-8").count("suspicious_support_value_ranges") >= 1


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
    ]
    assert result.trait_missing_values is not None
    assert result.trait_missing_values.missing_values == []
    assert "Bijux Dataset Report" in text


def test_render_phylo_inputs_report_writes_alignment_sections(tmp_path: Path) -> None:
    output = tmp_path / "phylo-inputs-report.html"
    result = render_phylo_inputs_report(
        tree_path=fixture("example_tree.nwk"),
        alignment_path=fixture("example_alignment.fasta"),
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
        "tree-validation",
        "tree-inspection",
        "tree-forensic",
        "alignment-summary",
        "alignment-quality",
        "alignment-forensic",
        "alignment-coding",
        "alignment-identity-matrix",
        "alignment-linkage",
    ]
    assert "Bijux Phylo Inputs Report" in text


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
    (bundle_root / "outputs" / "outputs" / "summary.txt").write_text("drift\n", encoding="utf-8")
    report = validate_bundle(bundle_root)
    assert report.valid is False
    assert report.mismatches[0].reason in {"checksum_mismatch", "size_mismatch"}


def test_bundle_file_paths_copies_explicit_input_and_output_files(tmp_path: Path) -> None:
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
    exit_code = main(["inspect", str(fixture("example_tree_zero_lengths.nwk")), "--json"])
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
    exit_code = main(["diagnose", "distances", str(fixture("example_tree.nwk")), "--out", str(output), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["outputs"] == [str(output)]
    assert output.read_text(encoding="utf-8") == (
        "tip\tdistance\n"
        "A\t0.3\n"
        "B\t0.3\n"
        "C\t0.3\n"
        "D\t0.3\n"
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
    exit_code = main(["diagnose", "ultrametric", str(fixture("example_tree_ladderized.nwk")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["data"]["ultrametric"] is False
    assert payload["metrics"]["max_deviation"] == 0.2


def test_cli_validate_writes_run_manifest(tmp_path: Path, capsys) -> None:
    manifest = tmp_path / "validate.manifest.json"
    exit_code = main(["validate", str(fixture("example_tree.nwk")), "--json", "--manifest", str(manifest)])
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
        "discrete-evolution",
        "diversification",
        "distance",
        "tree-set",
        "simulate",
        "benchmark",
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


def test_run_capability_demo_creates_expected_artifacts(tmp_path: Path) -> None:
    result = run_capability_demo(tmp_path / "demo")
    assert result.tree_report.exists()
    assert result.dataset_report.exists()
    assert result.phylo_inputs_report.exists()
    assert result.comparison_report.exists()
    assert result.evidence_bundle.exists()
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


def test_cli_demo_run_json_output_reports_generated_artifacts(tmp_path: Path, capsys) -> None:
    output = tmp_path / "demo"
    exit_code = main(["demo", "run", "--out", str(output), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["command"] == "demo"
    assert payload["metrics"]["artifact_count"] == 6
    assert payload["data"]["tree_report"] == str(output / "reports" / "tree-report.html")
    assert payload["data"]["evidence_bundle"] == str(output / "evidence-pack")


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
    assert payload["outputs"] == [str(output)]
    assert payload["data"]["report_kind"] == "tree"


def test_cli_report_dataset_json_output_uses_dataset_contract(tmp_path: Path, capsys) -> None:
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
    assert payload["data"]["report_kind"] == "dataset"
    assert payload["metrics"]["linked_taxa"] == 4


def test_cli_report_phylo_inputs_json_output_uses_alignment_contract(tmp_path: Path, capsys) -> None:
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
    assert payload["data"]["report_kind"] == "phylo-inputs"
    assert payload["metrics"]["alignment_length"] == 8
    assert payload["metrics"]["linked_taxa"] == 4


def test_cli_adapter_returns_typed_engine_error(capsys) -> None:
    exit_code = main(["adapter", "inspect", "iqtree", "--executable", "definitely-not-installed-engine", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == EngineUnavailableError.code
