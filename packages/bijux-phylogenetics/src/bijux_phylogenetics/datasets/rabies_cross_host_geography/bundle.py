# ruff: noqa: F401, F403, F405
from __future__ import annotations

import csv
from dataclasses import dataclass, replace
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from Bio import Phylo

from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.ancestral.tree_set import (
    DiscreteAncestralTreeSetReport,
    summarize_discrete_ancestral_tree_set,
)
from bijux_phylogenetics.biogeography import (
    BiogeographyReportPackageResult,
    build_biogeography_report_package,
)
from bijux_phylogenetics.trees import (
    CladeTableReport,
    CladeTableRow,
    extract_tree_clades,
    write_clade_table,
)
from bijux_phylogenetics.comparative.pgls import (
    PGLSResult,
    run_pgls,
    write_pgls_model_matrix_table,
)
from bijux_phylogenetics.comparative.pgls_categorical_contrasts import (
    PGLSCategoricalContrastReport,
    summarize_pgls_categorical_contrasts,
    write_pgls_categorical_contrast_table,
)
from bijux_phylogenetics.comparative.pgls_lambda_fit import (
    write_pgls_lambda_profile_table,
)
from bijux_phylogenetics.comparative.posterior_tree_pgls import (
    PosteriorTreePGLSReport,
    run_posterior_tree_pgls,
)
from bijux_phylogenetics.comparative.report_package import (
    ComparativeAnalysisSummaryRow,
    ComparativeCoefficientTableRow,
    ComparativeInterpretationRow,
    ComparativeResidualTableRow,
    ComparativeSignalTableRow,
    summarize_comparative_analysis,
    summarize_comparative_audit,
    summarize_comparative_coefficients,
    summarize_comparative_interpretation,
    summarize_comparative_residuals,
    summarize_comparative_signal,
    write_comparative_audit_table,
    write_comparative_coefficient_table,
    write_comparative_contrast_table,
    write_comparative_interpretation_table,
    write_comparative_model_comparison_table,
    write_comparative_residual_table,
    write_comparative_signal_table,
    write_comparative_summary_table,
)
from bijux_phylogenetics.comparative.reporting import (
    ComparativeMethodReport,
    build_comparative_method_report,
)
from bijux_phylogenetics.compare.reports import (
    ComparisonReportBuildResult,
    build_tree_comparison_report,
)
from bijux_phylogenetics.compare.topology import write_tree_comparison_table
from bijux_phylogenetics.core.alignment import (
    AlignmentQualityReport,
    SequenceQualityRankingReport,
)
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.phylo.topology import (
    TreeRootingReport,
    root_tree_on_outgroup,
    write_tree_rooting_report,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.datasets.rabies_method_sensitivity import (
    RabiesMethodSensitivityPanelWorkflowReport,
    run_rabies_method_sensitivity_panel_workflow,
)
from bijux_phylogenetics.diagnostics.conclusion_stability import (
    ConclusionStabilityReport,
    build_ancestral_state_stability_rows,
    build_comparative_coefficient_stability_rows,
    build_conclusion_stability_report,
    build_key_clade_stability_rows,
    build_support_value_stability_rows,
    write_ancestral_state_stability_table,
    write_comparative_coefficient_stability_table,
    write_conclusion_stability_report_html,
    write_conclusion_stability_summary_table,
    write_key_clade_stability_table,
    write_support_value_stability_table,
)
from bijux_phylogenetics.engines.fasta_to_tree import (
    FastaToTreeWorkflowReport,
    run_fasta_to_tree_workflow,
)
from bijux_phylogenetics.ecology import (
    HostSwitchingReport,
    summarize_host_switching,
    write_host_state_node_table,
    write_host_switch_branch_table,
    write_host_switch_count_table,
    write_host_switch_exclusion_table,
    write_host_switch_fit_table,
    write_host_switch_summary_table,
    write_unsupported_host_switch_claim_table,
)
from bijux_phylogenetics.io.biopython import tree_from_biophylo, tree_to_biophylo
from bijux_phylogenetics.io.fasta import load_permissive_fasta_records
from bijux_phylogenetics.io.fasta.quality import (
    build_alignment_quality_report,
    build_sequence_quality_ranking,
)
from bijux_phylogenetics.io.fasta.records import validate_fasta_input
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.simulation import write_tree_set
from bijux_phylogenetics.trees import (
    BootstrapTreeSetArtifactReport,
    BootstrapTreeSetSummaryReport,
    compute_clade_frequency_table,
    write_bootstrap_tree_set_artifacts,
)

from .models import *
from .shared import _checksum, _format_number, _html_list, _support_range_text, _table

def write_rabies_cross_host_geography_panel_workflow_bundle(
    output_root: Path,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
) -> RabiesCrossHostGeographyPanelWorkflowBundle:
    """Write the complete integrated workflow bundle for the packaged rabies panel."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    workflow = report.fasta_to_tree
    host_summary = report.host_switching.summary
    geography_summary = report.biogeography_report.state_report.summary
    migration_summary = report.biogeography_report.event_report.summary

    config_audit_path = _write_workflow_config_audit_table(
        output_root / "workflow-config-audit.tsv",
        report.config_audit_rows,
    )
    resolved_config_path = _write_resolved_workflow_config(
        output_root / "workflow-config.resolved.json",
        report.config,
    )
    input_validation_path = _write_input_validation_table(
        output_root / "input-validation.tsv",
        workflow=workflow,
    )
    alignment_quality_path = _write_alignment_quality_table(
        output_root / "alignment-quality.tsv",
        aligned=report.aligned_quality,
        trimmed=report.trimmed_quality,
    )
    alignment_sequence_ranking_path = _write_sequence_ranking_table(
        output_root / "alignment-sequence-ranking.tsv",
        report.trimmed_sequence_ranking,
    )

    alignment_path = _copy_output(
        workflow.output_paths["alignment"],
        output_root / workflow.output_paths["alignment"].name,
    )
    trimmed_alignment_path = _copy_output(
        workflow.output_paths["trimmed_alignment"],
        output_root / workflow.output_paths["trimmed_alignment"].name,
    )
    tree_path = _copy_output(
        report.rooted_tree_path,
        output_root / report.rooted_tree_path.name,
    )
    stable_rooting_report = replace(
        report.rooting_report, tree_path=Path(tree_path.name)
    )
    rooting_report_path = write_tree_rooting_report(
        output_root / f"{report.dataset.workflow_prefix}.rooting.tsv",
        stable_rooting_report,
    )
    model_table_path = _copy_output(
        workflow.output_paths["model_table"],
        output_root / workflow.output_paths["model_table"].name,
    )
    support_table_path = _copy_output(
        workflow.output_paths["support_table"],
        output_root / workflow.output_paths["support_table"].name,
    )
    log_path = _copy_output(
        workflow.output_paths["log"],
        output_root / workflow.output_paths["log"].name,
    )
    manifest_path = _copy_output(
        workflow.manifest_path,
        output_root / workflow.manifest_path.name,
    )
    engine_artifact_root = (
        output_root / "engine-artifacts" / report.dataset.workflow_prefix
    )
    shutil.copytree(workflow.engine_artifact_dir, engine_artifact_root)

    clade_report = extract_tree_clades(
        report.rooted_tree_path,
        metadata_path=report.dataset.metadata_path,
        taxon_column="taxon",
        metadata_columns=list(report.dataset.clade_metadata_columns),
    )
    stable_clade_report = _stabilize_clade_report(
        clade_report,
        stable_source_path=Path(tree_path.name),
    )
    clade_table_path = write_clade_table(
        output_root / "clade-table.tsv",
        stable_clade_report,
    )

    bootstrap_output_root = output_root / "bootstrap-review"
    bootstrap_artifacts = write_bootstrap_tree_set_artifacts(
        report.fasta_to_tree.bootstrap_workflow.output_paths["bootstrap_trees"],
        out_dir=bootstrap_output_root,
        prefix="bootstrap-review",
        consensus_threshold=report.config.bootstrap_consensus_threshold,
        robust_support_threshold=report.config.bootstrap_robust_support_threshold,
        max_tree_count=report.config.max_bootstrap_tree_count,
        memory_warning_threshold_bytes=report.config.memory_warning_threshold_bytes,
    )
    bootstrap_summary_path = _write_stable_bootstrap_summary_table(
        bootstrap_output_root / "bootstrap-review.summary.tsv",
        bootstrap_artifacts.summary_report,
    )
    bootstrap_tree_comparison_report = build_tree_comparison_report(
        tree_path,
        bootstrap_artifacts.output_paths["consensus_tree"],
        out_path=bootstrap_output_root
        / "rooted-tree-vs-bootstrap-consensus.report.html",
    )
    _stabilize_bundle_report_paths(
        bootstrap_tree_comparison_report.output_path,
        output_root=output_root,
    )
    bootstrap_tree_comparison_table_path = write_tree_comparison_table(
        bootstrap_output_root / "rooted-tree-vs-bootstrap-consensus.comparison.tsv",
        tree_path,
        bootstrap_artifacts.output_paths["consensus_tree"],
    )
    bootstrap_tree_comparison_summary_path = _write_bootstrap_tree_comparison_summary(
        bootstrap_output_root / "rooted-tree-vs-bootstrap-consensus.summary.tsv",
        bootstrap_tree_comparison_report,
    )

    host_switch_summary_path = write_host_switch_summary_table(
        output_root / "host-switch-summary.tsv",
        report.host_switching,
    )
    host_state_nodes_path = write_host_state_node_table(
        output_root / "host-state-nodes.tsv",
        report.host_switching,
    )
    host_switch_branches_path = write_host_switch_branch_table(
        output_root / "host-switch-branches.tsv",
        report.host_switching,
    )
    host_switch_counts_path = write_host_switch_count_table(
        output_root / "host-switch-counts.tsv",
        report.host_switching,
    )
    host_switch_fits_path = write_host_switch_fit_table(
        output_root / "host-switch-fits.tsv",
        report.host_switching,
    )
    host_switch_unsupported_path = write_unsupported_host_switch_claim_table(
        output_root / "host-switch-unsupported.tsv",
        report.host_switching,
    )
    host_switch_exclusions_path = write_host_switch_exclusion_table(
        output_root / "host-switch-exclusions.tsv",
        report.host_switching,
    )

    biogeography_output_root = output_root / "biogeography"
    shutil.copytree(report.biogeography_report.output_dir, biogeography_output_root)
    biogeography_report_path = biogeography_output_root / "biogeography-report.html"
    biogeography_tree_figure_path = (
        biogeography_output_root / "ancestral-region-tree.svg"
    )
    biogeography_map_path = biogeography_output_root / "geographic-region-map.html"

    comparative_traits_path = write_taxon_rows(
        output_root / "comparative-traits.tsv",
        columns=list(report.comparative_traits_rows[0].keys()),
        rows=report.comparative_traits_rows,
    )
    comparative_tree_path = _copy_output(
        report.comparative_tree_path,
        output_root / "comparative-tree.nwk",
    )
    comparative_repairs_path = _write_comparative_branch_repairs_table(
        output_root / "comparative-tree-adjustments.tsv",
        report.comparative_branch_repairs,
    )
    comparative_output_root = output_root / "comparative"
    comparative_output_root.mkdir(parents=True, exist_ok=True)
    comparative_report = report.comparative_report
    comparative_summary_row = summarize_comparative_analysis(comparative_report)
    comparative_coefficient_rows = summarize_comparative_coefficients(
        comparative_report
    )
    comparative_residual_rows = summarize_comparative_residuals(comparative_report)
    comparative_signal_row = summarize_comparative_signal(comparative_report)
    comparative_interpretation_rows = summarize_comparative_interpretation(
        comparative_report
    )
    comparative_audit_rows = summarize_comparative_audit(comparative_report)
    comparative_report_path = _write_comparative_report(
        comparative_output_root / "comparative-report.html",
        summary_row=comparative_summary_row,
        coefficient_rows=comparative_coefficient_rows,
        residual_rows=comparative_residual_rows,
        signal_row=comparative_signal_row,
        interpretation_rows=comparative_interpretation_rows,
        branch_repairs=report.comparative_branch_repairs,
    )
    comparative_summary_path = write_comparative_summary_table(
        comparative_output_root / "comparative-summary.tsv",
        comparative_summary_row,
    )
    comparative_coefficients_path = write_comparative_coefficient_table(
        comparative_output_root / "coefficient-table.tsv",
        comparative_coefficient_rows,
    )
    comparative_residuals_path = write_comparative_residual_table(
        comparative_output_root / "residual-summary.tsv",
        comparative_residual_rows,
    )
    comparative_signal_path = write_comparative_signal_table(
        comparative_output_root / "signal-summary.tsv",
        comparative_signal_row,
    )
    comparative_model_comparison_path = write_comparative_model_comparison_table(
        comparative_output_root / "model-comparison.tsv",
        comparative_report,
    )
    comparative_interpretation_path = write_comparative_interpretation_table(
        comparative_output_root / "interpretation-table.tsv",
        comparative_interpretation_rows,
    )
    comparative_audit_path = write_comparative_audit_table(
        comparative_output_root / "audit-table.tsv",
        comparative_audit_rows,
    )
    comparative_contrasts_path = write_comparative_contrast_table(
        comparative_output_root / "contrast-table.tsv",
        comparative_report,
    )
    comparative_model_matrix_path = comparative_output_root / "model-matrix.tsv"
    write_pgls_model_matrix_table(
        comparative_model_matrix_path,
        comparative_report.snapshot.pgls_inputs.model_matrix,
    )
    comparative_categorical_contrasts_path = write_pgls_categorical_contrast_table(
        comparative_output_root / "categorical-contrasts.tsv",
        report.comparative_categorical_contrasts,
    )
    comparative_lambda_profile_path = write_pgls_lambda_profile_table(
        comparative_output_root / "lambda-profile.tsv",
        comparative_report.snapshot.pgls_model.lambda_fit,
    )
    comparative_manifest_path = _write_comparative_manifest(
        comparative_output_root / "comparative.manifest.json",
        comparative_summary_row=comparative_summary_row,
        branch_repairs=report.comparative_branch_repairs,
        output_paths={
            "comparative_report": comparative_report_path,
            "summary_table": comparative_summary_path,
            "coefficient_table": comparative_coefficients_path,
            "residual_table": comparative_residuals_path,
            "signal_table": comparative_signal_path,
            "model_comparison_table": comparative_model_comparison_path,
            "interpretation_table": comparative_interpretation_path,
            "audit_table": comparative_audit_path,
            "contrast_table": comparative_contrasts_path,
            "model_matrix_table": comparative_model_matrix_path,
            "categorical_contrast_table": comparative_categorical_contrasts_path,
            "lambda_profile_table": comparative_lambda_profile_path,
        },
    )
    conclusion_stability_output_root = output_root / "conclusion-stability"
    conclusion_stability_output_root.mkdir(parents=True, exist_ok=True)
    conclusion_stability_summary_path = write_conclusion_stability_summary_table(
        conclusion_stability_output_root / "conclusion-stability-summary.tsv",
        report.conclusion_stability_report,
    )
    key_clade_stability_path = write_key_clade_stability_table(
        conclusion_stability_output_root / "key-clade-stability.tsv",
        report.conclusion_stability_report.key_clade_rows,
    )
    support_value_stability_path = write_support_value_stability_table(
        conclusion_stability_output_root / "support-value-stability.tsv",
        report.conclusion_stability_report.support_value_rows,
    )
    ancestral_state_stability_path = write_ancestral_state_stability_table(
        conclusion_stability_output_root / "ancestral-state-stability.tsv",
        report.conclusion_stability_report.ancestral_state_rows,
    )
    comparative_coefficient_stability_path = (
        write_comparative_coefficient_stability_table(
            conclusion_stability_output_root / "comparative-coefficient-stability.tsv",
            report.conclusion_stability_report.comparative_coefficient_rows,
        )
    )
    conclusion_stability_report_path = write_conclusion_stability_report_html(
        conclusion_stability_output_root / "conclusion-stability-report.html",
        report.conclusion_stability_report,
    )
    scientific_finding_rows = _build_scientific_finding_rows(
        report=report,
        bootstrap_tree_comparison_report=bootstrap_tree_comparison_report,
        comparative_summary_row=comparative_summary_row,
        comparative_interpretation_rows=comparative_interpretation_rows,
    )
    scientific_findings_path = _write_scientific_findings_table(
        output_root / "scientific-findings.tsv",
        scientific_finding_rows,
    )

    workflow_summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report=report,
        clade_row_count=len(stable_clade_report.rows),
        bootstrap_artifacts=bootstrap_artifacts,
        bootstrap_tree_comparison_report=bootstrap_tree_comparison_report,
        comparative_summary_row=comparative_summary_row,
        scientific_finding_count=len(scientific_finding_rows),
    )
    resource_observations_path = _write_resource_observation_table(
        output_root / "resource-observations.tsv",
        report=report,
        bootstrap_artifacts=bootstrap_artifacts,
    )
    final_report_path = _write_integrated_report(
        output_root / "rabies-cross-host-geography-report.html",
        report=report,
        workflow_summary_path=workflow_summary_path,
        bootstrap_artifacts=bootstrap_artifacts,
        bootstrap_tree_comparison_report=bootstrap_tree_comparison_report,
        clade_row_count=len(stable_clade_report.rows),
        comparative_summary_row=comparative_summary_row,
        comparative_interpretation_rows=comparative_interpretation_rows,
        comparative_branch_repair_count=len(report.comparative_branch_repairs),
        scientific_finding_rows=scientific_finding_rows,
        max_report_table_rows=report.config.max_report_table_rows,
    )
    final_manifest_path = _write_manifest(
        output_root / "rabies-cross-host-geography.manifest.json",
        report=report,
        comparative_summary_row=comparative_summary_row,
        bootstrap_artifacts=bootstrap_artifacts,
        bootstrap_tree_comparison_report=bootstrap_tree_comparison_report,
        clade_row_count=len(stable_clade_report.rows),
        scientific_finding_count=len(scientific_finding_rows),
        bundle_paths={
            "workflow_summary": workflow_summary_path,
            "config_audit": config_audit_path,
            "resolved_config": resolved_config_path,
            "input_validation": input_validation_path,
            "alignment_quality": alignment_quality_path,
            "alignment_sequence_ranking": alignment_sequence_ranking_path,
            "alignment": alignment_path,
            "trimmed_alignment": trimmed_alignment_path,
            "rooted_tree": tree_path,
            "rooting_report": rooting_report_path,
            "model_table": model_table_path,
            "support_table": support_table_path,
            "clade_table": clade_table_path,
            "bootstrap_summary": bootstrap_summary_path,
            "bootstrap_consensus_tree": bootstrap_artifacts.output_paths[
                "consensus_tree"
            ],
            "bootstrap_clade_frequencies": bootstrap_artifacts.output_paths[
                "clade_frequencies"
            ],
            "bootstrap_unstable_branches": bootstrap_artifacts.output_paths[
                "unstable_branches"
            ],
            "bootstrap_unstable_clades": bootstrap_artifacts.output_paths[
                "unstable_clades"
            ],
            "bootstrap_distance_matrix": bootstrap_artifacts.output_paths[
                "distance_matrix"
            ],
            "bootstrap_topology_clusters": bootstrap_artifacts.output_paths[
                "topology_clusters"
            ],
            "bootstrap_tree_comparison_summary": bootstrap_tree_comparison_summary_path,
            "bootstrap_tree_comparison_table": bootstrap_tree_comparison_table_path,
            "bootstrap_tree_comparison_report": (
                bootstrap_tree_comparison_report.output_path
            ),
            "host_switch_summary": host_switch_summary_path,
            "host_state_nodes": host_state_nodes_path,
            "host_switch_branches": host_switch_branches_path,
            "host_switch_counts": host_switch_counts_path,
            "host_switch_fits": host_switch_fits_path,
            "host_switch_unsupported": host_switch_unsupported_path,
            "host_switch_exclusions": host_switch_exclusions_path,
            "biogeography_report": biogeography_report_path,
            "biogeography_tree_figure": biogeography_tree_figure_path,
            "biogeography_map": biogeography_map_path,
            "comparative_traits": comparative_traits_path,
            "comparative_tree": comparative_tree_path,
            "comparative_repairs": comparative_repairs_path,
            "comparative_report": comparative_report_path,
            "comparative_summary": comparative_summary_path,
            "comparative_coefficients": comparative_coefficients_path,
            "comparative_residuals": comparative_residuals_path,
            "comparative_signal": comparative_signal_path,
            "comparative_model_comparison": comparative_model_comparison_path,
            "comparative_interpretation": comparative_interpretation_path,
            "comparative_audit": comparative_audit_path,
            "comparative_contrasts": comparative_contrasts_path,
            "comparative_model_matrix": comparative_model_matrix_path,
            "comparative_categorical_contrasts": comparative_categorical_contrasts_path,
            "comparative_lambda_profile": comparative_lambda_profile_path,
            "comparative_manifest": comparative_manifest_path,
            "conclusion_stability_summary": conclusion_stability_summary_path,
            "key_clade_stability": key_clade_stability_path,
            "support_value_stability": support_value_stability_path,
            "ancestral_state_stability": ancestral_state_stability_path,
            "comparative_coefficient_stability": (
                comparative_coefficient_stability_path
            ),
            "conclusion_stability_report": conclusion_stability_report_path,
            "scientific_findings": scientific_findings_path,
            "final_report": final_report_path,
        },
    )

    return RabiesCrossHostGeographyPanelWorkflowBundle(
        output_root=output_root,
        selected_model=workflow.selected_model,
        sequence_type=report.dataset.sequence_type,
        inferred_sequence_type=workflow.sequence_type,
        input_repair_applied=workflow.input_repair is not None,
        aligned_quality_score=report.aligned_quality.quality_score,
        trimmed_quality_score=report.trimmed_quality.quality_score,
        minimum_support=workflow.support_summary.minimum_support,
        maximum_support=workflow.support_summary.maximum_support,
        median_support=workflow.support_summary.median_support,
        weakly_supported_clade_count=workflow.support_summary.weakly_supported_clade_count,
        clade_row_count=len(stable_clade_report.rows),
        bootstrap_tree_count=bootstrap_artifacts.summary_report.tree_count,
        bootstrap_topology_count=(
            bootstrap_artifacts.summary_report.diversity.rooted_topology_count
        ),
        bootstrap_unstable_branch_count=(
            bootstrap_artifacts.summary_report.unstable_branch_count
        ),
        bootstrap_consensus_rooted_rf_distance=(
            bootstrap_tree_comparison_report.topology.rooted_robinson_foulds_distance
        ),
        bootstrap_consensus_same_unrooted_topology=(
            bootstrap_tree_comparison_report.topology.same_unrooted_topology
        ),
        bootstrap_consensus_high_support_conflict_count=len(
            [
                row
                for row in bootstrap_tree_comparison_report.support.conflicting_clades
                if row.conflict_classification == "high_support_conflict"
            ]
        ),
        bootstrap_consensus_branch_score_distance=(
            bootstrap_tree_comparison_report.branch_lengths.branch_score.branch_score_distance
        ),
        rooted_outgroup_taxa=tuple(report.rooting_report.rooted_outgroup_taxa),
        root_host=host_summary.root_host,
        root_host_confidence=host_summary.root_confidence,
        host_switch_count=host_summary.host_switch_count,
        certain_host_switch_count=host_summary.certain_host_switch_count,
        uncertain_host_switch_count=host_summary.uncertain_host_switch_count,
        root_region=geography_summary.root_region,
        root_region_probability=geography_summary.root_region_probability,
        changed_region_branch_count=geography_summary.changed_branch_count,
        migration_event_count=migration_summary.event_count,
        strongly_supported_migration_event_count=(
            migration_summary.strongly_supported_event_count
        ),
        comparative_selected_model=comparative_summary_row.selected_model,
        comparative_response=comparative_summary_row.response,
        comparative_formula=comparative_summary_row.formula,
        comparative_pgls_lambda=comparative_summary_row.pgls_lambda,
        comparative_pgls_r_squared=comparative_summary_row.pgls_r_squared,
        comparative_branch_repair_count=len(report.comparative_branch_repairs),
        conclusion_stable_count=report.conclusion_stability_report.summary.stable_count,
        conclusion_weak_count=report.conclusion_stability_report.summary.weak_count,
        conclusion_unstable_count=(
            report.conclusion_stability_report.summary.unstable_count
        ),
        timeout_seconds=report.config.timeout_seconds,
        max_bootstrap_tree_count=report.config.max_bootstrap_tree_count,
        max_report_table_rows=report.config.max_report_table_rows,
        memory_warning_threshold_bytes=report.config.memory_warning_threshold_bytes,
        workflow_runtime_seconds=report.fasta_to_tree.runtime_seconds,
        bootstrap_review_runtime_seconds=(
            bootstrap_artifacts.summary_report.processing.runtime_seconds
        ),
        bootstrap_review_peak_memory_bytes=(
            bootstrap_artifacts.summary_report.processing.peak_memory_bytes
        ),
        budget_warning_count=len(bootstrap_artifacts.budget_report.warning_messages),
        config_check_count=len(report.config_audit_rows),
        scientific_finding_count=len(scientific_finding_rows),
        workflow_summary_path=workflow_summary_path,
        resource_observations_path=resource_observations_path,
        config_audit_path=config_audit_path,
        resolved_config_path=resolved_config_path,
        input_validation_path=input_validation_path,
        alignment_quality_path=alignment_quality_path,
        alignment_sequence_ranking_path=alignment_sequence_ranking_path,
        alignment_path=alignment_path,
        trimmed_alignment_path=trimmed_alignment_path,
        tree_path=tree_path,
        rooting_report_path=rooting_report_path,
        model_table_path=model_table_path,
        support_table_path=support_table_path,
        log_path=log_path,
        manifest_path=manifest_path,
        engine_artifact_root=engine_artifact_root,
        clade_table_path=clade_table_path,
        bootstrap_output_root=bootstrap_output_root,
        bootstrap_summary_path=bootstrap_summary_path,
        bootstrap_consensus_tree_path=bootstrap_artifacts.output_paths[
            "consensus_tree"
        ],
        bootstrap_clade_frequencies_path=bootstrap_artifacts.output_paths[
            "clade_frequencies"
        ],
        bootstrap_unstable_branches_path=bootstrap_artifacts.output_paths[
            "unstable_branches"
        ],
        bootstrap_unstable_clades_path=bootstrap_artifacts.output_paths[
            "unstable_clades"
        ],
        bootstrap_distance_matrix_path=bootstrap_artifacts.output_paths[
            "distance_matrix"
        ],
        bootstrap_topology_clusters_path=bootstrap_artifacts.output_paths[
            "topology_clusters"
        ],
        bootstrap_tree_comparison_summary_path=(bootstrap_tree_comparison_summary_path),
        bootstrap_tree_comparison_table_path=bootstrap_tree_comparison_table_path,
        bootstrap_tree_comparison_report_path=(
            bootstrap_tree_comparison_report.output_path
        ),
        host_switch_summary_path=host_switch_summary_path,
        host_state_nodes_path=host_state_nodes_path,
        host_switch_branches_path=host_switch_branches_path,
        host_switch_counts_path=host_switch_counts_path,
        host_switch_fits_path=host_switch_fits_path,
        host_switch_unsupported_path=host_switch_unsupported_path,
        host_switch_exclusions_path=host_switch_exclusions_path,
        biogeography_output_root=biogeography_output_root,
        biogeography_report_path=biogeography_report_path,
        biogeography_tree_figure_path=biogeography_tree_figure_path,
        biogeography_map_path=biogeography_map_path,
        comparative_traits_path=comparative_traits_path,
        comparative_tree_path=comparative_tree_path,
        comparative_repairs_path=comparative_repairs_path,
        comparative_output_root=comparative_output_root,
        comparative_report_path=comparative_report_path,
        comparative_summary_path=comparative_summary_path,
        comparative_coefficients_path=comparative_coefficients_path,
        comparative_residuals_path=comparative_residuals_path,
        comparative_signal_path=comparative_signal_path,
        comparative_model_comparison_path=comparative_model_comparison_path,
        comparative_interpretation_path=comparative_interpretation_path,
        comparative_audit_path=comparative_audit_path,
        comparative_contrasts_path=comparative_contrasts_path,
        comparative_model_matrix_path=comparative_model_matrix_path,
        comparative_categorical_contrasts_path=comparative_categorical_contrasts_path,
        comparative_lambda_profile_path=comparative_lambda_profile_path,
        comparative_manifest_path=comparative_manifest_path,
        conclusion_stability_output_root=conclusion_stability_output_root,
        conclusion_stability_summary_path=conclusion_stability_summary_path,
        key_clade_stability_path=key_clade_stability_path,
        support_value_stability_path=support_value_stability_path,
        ancestral_state_stability_path=ancestral_state_stability_path,
        comparative_coefficient_stability_path=(comparative_coefficient_stability_path),
        conclusion_stability_report_path=conclusion_stability_report_path,
        scientific_findings_path=scientific_findings_path,
        final_report_path=final_report_path,
        final_manifest_path=final_manifest_path,
    )



def _copy_output(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    return Path(shutil.copy2(source, destination))


def _write_workflow_config_audit_table(
    path: Path,
    rows: list[RabiesWorkflowConfigAuditRow],
) -> Path:
    return write_taxon_rows(
        path,
        columns=["check_id", "status", "observed_value", "detail"],
        rows=[
            {
                "check_id": row.check_id,
                "status": row.status,
                "observed_value": row.observed_value,
                "detail": row.detail,
            }
            for row in rows
        ],
    )


def _write_resolved_workflow_config(
    path: Path,
    config: RabiesCrossHostGeographyPanelWorkflowConfig,
) -> Path:
    payload = {
        "report_kind": "rabies_cross_host_geography_workflow_config",
        "dataset_id": config.dataset_id,
        "label": config.label,
        "source_config": config.config_path.name,
        "input_files": {
            "sequences_path": {
                "path": config.sequences_path.name,
                "sha256": _checksum(config.sequences_path),
            },
            "metadata_path": {
                "path": config.metadata_path.name,
                "sha256": _checksum(config.metadata_path),
            },
            "centroids_path": {
                "path": config.centroids_path.name,
                "sha256": _checksum(config.centroids_path),
            },
        },
        "workflow": {
            "sequence_type": config.sequence_type,
            "workflow_prefix": config.workflow_prefix,
            "host_trait": config.host_trait,
            "geography_trait": config.geography_trait,
            "host_model": config.host_model,
            "geography_model": config.geography_model,
            "outgroup_taxa": list(config.outgroup_taxa),
            "iqtree_seed": config.iqtree_seed,
            "iqtree_threads": config.iqtree_threads,
            "bootstrap_replicates": config.bootstrap_replicates,
            "timeout_seconds": config.timeout_seconds,
            "max_bootstrap_tree_count": config.max_bootstrap_tree_count,
            "max_report_table_rows": config.max_report_table_rows,
            "memory_warning_threshold_bytes": (config.memory_warning_threshold_bytes),
            "alignment_mode": config.alignment_mode,
            "trimming_mode": config.trimming_mode,
            "trim_gap_threshold": config.trim_gap_threshold,
            "bootstrap_consensus_threshold": config.bootstrap_consensus_threshold,
            "bootstrap_robust_support_threshold": (
                config.bootstrap_robust_support_threshold
            ),
            "clade_metadata_columns": list(config.clade_metadata_columns),
            "comparative_formula": config.comparative_formula,
            "comparative_response": config.comparative_response,
            "comparative_branch_length_floor": (config.comparative_branch_length_floor),
        },
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _write_bootstrap_tree_comparison_summary(
    path: Path,
    comparison_report: ComparisonReportBuildResult,
) -> Path:
    high_support_conflict_count = len(
        [
            row
            for row in comparison_report.support.conflicting_clades
            if row.conflict_classification == "high_support_conflict"
        ]
    )
    row = {
        "left_tree": comparison_report.topology.left_path.name,
        "right_tree": comparison_report.topology.right_path.name,
        "shared_taxon_count": str(len(comparison_report.topology.shared_taxa)),
        "rooted_rf_distance": str(
            comparison_report.topology.rooted_robinson_foulds_distance
        ),
        "rooted_normalized_rf": _format_number(
            comparison_report.topology.rooted_normalized_robinson_foulds
        ),
        "topology_equal": (
            "true" if comparison_report.topology.topology_equal else "false"
        ),
        "same_unrooted_topology": (
            "true" if comparison_report.topology.same_unrooted_topology else "false"
        ),
        "same_taxa_different_rooting": (
            "true"
            if comparison_report.topology.same_taxa_different_rooting
            else "false"
        ),
        "same_topology_different_branch_lengths": (
            "true"
            if comparison_report.topology.same_topology_different_branch_lengths
            else "false"
        ),
        "support_conflict_count": str(
            len(comparison_report.support.conflicting_clades)
        ),
        "high_support_conflict_count": str(high_support_conflict_count),
        "branch_score_distance": _format_number(
            comparison_report.branch_lengths.branch_score.branch_score_distance
        ),
    }
    return write_taxon_rows(path, columns=list(row.keys()), rows=[row])


def _write_stable_bootstrap_summary_table(
    path: Path,
    report: BootstrapTreeSetSummaryReport,
) -> Path:
    row = {
        "tree_count": str(report.tree_count),
        "shared_taxon_count": str(len(report.shared_taxa)),
        "rooted_topology_count": str(report.diversity.rooted_topology_count),
        "dominant_topology_frequency": _format_number(
            report.diversity.dominant_topology_frequency
        ),
        "effective_topology_count": _format_number(
            report.diversity.effective_topology_count
        ),
        "mean_robinson_foulds_distance": _format_number(
            report.diversity.mean_robinson_foulds_distance
        ),
        "mean_normalized_robinson_foulds_distance": _format_number(
            report.diversity.mean_normalized_robinson_foulds_distance
        ),
        "consensus_threshold": _format_number(report.consensus_threshold),
        "robust_support_threshold": _format_number(report.robust_support_threshold),
        "unstable_branch_count": str(report.unstable_branch_count),
        "warning_count": str(len(report.warnings)),
        "consensus_newick": report.consensus.consensus_newick,
    }
    return write_taxon_rows(path, columns=list(row.keys()), rows=[row])


def _stabilize_bundle_report_paths(path: Path, *, output_root: Path) -> Path:
    text = path.read_text(encoding="utf-8")
    normalized_root = output_root.as_posix().rstrip("/") + "/"
    path.write_text(text.replace(normalized_root, ""), encoding="utf-8")
    return path


def _build_scientific_finding_rows(
    *,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    bootstrap_tree_comparison_report: ComparisonReportBuildResult,
    comparative_summary_row: ComparativeAnalysisSummaryRow,
    comparative_interpretation_rows: list[ComparativeInterpretationRow],
) -> list[RabiesScientificFindingRow]:
    host_summary = report.host_switching.summary
    geography_summary = report.biogeography_report.state_report.summary
    migration_summary = report.biogeography_report.event_report.summary
    bootstrap_question = (
        "Does the bootstrap consensus preserve the rooted ML conclusion?"
    )
    if bootstrap_tree_comparison_report.topology.topology_equal:
        bootstrap_claim = "The bootstrap consensus preserves the rooted ML topology on the shared taxon set."
    else:
        bootstrap_claim = "The bootstrap consensus differs from the rooted ML topology after support-driven summarization."
    comparative_claim = next(
        (
            row.claim
            for row in comparative_interpretation_rows
            if row.topic == "coefficient" and row.claim
        ),
        "The comparative layer did not expose one stable host-associated longitude shift.",
    )
    return [
        RabiesScientificFindingRow(
            finding_id="root_host_state",
            question="What host state anchors the rooted rabies panel?",
            claim=f"The rooted tree places the ancestral host state in {host_summary.root_host}.",
            evidence=(
                f"root host confidence {_format_number(host_summary.root_confidence)} "
                f"with outgroup {','.join(report.dataset.outgroup_taxa)}"
            ),
            caution=(
                "The panel is compact and grouped by broad host classes rather than species-level host states."
            ),
            source_artifact=report.dataset.workflow_prefix + ".rooting.tsv",
        ),
        RabiesScientificFindingRow(
            finding_id="root_region_state",
            question="What geographic regime anchors the rooted rabies panel?",
            claim=(
                f"The rooted tree places the ancestral region in {geography_summary.root_region}."
            ),
            evidence=(
                f"root region probability {_format_number(geography_summary.root_region_probability)} "
                f"across {geography_summary.changed_branch_count} changed branches"
            ),
            caution=(
                "Grouped macroregions simplify the raw locality labels so the result should be treated as regional rather than site-level history."
            ),
            source_artifact="biogeography/summary.tsv",
        ),
        RabiesScientificFindingRow(
            finding_id="host_switching",
            question="How much host-switching signal appears in the rooted tree?",
            claim=(
                f"The host reconstruction inferred {host_summary.host_switch_count} host-switch branch changes."
            ),
            evidence=(
                f"certain changes {host_summary.certain_host_switch_count}; "
                f"uncertain changes {host_summary.uncertain_host_switch_count}"
            ),
            caution=(
                "Branch-wise host changes depend on the grouped host coding and should not be over-read as one exhaustive host-jump catalogue."
            ),
            source_artifact="host-switch-summary.tsv",
        ),
        RabiesScientificFindingRow(
            finding_id="bootstrap_consensus",
            question=bootstrap_question,
            claim=bootstrap_claim,
            evidence=(
                f"rooted RF distance {bootstrap_tree_comparison_report.topology.rooted_robinson_foulds_distance}; "
                f"high-support conflicts "
                f"{len([row for row in bootstrap_tree_comparison_report.support.conflicting_clades if row.conflict_classification == 'high_support_conflict'])}"
            ),
            caution=(
                "Consensus trees can collapse low-support branches, so exact rooted agreement is stricter than shared major clades."
            ),
            source_artifact=(
                "bootstrap-review/rooted-tree-vs-bootstrap-consensus.summary.tsv"
            ),
        ),
        RabiesScientificFindingRow(
            finding_id="comparative_longitude",
            question=(
                "Do host-associated rabies lineages occupy one distinct longitudinal regime in this panel?"
            ),
            claim=comparative_claim,
            evidence=(
                f"selected model {comparative_summary_row.selected_model}; "
                f"PGLS lambda {_format_number(comparative_summary_row.pgls_lambda)}; "
                f"r-squared {_format_number(comparative_summary_row.pgls_r_squared)}"
            ),
            caution=(
                "The comparative claim is associational, uses only nine taxa, and retains residual-diagnostic cautions."
            ),
            source_artifact="comparative/interpretation-table.tsv",
        ),
        RabiesScientificFindingRow(
            finding_id="migration_events",
            question="How much regional movement is implied by the geographic reconstruction?",
            claim=(
                f"The biogeography layer inferred {migration_summary.event_count} migration events across the rooted tree."
            ),
            evidence=(
                f"strongly supported migration events {migration_summary.strongly_supported_event_count}"
            ),
            caution=(
                "Event counts summarize transitions over grouped regions and do not replace one dated dispersal analysis."
            ),
            source_artifact="biogeography/event-table.tsv",
        ),
    ]


def _write_scientific_findings_table(
    path: Path,
    rows: list[RabiesScientificFindingRow],
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "finding_id",
            "question",
            "claim",
            "evidence",
            "caution",
            "source_artifact",
        ],
        rows=[
            {
                "finding_id": row.finding_id,
                "question": row.question,
                "claim": row.claim,
                "evidence": row.evidence,
                "caution": row.caution,
                "source_artifact": row.source_artifact,
            }
            for row in rows
        ],
    )


def _write_input_validation_table(
    path: Path,
    *,
    workflow: FastaToTreeWorkflowReport,
) -> Path:
    validation = (
        workflow.input_validation
        if workflow.repaired_input_validation is None
        else workflow.repaired_input_validation
    )
    sequence_type_report = validation.sequence_type_report
    row = {
        "sequence_count": str(validation.summary.sequence_count),
        "detected_type": sequence_type_report.detected_type or "",
        "selected_type": sequence_type_report.selected_type or "",
        "confidence": sequence_type_report.confidence or "",
        "repair_required": "true"
        if (
            workflow.input_validation.duplicate_identifiers
            or workflow.input_validation.illegal_characters
            or workflow.input_validation.empty_sequences
        )
        else "false",
        "repair_applied": "true" if workflow.input_repair is not None else "false",
        "duplicate_identifier_count": str(
            len(workflow.input_validation.duplicate_identifiers)
        ),
        "illegal_character_count": str(
            len(workflow.input_validation.illegal_characters)
        ),
        "empty_sequence_count": str(len(workflow.input_validation.empty_sequences)),
        "warning_count": str(len(validation.warnings)),
        "warnings": " | ".join(validation.warnings),
    }
    return write_taxon_rows(path, columns=list(row.keys()), rows=[row])


def _write_alignment_quality_table(
    path: Path,
    *,
    aligned: AlignmentQualityReport,
    trimmed: AlignmentQualityReport,
) -> Path:
    rows = []
    for stage, report in (("aligned", aligned), ("trimmed", trimmed)):
        rows.append(
            {
                "stage": stage,
                "sequence_count": str(report.sequence_count),
                "alignment_length": str(report.alignment_length),
                "missing_data_fraction": _format_number(report.missing_data_fraction),
                "gap_fraction": _format_number(report.gap_fraction),
                "ambiguity_fraction": _format_number(report.ambiguity_fraction),
                "variable_site_count": str(report.variable_site_count),
                "parsimony_informative_site_count": str(
                    report.parsimony_informative_site_count
                ),
                "quality_score": _format_number(report.quality_score),
                "suspicious_alignment": (
                    "true" if report.suspicious_alignment else "false"
                ),
                "suspicious_reasons": " | ".join(report.suspicious_reasons),
            }
        )
    return write_taxon_rows(path, columns=list(rows[0].keys()), rows=rows)


def _write_sequence_ranking_table(
    path: Path,
    report: SequenceQualityRankingReport,
) -> Path:
    rows = [
        {
            "identifier": row.identifier,
            "rank": str(row.rank),
            "score": _format_number(row.score),
            "missing_fraction": _format_number(row.missing_fraction),
            "gap_fraction": _format_number(row.gap_fraction),
            "ambiguity_fraction": _format_number(row.ambiguity_fraction),
            "composition_outlier": "true" if row.composition_outlier else "false",
            "duplicate_status": row.duplicate_status,
            "note": row.note,
        }
        for row in report.rows
    ]
    return write_taxon_rows(path, columns=list(rows[0].keys()), rows=rows)


def _write_comparative_branch_repairs_table(
    path: Path,
    rows: list[RabiesComparativeBranchRepair],
) -> Path:
    if not rows:
        return write_taxon_rows(
            path,
            columns=[
                "node_label",
                "original_branch_length",
                "repaired_branch_length",
                "reason",
            ],
            rows=[],
        )
    return write_taxon_rows(
        path,
        columns=[
            "node_label",
            "original_branch_length",
            "repaired_branch_length",
            "reason",
        ],
        rows=[
            {
                "node_label": row.node_label,
                "original_branch_length": _format_number(row.original_branch_length),
                "repaired_branch_length": _format_number(row.repaired_branch_length),
                "reason": row.reason,
            }
            for row in rows
        ],
    )


def _write_workflow_summary_table(
    path: Path,
    *,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    clade_row_count: int,
    bootstrap_artifacts: BootstrapTreeSetArtifactReport,
    bootstrap_tree_comparison_report: ComparisonReportBuildResult,
    comparative_summary_row: ComparativeAnalysisSummaryRow,
    scientific_finding_count: int,
) -> Path:
    support = report.fasta_to_tree.support_summary
    host_summary = report.host_switching.summary
    geography_summary = report.biogeography_report.state_report.summary
    migration_summary = report.biogeography_report.event_report.summary
    bootstrap_summary = bootstrap_artifacts.summary_report
    row = {
        "dataset_id": report.dataset.dataset_id,
        "sequence_count": str(report.dataset.sequence_count),
        "sequence_type": report.dataset.sequence_type,
        "inferred_sequence_type": report.fasta_to_tree.sequence_type,
        "selected_model": report.fasta_to_tree.selected_model,
        "input_repair_applied": (
            "true" if report.fasta_to_tree.input_repair is not None else "false"
        ),
        "aligned_quality_score": _format_number(report.aligned_quality.quality_score),
        "trimmed_quality_score": _format_number(report.trimmed_quality.quality_score),
        "minimum_support": _format_number(support.minimum_support),
        "maximum_support": _format_number(support.maximum_support),
        "median_support": _format_number(support.median_support),
        "weakly_supported_clade_count": str(support.weakly_supported_clade_count),
        "clade_row_count": str(clade_row_count),
        "bootstrap_tree_count": str(bootstrap_summary.tree_count),
        "bootstrap_topology_count": str(
            bootstrap_summary.diversity.rooted_topology_count
        ),
        "bootstrap_unstable_branch_count": str(bootstrap_summary.unstable_branch_count),
        "bootstrap_consensus_rooted_rf_distance": str(
            bootstrap_tree_comparison_report.topology.rooted_robinson_foulds_distance
        ),
        "bootstrap_consensus_same_unrooted_topology": (
            "true"
            if bootstrap_tree_comparison_report.topology.same_unrooted_topology
            else "false"
        ),
        "bootstrap_consensus_high_support_conflict_count": str(
            len(
                [
                    row
                    for row in bootstrap_tree_comparison_report.support.conflicting_clades
                    if row.conflict_classification == "high_support_conflict"
                ]
            )
        ),
        "outgroup_taxa": ",".join(report.dataset.outgroup_taxa),
        "root_host": host_summary.root_host,
        "root_host_confidence": _format_number(host_summary.root_confidence),
        "host_switch_count": str(host_summary.host_switch_count),
        "certain_host_switch_count": str(host_summary.certain_host_switch_count),
        "uncertain_host_switch_count": str(host_summary.uncertain_host_switch_count),
        "root_region": geography_summary.root_region,
        "root_region_probability": _format_number(
            geography_summary.root_region_probability
        ),
        "changed_region_branch_count": str(geography_summary.changed_branch_count),
        "migration_event_count": str(migration_summary.event_count),
        "strongly_supported_migration_event_count": str(
            migration_summary.strongly_supported_event_count
        ),
        "comparative_response": comparative_summary_row.response,
        "comparative_formula": comparative_summary_row.formula,
        "comparative_selected_model": comparative_summary_row.selected_model,
        "comparative_pgls_lambda": _format_number(comparative_summary_row.pgls_lambda),
        "comparative_pgls_r_squared": _format_number(
            comparative_summary_row.pgls_r_squared
        ),
        "comparative_branch_repair_count": str(len(report.comparative_branch_repairs)),
        "conclusion_stable_count": str(
            report.conclusion_stability_report.summary.stable_count
        ),
        "conclusion_weak_count": str(
            report.conclusion_stability_report.summary.weak_count
        ),
        "conclusion_unstable_count": str(
            report.conclusion_stability_report.summary.unstable_count
        ),
        "timeout_seconds": _format_number(report.config.timeout_seconds),
        "max_bootstrap_tree_count": (
            ""
            if report.config.max_bootstrap_tree_count is None
            else str(report.config.max_bootstrap_tree_count)
        ),
        "max_report_table_rows": (
            ""
            if report.config.max_report_table_rows is None
            else str(report.config.max_report_table_rows)
        ),
        "memory_warning_threshold_bytes": (
            ""
            if report.config.memory_warning_threshold_bytes is None
            else str(report.config.memory_warning_threshold_bytes)
        ),
        "budget_warning_count": str(
            len(bootstrap_artifacts.budget_report.warning_messages)
        ),
        "config_check_count": str(len(report.config_audit_rows)),
        "scientific_finding_count": str(scientific_finding_count),
    }
    return write_taxon_rows(path, columns=list(row.keys()), rows=[row])


def _write_resource_observation_table(
    path: Path,
    *,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    bootstrap_artifacts: BootstrapTreeSetArtifactReport,
) -> Path:
    row = {
        "dataset_id": report.dataset.dataset_id,
        "timeout_seconds": _format_number(report.config.timeout_seconds),
        "workflow_runtime_seconds": _format_number(
            report.fasta_to_tree.runtime_seconds
        ),
        "bootstrap_review_runtime_seconds": _format_number(
            bootstrap_artifacts.summary_report.processing.runtime_seconds
        ),
        "bootstrap_review_peak_memory_bytes": str(
            bootstrap_artifacts.summary_report.processing.peak_memory_bytes
        ),
        "budget_warning_count": str(
            len(bootstrap_artifacts.budget_report.warning_messages)
        ),
        "budget_warnings": " | ".join(
            bootstrap_artifacts.budget_report.warning_messages
        ),
    }
    return write_taxon_rows(path, columns=list(row.keys()), rows=[row])



def _write_comparative_report(
    path: Path,
    *,
    summary_row: ComparativeAnalysisSummaryRow,
    coefficient_rows: list[ComparativeCoefficientTableRow],
    residual_rows: list[ComparativeResidualTableRow],
    signal_row: ComparativeSignalTableRow,
    interpretation_rows: list[ComparativeInterpretationRow],
    branch_repairs: list[RabiesComparativeBranchRepair],
) -> Path:
    key_claim = next(
        (
            row.claim
            for row in interpretation_rows
            if row.topic == "coefficient" and "nominally supported" in row.claim
        ),
        "no coefficient reached nominal support",
    )
    html = "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Rabies Comparative Report</title>",
            "  <style>",
            "    body { font-family: Georgia, 'Times New Roman', serif; margin: 0; background: linear-gradient(180deg, #f4f1ea 0%, #f3f7ef 100%); color: #173024; }",
            "    main { max-width: 1040px; margin: 0 auto; padding: 24px; }",
            "    h1, h2 { margin: 0 0 10px; }",
            "    p { line-height: 1.55; }",
            "    .panel { background: rgba(255,255,255,0.88); border: 1px solid rgba(23,48,36,0.12); border-radius: 18px; padding: 18px; margin-top: 18px; box-shadow: 0 14px 36px rgba(23,48,36,0.08); }",
            "    .cards { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 14px; margin-top: 18px; }",
            "    .card { background: rgba(255,255,255,0.88); border: 1px solid rgba(23,48,36,0.12); border-radius: 18px; padding: 16px; }",
            "    .label { color: #5f7469; font-size: 13px; text-transform: uppercase; letter-spacing: 0.04em; }",
            "    .value { display: block; font-size: 22px; margin-top: 6px; }",
            "    table { width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 14px; }",
            "    th, td { border-bottom: 1px solid rgba(23,48,36,0.10); padding: 8px 10px; text-align: left; vertical-align: top; }",
            "    th { color: #365443; }",
            "    ul { margin: 8px 0 0 18px; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Rabies Comparative Report</h1>",
            "  <p>This comparative section asks whether the host-associated lineages in the rabies demonstration tree are associated with a consistent eastward geographic placement when geography is summarized as regional longitude. The result is interpretive evidence, not causal proof, and it inherits the small-panel limits of this dataset.</p>",
            '  <section class="cards">',
            f'    <div class="card"><span class="label">formula</span><span class="value">{escape(summary_row.formula)}</span></div>',
            f'    <div class="card"><span class="label">analysis taxa</span><span class="value">{summary_row.analysis_taxa}</span></div>',
            f'    <div class="card"><span class="label">selected trait model</span><span class="value">{escape(summary_row.selected_model)}</span></div>',
            f'    <div class="card"><span class="label">pgls r-squared</span><span class="value">{_format_number(summary_row.pgls_r_squared)}</span></div>',
            "  </section>",
            '  <section class="panel">',
            "    <h2>Question and Answer</h2>",
            _html_list(
                [
                    "Question: does host association coincide with a consistent longitudinal shift in this rabies panel?",
                    f"Answer: {key_claim}.",
                    (
                        f"Phylogenetic signal remains strong for the response trait "
                        f"(Blomberg's K {_format_number(signal_row.blombergs_k)}, "
                        f"Pagel's lambda {_format_number(signal_row.pagels_lambda)})."
                    ),
                    (
                        "Interpret the coefficient evidence cautiously because the "
                        "residual diagnostics retain review warnings and the sample "
                        "is intentionally compact."
                    ),
                ]
            ),
            "  </section>",
            '  <section class="panel">',
            "    <h2>Coefficient Summary</h2>",
            _table(
                headers=[
                    "term",
                    "estimate",
                    "standard_error",
                    "p_value",
                    "significant",
                ],
                rows=[
                    [
                        row.term,
                        _format_number(row.estimate),
                        _format_number(row.standard_error),
                        _format_number(row.p_value),
                        "true" if row.significant else "false",
                    ]
                    for row in coefficient_rows
                ],
            ),
            "  </section>",
            '  <section class="panel">',
            "    <h2>Residual Diagnostics</h2>",
            _table(
                headers=[
                    "analysis",
                    "residual_variance",
                    "max_abs_standardized_residual",
                    "phylogenetic_residual_lambda",
                    "warnings",
                ],
                rows=[
                    [
                        row.analysis,
                        _format_number(row.residual_variance),
                        _format_number(row.max_abs_standardized_residual),
                        _format_number(row.phylogenetic_residual_lambda),
                        "; ".join(row.warnings),
                    ]
                    for row in residual_rows
                ],
            ),
            "  </section>",
            '  <section class="panel">',
            "    <h2>Comparative Tree Adjustments</h2>",
            _html_list(
                [
                    "The comparative fit uses the rooted demonstration tree after flooring any nonpositive nonroot branch lengths to a tiny positive value.",
                    f"Adjusted branch count: {len(branch_repairs)}",
                ]
            ),
            _table(
                headers=[
                    "node_label",
                    "original_branch_length",
                    "repaired_branch_length",
                    "reason",
                ],
                rows=[
                    [
                        row.node_label,
                        _format_number(row.original_branch_length),
                        _format_number(row.repaired_branch_length),
                        row.reason,
                    ]
                    for row in branch_repairs
                ]
                or [["", "", "", "no branch-length repair was required"]],
            ),
            "  </section>",
            '  <section class="panel">',
            "    <h2>Interpretation Ledger</h2>",
            _table(
                headers=["topic", "claim", "evidence", "caution"],
                rows=[
                    [row.topic, row.claim, row.evidence, row.caution]
                    for row in interpretation_rows
                ],
            ),
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )
    path.write_text(html + "\n", encoding="utf-8")
    return path


def _write_comparative_manifest(
    path: Path,
    *,
    comparative_summary_row: ComparativeAnalysisSummaryRow,
    branch_repairs: list[RabiesComparativeBranchRepair],
    output_paths: dict[str, Path],
) -> Path:
    payload = {
        "report_kind": "rabies_cross_host_geography_comparative_bundle",
        "metrics": {
            "response": comparative_summary_row.response,
            "formula": comparative_summary_row.formula,
            "analysis_taxa": comparative_summary_row.analysis_taxa,
            "selected_model": comparative_summary_row.selected_model,
            "pgls_lambda": comparative_summary_row.pgls_lambda,
            "pgls_r_squared": comparative_summary_row.pgls_r_squared,
            "branch_repair_count": len(branch_repairs),
        },
        "output_checksums": {
            key: _checksum(value) for key, value in output_paths.items()
        },
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _write_manifest(
    path: Path,
    *,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    comparative_summary_row: ComparativeAnalysisSummaryRow,
    bootstrap_artifacts: BootstrapTreeSetArtifactReport,
    bootstrap_tree_comparison_report: ComparisonReportBuildResult,
    clade_row_count: int,
    scientific_finding_count: int,
    bundle_paths: dict[str, Path],
) -> Path:
    manifest = {
        "report_kind": "rabies_cross_host_geography_workflow_bundle",
        "dataset_id": report.dataset.dataset_id,
        "input_checksums": {
            "workflow-config.json": _checksum(report.dataset.workflow_config_path),
            "sequences.fasta": _checksum(report.dataset.sequences_path),
            "metadata.csv": _checksum(report.dataset.metadata_path),
            "region-centroids.csv": _checksum(report.dataset.centroids_path),
        },
        "output_checksums": {
            key: _checksum(value) for key, value in bundle_paths.items()
        },
        "metrics": {
            "sequence_count": report.dataset.sequence_count,
            "selected_model": report.fasta_to_tree.selected_model,
            "minimum_support": report.fasta_to_tree.support_summary.minimum_support,
            "maximum_support": report.fasta_to_tree.support_summary.maximum_support,
            "host_switch_count": report.host_switching.summary.host_switch_count,
            "migration_event_count": report.biogeography_report.event_report.summary.event_count,
            "root_host": report.host_switching.summary.root_host,
            "root_region": report.biogeography_report.state_report.summary.root_region,
            "clade_row_count": clade_row_count,
            "bootstrap_tree_count": bootstrap_artifacts.summary_report.tree_count,
            "bootstrap_unstable_branch_count": (
                bootstrap_artifacts.summary_report.unstable_branch_count
            ),
            "bootstrap_consensus_rooted_rf_distance": (
                bootstrap_tree_comparison_report.topology.rooted_robinson_foulds_distance
            ),
            "bootstrap_consensus_same_unrooted_topology": (
                bootstrap_tree_comparison_report.topology.same_unrooted_topology
            ),
            "bootstrap_consensus_high_support_conflict_count": len(
                [
                    row
                    for row in bootstrap_tree_comparison_report.support.conflicting_clades
                    if row.conflict_classification == "high_support_conflict"
                ]
            ),
            "comparative_selected_model": comparative_summary_row.selected_model,
            "comparative_pgls_lambda": comparative_summary_row.pgls_lambda,
            "comparative_pgls_r_squared": comparative_summary_row.pgls_r_squared,
            "conclusion_stable_count": (
                report.conclusion_stability_report.summary.stable_count
            ),
            "conclusion_weak_count": (
                report.conclusion_stability_report.summary.weak_count
            ),
            "conclusion_unstable_count": (
                report.conclusion_stability_report.summary.unstable_count
            ),
            "config_check_count": len(report.config_audit_rows),
            "scientific_finding_count": scientific_finding_count,
        },
    }
    path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _write_integrated_report(
    path: Path,
    *,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    workflow_summary_path: Path,
    bootstrap_artifacts: BootstrapTreeSetArtifactReport,
    bootstrap_tree_comparison_report: ComparisonReportBuildResult,
    clade_row_count: int,
    comparative_summary_row: ComparativeAnalysisSummaryRow,
    comparative_interpretation_rows: list[ComparativeInterpretationRow],
    comparative_branch_repair_count: int,
    scientific_finding_rows: list[RabiesScientificFindingRow],
    max_report_table_rows: int | None,
) -> Path:
    support_summary = report.fasta_to_tree.support_summary
    host_summary = report.host_switching.summary
    geography_summary = report.biogeography_report.state_report.summary
    migration_summary = report.biogeography_report.event_report.summary
    bootstrap_summary = bootstrap_artifacts.summary_report
    core_question = (
        "Do the host-associated rabies lineages in this compact panel occupy one "
        "distinct geographic regime while retaining one coherent phylogenetic signal?"
    )
    core_answer = next(
        (
            row.claim
            for row in comparative_interpretation_rows
            if row.topic == "coefficient" and "nominally supported" in row.claim
        ),
        "the comparative layer did not recover a nominally supported host effect",
    )
    html = "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Bijux Rabies Host and Geography Workflow</title>",
            "  <style>",
            "    body { font-family: Georgia, 'Times New Roman', serif; margin: 0; background: linear-gradient(180deg, #f4f1ea 0%, #e7efe7 100%); color: #163222; }",
            "    main { max-width: 1360px; margin: 0 auto; padding: 24px; }",
            "    h1 { margin: 0 0 8px; font-size: 34px; }",
            "    h2 { margin: 0 0 10px; font-size: 22px; }",
            "    p { line-height: 1.55; }",
            "    .cards { display: grid; grid-template-columns: repeat(6, minmax(0, 1fr)); gap: 14px; margin: 18px 0 24px; }",
            "    .card, .panel { background: rgba(255,255,255,0.86); border: 1px solid rgba(22,50,34,0.12); border-radius: 18px; padding: 18px; box-shadow: 0 16px 42px rgba(22,50,34,0.08); }",
            "    .label { color: #5b7466; font-size: 13px; text-transform: uppercase; letter-spacing: 0.04em; }",
            "    .card strong { display: block; font-size: 21px; margin-top: 6px; }",
            "    .grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }",
            "    .full { grid-column: 1 / -1; }",
            "    .figure-shell { overflow: auto; }",
            "    .figure-shell img { width: 100%; height: auto; display: block; }",
            "    table { width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 14px; }",
            "    th, td { border-bottom: 1px solid rgba(22,50,34,0.10); padding: 8px 10px; text-align: left; vertical-align: top; }",
            "    th { color: #365443; }",
            "    ul { margin: 8px 0 0 18px; }",
            "    a { color: #16543a; }",
            "    iframe { width: 100%; min-height: 760px; border: 1px solid rgba(22,50,34,0.12); border-radius: 14px; background: white; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Bijux Rabies Host and Geography Workflow</h1>",
            "  <p>Complete end-to-end review for one real rabies nucleoprotein panel. The workflow starts from raw sequences plus combined host and geography metadata, validates the FASTA surface, aligns and trims the panel, infers a bootstrap-supported maximum-likelihood tree, roots that tree on one explicit outgroup, summarizes bootstrap topology uncertainty, extracts clades, reconstructs host and geographic histories, and fits one comparative model over a derived geographic trait.</p>",
            '  <section class="cards">',
            f'    <div class="card"><span class="label">sequences</span><strong>{report.dataset.sequence_count}</strong></div>',
            f'    <div class="card"><span class="label">selected model</span><strong>{escape(report.fasta_to_tree.selected_model)}</strong></div>',
            f'    <div class="card"><span class="label">aligned quality</span><strong>{_format_number(report.aligned_quality.quality_score)}</strong></div>',
            f'    <div class="card"><span class="label">trimmed quality</span><strong>{_format_number(report.trimmed_quality.quality_score)}</strong></div>',
            f'    <div class="card"><span class="label">root host</span><strong>{escape(host_summary.root_host)}</strong></div>',
            f'    <div class="card"><span class="label">root region</span><strong>{escape(geography_summary.root_region)}</strong></div>',
            "  </section>",
            '  <section class="panel">',
            "    <h2>Scientific Question</h2>",
            f"    <p>{escape(core_question)}</p>",
            f"    <p><strong>Working answer:</strong> {escape(core_answer)}. The comparative layer selects {escape(comparative_summary_row.selected_model)} as the better continuous-trait surface, but the residual diagnostics remain cautionary and the panel is intentionally small.</p>",
            _html_list(
                [
                    f"FASTA validation resolved the raw sequence type as {report.fasta_to_tree.sequence_type}.",
                    f"Bootstrap support spans {_support_range_text(support_summary.minimum_support, support_summary.maximum_support)} across the final rooted tree.",
                    f"Host reconstruction inferred {host_summary.host_switch_count} host-switch branches, with {host_summary.certain_host_switch_count} certain and {host_summary.uncertain_host_switch_count} uncertain changes.",
                    f"Geographic reconstruction inferred {migration_summary.event_count} migration events across {geography_summary.changed_branch_count} changed branches.",
                    f"Bootstrap replicate review retained {bootstrap_summary.tree_count} trees across {bootstrap_summary.diversity.rooted_topology_count} rooted topologies.",
                    (
                        "The rooted ML tree versus bootstrap consensus comparison "
                        f"returned rooted RF distance {bootstrap_tree_comparison_report.topology.rooted_robinson_foulds_distance}."
                    ),
                    f"The clade table contains {clade_row_count} node rows and the comparative tree required {comparative_branch_repair_count} explicit branch-length repair(s).",
                    (
                        "Bootstrap review emitted budget warnings: "
                        + "; ".join(bootstrap_artifacts.budget_report.warning_messages)
                    )
                    if bootstrap_artifacts.budget_report.warning_messages
                    else (
                        "Configured workflow budgets covered the bootstrap review "
                        "without tree-count failure or peak-memory warning."
                    ),
                ]
            ),
            "  </section>",
            '  <section class="grid" style="margin-top: 20px;">',
            '    <section class="panel">',
            "      <h2>Sequence-to-Tree Outputs</h2>",
            _html_list(
                [
                    "input validation: input-validation.tsv",
                    "alignment quality: alignment-quality.tsv",
                    "alignment sequence ranking: alignment-sequence-ranking.tsv",
                    "alignment: rabies-cross-host-geography-panel.aln",
                    "trimmed alignment: rabies-cross-host-geography-panel.trimmed.aln",
                    "rooted tree: rabies-cross-host-geography-panel.rooted.tree",
                    "support table: rabies-cross-host-geography-panel.support.tsv",
                    "workflow summary: workflow-summary.tsv",
                    "resource observations: resource-observations.tsv",
                ]
            ),
            _support_table(
                report.fasta_to_tree,
                max_rows=max_report_table_rows,
            ),
            "    </section>",
            '    <section class="panel">',
            "      <h2>Bootstrap and Clade Review</h2>",
            _html_list(
                [
                    f"bootstrap tree count: {bootstrap_summary.tree_count}",
                    f"rooted topology count: {bootstrap_summary.diversity.rooted_topology_count}",
                    f"unstable branch count: {bootstrap_summary.unstable_branch_count}",
                    f"clade row count: {clade_row_count}",
                    (
                        "see bootstrap-review/ for consensus, clade frequencies, "
                        "instability, distances, topology clusters, and rooted-tree comparison"
                    ),
                ]
            ),
            _table(
                headers=[
                    "tree_count",
                    "rooted_topology_count",
                    "dominant_topology_frequency",
                    "effective_topology_count",
                    "unstable_branch_count",
                ],
                rows=[
                    [
                        str(bootstrap_summary.tree_count),
                        str(bootstrap_summary.diversity.rooted_topology_count),
                        _format_number(
                            bootstrap_summary.diversity.dominant_topology_frequency
                        ),
                        _format_number(
                            bootstrap_summary.diversity.effective_topology_count
                        ),
                        str(bootstrap_summary.unstable_branch_count),
                    ]
                ],
            ),
            "    </section>",
            '    <section class="panel">',
            "      <h2>Host Switching</h2>",
            _html_list(
                [
                    f"workflow trait: {report.dataset.host_trait}",
                    f"root host confidence: {_format_number(host_summary.root_confidence)}",
                    f"host-switch rows: {len(report.host_switching.count_rows)}",
                    "see host-switch-summary.tsv, host-state-nodes.tsv, host-switch-branches.tsv, and host-switch-counts.tsv",
                ]
            ),
            _host_count_table(
                report.host_switching,
                max_rows=max_report_table_rows,
            ),
            "    </section>",
            '    <section class="panel">',
            "      <h2>Comparative Layer</h2>",
            _html_list(
                [
                    f"formula: {comparative_summary_row.formula}",
                    f"selected trait model: {comparative_summary_row.selected_model}",
                    f"PGLS lambda: {_format_number(comparative_summary_row.pgls_lambda)}",
                    f"PGLS r-squared: {_format_number(comparative_summary_row.pgls_r_squared)}",
                    "see comparative/ for coefficients, model comparison, diagnostics, signal summary, and interpretation tables",
                ]
            ),
            _table(
                headers=["topic", "claim", "evidence"],
                rows=[
                    [row.topic, row.claim, row.evidence]
                    for row in comparative_interpretation_rows[:5]
                ],
                max_rows=max_report_table_rows,
            ),
            "    </section>",
            '    <section class="panel full">',
            "      <h2>Biogeography</h2>",
            '      <p>The bundle includes the detailed biogeography package at <a href="biogeography/biogeography-report.html">biogeography/biogeography-report.html</a> together with the ancestral-region tree SVG and the self-contained geographic map.</p>',
            '      <div class="grid">',
            '        <div class="panel">',
            "          <h2>Ancestral-Region Tree</h2>",
            '          <div class="figure-shell">',
            '            <img src="biogeography/ancestral-region-tree.svg" alt="Ancestral region tree">',
            "          </div>",
            "        </div>",
            '        <div class="panel">',
            "          <h2>Geographic Map</h2>",
            '          <iframe src="biogeography/geographic-region-map.html" title="Geographic region map"></iframe>',
            "        </div>",
            "      </div>",
            _migration_event_table(
                report.biogeography_report,
                max_rows=max_report_table_rows,
            ),
            "    </section>",
            '    <section class="panel full">',
            "      <h2>Scientific Findings Ledger</h2>",
            _table(
                headers=[
                    "finding_id",
                    "question",
                    "claim",
                    "evidence",
                    "caution",
                    "source_artifact",
                ],
                rows=[
                    [
                        row.finding_id,
                        row.question,
                        row.claim,
                        row.evidence,
                        row.caution,
                        row.source_artifact,
                    ]
                    for row in scientific_finding_rows
                ],
                max_rows=max_report_table_rows,
            ),
            "    </section>",
            "  </section>",
            '  <section class="panel" style="margin-top: 20px;">',
            "    <h2>Key Files</h2>",
            _html_list(
                [
                    f'<a href="{workflow_summary_path.name}">{workflow_summary_path.name}</a>',
                    '<a href="resource-observations.tsv">resource-observations.tsv</a>',
                    '<a href="workflow-config-audit.tsv">workflow-config-audit.tsv</a>',
                    '<a href="clade-table.tsv">clade-table.tsv</a>',
                    '<a href="bootstrap-review/bootstrap-review.summary.tsv">bootstrap-review/bootstrap-review.summary.tsv</a>',
                    '<a href="bootstrap-review/rooted-tree-vs-bootstrap-consensus.summary.tsv">bootstrap-review/rooted-tree-vs-bootstrap-consensus.summary.tsv</a>',
                    '<a href="comparative/comparative-report.html">comparative/comparative-report.html</a>',
                    '<a href="comparative/interpretation-table.tsv">comparative/interpretation-table.tsv</a>',
                    '<a href="scientific-findings.tsv">scientific-findings.tsv</a>',
                    '<a href="host-switch-summary.tsv">host-switch-summary.tsv</a>',
                    '<a href="biogeography/event-table.tsv">biogeography/event-table.tsv</a>',
                    '<a href="rabies-cross-host-geography.manifest.json">rabies-cross-host-geography.manifest.json</a>',
                ]
            ),
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )
    path.write_text(html + "\n", encoding="utf-8")
    return path


def _support_table(
    report: FastaToTreeWorkflowReport,
    *,
    max_rows: int | None = None,
) -> str:
    return _table(
        headers=["node", "descendant_taxa", "support", "support_fraction"],
        rows=[
            [
                row.node,
                ", ".join(row.descendant_taxa),
                _format_number(row.support),
                _format_number(row.support_fraction),
            ]
            for row in report.support_rows
        ],
        max_rows=max_rows,
    )


def _host_count_table(
    report: HostSwitchingReport,
    *,
    max_rows: int | None = None,
) -> str:
    return _table(
        headers=[
            "transition",
            "certain_switch_count",
            "uncertain_switch_count",
            "total_switch_count",
        ],
        rows=[
            [
                row.transition,
                str(row.certain_switch_count),
                str(row.uncertain_switch_count),
                str(row.total_switch_count),
            ]
            for row in report.count_rows
        ],
        max_rows=max_rows,
    )


def _migration_event_table(
    report: BiogeographyReportPackageResult,
    *,
    max_rows: int | None = None,
) -> str:
    return _table(
        headers=[
            "branch_id",
            "source_region",
            "target_region",
            "support",
            "midpoint_depth",
        ],
        rows=[
            [
                row.branch_id,
                row.source_region,
                row.target_region,
                _format_number(row.support),
                _format_number(row.midpoint_depth),
            ]
            for row in report.event_report.event_rows
        ],
        max_rows=max_rows,
    )
