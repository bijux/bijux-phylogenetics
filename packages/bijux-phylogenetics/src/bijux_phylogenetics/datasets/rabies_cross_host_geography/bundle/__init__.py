# ruff: noqa: F401, F403, F405
from __future__ import annotations

import csv
from dataclasses import dataclass, replace
from pathlib import Path
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
from bijux_phylogenetics.comparative.pgls.categorical_contrasts import (
    PGLSCategoricalContrastReport,
    summarize_pgls_categorical_contrasts,
    write_pgls_categorical_contrast_table,
)
from bijux_phylogenetics.comparative.pgls.lambda_fit import (
    write_pgls_lambda_profile_table,
)
from bijux_phylogenetics.comparative.pgls.posterior_tree import (
    PosteriorTreePGLSReport,
    run_posterior_tree_pgls,
)
from bijux_phylogenetics.comparative.reporting.analysis_package import (
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
from bijux_phylogenetics.compare.presentation import (
    ComparisonReportBuildResult,
    build_tree_comparison_report,
)
from bijux_phylogenetics.compare.topology import write_tree_comparison_table
from bijux_phylogenetics.phylo.alignment import (
    AlignmentQualityReport,
    SequenceQualityRankingReport,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
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
from bijux_phylogenetics.engines.inference import (
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
    compute_clade_frequency_table,
    write_bootstrap_tree_set_artifacts,
)

from .bootstrap_review import (
    _stabilize_bundle_report_paths,
    _write_bootstrap_tree_comparison_summary,
    _write_stable_bootstrap_summary_table,
)
from .comparative_review import (
    _write_comparative_manifest,
    _write_comparative_report,
)
from .findings import (
    _build_scientific_finding_rows,
    _write_scientific_findings_table,
)
from .input_artifacts import (
    _copy_output,
    _write_alignment_quality_table,
    _write_comparative_branch_repairs_table,
    _write_input_validation_table,
    _write_resolved_workflow_config,
    _write_sequence_ranking_table,
    _write_workflow_config_audit_table,
)
from .package_ledger import (
    _write_manifest,
    _write_resource_observation_table,
    _write_workflow_summary_table,
)
from .integrated_report import _write_integrated_report
from ..models import *
from ..shared import _checksum, _format_number, _html_list, _support_range_text, _table

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


