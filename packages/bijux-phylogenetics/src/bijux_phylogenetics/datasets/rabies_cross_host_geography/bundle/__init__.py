from __future__ import annotations

from pathlib import Path
import shutil

from ..models import (
    RabiesCrossHostGeographyPanelWorkflowBundle,
    RabiesCrossHostGeographyPanelWorkflowReport,
)
from .bootstrap_bundle import _write_bootstrap_review_artifacts
from .comparative_bundle import _write_comparative_bundle_artifacts
from .finalization import _write_final_bundle_artifacts
from .host_geography import _write_host_geography_artifacts
from .stability_bundle import _write_conclusion_stability_artifacts
from .tree_inference import _write_tree_inference_artifacts


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

    tree_inference = _write_tree_inference_artifacts(output_root, report)
    bootstrap_review = _write_bootstrap_review_artifacts(
        output_root,
        report,
        rooted_tree_path=tree_inference.tree_path,
    )

    host_geography = _write_host_geography_artifacts(output_root, report)
    comparative_bundle = _write_comparative_bundle_artifacts(output_root, report)
    conclusion_stability = _write_conclusion_stability_artifacts(output_root, report)
    final_bundle = _write_final_bundle_artifacts(
        output_root,
        report=report,
        clade_row_count=len(tree_inference.clade_report.rows),
        bootstrap_artifacts=bootstrap_review.artifact_report,
        bootstrap_tree_comparison_report=bootstrap_review.comparison_report,
        comparative_summary_row=comparative_bundle.summary_row,
        comparative_interpretation_rows=comparative_bundle.interpretation_rows,
        comparative_branch_repair_count=len(report.comparative_branch_repairs),
        bundle_paths={
            "config_audit": tree_inference.config_audit_path,
            "resolved_config": tree_inference.resolved_config_path,
            "input_validation": tree_inference.input_validation_path,
            "alignment_quality": tree_inference.alignment_quality_path,
            "alignment_sequence_ranking": tree_inference.alignment_sequence_ranking_path,
            "alignment": tree_inference.alignment_path,
            "trimmed_alignment": tree_inference.trimmed_alignment_path,
            "rooted_tree": tree_inference.tree_path,
            "rooting_report": tree_inference.rooting_report_path,
            "model_table": tree_inference.model_table_path,
            "support_table": tree_inference.support_table_path,
            "clade_table": tree_inference.clade_table_path,
            "bootstrap_summary": bootstrap_review.summary_path,
            "bootstrap_consensus_tree": bootstrap_review.artifact_report.output_paths[
                "consensus_tree"
            ],
            "bootstrap_clade_frequencies": bootstrap_review.artifact_report.output_paths[
                "clade_frequencies"
            ],
            "bootstrap_unstable_branches": bootstrap_review.artifact_report.output_paths[
                "unstable_branches"
            ],
            "bootstrap_unstable_clades": bootstrap_review.artifact_report.output_paths[
                "unstable_clades"
            ],
            "bootstrap_distance_matrix": bootstrap_review.artifact_report.output_paths[
                "distance_matrix"
            ],
            "bootstrap_topology_clusters": bootstrap_review.artifact_report.output_paths[
                "topology_clusters"
            ],
            "bootstrap_tree_comparison_summary": bootstrap_review.comparison_summary_path,
            "bootstrap_tree_comparison_table": bootstrap_review.comparison_table_path,
            "bootstrap_tree_comparison_report": (
                bootstrap_review.comparison_report.output_path
            ),
            "host_switch_summary": host_geography.host_switch_summary_path,
            "host_state_nodes": host_geography.host_state_nodes_path,
            "host_switch_branches": host_geography.host_switch_branches_path,
            "host_switch_counts": host_geography.host_switch_counts_path,
            "host_switch_fits": host_geography.host_switch_fits_path,
            "host_switch_unsupported": host_geography.host_switch_unsupported_path,
            "host_switch_exclusions": host_geography.host_switch_exclusions_path,
            "biogeography_report": host_geography.biogeography_report_path,
            "biogeography_tree_figure": host_geography.biogeography_tree_figure_path,
            "biogeography_map": host_geography.biogeography_map_path,
            "comparative_traits": comparative_bundle.traits_path,
            "comparative_tree": comparative_bundle.tree_path,
            "comparative_repairs": tree_inference.comparative_repairs_path,
            "comparative_report": comparative_bundle.report_path,
            "comparative_summary": comparative_bundle.summary_path,
            "comparative_coefficients": comparative_bundle.coefficients_path,
            "comparative_residuals": comparative_bundle.residuals_path,
            "comparative_signal": comparative_bundle.signal_path,
            "comparative_model_comparison": comparative_bundle.model_comparison_path,
            "comparative_interpretation": comparative_bundle.interpretation_path,
            "comparative_audit": comparative_bundle.audit_path,
            "comparative_contrasts": comparative_bundle.contrasts_path,
            "comparative_model_matrix": comparative_bundle.model_matrix_path,
            "comparative_categorical_contrasts": comparative_bundle.categorical_contrasts_path,
            "comparative_lambda_profile": comparative_bundle.lambda_profile_path,
            "comparative_manifest": comparative_bundle.manifest_path,
            "conclusion_stability_summary": conclusion_stability.summary_path,
            "key_clade_stability": conclusion_stability.key_clade_stability_path,
            "support_value_stability": conclusion_stability.support_value_stability_path,
            "ancestral_state_stability": conclusion_stability.ancestral_state_stability_path,
            "comparative_coefficient_stability": (
                conclusion_stability.comparative_coefficient_stability_path
            ),
            "conclusion_stability_report": conclusion_stability.report_path,
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
        clade_row_count=len(tree_inference.clade_report.rows),
        bootstrap_tree_count=bootstrap_review.artifact_report.summary_report.tree_count,
        bootstrap_topology_count=(
            bootstrap_review.artifact_report.summary_report.diversity.rooted_topology_count
        ),
        bootstrap_unstable_branch_count=(
            bootstrap_review.artifact_report.summary_report.unstable_branch_count
        ),
        bootstrap_consensus_rooted_rf_distance=(
            bootstrap_review.comparison_report.topology.rooted_robinson_foulds_distance
        ),
        bootstrap_consensus_same_unrooted_topology=(
            bootstrap_review.comparison_report.topology.same_unrooted_topology
        ),
        bootstrap_consensus_high_support_conflict_count=len(
            [
                row
                for row in bootstrap_review.comparison_report.support.conflicting_clades
                if row.conflict_classification == "high_support_conflict"
            ]
        ),
        bootstrap_consensus_branch_score_distance=(
            bootstrap_review.comparison_report.branch_lengths.branch_score.branch_score_distance
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
        comparative_selected_model=comparative_bundle.summary_row.selected_model,
        comparative_response=comparative_bundle.summary_row.response,
        comparative_formula=comparative_bundle.summary_row.formula,
        comparative_pgls_lambda=comparative_bundle.summary_row.pgls_lambda,
        comparative_pgls_r_squared=comparative_bundle.summary_row.pgls_r_squared,
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
            bootstrap_review.artifact_report.summary_report.processing.runtime_seconds
        ),
        bootstrap_review_peak_memory_bytes=(
            bootstrap_review.artifact_report.summary_report.processing.peak_memory_bytes
        ),
        budget_warning_count=len(
            bootstrap_review.artifact_report.budget_report.warning_messages
        ),
        config_check_count=len(report.config_audit_rows),
        scientific_finding_count=len(final_bundle.scientific_finding_rows),
        workflow_summary_path=final_bundle.workflow_summary_path,
        resource_observations_path=final_bundle.resource_observations_path,
        config_audit_path=tree_inference.config_audit_path,
        resolved_config_path=tree_inference.resolved_config_path,
        input_validation_path=tree_inference.input_validation_path,
        alignment_quality_path=tree_inference.alignment_quality_path,
        alignment_sequence_ranking_path=tree_inference.alignment_sequence_ranking_path,
        alignment_path=tree_inference.alignment_path,
        trimmed_alignment_path=tree_inference.trimmed_alignment_path,
        tree_path=tree_inference.tree_path,
        rooting_report_path=tree_inference.rooting_report_path,
        model_table_path=tree_inference.model_table_path,
        support_table_path=tree_inference.support_table_path,
        log_path=tree_inference.log_path,
        manifest_path=tree_inference.manifest_path,
        engine_artifact_root=tree_inference.engine_artifact_root,
        clade_table_path=tree_inference.clade_table_path,
        bootstrap_output_root=bootstrap_review.output_root,
        bootstrap_summary_path=bootstrap_review.summary_path,
        bootstrap_consensus_tree_path=bootstrap_review.artifact_report.output_paths[
            "consensus_tree"
        ],
        bootstrap_clade_frequencies_path=bootstrap_review.artifact_report.output_paths[
            "clade_frequencies"
        ],
        bootstrap_unstable_branches_path=bootstrap_review.artifact_report.output_paths[
            "unstable_branches"
        ],
        bootstrap_unstable_clades_path=bootstrap_review.artifact_report.output_paths[
            "unstable_clades"
        ],
        bootstrap_distance_matrix_path=bootstrap_review.artifact_report.output_paths[
            "distance_matrix"
        ],
        bootstrap_topology_clusters_path=bootstrap_review.artifact_report.output_paths[
            "topology_clusters"
        ],
        bootstrap_tree_comparison_summary_path=(
            bootstrap_review.comparison_summary_path
        ),
        bootstrap_tree_comparison_table_path=bootstrap_review.comparison_table_path,
        bootstrap_tree_comparison_report_path=(
            bootstrap_review.comparison_report.output_path
        ),
        host_switch_summary_path=host_geography.host_switch_summary_path,
        host_state_nodes_path=host_geography.host_state_nodes_path,
        host_switch_branches_path=host_geography.host_switch_branches_path,
        host_switch_counts_path=host_geography.host_switch_counts_path,
        host_switch_fits_path=host_geography.host_switch_fits_path,
        host_switch_unsupported_path=host_geography.host_switch_unsupported_path,
        host_switch_exclusions_path=host_geography.host_switch_exclusions_path,
        biogeography_output_root=host_geography.biogeography_output_root,
        biogeography_report_path=host_geography.biogeography_report_path,
        biogeography_tree_figure_path=host_geography.biogeography_tree_figure_path,
        biogeography_map_path=host_geography.biogeography_map_path,
        comparative_traits_path=comparative_bundle.traits_path,
        comparative_tree_path=comparative_bundle.tree_path,
        comparative_repairs_path=tree_inference.comparative_repairs_path,
        comparative_output_root=comparative_bundle.output_root,
        comparative_report_path=comparative_bundle.report_path,
        comparative_summary_path=comparative_bundle.summary_path,
        comparative_coefficients_path=comparative_bundle.coefficients_path,
        comparative_residuals_path=comparative_bundle.residuals_path,
        comparative_signal_path=comparative_bundle.signal_path,
        comparative_model_comparison_path=comparative_bundle.model_comparison_path,
        comparative_interpretation_path=comparative_bundle.interpretation_path,
        comparative_audit_path=comparative_bundle.audit_path,
        comparative_contrasts_path=comparative_bundle.contrasts_path,
        comparative_model_matrix_path=comparative_bundle.model_matrix_path,
        comparative_categorical_contrasts_path=comparative_bundle.categorical_contrasts_path,
        comparative_lambda_profile_path=comparative_bundle.lambda_profile_path,
        comparative_manifest_path=comparative_bundle.manifest_path,
        conclusion_stability_output_root=conclusion_stability.output_root,
        conclusion_stability_summary_path=conclusion_stability.summary_path,
        key_clade_stability_path=conclusion_stability.key_clade_stability_path,
        support_value_stability_path=conclusion_stability.support_value_stability_path,
        ancestral_state_stability_path=conclusion_stability.ancestral_state_stability_path,
        comparative_coefficient_stability_path=(
            conclusion_stability.comparative_coefficient_stability_path
        ),
        conclusion_stability_report_path=conclusion_stability.report_path,
        scientific_findings_path=final_bundle.scientific_findings_path,
        final_report_path=final_bundle.final_report_path,
        final_manifest_path=final_bundle.final_manifest_path,
    )
