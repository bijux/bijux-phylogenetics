from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.demo.shared import (
    count_expected_output_entries,
    count_expected_output_files,
    emit_demo_result,
    resolve_demo_runner,
)
from bijux_phylogenetics.datasets import (
    run_rabies_cross_host_panel_demo,
    run_rabies_geographic_transition_panel_demo,
    run_rabies_method_sensitivity_panel_demo,
)


def add_rabies_demo_commands(demo_subparsers: Any) -> None:
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
        "--variant-id",
        action="append",
        dest="variant_ids",
        help="Restrict the workflow to one declared variant id. Repeat to preserve a specific subset order.",
    )
    demo_rabies_method_sensitivity.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_rabies_method_sensitivity)


def run_rabies_demo_command(args: Any) -> int | None:
    if args.demo_command == "rabies-cross-host-panel":
        result = run_rabies_cross_host_panel_demo(args.out)
        outputs = [
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
        ]
        return emit_demo_result(
            args,
            outputs=outputs,
            metrics={
                "artifact_count": len(outputs),
                "taxon_count": result.dataset.taxon_count,
                "workflow_trait": result.dataset.workflow_trait,
                "observed_host_group_count": (result.dataset.observed_host_group_count),
                "analysis_constraint_mode": (
                    result.workflow_bundle.analysis_constraint_mode
                ),
                "root_host": result.workflow_bundle.root_host,
                "root_confidence": result.workflow_bundle.root_confidence,
                "host_switch_count": result.workflow_bundle.host_switch_count,
                "certain_host_switch_count": (
                    result.workflow_bundle.certain_host_switch_count
                ),
                "uncertain_host_switch_count": (
                    result.workflow_bundle.uncertain_host_switch_count
                ),
                "reference_output_count": count_expected_output_entries(
                    result.dataset_export.expected_output_root
                ),
            },
            data=result,
            output_root=result.output_root,
        )

    if args.demo_command == "rabies-geographic-transition-panel":
        result = run_rabies_geographic_transition_panel_demo(args.out)
        outputs = [
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
        ]
        return emit_demo_result(
            args,
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
                "changed_branch_count": result.workflow_bundle.changed_branch_count,
                "strongly_supported_transition_count": (
                    result.workflow_bundle.strongly_supported_transition_count
                ),
                "migration_event_count": (result.workflow_bundle.migration_event_count),
                "strongly_supported_migration_event_count": (
                    result.workflow_bundle.strongly_supported_migration_event_count
                ),
                "reference_output_count": count_expected_output_entries(
                    result.dataset_export.expected_output_root
                ),
            },
            data=result,
            output_root=result.output_root,
        )

    if args.demo_command == "rabies-cross-host-geography-panel":
        result = resolve_demo_runner("run_rabies_cross_host_geography_panel_demo")(
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
        outputs = [
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
            result.artifact_inventory_path,
            result.reproducibility_checklist_path,
            result.package_manifest_path,
        ]
        return emit_demo_result(
            args,
            outputs=outputs,
            metrics={
                "artifact_count": len(outputs),
                "sequence_count": result.dataset.sequence_count,
                "config_path": str(result.dataset_export.workflow_config_path),
                "biological_question": (
                    "Do the host-associated rabies lineages in this compact panel occupy one distinct geographic regime while retaining one coherent phylogenetic signal?"
                ),
                "short_answer": (
                    "The rooted panel remains anchored in bat and north_asia, and `host_group[canid]` shows a nominally supported positive longitude association under the selected comparative model, but the inference remains cautionary because the panel is intentionally compact."
                ),
                "host_trait": result.dataset.host_trait,
                "geography_trait": result.dataset.geography_trait,
                "selected_model": result.workflow_bundle.selected_model,
                "aligned_quality_score": result.workflow_bundle.aligned_quality_score,
                "trimmed_quality_score": result.workflow_bundle.trimmed_quality_score,
                "minimum_support": result.workflow_bundle.minimum_support,
                "maximum_support": result.workflow_bundle.maximum_support,
                "root_host": result.workflow_bundle.root_host,
                "root_region": result.workflow_bundle.root_region,
                "host_switch_count": result.workflow_bundle.host_switch_count,
                "migration_event_count": result.workflow_bundle.migration_event_count,
                "clade_row_count": result.workflow_bundle.clade_row_count,
                "bootstrap_tree_count": result.workflow_bundle.bootstrap_tree_count,
                "timeout_seconds": result.workflow_bundle.timeout_seconds,
                "max_bootstrap_tree_count": (
                    result.workflow_bundle.max_bootstrap_tree_count
                ),
                "max_report_table_rows": (result.workflow_bundle.max_report_table_rows),
                "budget_warning_count": result.workflow_bundle.budget_warning_count,
                "bootstrap_review_runtime_seconds": (
                    result.workflow_bundle.bootstrap_review_runtime_seconds
                ),
                "bootstrap_review_peak_memory_bytes": (
                    result.workflow_bundle.bootstrap_review_peak_memory_bytes
                ),
                "bootstrap_consensus_rooted_rf_distance": (
                    result.workflow_bundle.bootstrap_consensus_rooted_rf_distance
                ),
                "comparative_formula": result.workflow_bundle.comparative_formula,
                "comparative_selected_model": (
                    result.workflow_bundle.comparative_selected_model
                ),
                "conclusion_stable_count": (
                    result.workflow_bundle.conclusion_stable_count
                ),
                "conclusion_weak_count": result.workflow_bundle.conclusion_weak_count,
                "conclusion_unstable_count": (
                    result.workflow_bundle.conclusion_unstable_count
                ),
                "config_check_count": result.workflow_bundle.config_check_count,
                "scientific_finding_count": (
                    result.workflow_bundle.scientific_finding_count
                ),
                "package_artifact_count": (
                    sum(
                        1
                        for _ in result.artifact_inventory_path.open(
                            "r", encoding="utf-8"
                        )
                    )
                    - 1
                ),
                "package_checklist_item_count": (
                    sum(
                        1
                        for _ in result.reproducibility_checklist_path.open(
                            "r", encoding="utf-8"
                        )
                    )
                    - 1
                ),
                "reference_output_count": count_expected_output_files(
                    result.dataset_export.expected_output_root
                ),
            },
            data=result,
            output_root=result.output_root,
        )

    if args.demo_command != "rabies-method-sensitivity-panel":
        return None

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
        variant_ids=(
            tuple(args.variant_ids) if getattr(args, "variant_ids", None) else None
        ),
    )
    outputs = [
        result.dataset_export.readme_path,
        result.dataset_export.config_path,
        result.dataset_export.sequences_path,
        result.dataset_export.metadata_path,
        result.workflow_bundle.workflow_summary_path,
        result.workflow_bundle.variant_summary_path,
        result.workflow_bundle.parallel_summary_path,
        result.workflow_bundle.execution_record_path,
        result.workflow_bundle.preprocessing_comparison_path,
        result.workflow_bundle.stable_clades_path,
        result.workflow_bundle.changed_clades_path,
        result.workflow_bundle.conclusion_summary_path,
        result.workflow_bundle.config_path,
        result.workflow_bundle.manifest_path,
        result.workflow_bundle.report_manifest_path,
        result.workflow_bundle.slurm_job_plan_path,
        result.workflow_bundle.slurm_assumptions_path,
        result.workflow_bundle.slurm_summary_path,
        result.workflow_bundle.slurm_array_partitions_path,
        result.workflow_bundle.slurm_array_members_path,
        result.workflow_bundle.slurm_array_strategy_path,
        result.workflow_bundle.slurm_job_evidence_index_path,
        result.workflow_bundle.slurm_job_evidence_summary_path,
        result.workflow_bundle.slurm_storage_categories_path,
        result.workflow_bundle.slurm_storage_variants_path,
        result.workflow_bundle.slurm_storage_summary_path,
        result.workflow_bundle.slurm_storage_report_path,
        result.workflow_bundle.slurm_output_explosion_checks_path,
        result.workflow_bundle.slurm_output_explosion_variants_path,
        result.workflow_bundle.slurm_output_explosion_summary_path,
        result.workflow_bundle.slurm_output_explosion_report_path,
        result.workflow_bundle.slurm_tree_retention_checks_path,
        result.workflow_bundle.slurm_tree_retention_files_path,
        result.workflow_bundle.slurm_tree_retention_summary_path,
        result.workflow_bundle.slurm_tree_retention_report_path,
        result.workflow_bundle.slurm_merge_checks_path,
        result.workflow_bundle.slurm_merge_variants_path,
        result.workflow_bundle.slurm_merge_summary_path,
        result.workflow_bundle.slurm_merge_report_path,
        result.workflow_bundle.slurm_output_freshness_path,
        result.workflow_bundle.slurm_output_freshness_checks_path,
        result.workflow_bundle.slurm_output_freshness_summary_path,
        result.workflow_bundle.slurm_job_status_path,
        result.workflow_bundle.slurm_partition_status_path,
        result.workflow_bundle.slurm_workflow_status_path,
        result.workflow_bundle.slurm_failure_recovery_jobs_path,
        result.workflow_bundle.slurm_failure_recovery_partitions_path,
        result.workflow_bundle.slurm_failure_recovery_summary_path,
        result.workflow_bundle.slurm_failure_recovery_report_path,
        result.workflow_bundle.reproducibility_checks_path,
        result.workflow_bundle.reproducibility_variant_audit_path,
        result.workflow_bundle.reproducibility_audit_path,
        result.workflow_bundle.report_path,
        result.overview_path,
    ]
    return emit_demo_result(
        args,
        outputs=outputs,
        metrics={
            "artifact_count": len(outputs),
            "taxon_count": result.dataset.taxon_count,
            "variant_count": result.workflow_bundle.variant_count,
            "parallel_workers": result.workflow_bundle.parallel_workers,
            "execution_mode": result.workflow_bundle.execution_mode,
            "stable_clade_count": result.workflow_bundle.stable_clade_count,
            "changed_clade_count": result.workflow_bundle.changed_clade_count,
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
            "report_html_size_bytes": result.workflow_bundle.report_html_size_bytes,
            "report_linked_artifact_bytes": (
                result.workflow_bundle.report_linked_artifact_bytes
            ),
            "report_total_output_bytes": (
                result.workflow_bundle.report_total_output_bytes
            ),
            "slurm_job_count": result.workflow_bundle.slurm_job_count,
            "slurm_total_estimated_core_hours": (
                result.workflow_bundle.slurm_total_estimated_core_hours
            ),
            "slurm_maximum_estimated_memory_mib": (
                result.workflow_bundle.slurm_maximum_estimated_memory_mib
            ),
            "slurm_maximum_estimated_wallclock_minutes": (
                result.workflow_bundle.slurm_maximum_estimated_wallclock_minutes
            ),
            "slurm_total_estimated_scratch_mib": (
                result.workflow_bundle.slurm_total_estimated_scratch_mib
            ),
            "slurm_total_estimated_output_mib": (
                result.workflow_bundle.slurm_total_estimated_output_mib
            ),
            "slurm_array_partition_count": (
                result.workflow_bundle.slurm_array_partition_count
            ),
            "slurm_array_script_count": (
                result.workflow_bundle.slurm_array_script_count
            ),
            "slurm_array_largest_partition_size": (
                result.workflow_bundle.slurm_array_largest_partition_size
            ),
            "slurm_job_evidence_file_count": (
                result.workflow_bundle.slurm_job_evidence_file_count
            ),
            "slurm_job_evidence_total_runtime_seconds": (
                result.workflow_bundle.slurm_job_evidence_total_runtime_seconds
            ),
            "slurm_job_evidence_total_output_byte_count": (
                result.workflow_bundle.slurm_job_evidence_total_output_byte_count
            ),
            "slurm_storage_total_estimated_mib": (
                result.workflow_bundle.slurm_storage_total_estimated_mib
            ),
            "slurm_storage_output_byte_count": (
                result.workflow_bundle.slurm_storage_output_byte_count
            ),
            "slurm_storage_log_byte_count": (
                result.workflow_bundle.slurm_storage_log_byte_count
            ),
            "slurm_storage_tree_byte_count": (
                result.workflow_bundle.slurm_storage_tree_byte_count
            ),
            "slurm_storage_posterior_sample_byte_count": (
                result.workflow_bundle.slurm_storage_posterior_sample_byte_count
            ),
            "slurm_storage_report_byte_count": (
                result.workflow_bundle.slurm_storage_report_byte_count
            ),
            "slurm_storage_largest_variant_id": (
                result.workflow_bundle.slurm_storage_largest_variant_id
            ),
            "slurm_output_explosion_status": (
                result.workflow_bundle.slurm_output_explosion_status
            ),
            "slurm_output_explosion_global_issue_count": (
                result.workflow_bundle.slurm_output_explosion_global_issue_count
            ),
            "slurm_output_explosion_warning_variant_count": (
                result.workflow_bundle.slurm_output_explosion_warning_variant_count
            ),
            "slurm_output_explosion_high_risk_variant_count": (
                result.workflow_bundle.slurm_output_explosion_high_risk_variant_count
            ),
            "slurm_tree_retention_status": (
                result.workflow_bundle.slurm_tree_retention_status
            ),
            "slurm_tree_set_file_count": (
                result.workflow_bundle.slurm_tree_set_file_count
            ),
            "slurm_tree_posterior_sample_file_count": (
                result.workflow_bundle.slurm_tree_posterior_sample_file_count
            ),
            "slurm_tree_thinning_recommended_file_count": (
                result.workflow_bundle.slurm_tree_thinning_recommended_file_count
            ),
            "slurm_tree_thinning_required_file_count": (
                result.workflow_bundle.slurm_tree_thinning_required_file_count
            ),
            "slurm_tree_compression_recommended_file_count": (
                result.workflow_bundle.slurm_tree_compression_recommended_file_count
            ),
            "slurm_tree_compression_required_file_count": (
                result.workflow_bundle.slurm_tree_compression_required_file_count
            ),
            "slurm_merge_status": result.workflow_bundle.slurm_merge_status,
            "slurm_merge_ready": result.workflow_bundle.slurm_merge_ready,
            "slurm_mergeable_variant_count": (
                result.workflow_bundle.slurm_mergeable_variant_count
            ),
            "slurm_merge_failed_check_count": (
                result.workflow_bundle.slurm_merge_failed_check_count
            ),
            "slurm_output_freshness_check_count": (
                result.workflow_bundle.slurm_output_freshness_check_count
            ),
            "slurm_output_freshness_failed_check_count": (
                result.workflow_bundle.slurm_output_freshness_failed_check_count
            ),
            "slurm_fresh_output_job_count": (
                result.workflow_bundle.slurm_fresh_output_job_count
            ),
            "slurm_stale_output_job_count": (
                result.workflow_bundle.slurm_stale_output_job_count
            ),
            "slurm_completed_job_count": result.workflow_bundle.slurm_completed_job_count,
            "slurm_failed_job_count": result.workflow_bundle.slurm_failed_job_count,
            "slurm_pending_job_count": result.workflow_bundle.slurm_pending_job_count,
            "slurm_stale_job_count": result.workflow_bundle.slurm_stale_job_count,
            "slurm_failure_recovery_status": (
                result.workflow_bundle.slurm_failure_recovery_status
            ),
            "slurm_failure_recovery_rerunnable_job_count": (
                result.workflow_bundle.slurm_failure_recovery_rerunnable_job_count
            ),
            "slurm_failure_recovery_blocked_job_count": (
                result.workflow_bundle.slurm_failure_recovery_blocked_job_count
            ),
            "slurm_failure_recovery_partition_count": (
                result.workflow_bundle.slurm_failure_recovery_partition_count
            ),
            "reproducibility_passed": result.workflow_bundle.reproducibility_passed,
            "reproducibility_check_count": (
                result.workflow_bundle.reproducibility_check_count
            ),
            "reproducibility_failed_check_count": (
                result.workflow_bundle.reproducibility_failed_check_count
            ),
            "reproducibility_failed_variant_count": (
                result.workflow_bundle.reproducibility_failed_variant_count
            ),
            "reference_output_count": count_expected_output_files(
                result.dataset_export.expected_output_root
            ),
        },
        data=result,
        output_root=result.output_root,
    )
