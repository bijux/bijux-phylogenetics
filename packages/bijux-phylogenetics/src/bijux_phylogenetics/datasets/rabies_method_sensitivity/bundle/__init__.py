from __future__ import annotations

from pathlib import Path
import shutil

from ..audit import (
    audit_rabies_method_sensitivity_workflow_bundle,
    write_rabies_method_sensitivity_reproducibility_audit_json,
    write_rabies_method_sensitivity_reproducibility_checks_table,
    write_rabies_method_sensitivity_variant_audit_table,
)
from ..models import (
    RabiesMethodSensitivityPanelWorkflowBundle,
    RabiesMethodSensitivityPanelWorkflowReport,
)
from .package_ledger import (
    _write_clade_table,
    _write_conclusion_summary_table,
    _write_parallel_execution_summary_table,
    _write_preprocessing_comparison_table,
    _write_variant_summary_table,
    _write_workflow_summary_table,
)
from .package_manifest import (
    _sha256,
    _write_manifest,
    _write_resolved_config,
)
from .shared import (
    _format_float,
    _format_optional_bool,
    _format_optional_float,
    _write_tsv,
)
from .slurm_artifacts import _write_slurm_bundle_artifacts
from .variant_artifacts import (
    _copy_output,
    _copy_task_logs,
    _write_rooting_summary_table,
    _write_variant_outputs,
)
from ..reporting import _write_report, _write_report_manifest


def write_rabies_method_sensitivity_panel_workflow_bundle(
    output_root: Path,
    report: RabiesMethodSensitivityPanelWorkflowReport,
) -> RabiesMethodSensitivityPanelWorkflowBundle:
    """Write the governed reviewer-facing bundle for the method-sensitivity workflow."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    workflow_summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report,
    )
    variant_summary_path = _write_variant_summary_table(
        output_root / "variant-summary.tsv",
        report,
    )
    parallel_summary_path = _write_parallel_execution_summary_table(
        output_root / "parallel-execution-summary.tsv",
        report,
    )
    preprocessing_comparison_path = _write_preprocessing_comparison_table(
        output_root / "preprocessing-rooted-comparisons.tsv",
        report.preprocessing_comparison_rows,
    )
    stable_clades_path = _write_clade_table(
        output_root / "stable-clades.tsv",
        report.stable_clade_rows,
    )
    changed_clades_path = _write_clade_table(
        output_root / "changed-clades.tsv",
        report.changed_clade_rows,
    )
    conclusion_summary_path = _write_conclusion_summary_table(
        output_root / "method-conclusion-summary.tsv",
        report.conclusion_rows,
    )
    config_path = _write_resolved_config(
        output_root / "workflow-config.resolved.json",
        report,
    )
    execution_record_path = _copy_output(
        report.execution_record_path,
        output_root / report.execution_record_path.name,
    )
    task_logs_root = _copy_task_logs(output_root / "parallel-logs", report.task_records)
    variants_root = _write_variant_outputs(
        output_root / "variants", report.variant_runs
    )
    slurm_artifacts = _write_slurm_bundle_artifacts(
        output_root,
        report,
        execution_record_path=execution_record_path,
    )
    manifest_path = _write_manifest(
        output_root / "rabies-method-sensitivity.manifest.json",
        report=report,
        bundle_paths={
            "workflow_summary": workflow_summary_path,
            "variant_summary": variant_summary_path,
            "parallel_summary": parallel_summary_path,
            "preprocessing_comparison": preprocessing_comparison_path,
            "stable_clades": stable_clades_path,
            "changed_clades": changed_clades_path,
            "conclusion_summary": conclusion_summary_path,
            "config": config_path,
            "execution_record": execution_record_path,
            "task_logs_root": task_logs_root,
            "variants_root": variants_root,
            "slurm_job_plan": slurm_artifacts.job_plan_path,
            "slurm_assumptions": slurm_artifacts.assumptions_path,
            "slurm_summary": slurm_artifacts.summary_path,
            "slurm_array_partitions": slurm_artifacts.array_partitions_path,
            "slurm_array_members": slurm_artifacts.array_members_path,
            "slurm_array_strategy": slurm_artifacts.array_strategy_path,
            "slurm_array_scripts_root": slurm_artifacts.array_scripts_root,
            "slurm_job_evidence_root": slurm_artifacts.job_evidence_report.evidence_root,
            "slurm_job_evidence_index": slurm_artifacts.job_evidence_report.index_path,
            "slurm_job_evidence_summary": slurm_artifacts.job_evidence_report.summary_path,
            "slurm_merge_checks": slurm_artifacts.merge_checks_path,
            "slurm_merge_variants": slurm_artifacts.merge_variants_path,
            "slurm_merge_summary": slurm_artifacts.merge_summary_path,
            "slurm_merge_report": slurm_artifacts.merge_report_path,
            "slurm_output_freshness": slurm_artifacts.output_freshness_path,
            "slurm_output_freshness_checks": slurm_artifacts.output_freshness_checks_path,
            "slurm_output_freshness_summary": slurm_artifacts.output_freshness_summary_path,
            "slurm_job_status": slurm_artifacts.job_status_path,
            "slurm_partition_status": slurm_artifacts.partition_status_path,
            "slurm_workflow_status": slurm_artifacts.workflow_status_path,
            "slurm_failure_recovery_jobs": slurm_artifacts.failure_recovery_jobs_path,
            "slurm_failure_recovery_partitions": slurm_artifacts.failure_recovery_partitions_path,
            "slurm_failure_recovery_summary": slurm_artifacts.failure_recovery_summary_path,
            "slurm_failure_recovery_report": slurm_artifacts.failure_recovery_report_path,
        },
    )
    report_manifest_path = _write_report_manifest(
        output_root
        / "report-artifacts"
        / "rabies-method-sensitivity-report.manifest.json",
        report=report,
        bundle_paths={
            "workflow_summary": workflow_summary_path,
            "variant_summary": variant_summary_path,
            "parallel_summary": parallel_summary_path,
            "preprocessing_comparison": preprocessing_comparison_path,
            "stable_clades": stable_clades_path,
            "changed_clades": changed_clades_path,
            "conclusion_summary": conclusion_summary_path,
            "config": config_path,
            "execution_record": execution_record_path,
            "workflow_manifest": manifest_path,
            "slurm_job_plan": slurm_artifacts.job_plan_path,
            "slurm_assumptions": slurm_artifacts.assumptions_path,
            "slurm_summary": slurm_artifacts.summary_path,
            "slurm_array_partitions": slurm_artifacts.array_partitions_path,
            "slurm_array_members": slurm_artifacts.array_members_path,
            "slurm_array_strategy": slurm_artifacts.array_strategy_path,
            "slurm_job_evidence_index": slurm_artifacts.job_evidence_report.index_path,
            "slurm_job_evidence_summary": slurm_artifacts.job_evidence_report.summary_path,
            "slurm_merge_checks": slurm_artifacts.merge_checks_path,
            "slurm_merge_variants": slurm_artifacts.merge_variants_path,
            "slurm_merge_summary": slurm_artifacts.merge_summary_path,
            "slurm_merge_report": slurm_artifacts.merge_report_path,
            "slurm_output_freshness": slurm_artifacts.output_freshness_path,
            "slurm_output_freshness_checks": slurm_artifacts.output_freshness_checks_path,
            "slurm_output_freshness_summary": slurm_artifacts.output_freshness_summary_path,
            "slurm_job_status": slurm_artifacts.job_status_path,
            "slurm_partition_status": slurm_artifacts.partition_status_path,
            "slurm_workflow_status": slurm_artifacts.workflow_status_path,
            "slurm_failure_recovery_jobs": slurm_artifacts.failure_recovery_jobs_path,
            "slurm_failure_recovery_partitions": slurm_artifacts.failure_recovery_partitions_path,
            "slurm_failure_recovery_summary": slurm_artifacts.failure_recovery_summary_path,
            "slurm_failure_recovery_report": slurm_artifacts.failure_recovery_report_path,
        },
        sha256=_sha256,
    )
    reproducibility_report = audit_rabies_method_sensitivity_workflow_bundle(
        output_root,
        sequences_path=report.dataset.sequences_path,
        metadata_path=report.dataset.metadata_path,
    )
    reproducibility_checks_path = (
        write_rabies_method_sensitivity_reproducibility_checks_table(
            output_root / "reproducibility-checks.tsv",
            reproducibility_report,
        )
    )
    reproducibility_variant_audit_path = (
        write_rabies_method_sensitivity_variant_audit_table(
            output_root / "reproducibility-variants.tsv",
            reproducibility_report,
        )
    )
    reproducibility_audit_path = (
        write_rabies_method_sensitivity_reproducibility_audit_json(
            output_root / "reproducibility-audit.json",
            reproducibility_report,
        )
    )
    report_linked_files = (
        workflow_summary_path,
        variant_summary_path,
        parallel_summary_path,
        preprocessing_comparison_path,
        stable_clades_path,
        changed_clades_path,
        conclusion_summary_path,
        config_path,
        execution_record_path,
        manifest_path,
        report_manifest_path,
        slurm_artifacts.job_plan_path,
        slurm_artifacts.assumptions_path,
        slurm_artifacts.summary_path,
        slurm_artifacts.array_partitions_path,
        slurm_artifacts.array_members_path,
        slurm_artifacts.array_strategy_path,
        slurm_artifacts.job_evidence_report.index_path,
        slurm_artifacts.job_evidence_report.summary_path,
        slurm_artifacts.storage_categories_path,
        slurm_artifacts.storage_variants_path,
        slurm_artifacts.storage_summary_path,
        slurm_artifacts.storage_report_path,
        slurm_artifacts.output_explosion_checks_path,
        slurm_artifacts.output_explosion_variants_path,
        slurm_artifacts.output_explosion_summary_path,
        slurm_artifacts.output_explosion_report_path,
        slurm_artifacts.tree_retention_checks_path,
        slurm_artifacts.tree_retention_files_path,
        slurm_artifacts.tree_retention_summary_path,
        slurm_artifacts.tree_retention_report_path,
        slurm_artifacts.merge_checks_path,
        slurm_artifacts.merge_variants_path,
        slurm_artifacts.merge_summary_path,
        slurm_artifacts.merge_report_path,
        slurm_artifacts.output_freshness_path,
        slurm_artifacts.output_freshness_checks_path,
        slurm_artifacts.output_freshness_summary_path,
        slurm_artifacts.job_status_path,
        slurm_artifacts.partition_status_path,
        slurm_artifacts.workflow_status_path,
        slurm_artifacts.failure_recovery_jobs_path,
        slurm_artifacts.failure_recovery_partitions_path,
        slurm_artifacts.failure_recovery_summary_path,
        slurm_artifacts.failure_recovery_report_path,
    )
    report_linked_artifact_count = len(report_linked_files)
    report_path = _write_report(
        output_root / "rabies-method-sensitivity-report.html",
        report=report,
        bundle_paths={
            "workflow_summary": workflow_summary_path,
            "variant_summary": variant_summary_path,
            "parallel_summary": parallel_summary_path,
            "preprocessing_comparison": preprocessing_comparison_path,
            "stable_clades": stable_clades_path,
            "changed_clades": changed_clades_path,
            "conclusion_summary": conclusion_summary_path,
            "config": config_path,
            "execution_record": execution_record_path,
            "workflow_manifest": manifest_path,
            "slurm_job_plan": slurm_artifacts.job_plan_path,
            "slurm_assumptions": slurm_artifacts.assumptions_path,
            "slurm_summary": slurm_artifacts.summary_path,
            "slurm_array_partitions": slurm_artifacts.array_partitions_path,
            "slurm_array_members": slurm_artifacts.array_members_path,
            "slurm_array_strategy": slurm_artifacts.array_strategy_path,
            "slurm_job_evidence_index": slurm_artifacts.job_evidence_report.index_path,
            "slurm_job_evidence_summary": slurm_artifacts.job_evidence_report.summary_path,
            "slurm_storage_categories": slurm_artifacts.storage_categories_path,
            "slurm_storage_variants": slurm_artifacts.storage_variants_path,
            "slurm_storage_summary": slurm_artifacts.storage_summary_path,
            "slurm_storage_report": slurm_artifacts.storage_report_path,
            "slurm_output_explosion_checks": slurm_artifacts.output_explosion_checks_path,
            "slurm_output_explosion_variants": slurm_artifacts.output_explosion_variants_path,
            "slurm_output_explosion_summary": slurm_artifacts.output_explosion_summary_path,
            "slurm_output_explosion_report": slurm_artifacts.output_explosion_report_path,
            "slurm_tree_retention_checks": slurm_artifacts.tree_retention_checks_path,
            "slurm_tree_retention_files": slurm_artifacts.tree_retention_files_path,
            "slurm_tree_retention_summary": slurm_artifacts.tree_retention_summary_path,
            "slurm_tree_retention_report": slurm_artifacts.tree_retention_report_path,
            "slurm_merge_checks": slurm_artifacts.merge_checks_path,
            "slurm_merge_variants": slurm_artifacts.merge_variants_path,
            "slurm_merge_summary": slurm_artifacts.merge_summary_path,
            "slurm_merge_report": slurm_artifacts.merge_report_path,
            "slurm_output_freshness": slurm_artifacts.output_freshness_path,
            "slurm_output_freshness_checks": slurm_artifacts.output_freshness_checks_path,
            "slurm_output_freshness_summary": slurm_artifacts.output_freshness_summary_path,
            "slurm_job_status": slurm_artifacts.job_status_path,
            "slurm_partition_status": slurm_artifacts.partition_status_path,
            "slurm_workflow_status": slurm_artifacts.workflow_status_path,
            "slurm_failure_recovery_jobs": slurm_artifacts.failure_recovery_jobs_path,
            "slurm_failure_recovery_partitions": slurm_artifacts.failure_recovery_partitions_path,
            "slurm_failure_recovery_summary": slurm_artifacts.failure_recovery_summary_path,
            "slurm_failure_recovery_report": slurm_artifacts.failure_recovery_report_path,
            "reproducibility_checks": reproducibility_checks_path,
            "reproducibility_variant_audit": reproducibility_variant_audit_path,
            "reproducibility_audit": reproducibility_audit_path,
        },
        report_manifest_path=report_manifest_path,
        reproducibility_report=reproducibility_report,
        slurm_planning_report=slurm_artifacts.planning_report,
        slurm_array_strategy_report=slurm_artifacts.array_strategy_report,
        slurm_job_evidence_report=slurm_artifacts.job_evidence_report,
        slurm_storage_report=slurm_artifacts.storage_report,
        slurm_output_explosion_report=slurm_artifacts.output_explosion_report,
        slurm_tree_retention_report=slurm_artifacts.tree_retention_report,
        slurm_merge_report=slurm_artifacts.merge_report,
        slurm_output_freshness_report=slurm_artifacts.output_freshness_report,
        slurm_status_report=slurm_artifacts.status_report,
        slurm_failure_recovery_report=slurm_artifacts.failure_recovery_report,
    )
    report_html_size_bytes = report_path.stat().st_size
    report_linked_artifact_bytes = sum(
        path.stat().st_size for path in report_linked_files
    )
    report_total_output_bytes = report_html_size_bytes + report_linked_artifact_bytes
    preprocessing_change_pair_count = sum(
        1
        for row in report.preprocessing_comparison_rows
        if row.robinson_foulds_distance > 0 or row.same_taxa_different_rooting
    )
    rooted_engine_change_variant_count = sum(
        1
        for variant in report.variant_runs
        if variant.rooted_engine_comparison.robinson_foulds_distance > 0
        or variant.rooted_engine_comparison.same_taxa_different_rooting
    )
    serious_conflict_variant_count = sum(
        1
        for variant in report.variant_runs
        if variant.inference_comparison.conclusion_summary.serious_conflict_count > 0
    )
    return RabiesMethodSensitivityPanelWorkflowBundle(
        output_root=output_root,
        variant_count=len(report.variant_runs),
        stable_clade_count=len(report.stable_clade_rows),
        changed_clade_count=len(report.changed_clade_rows),
        preprocessing_change_pair_count=preprocessing_change_pair_count,
        rooted_engine_change_variant_count=rooted_engine_change_variant_count,
        serious_conflict_variant_count=serious_conflict_variant_count,
        execution_record_path=execution_record_path,
        parallel_workers=report.parallel_workers,
        execution_mode=report.execution_mode,
        workflow_summary_path=workflow_summary_path,
        variant_summary_path=variant_summary_path,
        parallel_summary_path=parallel_summary_path,
        preprocessing_comparison_path=preprocessing_comparison_path,
        stable_clades_path=stable_clades_path,
        changed_clades_path=changed_clades_path,
        conclusion_summary_path=conclusion_summary_path,
        config_path=config_path,
        manifest_path=manifest_path,
        report_manifest_path=report_manifest_path,
        slurm_job_plan_path=slurm_artifacts.job_plan_path,
        slurm_assumptions_path=slurm_artifacts.assumptions_path,
        slurm_summary_path=slurm_artifacts.summary_path,
        slurm_array_partitions_path=slurm_artifacts.array_partitions_path,
        slurm_array_members_path=slurm_artifacts.array_members_path,
        slurm_array_strategy_path=slurm_artifacts.array_strategy_path,
        slurm_array_scripts_root=slurm_artifacts.array_scripts_root,
        slurm_job_evidence_root=slurm_artifacts.job_evidence_report.evidence_root,
        slurm_job_evidence_index_path=slurm_artifacts.job_evidence_report.index_path,
        slurm_job_evidence_summary_path=slurm_artifacts.job_evidence_report.summary_path,
        slurm_storage_categories_path=slurm_artifacts.storage_categories_path,
        slurm_storage_variants_path=slurm_artifacts.storage_variants_path,
        slurm_storage_summary_path=slurm_artifacts.storage_summary_path,
        slurm_storage_report_path=slurm_artifacts.storage_report_path,
        slurm_output_explosion_checks_path=slurm_artifacts.output_explosion_checks_path,
        slurm_output_explosion_variants_path=slurm_artifacts.output_explosion_variants_path,
        slurm_output_explosion_summary_path=slurm_artifacts.output_explosion_summary_path,
        slurm_output_explosion_report_path=slurm_artifacts.output_explosion_report_path,
        slurm_tree_retention_checks_path=slurm_artifacts.tree_retention_checks_path,
        slurm_tree_retention_files_path=slurm_artifacts.tree_retention_files_path,
        slurm_tree_retention_summary_path=slurm_artifacts.tree_retention_summary_path,
        slurm_tree_retention_report_path=slurm_artifacts.tree_retention_report_path,
        slurm_merge_checks_path=slurm_artifacts.merge_checks_path,
        slurm_merge_variants_path=slurm_artifacts.merge_variants_path,
        slurm_merge_summary_path=slurm_artifacts.merge_summary_path,
        slurm_merge_report_path=slurm_artifacts.merge_report_path,
        slurm_job_count=slurm_artifacts.planning_report.job_count,
        slurm_total_estimated_core_hours=(
            slurm_artifacts.planning_report.total_estimated_core_hours
        ),
        slurm_maximum_estimated_memory_mib=(
            slurm_artifacts.planning_report.maximum_estimated_memory_mib
        ),
        slurm_maximum_estimated_wallclock_minutes=(
            slurm_artifacts.planning_report.maximum_estimated_wallclock_minutes
        ),
        slurm_total_estimated_scratch_mib=(
            slurm_artifacts.planning_report.total_estimated_scratch_mib
        ),
        slurm_total_estimated_output_mib=(
            slurm_artifacts.planning_report.total_estimated_output_mib
        ),
        slurm_array_partition_count=slurm_artifacts.array_strategy_report.partition_count,
        slurm_array_script_count=slurm_artifacts.array_strategy_report.script_count,
        slurm_array_largest_partition_size=(
            slurm_artifacts.array_strategy_report.largest_partition_size
        ),
        slurm_job_evidence_file_count=(
            slurm_artifacts.job_evidence_report.total_artifact_file_count
        ),
        slurm_job_evidence_total_runtime_seconds=(
            slurm_artifacts.job_evidence_report.total_runtime_seconds
        ),
        slurm_job_evidence_total_output_byte_count=(
            slurm_artifacts.job_evidence_report.total_output_byte_count
        ),
        slurm_storage_total_estimated_mib=(
            slurm_artifacts.storage_report.total_estimated_storage_mib
        ),
        slurm_storage_output_byte_count=slurm_artifacts.storage_report.output_byte_count,
        slurm_storage_log_byte_count=slurm_artifacts.storage_report.log_byte_count,
        slurm_storage_tree_byte_count=slurm_artifacts.storage_report.tree_byte_count,
        slurm_storage_posterior_sample_byte_count=(
            slurm_artifacts.storage_report.posterior_sample_byte_count
        ),
        slurm_storage_report_byte_count=slurm_artifacts.storage_report.report_byte_count,
        slurm_storage_largest_variant_id=slurm_artifacts.storage_report.largest_variant_id,
        slurm_output_explosion_status=(
            slurm_artifacts.output_explosion_report.overall_risk_status
        ),
        slurm_output_explosion_global_issue_count=(
            slurm_artifacts.output_explosion_report.global_issue_count
        ),
        slurm_output_explosion_warning_variant_count=(
            slurm_artifacts.output_explosion_report.warning_variant_count
        ),
        slurm_output_explosion_high_risk_variant_count=(
            slurm_artifacts.output_explosion_report.high_risk_variant_count
        ),
        slurm_tree_retention_status=(
            slurm_artifacts.tree_retention_report.overall_policy_status
        ),
        slurm_tree_set_file_count=slurm_artifacts.tree_retention_report.tree_set_file_count,
        slurm_tree_posterior_sample_file_count=(
            slurm_artifacts.tree_retention_report.posterior_sample_file_count
        ),
        slurm_tree_thinning_recommended_file_count=(
            slurm_artifacts.tree_retention_report.thinning_recommended_file_count
        ),
        slurm_tree_thinning_required_file_count=(
            slurm_artifacts.tree_retention_report.thinning_required_file_count
        ),
        slurm_tree_compression_recommended_file_count=(
            slurm_artifacts.tree_retention_report.compression_recommended_file_count
        ),
        slurm_tree_compression_required_file_count=(
            slurm_artifacts.tree_retention_report.compression_required_file_count
        ),
        slurm_merge_status=slurm_artifacts.merge_report.merge_status,
        slurm_merge_ready=slurm_artifacts.merge_report.merge_ready,
        slurm_mergeable_variant_count=slurm_artifacts.merge_report.mergeable_variant_count,
        slurm_merge_failed_check_count=slurm_artifacts.merge_report.failed_check_count,
        slurm_output_freshness_path=slurm_artifacts.output_freshness_path,
        slurm_output_freshness_checks_path=slurm_artifacts.output_freshness_checks_path,
        slurm_output_freshness_summary_path=slurm_artifacts.output_freshness_summary_path,
        slurm_job_status_path=slurm_artifacts.job_status_path,
        slurm_partition_status_path=slurm_artifacts.partition_status_path,
        slurm_workflow_status_path=slurm_artifacts.workflow_status_path,
        slurm_failure_recovery_jobs_path=slurm_artifacts.failure_recovery_jobs_path,
        slurm_failure_recovery_partitions_path=slurm_artifacts.failure_recovery_partitions_path,
        slurm_failure_recovery_summary_path=slurm_artifacts.failure_recovery_summary_path,
        slurm_failure_recovery_report_path=slurm_artifacts.failure_recovery_report_path,
        slurm_output_freshness_check_count=slurm_artifacts.output_freshness_report.check_count,
        slurm_output_freshness_failed_check_count=(
            slurm_artifacts.output_freshness_report.failed_check_count
        ),
        slurm_fresh_output_job_count=slurm_artifacts.output_freshness_report.fresh_job_count,
        slurm_stale_output_job_count=slurm_artifacts.output_freshness_report.stale_job_count,
        slurm_completed_job_count=slurm_artifacts.status_report.completed_job_count,
        slurm_failed_job_count=slurm_artifacts.status_report.failed_job_count,
        slurm_pending_job_count=slurm_artifacts.status_report.pending_job_count,
        slurm_stale_job_count=slurm_artifacts.status_report.stale_job_count,
        slurm_failure_recovery_status=(
            slurm_artifacts.failure_recovery_report.overall_recovery_status
        ),
        slurm_failure_recovery_rerunnable_job_count=(
            slurm_artifacts.failure_recovery_report.rerunnable_job_count
        ),
        slurm_failure_recovery_blocked_job_count=(
            slurm_artifacts.failure_recovery_report.blocked_job_count
        ),
        slurm_failure_recovery_partition_count=(
            slurm_artifacts.failure_recovery_report.recovery_partition_count
        ),
        reproducibility_checks_path=reproducibility_checks_path,
        reproducibility_variant_audit_path=reproducibility_variant_audit_path,
        reproducibility_audit_path=reproducibility_audit_path,
        reproducibility_passed=reproducibility_report.all_passed,
        reproducibility_check_count=reproducibility_report.check_count,
        reproducibility_failed_check_count=reproducibility_report.failed_check_count,
        reproducibility_failed_variant_count=reproducibility_report.failed_variant_count,
        report_path=report_path,
        report_linked_artifact_count=report_linked_artifact_count,
        report_html_size_bytes=report_html_size_bytes,
        report_linked_artifact_bytes=report_linked_artifact_bytes,
        report_total_output_bytes=report_total_output_bytes,
        task_logs_root=task_logs_root,
        variants_root=variants_root,
    )

