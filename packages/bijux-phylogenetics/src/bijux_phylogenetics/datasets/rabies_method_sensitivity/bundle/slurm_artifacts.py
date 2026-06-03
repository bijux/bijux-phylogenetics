from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..models import RabiesMethodSensitivityPanelWorkflowReport
from ..slurm import (
    build_rabies_method_sensitivity_slurm_array_strategy_report,
    build_rabies_method_sensitivity_slurm_failure_recovery_report,
    build_rabies_method_sensitivity_slurm_merge_report,
    build_rabies_method_sensitivity_slurm_output_explosion_report,
    build_rabies_method_sensitivity_slurm_output_freshness_report,
    build_rabies_method_sensitivity_slurm_planning_report,
    build_rabies_method_sensitivity_slurm_status_report,
    build_rabies_method_sensitivity_slurm_storage_report,
    build_rabies_method_sensitivity_slurm_tree_retention_report,
    write_rabies_method_sensitivity_slurm_array_members_table,
    write_rabies_method_sensitivity_slurm_array_partition_scripts,
    write_rabies_method_sensitivity_slurm_array_partitions_table,
    write_rabies_method_sensitivity_slurm_array_strategy_json,
    write_rabies_method_sensitivity_slurm_assumptions_table,
    write_rabies_method_sensitivity_slurm_failure_recovery_html_report,
    write_rabies_method_sensitivity_slurm_failure_recovery_jobs_table,
    write_rabies_method_sensitivity_slurm_failure_recovery_partitions_table,
    write_rabies_method_sensitivity_slurm_failure_recovery_summary_json,
    write_rabies_method_sensitivity_slurm_job_evidence_bundle,
    write_rabies_method_sensitivity_slurm_job_plan_table,
    write_rabies_method_sensitivity_slurm_job_status_table,
    write_rabies_method_sensitivity_slurm_merge_checks_table,
    write_rabies_method_sensitivity_slurm_merge_html_report,
    write_rabies_method_sensitivity_slurm_merge_summary_json,
    write_rabies_method_sensitivity_slurm_merge_variants_table,
    write_rabies_method_sensitivity_slurm_output_explosion_checks_table,
    write_rabies_method_sensitivity_slurm_output_explosion_html_report,
    write_rabies_method_sensitivity_slurm_output_explosion_summary_json,
    write_rabies_method_sensitivity_slurm_output_explosion_variants_table,
    write_rabies_method_sensitivity_slurm_output_freshness_checks_table,
    write_rabies_method_sensitivity_slurm_output_freshness_json,
    write_rabies_method_sensitivity_slurm_output_freshness_table,
    write_rabies_method_sensitivity_slurm_partition_status_table,
    write_rabies_method_sensitivity_slurm_status_json,
    write_rabies_method_sensitivity_slurm_storage_categories_table,
    write_rabies_method_sensitivity_slurm_storage_html_report,
    write_rabies_method_sensitivity_slurm_storage_summary_json,
    write_rabies_method_sensitivity_slurm_storage_variants_table,
    write_rabies_method_sensitivity_slurm_summary_json,
    write_rabies_method_sensitivity_slurm_tree_retention_checks_table,
    write_rabies_method_sensitivity_slurm_tree_retention_files_table,
    write_rabies_method_sensitivity_slurm_tree_retention_html_report,
    write_rabies_method_sensitivity_slurm_tree_retention_summary_json,
)


@dataclass(slots=True)
class RabiesMethodSensitivitySlurmBundleArtifacts:
    """Owned SLURM review artifacts written for the workflow bundle."""

    planning_report: object
    array_strategy_report: object
    output_freshness_report: object
    status_report: object
    failure_recovery_report: object
    job_evidence_report: object
    merge_report: object
    storage_report: object
    output_explosion_report: object
    tree_retention_report: object
    job_plan_path: Path
    assumptions_path: Path
    summary_path: Path
    array_partitions_path: Path
    array_members_path: Path
    array_strategy_path: Path
    array_scripts_root: Path
    output_freshness_path: Path
    output_freshness_checks_path: Path
    output_freshness_summary_path: Path
    job_status_path: Path
    partition_status_path: Path
    workflow_status_path: Path
    failure_recovery_jobs_path: Path
    failure_recovery_partitions_path: Path
    failure_recovery_summary_path: Path
    failure_recovery_report_path: Path
    merge_checks_path: Path
    merge_variants_path: Path
    merge_summary_path: Path
    merge_report_path: Path
    storage_categories_path: Path
    storage_variants_path: Path
    storage_summary_path: Path
    storage_report_path: Path
    output_explosion_checks_path: Path
    output_explosion_variants_path: Path
    output_explosion_summary_path: Path
    output_explosion_report_path: Path
    tree_retention_checks_path: Path
    tree_retention_files_path: Path
    tree_retention_summary_path: Path
    tree_retention_report_path: Path


def _write_slurm_bundle_artifacts(
    output_root: Path,
    report: RabiesMethodSensitivityPanelWorkflowReport,
    *,
    execution_record_path: Path,
) -> RabiesMethodSensitivitySlurmBundleArtifacts:
    planning_report = build_rabies_method_sensitivity_slurm_planning_report(report)
    job_plan_path = write_rabies_method_sensitivity_slurm_job_plan_table(
        output_root / "slurm-job-plan.tsv",
        planning_report,
    )
    assumptions_path = write_rabies_method_sensitivity_slurm_assumptions_table(
        output_root / "slurm-estimation-assumptions.tsv",
        planning_report,
    )
    summary_path = write_rabies_method_sensitivity_slurm_summary_json(
        output_root / "slurm-planning-summary.json",
        planning_report,
    )
    array_strategy_report = build_rabies_method_sensitivity_slurm_array_strategy_report(
        planning_report
    )
    array_partitions_path = (
        write_rabies_method_sensitivity_slurm_array_partitions_table(
            output_root / "slurm-array-partitions.tsv",
            array_strategy_report,
        )
    )
    array_members_path = write_rabies_method_sensitivity_slurm_array_members_table(
        output_root / "slurm-array-members.tsv",
        array_strategy_report,
    )
    array_strategy_path = write_rabies_method_sensitivity_slurm_array_strategy_json(
        output_root / "slurm-array-strategy.json",
        array_strategy_report,
    )
    array_scripts_root = write_rabies_method_sensitivity_slurm_array_partition_scripts(
        output_root / "slurm-arrays",
        array_strategy_report,
    )
    output_freshness_report = (
        build_rabies_method_sensitivity_slurm_output_freshness_report(
            output_root,
            dataset=report.dataset,
        )
    )
    output_freshness_path = (
        write_rabies_method_sensitivity_slurm_output_freshness_table(
            output_root / "slurm-output-freshness.tsv",
            output_freshness_report,
        )
    )
    output_freshness_checks_path = (
        write_rabies_method_sensitivity_slurm_output_freshness_checks_table(
            output_root / "slurm-output-freshness-checks.tsv",
            output_freshness_report,
        )
    )
    output_freshness_summary_path = (
        write_rabies_method_sensitivity_slurm_output_freshness_json(
            output_root / "slurm-output-freshness.json",
            output_freshness_report,
        )
    )
    status_report = build_rabies_method_sensitivity_slurm_status_report(
        output_root,
        dataset=report.dataset,
    )
    job_status_path = write_rabies_method_sensitivity_slurm_job_status_table(
        output_root / "slurm-job-status.tsv",
        status_report,
    )
    partition_status_path = (
        write_rabies_method_sensitivity_slurm_partition_status_table(
            output_root / "slurm-partition-status.tsv",
            status_report,
        )
    )
    workflow_status_path = write_rabies_method_sensitivity_slurm_status_json(
        output_root / "slurm-workflow-status.json",
        status_report,
    )
    failure_recovery_report = (
        build_rabies_method_sensitivity_slurm_failure_recovery_report(output_root)
    )
    failure_recovery_jobs_path = (
        write_rabies_method_sensitivity_slurm_failure_recovery_jobs_table(
            output_root / "slurm-failure-recovery-jobs.tsv",
            failure_recovery_report,
        )
    )
    failure_recovery_partitions_path = (
        write_rabies_method_sensitivity_slurm_failure_recovery_partitions_table(
            output_root / "slurm-failure-recovery-partitions.tsv",
            failure_recovery_report,
        )
    )
    failure_recovery_summary_path = (
        write_rabies_method_sensitivity_slurm_failure_recovery_summary_json(
            output_root / "slurm-failure-recovery-report.json",
            failure_recovery_report,
        )
    )
    failure_recovery_report_path = (
        write_rabies_method_sensitivity_slurm_failure_recovery_html_report(
            output_root / "slurm-failure-recovery-report.html",
            failure_recovery_report,
        )
    )
    job_evidence_report = write_rabies_method_sensitivity_slurm_job_evidence_bundle(
        output_root / "slurm-job-evidence",
        bundle_root=output_root,
        dataset_id=report.dataset.dataset_id,
        workflow_prefix=report.dataset.workflow_prefix,
        execution_mode=report.execution_mode,
        parallel_workers=report.parallel_workers,
        task_records=report.task_records,
        variant_runs=report.variant_runs,
        array_strategy_report=array_strategy_report,
        execution_record_path=execution_record_path,
        workflow_manifest_path=output_root / "rabies-method-sensitivity.manifest.json",
    )
    merge_report = build_rabies_method_sensitivity_slurm_merge_report(output_root)
    merge_checks_path = write_rabies_method_sensitivity_slurm_merge_checks_table(
        output_root / "slurm-merge-checks.tsv",
        merge_report,
    )
    merge_variants_path = write_rabies_method_sensitivity_slurm_merge_variants_table(
        output_root / "slurm-merge-variants.tsv",
        merge_report,
    )
    merge_summary_path = write_rabies_method_sensitivity_slurm_merge_summary_json(
        output_root / "slurm-merge-report.json",
        merge_report,
    )
    merge_report_path = write_rabies_method_sensitivity_slurm_merge_html_report(
        output_root / "slurm-merge-report.html",
        merge_report,
    )
    storage_report = build_rabies_method_sensitivity_slurm_storage_report(output_root)
    storage_categories_path = (
        write_rabies_method_sensitivity_slurm_storage_categories_table(
            output_root / "slurm-storage-categories.tsv",
            storage_report,
        )
    )
    storage_variants_path = (
        write_rabies_method_sensitivity_slurm_storage_variants_table(
            output_root / "slurm-storage-variants.tsv",
            storage_report,
        )
    )
    storage_summary_path = write_rabies_method_sensitivity_slurm_storage_summary_json(
        output_root / "slurm-storage-report.json",
        storage_report,
    )
    storage_report_path = write_rabies_method_sensitivity_slurm_storage_html_report(
        output_root / "slurm-storage-report.html",
        storage_report,
    )
    output_explosion_report = (
        build_rabies_method_sensitivity_slurm_output_explosion_report(output_root)
    )
    output_explosion_checks_path = (
        write_rabies_method_sensitivity_slurm_output_explosion_checks_table(
            output_root / "slurm-output-explosion-checks.tsv",
            output_explosion_report,
        )
    )
    output_explosion_variants_path = (
        write_rabies_method_sensitivity_slurm_output_explosion_variants_table(
            output_root / "slurm-output-explosion-variants.tsv",
            output_explosion_report,
        )
    )
    output_explosion_summary_path = (
        write_rabies_method_sensitivity_slurm_output_explosion_summary_json(
            output_root / "slurm-output-explosion-report.json",
            output_explosion_report,
        )
    )
    output_explosion_report_path = (
        write_rabies_method_sensitivity_slurm_output_explosion_html_report(
            output_root / "slurm-output-explosion-report.html",
            output_explosion_report,
        )
    )
    tree_retention_report = build_rabies_method_sensitivity_slurm_tree_retention_report(
        output_root
    )
    tree_retention_checks_path = (
        write_rabies_method_sensitivity_slurm_tree_retention_checks_table(
            output_root / "slurm-tree-retention-checks.tsv",
            tree_retention_report,
        )
    )
    tree_retention_files_path = (
        write_rabies_method_sensitivity_slurm_tree_retention_files_table(
            output_root / "slurm-tree-retention-files.tsv",
            tree_retention_report,
        )
    )
    tree_retention_summary_path = (
        write_rabies_method_sensitivity_slurm_tree_retention_summary_json(
            output_root / "slurm-tree-retention-policy.json",
            tree_retention_report,
        )
    )
    tree_retention_report_path = (
        write_rabies_method_sensitivity_slurm_tree_retention_html_report(
            output_root / "slurm-tree-retention-policy.html",
            tree_retention_report,
        )
    )
    return RabiesMethodSensitivitySlurmBundleArtifacts(
        planning_report=planning_report,
        array_strategy_report=array_strategy_report,
        output_freshness_report=output_freshness_report,
        status_report=status_report,
        failure_recovery_report=failure_recovery_report,
        job_evidence_report=job_evidence_report,
        merge_report=merge_report,
        storage_report=storage_report,
        output_explosion_report=output_explosion_report,
        tree_retention_report=tree_retention_report,
        job_plan_path=job_plan_path,
        assumptions_path=assumptions_path,
        summary_path=summary_path,
        array_partitions_path=array_partitions_path,
        array_members_path=array_members_path,
        array_strategy_path=array_strategy_path,
        array_scripts_root=array_scripts_root,
        output_freshness_path=output_freshness_path,
        output_freshness_checks_path=output_freshness_checks_path,
        output_freshness_summary_path=output_freshness_summary_path,
        job_status_path=job_status_path,
        partition_status_path=partition_status_path,
        workflow_status_path=workflow_status_path,
        failure_recovery_jobs_path=failure_recovery_jobs_path,
        failure_recovery_partitions_path=failure_recovery_partitions_path,
        failure_recovery_summary_path=failure_recovery_summary_path,
        failure_recovery_report_path=failure_recovery_report_path,
        merge_checks_path=merge_checks_path,
        merge_variants_path=merge_variants_path,
        merge_summary_path=merge_summary_path,
        merge_report_path=merge_report_path,
        storage_categories_path=storage_categories_path,
        storage_variants_path=storage_variants_path,
        storage_summary_path=storage_summary_path,
        storage_report_path=storage_report_path,
        output_explosion_checks_path=output_explosion_checks_path,
        output_explosion_variants_path=output_explosion_variants_path,
        output_explosion_summary_path=output_explosion_summary_path,
        output_explosion_report_path=output_explosion_report_path,
        tree_retention_checks_path=tree_retention_checks_path,
        tree_retention_files_path=tree_retention_files_path,
        tree_retention_summary_path=tree_retention_summary_path,
        tree_retention_report_path=tree_retention_report_path,
    )
