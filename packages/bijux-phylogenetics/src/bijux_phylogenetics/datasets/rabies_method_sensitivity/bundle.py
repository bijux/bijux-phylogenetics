from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil

from .audit import (
    RabiesMethodSensitivityReproducibilityAuditReport,
    audit_rabies_method_sensitivity_workflow_bundle,
    write_rabies_method_sensitivity_reproducibility_audit_json,
    write_rabies_method_sensitivity_reproducibility_checks_table,
    write_rabies_method_sensitivity_variant_audit_table,
)
from .slurm import (
    RabiesMethodSensitivitySlurmPlanningReport,
    build_rabies_method_sensitivity_slurm_planning_report,
    write_rabies_method_sensitivity_slurm_assumptions_table,
    write_rabies_method_sensitivity_slurm_job_plan_table,
    write_rabies_method_sensitivity_slurm_summary_json,
)
from .slurm import (
    RabiesMethodSensitivitySlurmArrayStrategyReport,
    build_rabies_method_sensitivity_slurm_array_strategy_report,
    write_rabies_method_sensitivity_slurm_array_members_table,
    write_rabies_method_sensitivity_slurm_array_partition_scripts,
    write_rabies_method_sensitivity_slurm_array_partitions_table,
    write_rabies_method_sensitivity_slurm_array_strategy_json,
)
from .slurm import (
    RabiesMethodSensitivitySlurmOutputFreshnessReport,
    build_rabies_method_sensitivity_slurm_output_freshness_report,
    write_rabies_method_sensitivity_slurm_output_freshness_checks_table,
    write_rabies_method_sensitivity_slurm_output_freshness_json,
    write_rabies_method_sensitivity_slurm_output_freshness_table,
)
from .slurm import (
    RabiesMethodSensitivitySlurmFailureRecoveryReport,
    build_rabies_method_sensitivity_slurm_failure_recovery_report,
    write_rabies_method_sensitivity_slurm_failure_recovery_html_report,
    write_rabies_method_sensitivity_slurm_failure_recovery_jobs_table,
    write_rabies_method_sensitivity_slurm_failure_recovery_partitions_table,
    write_rabies_method_sensitivity_slurm_failure_recovery_summary_json,
)
from .slurm import (
    RabiesMethodSensitivitySlurmJobEvidenceReport,
    write_rabies_method_sensitivity_slurm_job_evidence_bundle,
)
from .slurm import (
    RabiesMethodSensitivitySlurmMergeReport,
    build_rabies_method_sensitivity_slurm_merge_report,
    write_rabies_method_sensitivity_slurm_merge_checks_table,
    write_rabies_method_sensitivity_slurm_merge_html_report,
    write_rabies_method_sensitivity_slurm_merge_summary_json,
    write_rabies_method_sensitivity_slurm_merge_variants_table,
)
from .slurm import (
    RabiesMethodSensitivitySlurmStatusReport,
    build_rabies_method_sensitivity_slurm_status_report,
    write_rabies_method_sensitivity_slurm_job_status_table,
    write_rabies_method_sensitivity_slurm_partition_status_table,
    write_rabies_method_sensitivity_slurm_status_json,
)
from .slurm import (
    RabiesMethodSensitivitySlurmOutputExplosionReport,
    build_rabies_method_sensitivity_slurm_output_explosion_report,
    write_rabies_method_sensitivity_slurm_output_explosion_checks_table,
    write_rabies_method_sensitivity_slurm_output_explosion_html_report,
    write_rabies_method_sensitivity_slurm_output_explosion_summary_json,
    write_rabies_method_sensitivity_slurm_output_explosion_variants_table,
)
from .slurm import (
    RabiesMethodSensitivitySlurmTreeRetentionReport,
    build_rabies_method_sensitivity_slurm_tree_retention_report,
    write_rabies_method_sensitivity_slurm_tree_retention_checks_table,
    write_rabies_method_sensitivity_slurm_tree_retention_files_table,
    write_rabies_method_sensitivity_slurm_tree_retention_html_report,
    write_rabies_method_sensitivity_slurm_tree_retention_summary_json,
)
from .slurm import (
    RabiesMethodSensitivitySlurmStorageReport,
    build_rabies_method_sensitivity_slurm_storage_report,
    write_rabies_method_sensitivity_slurm_storage_categories_table,
    write_rabies_method_sensitivity_slurm_storage_html_report,
    write_rabies_method_sensitivity_slurm_storage_summary_json,
    write_rabies_method_sensitivity_slurm_storage_variants_table,
)
from bijux_phylogenetics.render.html import write_html_report

from .models import (
    RabiesMethodSensitivityCladeRow,
    RabiesMethodSensitivityConclusionRow,
    RabiesMethodSensitivityPanelWorkflowBundle,
    RabiesMethodSensitivityPanelWorkflowReport,
    RabiesMethodSensitivityPreprocessingComparisonRow,
    RabiesMethodSensitivityTaskRecord,
    RabiesMethodSensitivityVariant,
    RabiesMethodSensitivityVariantRun,
)


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
    slurm_planning_report = build_rabies_method_sensitivity_slurm_planning_report(
        report
    )
    slurm_job_plan_path = write_rabies_method_sensitivity_slurm_job_plan_table(
        output_root / "slurm-job-plan.tsv",
        slurm_planning_report,
    )
    slurm_assumptions_path = (
        write_rabies_method_sensitivity_slurm_assumptions_table(
            output_root / "slurm-estimation-assumptions.tsv",
            slurm_planning_report,
        )
    )
    slurm_summary_path = write_rabies_method_sensitivity_slurm_summary_json(
        output_root / "slurm-planning-summary.json",
        slurm_planning_report,
    )
    slurm_array_strategy_report = (
        build_rabies_method_sensitivity_slurm_array_strategy_report(
            slurm_planning_report
        )
    )
    slurm_array_partitions_path = (
        write_rabies_method_sensitivity_slurm_array_partitions_table(
            output_root / "slurm-array-partitions.tsv",
            slurm_array_strategy_report,
        )
    )
    slurm_array_members_path = (
        write_rabies_method_sensitivity_slurm_array_members_table(
            output_root / "slurm-array-members.tsv",
            slurm_array_strategy_report,
        )
    )
    slurm_array_strategy_path = (
        write_rabies_method_sensitivity_slurm_array_strategy_json(
            output_root / "slurm-array-strategy.json",
            slurm_array_strategy_report,
        )
    )
    slurm_array_scripts_root = (
        write_rabies_method_sensitivity_slurm_array_partition_scripts(
            output_root / "slurm-arrays",
            slurm_array_strategy_report,
        )
    )
    slurm_output_freshness_report = (
        build_rabies_method_sensitivity_slurm_output_freshness_report(
            output_root,
            dataset=report.dataset,
        )
    )
    slurm_output_freshness_path = (
        write_rabies_method_sensitivity_slurm_output_freshness_table(
            output_root / "slurm-output-freshness.tsv",
            slurm_output_freshness_report,
        )
    )
    slurm_output_freshness_checks_path = (
        write_rabies_method_sensitivity_slurm_output_freshness_checks_table(
            output_root / "slurm-output-freshness-checks.tsv",
            slurm_output_freshness_report,
        )
    )
    slurm_output_freshness_summary_path = (
        write_rabies_method_sensitivity_slurm_output_freshness_json(
            output_root / "slurm-output-freshness.json",
            slurm_output_freshness_report,
        )
    )
    slurm_status_report = build_rabies_method_sensitivity_slurm_status_report(
        output_root,
        dataset=report.dataset,
    )
    slurm_job_status_path = write_rabies_method_sensitivity_slurm_job_status_table(
        output_root / "slurm-job-status.tsv",
        slurm_status_report,
    )
    slurm_partition_status_path = (
        write_rabies_method_sensitivity_slurm_partition_status_table(
            output_root / "slurm-partition-status.tsv",
            slurm_status_report,
        )
    )
    slurm_workflow_status_path = write_rabies_method_sensitivity_slurm_status_json(
        output_root / "slurm-workflow-status.json",
        slurm_status_report,
    )
    slurm_failure_recovery_report = (
        build_rabies_method_sensitivity_slurm_failure_recovery_report(output_root)
    )
    slurm_failure_recovery_jobs_path = (
        write_rabies_method_sensitivity_slurm_failure_recovery_jobs_table(
            output_root / "slurm-failure-recovery-jobs.tsv",
            slurm_failure_recovery_report,
        )
    )
    slurm_failure_recovery_partitions_path = (
        write_rabies_method_sensitivity_slurm_failure_recovery_partitions_table(
            output_root / "slurm-failure-recovery-partitions.tsv",
            slurm_failure_recovery_report,
        )
    )
    slurm_failure_recovery_summary_path = (
        write_rabies_method_sensitivity_slurm_failure_recovery_summary_json(
            output_root / "slurm-failure-recovery-report.json",
            slurm_failure_recovery_report,
        )
    )
    slurm_failure_recovery_report_path = (
        write_rabies_method_sensitivity_slurm_failure_recovery_html_report(
            output_root / "slurm-failure-recovery-report.html",
            slurm_failure_recovery_report,
        )
    )
    slurm_job_evidence_report = (
        write_rabies_method_sensitivity_slurm_job_evidence_bundle(
            output_root / "slurm-job-evidence",
            bundle_root=output_root,
            dataset_id=report.dataset.dataset_id,
            workflow_prefix=report.dataset.workflow_prefix,
            execution_mode=report.execution_mode,
            parallel_workers=report.parallel_workers,
            task_records=report.task_records,
            variant_runs=report.variant_runs,
            array_strategy_report=slurm_array_strategy_report,
            execution_record_path=execution_record_path,
            workflow_manifest_path=output_root
            / "rabies-method-sensitivity.manifest.json",
        )
    )
    slurm_merge_report = build_rabies_method_sensitivity_slurm_merge_report(output_root)
    slurm_merge_checks_path = write_rabies_method_sensitivity_slurm_merge_checks_table(
        output_root / "slurm-merge-checks.tsv",
        slurm_merge_report,
    )
    slurm_merge_variants_path = (
        write_rabies_method_sensitivity_slurm_merge_variants_table(
            output_root / "slurm-merge-variants.tsv",
            slurm_merge_report,
        )
    )
    slurm_merge_summary_path = write_rabies_method_sensitivity_slurm_merge_summary_json(
        output_root / "slurm-merge-report.json",
        slurm_merge_report,
    )
    slurm_merge_report_path = write_rabies_method_sensitivity_slurm_merge_html_report(
        output_root / "slurm-merge-report.html",
        slurm_merge_report,
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
            "slurm_job_plan": slurm_job_plan_path,
            "slurm_assumptions": slurm_assumptions_path,
            "slurm_summary": slurm_summary_path,
            "slurm_array_partitions": slurm_array_partitions_path,
            "slurm_array_members": slurm_array_members_path,
            "slurm_array_strategy": slurm_array_strategy_path,
            "slurm_array_scripts_root": slurm_array_scripts_root,
            "slurm_job_evidence_root": slurm_job_evidence_report.evidence_root,
            "slurm_job_evidence_index": slurm_job_evidence_report.index_path,
            "slurm_job_evidence_summary": slurm_job_evidence_report.summary_path,
            "slurm_merge_checks": slurm_merge_checks_path,
            "slurm_merge_variants": slurm_merge_variants_path,
            "slurm_merge_summary": slurm_merge_summary_path,
            "slurm_merge_report": slurm_merge_report_path,
            "slurm_output_freshness": slurm_output_freshness_path,
            "slurm_output_freshness_checks": slurm_output_freshness_checks_path,
            "slurm_output_freshness_summary": slurm_output_freshness_summary_path,
            "slurm_job_status": slurm_job_status_path,
            "slurm_partition_status": slurm_partition_status_path,
            "slurm_workflow_status": slurm_workflow_status_path,
            "slurm_failure_recovery_jobs": slurm_failure_recovery_jobs_path,
            "slurm_failure_recovery_partitions": slurm_failure_recovery_partitions_path,
            "slurm_failure_recovery_summary": slurm_failure_recovery_summary_path,
            "slurm_failure_recovery_report": slurm_failure_recovery_report_path,
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
            "slurm_job_plan": slurm_job_plan_path,
            "slurm_assumptions": slurm_assumptions_path,
            "slurm_summary": slurm_summary_path,
            "slurm_array_partitions": slurm_array_partitions_path,
            "slurm_array_members": slurm_array_members_path,
            "slurm_array_strategy": slurm_array_strategy_path,
            "slurm_job_evidence_index": slurm_job_evidence_report.index_path,
            "slurm_job_evidence_summary": slurm_job_evidence_report.summary_path,
            "slurm_merge_checks": slurm_merge_checks_path,
            "slurm_merge_variants": slurm_merge_variants_path,
            "slurm_merge_summary": slurm_merge_summary_path,
            "slurm_merge_report": slurm_merge_report_path,
            "slurm_output_freshness": slurm_output_freshness_path,
            "slurm_output_freshness_checks": slurm_output_freshness_checks_path,
            "slurm_output_freshness_summary": slurm_output_freshness_summary_path,
            "slurm_job_status": slurm_job_status_path,
            "slurm_partition_status": slurm_partition_status_path,
            "slurm_workflow_status": slurm_workflow_status_path,
            "slurm_failure_recovery_jobs": slurm_failure_recovery_jobs_path,
            "slurm_failure_recovery_partitions": slurm_failure_recovery_partitions_path,
            "slurm_failure_recovery_summary": slurm_failure_recovery_summary_path,
            "slurm_failure_recovery_report": slurm_failure_recovery_report_path,
        },
    )
    slurm_storage_report = build_rabies_method_sensitivity_slurm_storage_report(
        output_root
    )
    slurm_storage_categories_path = (
        write_rabies_method_sensitivity_slurm_storage_categories_table(
            output_root / "slurm-storage-categories.tsv",
            slurm_storage_report,
        )
    )
    slurm_storage_variants_path = (
        write_rabies_method_sensitivity_slurm_storage_variants_table(
            output_root / "slurm-storage-variants.tsv",
            slurm_storage_report,
        )
    )
    slurm_storage_summary_path = (
        write_rabies_method_sensitivity_slurm_storage_summary_json(
            output_root / "slurm-storage-report.json",
            slurm_storage_report,
        )
    )
    slurm_storage_report_path = (
        write_rabies_method_sensitivity_slurm_storage_html_report(
            output_root / "slurm-storage-report.html",
            slurm_storage_report,
        )
    )
    slurm_output_explosion_report = (
        build_rabies_method_sensitivity_slurm_output_explosion_report(output_root)
    )
    slurm_output_explosion_checks_path = (
        write_rabies_method_sensitivity_slurm_output_explosion_checks_table(
            output_root / "slurm-output-explosion-checks.tsv",
            slurm_output_explosion_report,
        )
    )
    slurm_output_explosion_variants_path = (
        write_rabies_method_sensitivity_slurm_output_explosion_variants_table(
            output_root / "slurm-output-explosion-variants.tsv",
            slurm_output_explosion_report,
        )
    )
    slurm_output_explosion_summary_path = (
        write_rabies_method_sensitivity_slurm_output_explosion_summary_json(
            output_root / "slurm-output-explosion-report.json",
            slurm_output_explosion_report,
        )
    )
    slurm_output_explosion_report_path = (
        write_rabies_method_sensitivity_slurm_output_explosion_html_report(
            output_root / "slurm-output-explosion-report.html",
            slurm_output_explosion_report,
        )
    )
    slurm_tree_retention_report = (
        build_rabies_method_sensitivity_slurm_tree_retention_report(output_root)
    )
    slurm_tree_retention_checks_path = (
        write_rabies_method_sensitivity_slurm_tree_retention_checks_table(
            output_root / "slurm-tree-retention-checks.tsv",
            slurm_tree_retention_report,
        )
    )
    slurm_tree_retention_files_path = (
        write_rabies_method_sensitivity_slurm_tree_retention_files_table(
            output_root / "slurm-tree-retention-files.tsv",
            slurm_tree_retention_report,
        )
    )
    slurm_tree_retention_summary_path = (
        write_rabies_method_sensitivity_slurm_tree_retention_summary_json(
            output_root / "slurm-tree-retention-policy.json",
            slurm_tree_retention_report,
        )
    )
    slurm_tree_retention_report_path = (
        write_rabies_method_sensitivity_slurm_tree_retention_html_report(
            output_root / "slurm-tree-retention-policy.html",
            slurm_tree_retention_report,
        )
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
        slurm_job_plan_path,
        slurm_assumptions_path,
        slurm_summary_path,
        slurm_array_partitions_path,
        slurm_array_members_path,
        slurm_array_strategy_path,
        slurm_job_evidence_report.index_path,
        slurm_job_evidence_report.summary_path,
        slurm_storage_categories_path,
        slurm_storage_variants_path,
        slurm_storage_summary_path,
        slurm_storage_report_path,
        slurm_output_explosion_checks_path,
        slurm_output_explosion_variants_path,
        slurm_output_explosion_summary_path,
        slurm_output_explosion_report_path,
        slurm_tree_retention_checks_path,
        slurm_tree_retention_files_path,
        slurm_tree_retention_summary_path,
        slurm_tree_retention_report_path,
        slurm_merge_checks_path,
        slurm_merge_variants_path,
        slurm_merge_summary_path,
        slurm_merge_report_path,
        slurm_output_freshness_path,
        slurm_output_freshness_checks_path,
        slurm_output_freshness_summary_path,
        slurm_job_status_path,
        slurm_partition_status_path,
        slurm_workflow_status_path,
        slurm_failure_recovery_jobs_path,
        slurm_failure_recovery_partitions_path,
        slurm_failure_recovery_summary_path,
        slurm_failure_recovery_report_path,
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
            "slurm_job_plan": slurm_job_plan_path,
            "slurm_assumptions": slurm_assumptions_path,
            "slurm_summary": slurm_summary_path,
            "slurm_array_partitions": slurm_array_partitions_path,
            "slurm_array_members": slurm_array_members_path,
            "slurm_array_strategy": slurm_array_strategy_path,
            "slurm_job_evidence_index": slurm_job_evidence_report.index_path,
            "slurm_job_evidence_summary": slurm_job_evidence_report.summary_path,
            "slurm_storage_categories": slurm_storage_categories_path,
            "slurm_storage_variants": slurm_storage_variants_path,
            "slurm_storage_summary": slurm_storage_summary_path,
            "slurm_storage_report": slurm_storage_report_path,
            "slurm_output_explosion_checks": slurm_output_explosion_checks_path,
            "slurm_output_explosion_variants": slurm_output_explosion_variants_path,
            "slurm_output_explosion_summary": slurm_output_explosion_summary_path,
            "slurm_output_explosion_report": slurm_output_explosion_report_path,
            "slurm_tree_retention_checks": slurm_tree_retention_checks_path,
            "slurm_tree_retention_files": slurm_tree_retention_files_path,
            "slurm_tree_retention_summary": slurm_tree_retention_summary_path,
            "slurm_tree_retention_report": slurm_tree_retention_report_path,
            "slurm_merge_checks": slurm_merge_checks_path,
            "slurm_merge_variants": slurm_merge_variants_path,
            "slurm_merge_summary": slurm_merge_summary_path,
            "slurm_merge_report": slurm_merge_report_path,
            "slurm_output_freshness": slurm_output_freshness_path,
            "slurm_output_freshness_checks": slurm_output_freshness_checks_path,
            "slurm_output_freshness_summary": slurm_output_freshness_summary_path,
            "slurm_job_status": slurm_job_status_path,
            "slurm_partition_status": slurm_partition_status_path,
            "slurm_workflow_status": slurm_workflow_status_path,
            "slurm_failure_recovery_jobs": slurm_failure_recovery_jobs_path,
            "slurm_failure_recovery_partitions": slurm_failure_recovery_partitions_path,
            "slurm_failure_recovery_summary": slurm_failure_recovery_summary_path,
            "slurm_failure_recovery_report": slurm_failure_recovery_report_path,
            "reproducibility_checks": reproducibility_checks_path,
            "reproducibility_variant_audit": reproducibility_variant_audit_path,
            "reproducibility_audit": reproducibility_audit_path,
        },
        report_manifest_path=report_manifest_path,
        reproducibility_report=reproducibility_report,
        slurm_planning_report=slurm_planning_report,
        slurm_array_strategy_report=slurm_array_strategy_report,
        slurm_job_evidence_report=slurm_job_evidence_report,
        slurm_storage_report=slurm_storage_report,
        slurm_output_explosion_report=slurm_output_explosion_report,
        slurm_tree_retention_report=slurm_tree_retention_report,
        slurm_merge_report=slurm_merge_report,
        slurm_output_freshness_report=slurm_output_freshness_report,
        slurm_status_report=slurm_status_report,
        slurm_failure_recovery_report=slurm_failure_recovery_report,
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
        slurm_job_plan_path=slurm_job_plan_path,
        slurm_assumptions_path=slurm_assumptions_path,
        slurm_summary_path=slurm_summary_path,
        slurm_array_partitions_path=slurm_array_partitions_path,
        slurm_array_members_path=slurm_array_members_path,
        slurm_array_strategy_path=slurm_array_strategy_path,
        slurm_array_scripts_root=slurm_array_scripts_root,
        slurm_job_evidence_root=slurm_job_evidence_report.evidence_root,
        slurm_job_evidence_index_path=slurm_job_evidence_report.index_path,
        slurm_job_evidence_summary_path=slurm_job_evidence_report.summary_path,
        slurm_storage_categories_path=slurm_storage_categories_path,
        slurm_storage_variants_path=slurm_storage_variants_path,
        slurm_storage_summary_path=slurm_storage_summary_path,
        slurm_storage_report_path=slurm_storage_report_path,
        slurm_output_explosion_checks_path=slurm_output_explosion_checks_path,
        slurm_output_explosion_variants_path=slurm_output_explosion_variants_path,
        slurm_output_explosion_summary_path=slurm_output_explosion_summary_path,
        slurm_output_explosion_report_path=slurm_output_explosion_report_path,
        slurm_tree_retention_checks_path=slurm_tree_retention_checks_path,
        slurm_tree_retention_files_path=slurm_tree_retention_files_path,
        slurm_tree_retention_summary_path=slurm_tree_retention_summary_path,
        slurm_tree_retention_report_path=slurm_tree_retention_report_path,
        slurm_merge_checks_path=slurm_merge_checks_path,
        slurm_merge_variants_path=slurm_merge_variants_path,
        slurm_merge_summary_path=slurm_merge_summary_path,
        slurm_merge_report_path=slurm_merge_report_path,
        slurm_job_count=slurm_planning_report.job_count,
        slurm_total_estimated_core_hours=(
            slurm_planning_report.total_estimated_core_hours
        ),
        slurm_maximum_estimated_memory_mib=(
            slurm_planning_report.maximum_estimated_memory_mib
        ),
        slurm_maximum_estimated_wallclock_minutes=(
            slurm_planning_report.maximum_estimated_wallclock_minutes
        ),
        slurm_total_estimated_scratch_mib=(
            slurm_planning_report.total_estimated_scratch_mib
        ),
        slurm_total_estimated_output_mib=(
            slurm_planning_report.total_estimated_output_mib
        ),
        slurm_array_partition_count=slurm_array_strategy_report.partition_count,
        slurm_array_script_count=slurm_array_strategy_report.script_count,
        slurm_array_largest_partition_size=(
            slurm_array_strategy_report.largest_partition_size
        ),
        slurm_job_evidence_file_count=(
            slurm_job_evidence_report.total_artifact_file_count
        ),
        slurm_job_evidence_total_runtime_seconds=(
            slurm_job_evidence_report.total_runtime_seconds
        ),
        slurm_job_evidence_total_output_byte_count=(
            slurm_job_evidence_report.total_output_byte_count
        ),
        slurm_storage_total_estimated_mib=(
            slurm_storage_report.total_estimated_storage_mib
        ),
        slurm_storage_output_byte_count=slurm_storage_report.output_byte_count,
        slurm_storage_log_byte_count=slurm_storage_report.log_byte_count,
        slurm_storage_tree_byte_count=slurm_storage_report.tree_byte_count,
        slurm_storage_posterior_sample_byte_count=(
            slurm_storage_report.posterior_sample_byte_count
        ),
        slurm_storage_report_byte_count=slurm_storage_report.report_byte_count,
        slurm_storage_largest_variant_id=slurm_storage_report.largest_variant_id,
        slurm_output_explosion_status=(
            slurm_output_explosion_report.overall_risk_status
        ),
        slurm_output_explosion_global_issue_count=(
            slurm_output_explosion_report.global_issue_count
        ),
        slurm_output_explosion_warning_variant_count=(
            slurm_output_explosion_report.warning_variant_count
        ),
        slurm_output_explosion_high_risk_variant_count=(
            slurm_output_explosion_report.high_risk_variant_count
        ),
        slurm_tree_retention_status=(
            slurm_tree_retention_report.overall_policy_status
        ),
        slurm_tree_set_file_count=slurm_tree_retention_report.tree_set_file_count,
        slurm_tree_posterior_sample_file_count=(
            slurm_tree_retention_report.posterior_sample_file_count
        ),
        slurm_tree_thinning_recommended_file_count=(
            slurm_tree_retention_report.thinning_recommended_file_count
        ),
        slurm_tree_thinning_required_file_count=(
            slurm_tree_retention_report.thinning_required_file_count
        ),
        slurm_tree_compression_recommended_file_count=(
            slurm_tree_retention_report.compression_recommended_file_count
        ),
        slurm_tree_compression_required_file_count=(
            slurm_tree_retention_report.compression_required_file_count
        ),
        slurm_merge_status=slurm_merge_report.merge_status,
        slurm_merge_ready=slurm_merge_report.merge_ready,
        slurm_mergeable_variant_count=slurm_merge_report.mergeable_variant_count,
        slurm_merge_failed_check_count=slurm_merge_report.failed_check_count,
        slurm_output_freshness_path=slurm_output_freshness_path,
        slurm_output_freshness_checks_path=slurm_output_freshness_checks_path,
        slurm_output_freshness_summary_path=slurm_output_freshness_summary_path,
        slurm_job_status_path=slurm_job_status_path,
        slurm_partition_status_path=slurm_partition_status_path,
        slurm_workflow_status_path=slurm_workflow_status_path,
        slurm_failure_recovery_jobs_path=slurm_failure_recovery_jobs_path,
        slurm_failure_recovery_partitions_path=slurm_failure_recovery_partitions_path,
        slurm_failure_recovery_summary_path=slurm_failure_recovery_summary_path,
        slurm_failure_recovery_report_path=slurm_failure_recovery_report_path,
        slurm_output_freshness_check_count=slurm_output_freshness_report.check_count,
        slurm_output_freshness_failed_check_count=(
            slurm_output_freshness_report.failed_check_count
        ),
        slurm_fresh_output_job_count=slurm_output_freshness_report.fresh_job_count,
        slurm_stale_output_job_count=slurm_output_freshness_report.stale_job_count,
        slurm_completed_job_count=slurm_status_report.completed_job_count,
        slurm_failed_job_count=slurm_status_report.failed_job_count,
        slurm_pending_job_count=slurm_status_report.pending_job_count,
        slurm_stale_job_count=slurm_status_report.stale_job_count,
        slurm_failure_recovery_status=(
            slurm_failure_recovery_report.overall_recovery_status
        ),
        slurm_failure_recovery_rerunnable_job_count=(
            slurm_failure_recovery_report.rerunnable_job_count
        ),
        slurm_failure_recovery_blocked_job_count=(
            slurm_failure_recovery_report.blocked_job_count
        ),
        slurm_failure_recovery_partition_count=(
            slurm_failure_recovery_report.recovery_partition_count
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

def _write_workflow_summary_table(
    path: Path, report: RabiesMethodSensitivityPanelWorkflowReport
) -> Path:
    serious_conflicts = [
        variant.inference_comparison.conclusion_summary.serious_conflict_count
        for variant in report.variant_runs
    ]
    rows = [
        [
            "dataset_id",
            "variant_count",
            "stable_clade_count",
            "changed_clade_count",
            "preprocessing_change_pair_count",
            "rooted_engine_change_variant_count",
            "serious_conflict_variant_count",
            "maximum_serious_conflict_count",
        ],
        [
            report.dataset.dataset_id,
            str(len(report.variant_runs)),
            str(len(report.stable_clade_rows)),
            str(len(report.changed_clade_rows)),
            str(
                sum(
                    1
                    for row in report.preprocessing_comparison_rows
                    if row.robinson_foulds_distance > 0
                    or row.same_taxa_different_rooting
                )
            ),
            str(
                sum(
                    1
                    for variant in report.variant_runs
                    if variant.rooted_engine_comparison.robinson_foulds_distance > 0
                    or variant.rooted_engine_comparison.same_taxa_different_rooting
                )
            ),
            str(sum(1 for value in serious_conflicts if value > 0)),
            str(max(serious_conflicts)),
        ],
    ]
    return _write_tsv(path, rows)


def _write_variant_summary_table(
    path: Path, report: RabiesMethodSensitivityPanelWorkflowReport
) -> Path:
    rows = [
        [
            "variant_id",
            "alignment_mode",
            "trimming_mode",
            "trim_gap_threshold",
            "alignment_length",
            "trimmed_alignment_length",
            "selected_model",
            "minimum_support",
            "maximum_support",
            "stable_clade_count",
            "engine_specific_clade_count",
            "serious_conflict_count",
            "rooted_engine_rf_distance",
            "rooted_engine_same_taxa_different_rooting",
        ]
    ]
    for variant in report.variant_runs:
        summary = variant.inference_comparison.conclusion_summary
        rows.append(
            [
                variant.config.variant_id,
                variant.config.alignment_mode,
                variant.config.trimming_mode,
                _format_float(variant.config.trim_gap_threshold),
                str(variant.alignment_length),
                str(variant.trimmed_alignment_length),
                variant.inference_comparison.selected_model,
                _format_optional_float(
                    variant.inference_comparison.iqtree_support_workflow.bootstrap_support_summary.minimum_support
                    if variant.inference_comparison.iqtree_support_workflow.bootstrap_support_summary
                    is not None
                    else None
                ),
                _format_optional_float(
                    variant.inference_comparison.iqtree_support_workflow.bootstrap_support_summary.maximum_support
                    if variant.inference_comparison.iqtree_support_workflow.bootstrap_support_summary
                    is not None
                    else None
                ),
                str(summary.stable_clade_count),
                str(summary.engine_specific_clade_count),
                str(summary.serious_conflict_count),
                str(variant.rooted_engine_comparison.robinson_foulds_distance),
                str(
                    variant.rooted_engine_comparison.same_taxa_different_rooting
                ).lower(),
            ]
        )
    return _write_tsv(path, rows)


def _write_preprocessing_comparison_table(
    path: Path,
    rows: tuple[RabiesMethodSensitivityPreprocessingComparisonRow, ...],
) -> Path:
    rendered = [
        [
            "left_variant_id",
            "right_variant_id",
            "comparison_axis",
            "robinson_foulds_distance",
            "normalized_robinson_foulds",
            "same_taxa_different_rooting",
        ]
    ]
    for row in rows:
        rendered.append(
            [
                row.left_variant_id,
                row.right_variant_id,
                row.comparison_axis,
                str(row.robinson_foulds_distance),
                _format_float(row.normalized_robinson_foulds),
                str(row.same_taxa_different_rooting).lower(),
            ]
        )
    return _write_tsv(path, rendered)


def _write_clade_table(
    path: Path, rows: tuple[RabiesMethodSensitivityCladeRow, ...]
) -> Path:
    rendered = [
        [
            "split_id",
            "conclusion_class",
            "evidence_class",
            "occurrence_count",
            "variant_count",
            "detail",
        ]
    ]
    for row in rows:
        rendered.append(
            [
                row.split_id,
                row.conclusion_class,
                row.evidence_class,
                str(row.occurrence_count),
                str(row.variant_count),
                row.detail,
            ]
        )
    return _write_tsv(path, rendered)


def _write_conclusion_summary_table(
    path: Path, rows: tuple[RabiesMethodSensitivityConclusionRow, ...]
) -> Path:
    rendered = [
        [
            "conclusion_id",
            "method_axis",
            "stability_status",
            "claim",
            "evidence",
            "caution",
        ]
    ]
    for row in rows:
        rendered.append(
            [
                row.conclusion_id,
                row.method_axis,
                row.stability_status,
                row.claim,
                row.evidence,
                row.caution,
            ]
        )
    return _write_tsv(path, rendered)


def _write_parallel_execution_summary_table(
    path: Path, report: RabiesMethodSensitivityPanelWorkflowReport
) -> Path:
    rows = [
        [
            "variant_id",
            "label",
            "execution_mode",
            "status",
            "log_path",
            "error_code",
        ]
    ]
    for task in report.task_records:
        rows.append(
            [
                task.variant_id,
                task.label,
                task.execution_mode,
                task.status,
                Path("parallel-logs", task.log_path.name).as_posix(),
                "" if task.error_code is None else task.error_code,
            ]
        )
    return _write_tsv(path, rows)


def _write_resolved_config(
    path: Path, report: RabiesMethodSensitivityPanelWorkflowReport
) -> Path:
    payload = {
        "dataset_id": report.dataset.dataset_id,
        "label": report.dataset.label,
        "sequence_type": report.dataset.sequence_type,
        "workflow_prefix": report.dataset.workflow_prefix,
        "outgroup_taxa": list(report.dataset.outgroup_taxa),
        "iqtree_seed": report.iqtree_seed,
        "iqtree_threads": report.iqtree_threads,
        "bootstrap_replicates": report.bootstrap_replicates,
        "parallel_workers": report.parallel_workers,
        "execution_mode": report.execution_mode,
        "selected_variant_ids": [
            variant.variant_id for variant in report.dataset.variants
        ],
        "input_checksums": {
            "sequences.fasta": _sha256(report.dataset.sequences_path),
            "metadata.csv": _sha256(report.dataset.metadata_path),
        },
        "variants": [
            {
                "variant_id": variant.variant_id,
                "label": variant.label,
                "alignment_mode": variant.alignment_mode,
                "trimming_mode": variant.trimming_mode,
                "trim_gap_threshold": variant.trim_gap_threshold,
            }
            for variant in report.dataset.variants
        ],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _copy_task_logs(
    output_root: Path,
    task_records: tuple[RabiesMethodSensitivityTaskRecord, ...],
) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    for task in task_records:
        _copy_output(task.log_path, output_root / task.log_path.name)
    return output_root


def _write_variant_outputs(
    output_root: Path, variant_runs: tuple[RabiesMethodSensitivityVariantRun, ...]
) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    for variant in variant_runs:
        variant_root = output_root / variant.config.variant_id
        variant_root.mkdir(parents=True, exist_ok=True)
        _copy_output(
            variant.alignment_workflow.output_paths["alignment"],
            variant_root / f"{variant.config.variant_id}.aln",
        )
        _copy_output(
            variant.trimming_workflow.output_paths["trimmed_alignment"],
            variant_root / f"{variant.config.variant_id}.trimmed.aln",
        )
        _copy_output(
            variant.inference_comparison.output_paths["fasttree_tree"],
            variant_root / "fasttree.nwk",
        )
        _copy_output(
            variant.inference_comparison.output_paths["iqtree_support_tree"],
            variant_root / "iqtree-support.nwk",
        )
        _copy_output(
            variant.rooted_fasttree_path,
            variant_root / "rooted-fasttree.nwk",
        )
        _copy_output(
            variant.rooted_iqtree_path,
            variant_root / "rooted-iqtree-support.nwk",
        )
        _write_rooting_summary_table(
            variant_root / "rooting-summary.tsv",
            variant,
        )
        _copy_output(
            variant.inference_comparison.output_paths["stability_summary"],
            variant_root / "unrooted-stability-summary.tsv",
        )
        _copy_output(
            variant.inference_comparison.output_paths["conclusion_table"],
            variant_root / "unrooted-conclusions.tsv",
        )
        _copy_output(
            variant.inference_comparison.output_paths["support_weighted_conflicts"],
            variant_root / "unrooted-support-weighted-conflicts.tsv",
        )
        _copy_output(
            variant.inference_comparison.output_paths["shared_clades"],
            variant_root / "unrooted-shared-clades.tsv",
        )
        _copy_output(
            variant.inference_comparison.output_paths["conflicting_clades"],
            variant_root / "unrooted-conflicting-clades.tsv",
        )
        _copy_output(
            variant.inference_comparison.output_paths["comparison_table"],
            variant_root / "unrooted-comparison.tsv",
        )
        _copy_output(
            variant.rooted_engine_comparison_table_path,
            variant_root / "rooted-engine-comparison.tsv",
        )
    return output_root


def _write_rooting_summary_table(
    path: Path, variant: RabiesMethodSensitivityVariantRun
) -> Path:
    rows = [
        [
            "engine_name",
            "requested_taxa",
            "matched_taxa",
            "outgroup_monophyletic",
            "rooted_outgroup_taxa",
            "warning_count",
        ],
        [
            "fasttree",
            ",".join(variant.fasttree_rooting.requested_taxa),
            ",".join(variant.fasttree_rooting.matched_taxa),
            _format_optional_bool(variant.fasttree_rooting.outgroup_monophyletic),
            ",".join(variant.fasttree_rooting.rooted_outgroup_taxa),
            str(len(variant.fasttree_rooting.warnings)),
        ],
        [
            "iqtree",
            ",".join(variant.iqtree_rooting.requested_taxa),
            ",".join(variant.iqtree_rooting.matched_taxa),
            _format_optional_bool(variant.iqtree_rooting.outgroup_monophyletic),
            ",".join(variant.iqtree_rooting.rooted_outgroup_taxa),
            str(len(variant.iqtree_rooting.warnings)),
        ],
    ]
    return _write_tsv(path, rows)


def _write_manifest(
    path: Path,
    *,
    report: RabiesMethodSensitivityPanelWorkflowReport,
    bundle_paths: dict[str, Path],
) -> Path:
    payload = {
        "dataset_id": report.dataset.dataset_id,
        "label": report.dataset.label,
        "report_kind": "rabies_method_sensitivity_workflow_bundle",
        "variant_count": len(report.variant_runs),
        "parallel_execution": {
            "execution_mode": report.execution_mode,
            "parallel_workers": report.parallel_workers,
            "requested_task_count": len(report.task_records),
            "completed_task_count": len(
                [task for task in report.task_records if task.status == "succeeded"]
            ),
            "failed_task_count": len(
                [task for task in report.task_records if task.status != "succeeded"]
            ),
        },
        "task_records": [
            {
                "variant_id": task.variant_id,
                "label": task.label,
                "status": task.status,
                "execution_mode": task.execution_mode,
                "log_path": Path("parallel-logs", task.log_path.name).as_posix(),
                "output_root": Path("variants", task.variant_id).as_posix(),
                "error_code": task.error_code,
                "error_message": task.error_message,
            }
            for task in report.task_records
        ],
        "output_paths": {
            key: value.name
            if value.parent == path.parent
            else value.relative_to(path.parent).as_posix()
            for key, value in bundle_paths.items()
        },
        "output_checksums": {
            key: _sha256(value)
            for key, value in bundle_paths.items()
            if value.is_file()
        },
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _relative_bundle_path(base_path: Path, value: Path) -> str:
    return Path(os.path.relpath(value, start=base_path.parent)).as_posix()


def _write_report_manifest(
    path: Path,
    *,
    report: RabiesMethodSensitivityPanelWorkflowReport,
    bundle_paths: dict[str, Path],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    linked_files = {
        key: value for key, value in bundle_paths.items() if value.is_file()
    }
    payload = {
        "dataset_id": report.dataset.dataset_id,
        "report_kind": "rabies_method_sensitivity_html_report",
        "variant_count": len(report.variant_runs),
        "parallel_workers": report.parallel_workers,
        "execution_mode": report.execution_mode,
        "stable_clade_count": len(report.stable_clade_rows),
        "changed_clade_count": len(report.changed_clade_rows),
        "linked_artifact_count": len(linked_files),
        "linked_artifacts": {
            key: {
                "path": _relative_bundle_path(path, value),
                "byte_count": value.stat().st_size,
                "sha256": _sha256(value),
            }
            for key, value in linked_files.items()
        },
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _write_report(
    path: Path,
    *,
    report: RabiesMethodSensitivityPanelWorkflowReport,
    bundle_paths: dict[str, Path],
    report_manifest_path: Path,
    reproducibility_report: RabiesMethodSensitivityReproducibilityAuditReport,
    slurm_planning_report: RabiesMethodSensitivitySlurmPlanningReport,
    slurm_array_strategy_report: RabiesMethodSensitivitySlurmArrayStrategyReport,
    slurm_job_evidence_report: RabiesMethodSensitivitySlurmJobEvidenceReport,
    slurm_storage_report: RabiesMethodSensitivitySlurmStorageReport,
    slurm_output_explosion_report: RabiesMethodSensitivitySlurmOutputExplosionReport,
    slurm_tree_retention_report: RabiesMethodSensitivitySlurmTreeRetentionReport,
    slurm_merge_report: RabiesMethodSensitivitySlurmMergeReport,
    slurm_output_freshness_report: RabiesMethodSensitivitySlurmOutputFreshnessReport,
    slurm_status_report: RabiesMethodSensitivitySlurmStatusReport,
    slurm_failure_recovery_report: RabiesMethodSensitivitySlurmFailureRecoveryReport,
) -> Path:
    variant_lines = [
        (
            f"{variant.config.variant_id}: model {variant.inference_comparison.selected_model}; "
            f"unrooted serious conflicts {variant.inference_comparison.conclusion_summary.serious_conflict_count}; "
            f"rooted engine RF {variant.rooted_engine_comparison.robinson_foulds_distance}"
        )
        for variant in report.variant_runs
    ]
    conclusion_lines: list[str] = []
    for row in report.conclusion_rows:
        conclusion_lines.extend(
            [
                f"{row.conclusion_id}: {row.stability_status}",
                f"claim: {row.claim}",
                f"evidence: {row.evidence}",
                f"caution: {row.caution}",
                "",
            ]
        )
    sections = [
        (
            "workflow-summary",
            "\n".join(
                [
                    f"dataset: {report.dataset.label}",
                    f"variants: {len(report.variant_runs)}",
                    f"execution mode: {report.execution_mode}",
                    f"parallel workers: {report.parallel_workers}",
                    f"stable clades across variants: {len(report.stable_clade_rows)}",
                    f"changed clades across variants: {len(report.changed_clade_rows)}",
                    *variant_lines,
                ]
            ),
        ),
        (
            "conclusions",
            "\n".join(conclusion_lines).strip(),
        ),
        (
            "reproducibility-audit",
            "\n".join(
                [
                    f"all passed: {str(getattr(reproducibility_report, 'all_passed')).lower()}",
                    f"top-level checks: {getattr(reproducibility_report, 'check_count')}",
                    (
                        "failed top-level checks: "
                        f"{getattr(reproducibility_report, 'failed_check_count')}"
                    ),
                    (
                        "failed variants: "
                        f"{getattr(reproducibility_report, 'failed_variant_count')}"
                    ),
                ]
            ),
        ),
        (
            "slurm-job-planning",
            "\n".join(
                [
                    f"planned jobs: {slurm_planning_report.job_count}",
                    (
                        "maximum estimated memory MiB: "
                        f"{slurm_planning_report.maximum_estimated_memory_mib}"
                    ),
                    (
                        "maximum estimated wallclock minutes: "
                        f"{slurm_planning_report.maximum_estimated_wallclock_minutes}"
                    ),
                    (
                        "total estimated core hours: "
                        f"{_format_float(slurm_planning_report.total_estimated_core_hours)}"
                    ),
                    (
                        "total estimated scratch MiB: "
                        f"{slurm_planning_report.total_estimated_scratch_mib}"
                    ),
                    (
                        "total estimated output MiB: "
                        f"{slurm_planning_report.total_estimated_output_mib}"
                    ),
                ]
            ),
        ),
        (
            "slurm-array-partitioning",
            "\n".join(
                [
                    (
                        "array partitions: "
                        f"{slurm_array_strategy_report.partition_count}"
                    ),
                    (
                        "largest partition size: "
                        f"{slurm_array_strategy_report.largest_partition_size}"
                    ),
                    (
                        "partition scripts: "
                        f"{slurm_array_strategy_report.script_count}"
                    ),
                    (
                        "total array jobs: "
                        f"{slurm_array_strategy_report.total_job_count}"
                    ),
                ]
            ),
        ),
        (
            "slurm-job-evidence",
            "\n".join(
                [
                    f"job evidence packages: {slurm_job_evidence_report.job_count}",
                    (
                        "job evidence artifact files: "
                        f"{slurm_job_evidence_report.total_artifact_file_count}"
                    ),
                    (
                        "job evidence total runtime seconds: "
                        f"{_format_float(slurm_job_evidence_report.total_runtime_seconds)}"
                    ),
                    (
                        "job evidence total output bytes: "
                        f"{slurm_job_evidence_report.total_output_byte_count}"
                    ),
                ]
            ),
        ),
        (
            "slurm-storage-estimate",
            "\n".join(
                [
                    (
                        "estimated retained storage MiB: "
                        f"{slurm_storage_report.total_estimated_storage_mib}"
                    ),
                    (
                        "workflow outputs bytes: "
                        f"{slurm_storage_report.output_byte_count}"
                    ),
                    f"log bytes: {slurm_storage_report.log_byte_count}",
                    f"tree bytes: {slurm_storage_report.tree_byte_count}",
                    (
                        "posterior sample bytes: "
                        f"{slurm_storage_report.posterior_sample_byte_count}"
                    ),
                    f"review artifact bytes: {slurm_storage_report.report_byte_count}",
                    (
                        "largest variant by retained bytes: "
                        f"{slurm_storage_report.largest_variant_id}"
                    ),
                ]
            ),
        ),
        (
            "slurm-output-explosion-risk",
            "\n".join(
                [
                    (
                        "overall risk status: "
                        f"{slurm_output_explosion_report.overall_risk_status}"
                    ),
                    (
                        "global issues: "
                        f"{slurm_output_explosion_report.global_issue_count}"
                    ),
                    (
                        "warning variants: "
                        f"{slurm_output_explosion_report.warning_variant_count}"
                    ),
                    (
                        "high-risk variants: "
                        f"{slurm_output_explosion_report.high_risk_variant_count}"
                    ),
                    (
                        "largest variant output share: "
                        f"{_format_float(slurm_output_explosion_report.largest_variant_output_share)}"
                    ),
                    (
                        "posterior sample bytes: "
                        f"{slurm_output_explosion_report.total_posterior_sample_byte_count}"
                    ),
                ]
            ),
        ),
        (
            "slurm-tree-retention-policy",
            "\n".join(
                [
                    (
                        "overall policy status: "
                        f"{slurm_tree_retention_report.overall_policy_status}"
                    ),
                    (
                        "tree-set files: "
                        f"{slurm_tree_retention_report.tree_set_file_count}"
                    ),
                    (
                        "posterior sample files: "
                        f"{slurm_tree_retention_report.posterior_sample_file_count}"
                    ),
                    (
                        "thinning required files: "
                        f"{slurm_tree_retention_report.thinning_required_file_count}"
                    ),
                    (
                        "compression required files: "
                        f"{slurm_tree_retention_report.compression_required_file_count}"
                    ),
                    (
                        "largest tree set: "
                        f"{slurm_tree_retention_report.largest_tree_set_path}"
                    ),
                ]
            ),
        ),
        (
            "slurm-merge-report",
            "\n".join(
                [
                    f"merge status: {slurm_merge_report.merge_status}",
                    f"merge ready: {str(slurm_merge_report.merge_ready).lower()}",
                    (
                        "mergeable variants: "
                        f"{slurm_merge_report.mergeable_variant_count}"
                    ),
                    (
                        "failed merge checks: "
                        f"{slurm_merge_report.failed_check_count}"
                    ),
                    (
                        "merged stable clades: "
                        f"{slurm_merge_report.stable_clade_count}"
                    ),
                    (
                        "merged changed clades: "
                        f"{slurm_merge_report.changed_clade_count}"
                    ),
                ]
            ),
        ),
        (
            "slurm-output-freshness",
            "\n".join(
                [
                    (
                        "all outputs fresh: "
                        f"{str(slurm_output_freshness_report.all_outputs_fresh).lower()}"
                    ),
                    (
                        "fresh jobs: "
                        f"{slurm_output_freshness_report.fresh_job_count}"
                    ),
                    (
                        "stale jobs: "
                        f"{slurm_output_freshness_report.stale_job_count}"
                    ),
                    (
                        "freshness checks: "
                        f"{slurm_output_freshness_report.check_count}"
                    ),
                    (
                        "failed freshness checks: "
                        f"{slurm_output_freshness_report.failed_check_count}"
                    ),
                ]
            ),
        ),
        (
            "slurm-workflow-status",
            "\n".join(
                [
                    f"workflow status: {slurm_status_report.workflow_status}",
                    f"active run state: {slurm_status_report.active_run_state}",
                    f"completed jobs: {slurm_status_report.completed_job_count}",
                    f"failed jobs: {slurm_status_report.failed_job_count}",
                    f"pending jobs: {slurm_status_report.pending_job_count}",
                    f"stale jobs: {slurm_status_report.stale_job_count}",
                ]
            ),
        ),
        (
            "slurm-failure-recovery",
            "\n".join(
                [
                    (
                        "overall recovery status: "
                        f"{slurm_failure_recovery_report.overall_recovery_status}"
                    ),
                    (
                        "rerunnable jobs: "
                        f"{slurm_failure_recovery_report.rerunnable_job_count}"
                    ),
                    (
                        "blocked jobs: "
                        f"{slurm_failure_recovery_report.blocked_job_count}"
                    ),
                    (
                        "recovery partitions: "
                        f"{slurm_failure_recovery_report.recovery_partition_count}"
                    ),
                    (
                        "workflow state: "
                        f"{slurm_failure_recovery_report.workflow_status}"
                    ),
                    (
                        "active run state: "
                        f"{slurm_failure_recovery_report.active_run_state}"
                    ),
                ]
            ),
        ),
        (
            "artifacts",
            "\n".join(
                [
                    f"workflow summary: {bundle_paths['workflow_summary'].name}",
                    f"variant summary: {bundle_paths['variant_summary'].name}",
                    f"parallel execution summary: {bundle_paths['parallel_summary'].name}",
                    f"preprocessing rooted comparisons: {bundle_paths['preprocessing_comparison'].name}",
                    f"stable clades: {bundle_paths['stable_clades'].name}",
                    f"changed clades: {bundle_paths['changed_clades'].name}",
                    f"method conclusions: {bundle_paths['conclusion_summary'].name}",
                    f"resolved config: {bundle_paths['config'].name}",
                    f"workflow manifest: {bundle_paths['workflow_manifest'].name}",
                    f"slurm job plan: {bundle_paths['slurm_job_plan'].name}",
                    (
                        "slurm estimation assumptions: "
                        f"{bundle_paths['slurm_assumptions'].name}"
                    ),
                    f"slurm planning summary: {bundle_paths['slurm_summary'].name}",
                    (
                        "slurm array partitions: "
                        f"{bundle_paths['slurm_array_partitions'].name}"
                    ),
                    (
                        "slurm array members: "
                        f"{bundle_paths['slurm_array_members'].name}"
                    ),
                    (
                        "slurm array strategy: "
                        f"{bundle_paths['slurm_array_strategy'].name}"
                    ),
                    (
                        "slurm job evidence index: "
                        f"{bundle_paths['slurm_job_evidence_index'].name}"
                    ),
                    (
                        "slurm job evidence summary: "
                        f"{bundle_paths['slurm_job_evidence_summary'].name}"
                    ),
                    (
                        "slurm storage categories: "
                        f"{bundle_paths['slurm_storage_categories'].name}"
                    ),
                    (
                        "slurm storage variants: "
                        f"{bundle_paths['slurm_storage_variants'].name}"
                    ),
                    (
                        "slurm storage summary: "
                        f"{bundle_paths['slurm_storage_summary'].name}"
                    ),
                    (
                        "slurm storage report: "
                        f"{bundle_paths['slurm_storage_report'].name}"
                    ),
                    (
                        "slurm output explosion checks: "
                        f"{bundle_paths['slurm_output_explosion_checks'].name}"
                    ),
                    (
                        "slurm output explosion variants: "
                        f"{bundle_paths['slurm_output_explosion_variants'].name}"
                    ),
                    (
                        "slurm output explosion summary: "
                        f"{bundle_paths['slurm_output_explosion_summary'].name}"
                    ),
                    (
                        "slurm output explosion report: "
                        f"{bundle_paths['slurm_output_explosion_report'].name}"
                    ),
                    (
                        "slurm tree retention checks: "
                        f"{bundle_paths['slurm_tree_retention_checks'].name}"
                    ),
                    (
                        "slurm tree retention files: "
                        f"{bundle_paths['slurm_tree_retention_files'].name}"
                    ),
                    (
                        "slurm tree retention summary: "
                        f"{bundle_paths['slurm_tree_retention_summary'].name}"
                    ),
                    (
                        "slurm tree retention report: "
                        f"{bundle_paths['slurm_tree_retention_report'].name}"
                    ),
                    (
                        "slurm merge checks: "
                        f"{bundle_paths['slurm_merge_checks'].name}"
                    ),
                    (
                        "slurm merge variants: "
                        f"{bundle_paths['slurm_merge_variants'].name}"
                    ),
                    (
                        "slurm merge summary: "
                        f"{bundle_paths['slurm_merge_summary'].name}"
                    ),
                    (
                        "slurm merge report: "
                        f"{bundle_paths['slurm_merge_report'].name}"
                    ),
                    (
                        "slurm output freshness: "
                        f"{bundle_paths['slurm_output_freshness'].name}"
                    ),
                    (
                        "slurm output freshness checks: "
                        f"{bundle_paths['slurm_output_freshness_checks'].name}"
                    ),
                    (
                        "slurm output freshness summary: "
                        f"{bundle_paths['slurm_output_freshness_summary'].name}"
                    ),
                    f"slurm job status: {bundle_paths['slurm_job_status'].name}",
                    (
                        "slurm partition status: "
                        f"{bundle_paths['slurm_partition_status'].name}"
                    ),
                    (
                        "slurm workflow status: "
                        f"{bundle_paths['slurm_workflow_status'].name}"
                    ),
                    (
                        "slurm failure recovery jobs: "
                        f"{bundle_paths['slurm_failure_recovery_jobs'].name}"
                    ),
                    (
                        "slurm failure recovery partitions: "
                        f"{bundle_paths['slurm_failure_recovery_partitions'].name}"
                    ),
                    (
                        "slurm failure recovery summary: "
                        f"{bundle_paths['slurm_failure_recovery_summary'].name}"
                    ),
                    (
                        "slurm failure recovery report: "
                        f"{bundle_paths['slurm_failure_recovery_report'].name}"
                    ),
                    f"reproducibility checks: {bundle_paths['reproducibility_checks'].name}",
                    (
                        "reproducibility variant audit: "
                        f"{bundle_paths['reproducibility_variant_audit'].name}"
                    ),
                    (
                        "reproducibility audit: "
                        f"{bundle_paths['reproducibility_audit'].name}"
                    ),
                ]
            ),
        ),
    ]
    report_manifest = json.loads(report_manifest_path.read_text(encoding="utf-8"))
    artifact_links = [
        (
            key.replace("_", "-"),
            _relative_bundle_path(path, value),
            f"{value.stat().st_size} bytes",
        )
        for key, value in bundle_paths.items()
        if value.is_file()
    ]
    return write_html_report(
        title="Rabies Method-Sensitivity Report",
        sections=sections,
        out_path=path,
        embedded_json={
            "dataset_id": report.dataset.dataset_id,
            "variant_count": len(report.variant_runs),
            "parallel_workers": report.parallel_workers,
            "execution_mode": report.execution_mode,
            "stable_clade_count": len(report.stable_clade_rows),
            "changed_clade_count": len(report.changed_clade_rows),
            "report_manifest_path": _relative_bundle_path(path, report_manifest_path),
            "reproducibility_passed": getattr(reproducibility_report, "all_passed"),
            "reproducibility_check_count": getattr(
                reproducibility_report, "check_count"
            ),
            "reproducibility_failed_check_count": getattr(
                reproducibility_report, "failed_check_count"
            ),
            "reproducibility_failed_variant_count": getattr(
                reproducibility_report, "failed_variant_count"
            ),
            "slurm_job_count": slurm_planning_report.job_count,
            "slurm_total_estimated_core_hours": (
                slurm_planning_report.total_estimated_core_hours
            ),
            "slurm_maximum_estimated_memory_mib": (
                slurm_planning_report.maximum_estimated_memory_mib
            ),
            "slurm_maximum_estimated_wallclock_minutes": (
                slurm_planning_report.maximum_estimated_wallclock_minutes
            ),
            "slurm_array_partition_count": (
                slurm_array_strategy_report.partition_count
            ),
            "slurm_array_script_count": slurm_array_strategy_report.script_count,
            "slurm_array_largest_partition_size": (
                slurm_array_strategy_report.largest_partition_size
            ),
            "slurm_job_evidence_file_count": (
                slurm_job_evidence_report.total_artifact_file_count
            ),
            "slurm_job_evidence_total_runtime_seconds": (
                slurm_job_evidence_report.total_runtime_seconds
            ),
            "slurm_job_evidence_total_output_byte_count": (
                slurm_job_evidence_report.total_output_byte_count
            ),
            "slurm_storage_total_estimated_mib": (
                slurm_storage_report.total_estimated_storage_mib
            ),
            "slurm_storage_output_byte_count": (
                slurm_storage_report.output_byte_count
            ),
            "slurm_storage_log_byte_count": slurm_storage_report.log_byte_count,
            "slurm_storage_tree_byte_count": slurm_storage_report.tree_byte_count,
            "slurm_storage_posterior_sample_byte_count": (
                slurm_storage_report.posterior_sample_byte_count
            ),
            "slurm_storage_report_byte_count": (
                slurm_storage_report.report_byte_count
            ),
            "slurm_storage_largest_variant_id": (
                slurm_storage_report.largest_variant_id
            ),
            "slurm_output_explosion_status": (
                slurm_output_explosion_report.overall_risk_status
            ),
            "slurm_output_explosion_global_issue_count": (
                slurm_output_explosion_report.global_issue_count
            ),
            "slurm_output_explosion_warning_variant_count": (
                slurm_output_explosion_report.warning_variant_count
            ),
            "slurm_output_explosion_high_risk_variant_count": (
                slurm_output_explosion_report.high_risk_variant_count
            ),
            "slurm_tree_retention_status": (
                slurm_tree_retention_report.overall_policy_status
            ),
            "slurm_tree_set_file_count": (
                slurm_tree_retention_report.tree_set_file_count
            ),
            "slurm_tree_posterior_sample_file_count": (
                slurm_tree_retention_report.posterior_sample_file_count
            ),
            "slurm_tree_thinning_recommended_file_count": (
                slurm_tree_retention_report.thinning_recommended_file_count
            ),
            "slurm_tree_thinning_required_file_count": (
                slurm_tree_retention_report.thinning_required_file_count
            ),
            "slurm_tree_compression_recommended_file_count": (
                slurm_tree_retention_report.compression_recommended_file_count
            ),
            "slurm_tree_compression_required_file_count": (
                slurm_tree_retention_report.compression_required_file_count
            ),
            "slurm_merge_status": slurm_merge_report.merge_status,
            "slurm_merge_ready": slurm_merge_report.merge_ready,
            "slurm_mergeable_variant_count": (
                slurm_merge_report.mergeable_variant_count
            ),
            "slurm_merge_failed_check_count": (
                slurm_merge_report.failed_check_count
            ),
            "slurm_output_freshness_check_count": (
                slurm_output_freshness_report.check_count
            ),
            "slurm_output_freshness_failed_check_count": (
                slurm_output_freshness_report.failed_check_count
            ),
            "slurm_fresh_output_job_count": (
                slurm_output_freshness_report.fresh_job_count
            ),
            "slurm_stale_output_job_count": (
                slurm_output_freshness_report.stale_job_count
            ),
            "slurm_completed_job_count": slurm_status_report.completed_job_count,
            "slurm_failed_job_count": slurm_status_report.failed_job_count,
            "slurm_pending_job_count": slurm_status_report.pending_job_count,
            "slurm_stale_job_count": slurm_status_report.stale_job_count,
            "slurm_active_run_state": slurm_status_report.active_run_state,
            "slurm_failure_recovery_status": (
                slurm_failure_recovery_report.overall_recovery_status
            ),
            "slurm_failure_recovery_rerunnable_job_count": (
                slurm_failure_recovery_report.rerunnable_job_count
            ),
            "slurm_failure_recovery_blocked_job_count": (
                slurm_failure_recovery_report.blocked_job_count
            ),
            "slurm_failure_recovery_partition_count": (
                slurm_failure_recovery_report.recovery_partition_count
            ),
        },
        summary_metrics=[
            ("variants", len(report.variant_runs)),
            ("execution mode", report.execution_mode),
            ("parallel workers", report.parallel_workers),
            ("stable clades", len(report.stable_clade_rows)),
            ("changed clades", len(report.changed_clade_rows)),
            (
                "slurm planned jobs",
                slurm_planning_report.job_count,
            ),
            (
                "slurm max memory MiB",
                slurm_planning_report.maximum_estimated_memory_mib,
            ),
            (
                "slurm array partitions",
                slurm_array_strategy_report.partition_count,
            ),
            (
                "slurm job evidence files",
                slurm_job_evidence_report.total_artifact_file_count,
            ),
            (
                "slurm storage MiB",
                slurm_storage_report.total_estimated_storage_mib,
            ),
            (
                "slurm largest storage variant",
                slurm_storage_report.largest_variant_id,
            ),
            (
                "slurm output explosion risk",
                slurm_output_explosion_report.overall_risk_status,
            ),
            (
                "slurm output explosion high-risk variants",
                slurm_output_explosion_report.high_risk_variant_count,
            ),
            (
                "slurm tree retention status",
                slurm_tree_retention_report.overall_policy_status,
            ),
            (
                "slurm tree-set files",
                slurm_tree_retention_report.tree_set_file_count,
            ),
            (
                "slurm merge ready",
                str(slurm_merge_report.merge_ready).lower(),
            ),
            (
                "slurm mergeable variants",
                slurm_merge_report.mergeable_variant_count,
            ),
            (
                "slurm fresh output jobs",
                slurm_output_freshness_report.fresh_job_count,
            ),
            (
                "slurm stale output jobs",
                slurm_output_freshness_report.stale_job_count,
            ),
            (
                "slurm completed jobs",
                slurm_status_report.completed_job_count,
            ),
            (
                "slurm stale jobs",
                slurm_status_report.stale_job_count,
            ),
            (
                "slurm recovery status",
                slurm_failure_recovery_report.overall_recovery_status,
            ),
            (
                "slurm rerunnable jobs",
                slurm_failure_recovery_report.rerunnable_job_count,
            ),
            (
                "reproducibility passed",
                str(getattr(reproducibility_report, "all_passed")).lower(),
            ),
            (
                "reproducibility checks",
                getattr(reproducibility_report, "check_count"),
            ),
            ("linked artifacts", report_manifest["linked_artifact_count"]),
        ],
        artifact_links=[
            *artifact_links,
            (
                "report-manifest",
                _relative_bundle_path(path, report_manifest_path),
                f"{report_manifest_path.stat().st_size} bytes",
            ),
        ],
    )

def _copy_output(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    return Path(shutil.copy2(source, destination))


def _write_tsv(path: Path, rows: list[list[str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join("\t".join(row) for row in rows) + "\n",
        encoding="utf-8",
    )
    return path


def _format_float(value: float) -> str:
    return format(value, ".12g")


def _format_optional_float(value: float | None) -> str:
    return "" if value is None else format(value, ".12g")


def _format_optional_bool(value: bool | None) -> str:
    if value is None:
        return ""
    return str(value).lower()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()
