from __future__ import annotations

import json
from pathlib import Path

from ..audit import RabiesMethodSensitivityReproducibilityAuditReport
from ..models import RabiesMethodSensitivityPanelWorkflowReport
from ..slurm import (
    RabiesMethodSensitivitySlurmArrayStrategyReport,
    RabiesMethodSensitivitySlurmFailureRecoveryReport,
    RabiesMethodSensitivitySlurmJobEvidenceReport,
    RabiesMethodSensitivitySlurmMergeReport,
    RabiesMethodSensitivitySlurmOutputExplosionReport,
    RabiesMethodSensitivitySlurmOutputFreshnessReport,
    RabiesMethodSensitivitySlurmPlanningReport,
    RabiesMethodSensitivitySlurmStatusReport,
    RabiesMethodSensitivitySlurmStorageReport,
    RabiesMethodSensitivitySlurmTreeRetentionReport,
)
from .shared import _format_float, _relative_bundle_path


def _build_report_sections(
    *,
    report: RabiesMethodSensitivityPanelWorkflowReport,
    bundle_paths: dict[str, Path],
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
) -> list[tuple[str, str]]:
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
    return [
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
        ("conclusions", "\n".join(conclusion_lines).strip()),
        (
            "reproducibility-audit",
            "\n".join(
                [
                    f"all passed: {str(reproducibility_report.all_passed).lower()}",
                    f"top-level checks: {reproducibility_report.check_count}",
                    (
                        "failed top-level checks: "
                        f"{reproducibility_report.failed_check_count}"
                    ),
                    (f"failed variants: {reproducibility_report.failed_variant_count}"),
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
                    (f"partition scripts: {slurm_array_strategy_report.script_count}"),
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
                    (f"failed merge checks: {slurm_merge_report.failed_check_count}"),
                    (f"merged stable clades: {slurm_merge_report.stable_clade_count}"),
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
                    (f"fresh jobs: {slurm_output_freshness_report.fresh_job_count}"),
                    (f"stale jobs: {slurm_output_freshness_report.stale_job_count}"),
                    (f"freshness checks: {slurm_output_freshness_report.check_count}"),
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
                    f"slurm estimation assumptions: {bundle_paths['slurm_assumptions'].name}",
                    f"slurm planning summary: {bundle_paths['slurm_summary'].name}",
                    f"slurm array partitions: {bundle_paths['slurm_array_partitions'].name}",
                    f"slurm array members: {bundle_paths['slurm_array_members'].name}",
                    f"slurm array strategy: {bundle_paths['slurm_array_strategy'].name}",
                    f"slurm job evidence index: {bundle_paths['slurm_job_evidence_index'].name}",
                    f"slurm job evidence summary: {bundle_paths['slurm_job_evidence_summary'].name}",
                    f"slurm storage categories: {bundle_paths['slurm_storage_categories'].name}",
                    f"slurm storage variants: {bundle_paths['slurm_storage_variants'].name}",
                    f"slurm storage summary: {bundle_paths['slurm_storage_summary'].name}",
                    f"slurm storage report: {bundle_paths['slurm_storage_report'].name}",
                    f"slurm output explosion checks: {bundle_paths['slurm_output_explosion_checks'].name}",
                    f"slurm output explosion variants: {bundle_paths['slurm_output_explosion_variants'].name}",
                    f"slurm output explosion summary: {bundle_paths['slurm_output_explosion_summary'].name}",
                    f"slurm output explosion report: {bundle_paths['slurm_output_explosion_report'].name}",
                    f"slurm tree retention checks: {bundle_paths['slurm_tree_retention_checks'].name}",
                    f"slurm tree retention files: {bundle_paths['slurm_tree_retention_files'].name}",
                    f"slurm tree retention summary: {bundle_paths['slurm_tree_retention_summary'].name}",
                    f"slurm tree retention report: {bundle_paths['slurm_tree_retention_report'].name}",
                    f"slurm merge checks: {bundle_paths['slurm_merge_checks'].name}",
                    f"slurm merge variants: {bundle_paths['slurm_merge_variants'].name}",
                    f"slurm merge summary: {bundle_paths['slurm_merge_summary'].name}",
                    f"slurm merge report: {bundle_paths['slurm_merge_report'].name}",
                    f"slurm output freshness: {bundle_paths['slurm_output_freshness'].name}",
                    f"slurm output freshness checks: {bundle_paths['slurm_output_freshness_checks'].name}",
                    f"slurm output freshness summary: {bundle_paths['slurm_output_freshness_summary'].name}",
                    f"slurm job status: {bundle_paths['slurm_job_status'].name}",
                    f"slurm partition status: {bundle_paths['slurm_partition_status'].name}",
                    f"slurm workflow status: {bundle_paths['slurm_workflow_status'].name}",
                    f"slurm failure recovery jobs: {bundle_paths['slurm_failure_recovery_jobs'].name}",
                    f"slurm failure recovery partitions: {bundle_paths['slurm_failure_recovery_partitions'].name}",
                    f"slurm failure recovery summary: {bundle_paths['slurm_failure_recovery_summary'].name}",
                    f"slurm failure recovery report: {bundle_paths['slurm_failure_recovery_report'].name}",
                    f"reproducibility checks: {bundle_paths['reproducibility_checks'].name}",
                    f"reproducibility variant audit: {bundle_paths['reproducibility_variant_audit'].name}",
                    f"reproducibility audit: {bundle_paths['reproducibility_audit'].name}",
                ]
            ),
        ),
    ]


def _build_report_artifact_links(
    path: Path,
    bundle_paths: dict[str, Path],
) -> list[tuple[str, str, str]]:
    return [
        (
            key.replace("_", "-"),
            _relative_bundle_path(path, value),
            f"{value.stat().st_size} bytes",
        )
        for key, value in bundle_paths.items()
        if value.is_file()
    ]


def _build_report_embedded_json(
    *,
    path: Path,
    report: RabiesMethodSensitivityPanelWorkflowReport,
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
) -> dict[str, object]:
    return {
        "dataset_id": report.dataset.dataset_id,
        "variant_count": len(report.variant_runs),
        "parallel_workers": report.parallel_workers,
        "execution_mode": report.execution_mode,
        "stable_clade_count": len(report.stable_clade_rows),
        "changed_clade_count": len(report.changed_clade_rows),
        "report_manifest_path": _relative_bundle_path(path, report_manifest_path),
        "reproducibility_passed": reproducibility_report.all_passed,
        "reproducibility_check_count": reproducibility_report.check_count,
        "reproducibility_failed_check_count": reproducibility_report.failed_check_count,
        "reproducibility_failed_variant_count": reproducibility_report.failed_variant_count,
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
        "slurm_array_partition_count": slurm_array_strategy_report.partition_count,
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
        "slurm_storage_output_byte_count": slurm_storage_report.output_byte_count,
        "slurm_storage_log_byte_count": slurm_storage_report.log_byte_count,
        "slurm_storage_tree_byte_count": slurm_storage_report.tree_byte_count,
        "slurm_storage_posterior_sample_byte_count": (
            slurm_storage_report.posterior_sample_byte_count
        ),
        "slurm_storage_report_byte_count": slurm_storage_report.report_byte_count,
        "slurm_storage_largest_variant_id": slurm_storage_report.largest_variant_id,
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
        "slurm_tree_retention_status": slurm_tree_retention_report.overall_policy_status,
        "slurm_tree_set_file_count": slurm_tree_retention_report.tree_set_file_count,
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
        "slurm_mergeable_variant_count": slurm_merge_report.mergeable_variant_count,
        "slurm_merge_failed_check_count": slurm_merge_report.failed_check_count,
        "slurm_output_freshness_check_count": (
            slurm_output_freshness_report.check_count
        ),
        "slurm_output_freshness_failed_check_count": (
            slurm_output_freshness_report.failed_check_count
        ),
        "slurm_fresh_output_job_count": slurm_output_freshness_report.fresh_job_count,
        "slurm_stale_output_job_count": slurm_output_freshness_report.stale_job_count,
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
    }


def _build_report_summary_metrics(
    *,
    report: RabiesMethodSensitivityPanelWorkflowReport,
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
) -> list[tuple[str, object]]:
    report_manifest = json.loads(report_manifest_path.read_text(encoding="utf-8"))
    return [
        ("variants", len(report.variant_runs)),
        ("execution mode", report.execution_mode),
        ("parallel workers", report.parallel_workers),
        ("stable clades", len(report.stable_clade_rows)),
        ("changed clades", len(report.changed_clade_rows)),
        ("slurm planned jobs", slurm_planning_report.job_count),
        ("slurm max memory MiB", slurm_planning_report.maximum_estimated_memory_mib),
        ("slurm array partitions", slurm_array_strategy_report.partition_count),
        (
            "slurm job evidence files",
            slurm_job_evidence_report.total_artifact_file_count,
        ),
        ("slurm storage MiB", slurm_storage_report.total_estimated_storage_mib),
        ("slurm largest storage variant", slurm_storage_report.largest_variant_id),
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
        ("slurm tree-set files", slurm_tree_retention_report.tree_set_file_count),
        ("slurm merge ready", str(slurm_merge_report.merge_ready).lower()),
        ("slurm mergeable variants", slurm_merge_report.mergeable_variant_count),
        ("slurm fresh output jobs", slurm_output_freshness_report.fresh_job_count),
        ("slurm stale output jobs", slurm_output_freshness_report.stale_job_count),
        ("slurm completed jobs", slurm_status_report.completed_job_count),
        ("slurm stale jobs", slurm_status_report.stale_job_count),
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
            str(reproducibility_report.all_passed).lower(),
        ),
        ("reproducibility checks", reproducibility_report.check_count),
        ("linked artifacts", report_manifest["linked_artifact_count"]),
    ]
