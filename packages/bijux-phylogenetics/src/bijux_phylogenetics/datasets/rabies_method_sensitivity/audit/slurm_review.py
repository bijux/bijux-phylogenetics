from __future__ import annotations

from pathlib import Path

from .contracts import RabiesMethodSensitivityReproducibilityCheckRow
from .inventory import RabiesMethodSensitivityAuditSnapshot
from .review_context import RabiesMethodSensitivityAuditReviewContext


def record_rabies_method_sensitivity_slurm_checks(
    *,
    checks: list[RabiesMethodSensitivityReproducibilityCheckRow],
    snapshot: RabiesMethodSensitivityAuditSnapshot,
    review_context: RabiesMethodSensitivityAuditReviewContext,
) -> None:
    def add_check(
        check_id: str,
        *,
        surface: str,
        condition: bool,
        expected: object,
        observed: object,
        detail: str,
    ) -> None:
        checks.append(
            RabiesMethodSensitivityReproducibilityCheckRow(
                check_id=check_id,
                surface=surface,
                status="passed" if condition else "failed",
                expected="" if expected is None else str(expected),
                observed="" if observed is None else str(observed),
                detail=detail,
            )
        )

    slurm_partition_ids = review_context.slurm_partition_ids
    slurm_script_paths = review_context.slurm_script_paths
    member_partition_ids = review_context.member_partition_ids
    member_variant_ids = review_context.member_variant_ids
    config_variant_ids = review_context.config_variant_ids
    add_check(
        "slurm-arrays:partition-coverage",
        surface="slurm-arrays",
        condition=bool(slurm_partition_ids)
        and slurm_partition_ids == member_partition_ids,
        expected=sorted(slurm_partition_ids),
        observed=sorted(member_partition_ids),
        detail="array member rows cover the same partition ids as the partition summary",
    )
    add_check(
        "slurm-arrays:member-coverage",
        surface="slurm-arrays",
        condition=config_variant_ids == member_variant_ids,
        expected=config_variant_ids,
        observed=member_variant_ids,
        detail="array member rows cover the configured variant ids",
    )
    slurm_job_evidence_variant_ids = review_context.slurm_job_evidence_variant_ids
    add_check(
        "slurm-job-evidence:job-coverage",
        surface="slurm-job-evidence",
        condition=config_variant_ids == slurm_job_evidence_variant_ids,
        expected=config_variant_ids,
        observed=slurm_job_evidence_variant_ids,
        detail="job-evidence rows cover the configured variant ids",
    )
    add_check(
        "slurm-job-evidence:job-count",
        surface="slurm-job-evidence",
        condition=int(snapshot.slurm_job_evidence_summary["job_count"])
        == len(snapshot.slurm_job_evidence_rows),
        expected=snapshot.slurm_job_evidence_summary["job_count"],
        observed=len(snapshot.slurm_job_evidence_rows),
        detail="job-evidence summary job_count matches the written job-evidence rows",
    )
    job_evidence_artifact_file_count = len(
        [path for path in snapshot.slurm_job_evidence_root.rglob("*") if path.is_file()]
    )
    add_check(
        "slurm-job-evidence:artifact-file-count",
        surface="slurm-job-evidence",
        condition=int(snapshot.slurm_job_evidence_summary["total_artifact_file_count"])
        == job_evidence_artifact_file_count,
        expected=snapshot.slurm_job_evidence_summary["total_artifact_file_count"],
        observed=job_evidence_artifact_file_count,
        detail="job-evidence summary total_artifact_file_count matches the written evidence files",
    )
    slurm_storage_category_ids = review_context.slurm_storage_category_ids
    slurm_storage_variant_ids = review_context.slurm_storage_variant_ids
    add_check(
        "slurm-storage:category-coverage",
        surface="slurm-storage",
        condition=slurm_storage_category_ids
        == ["logs", "outputs", "posterior_samples", "reports", "trees"],
        expected=["logs", "outputs", "posterior_samples", "reports", "trees"],
        observed=slurm_storage_category_ids,
        detail="storage category rows cover the explicit retained-storage categories",
    )
    add_check(
        "slurm-storage:variant-coverage",
        surface="slurm-storage",
        condition=config_variant_ids == slurm_storage_variant_ids,
        expected=config_variant_ids,
        observed=slurm_storage_variant_ids,
        detail="storage variant rows cover the configured variant ids",
    )
    storage_total_file_count = sum(
        int(row["total_file_count"]) for row in snapshot.slurm_storage_category_rows
    )
    storage_total_byte_count = sum(
        int(row["total_byte_count"]) for row in snapshot.slurm_storage_category_rows
    )
    add_check(
        "slurm-storage:summary-totals",
        surface="slurm-storage",
        condition=int(snapshot.slurm_storage_summary["total_file_count"])
        == storage_total_file_count
        and int(snapshot.slurm_storage_summary["total_byte_count"])
        == storage_total_byte_count,
        expected=(
            int(snapshot.slurm_storage_summary["total_file_count"]),
            int(snapshot.slurm_storage_summary["total_byte_count"]),
        ),
        observed=(storage_total_file_count, storage_total_byte_count),
        detail="storage summary totals match the written category ledger",
    )
    largest_storage_variant_row = max(
        snapshot.slurm_storage_variant_rows,
        key=lambda row: (int(row["total_byte_count"]), str(row["variant_id"])),
    )
    add_check(
        "slurm-storage:largest-variant",
        surface="slurm-storage",
        condition=str(snapshot.slurm_storage_summary["largest_variant_id"])
        == str(largest_storage_variant_row["variant_id"])
        and int(snapshot.slurm_storage_summary["largest_variant_total_byte_count"])
        == int(largest_storage_variant_row["total_byte_count"]),
        expected=(
            str(snapshot.slurm_storage_summary["largest_variant_id"]),
            int(snapshot.slurm_storage_summary["largest_variant_total_byte_count"]),
        ),
        observed=(
            str(largest_storage_variant_row["variant_id"]),
            int(largest_storage_variant_row["total_byte_count"]),
        ),
        detail="storage summary largest-variant fields match the written per-variant storage ledger",
    )
    slurm_output_explosion_variant_ids = (
        review_context.slurm_output_explosion_variant_ids
    )
    add_check(
        "slurm-output-explosion:variant-coverage",
        surface="slurm-output-explosion",
        condition=config_variant_ids == slurm_output_explosion_variant_ids,
        expected=config_variant_ids,
        observed=slurm_output_explosion_variant_ids,
        detail="output-explosion rows cover the configured variant ids",
    )
    add_check(
        "slurm-output-explosion:check-count",
        surface="slurm-output-explosion",
        condition=int(snapshot.slurm_output_explosion_summary["check_count"])
        == len(snapshot.slurm_output_explosion_check_rows),
        expected=snapshot.slurm_output_explosion_summary["check_count"],
        observed=len(snapshot.slurm_output_explosion_check_rows),
        detail="output-explosion summary check_count matches the written check rows",
    )
    observed_output_explosion_variant_total = (
        sum(
            1
            for row in snapshot.slurm_output_explosion_variant_rows
            if str(row["risk_status"]) == "low"
        )
        + sum(
            1
            for row in snapshot.slurm_output_explosion_variant_rows
            if str(row["risk_status"]) == "warning"
        )
        + sum(
            1
            for row in snapshot.slurm_output_explosion_variant_rows
            if str(row["risk_status"]) == "high"
        )
    )
    add_check(
        "slurm-output-explosion:variant-counts",
        surface="slurm-output-explosion",
        condition=int(snapshot.slurm_output_explosion_summary["variant_count"])
        == observed_output_explosion_variant_total,
        expected=snapshot.slurm_output_explosion_summary["variant_count"],
        observed=observed_output_explosion_variant_total,
        detail="output-explosion summary variant_count matches the risk-classified variant rows",
    )
    observed_overall_output_explosion_status = "high"
    if (
        int(snapshot.slurm_output_explosion_summary["failed_check_count"]) == 0
        and int(snapshot.slurm_output_explosion_summary["high_risk_variant_count"]) == 0
        and int(snapshot.slurm_output_explosion_summary["global_issue_count"]) == 0
        and int(snapshot.slurm_output_explosion_summary["warning_variant_count"]) == 0
    ):
        observed_overall_output_explosion_status = "low"
    elif (
        int(snapshot.slurm_output_explosion_summary["failed_check_count"]) == 0
        and int(snapshot.slurm_output_explosion_summary["high_risk_variant_count"]) == 0
    ):
        observed_overall_output_explosion_status = "warning"
    add_check(
        "slurm-output-explosion:overall-risk",
        surface="slurm-output-explosion",
        condition=str(snapshot.slurm_output_explosion_summary["overall_risk_status"])
        == observed_overall_output_explosion_status,
        expected=snapshot.slurm_output_explosion_summary["overall_risk_status"],
        observed=observed_overall_output_explosion_status,
        detail="output-explosion summary overall_risk_status matches the written risk counts",
    )
    slurm_tree_retention_variant_ids = review_context.slurm_tree_retention_variant_ids
    add_check(
        "slurm-tree-retention:variant-coverage",
        surface="slurm-tree-retention",
        condition=config_variant_ids == slurm_tree_retention_variant_ids,
        expected=config_variant_ids,
        observed=slurm_tree_retention_variant_ids,
        detail="tree-retention file rows cover the configured variant ids",
    )
    add_check(
        "slurm-tree-retention:check-count",
        surface="slurm-tree-retention",
        condition=int(snapshot.slurm_tree_retention_summary["check_count"])
        == len(snapshot.slurm_tree_retention_check_rows),
        expected=snapshot.slurm_tree_retention_summary["check_count"],
        observed=len(snapshot.slurm_tree_retention_check_rows),
        detail="tree-retention summary check_count matches the written check rows",
    )
    observed_tree_set_file_count = sum(
        1
        for row in snapshot.slurm_tree_retention_file_rows
        if int(row["tree_count"]) > 1
    )
    add_check(
        "slurm-tree-retention:tree-set-count",
        surface="slurm-tree-retention",
        condition=int(snapshot.slurm_tree_retention_summary["tree_set_file_count"])
        == observed_tree_set_file_count,
        expected=snapshot.slurm_tree_retention_summary["tree_set_file_count"],
        observed=observed_tree_set_file_count,
        detail="tree-retention summary tree_set_file_count matches the inspected file rows",
    )
    observed_tree_retention_status = "required"
    if (
        int(snapshot.slurm_tree_retention_summary["failed_check_count"]) == 0
        and int(snapshot.slurm_tree_retention_summary["thinning_required_file_count"])
        == 0
        and int(
            snapshot.slurm_tree_retention_summary["compression_required_file_count"]
        )
        == 0
        and int(
            snapshot.slurm_tree_retention_summary["thinning_recommended_file_count"]
        )
        == 0
        and int(
            snapshot.slurm_tree_retention_summary["compression_recommended_file_count"]
        )
        == 0
    ):
        observed_tree_retention_status = "no_action"
    elif int(snapshot.slurm_tree_retention_summary["failed_check_count"]) == 0 and (
        int(snapshot.slurm_tree_retention_summary["thinning_required_file_count"]) == 0
        and int(
            snapshot.slurm_tree_retention_summary["compression_required_file_count"]
        )
        == 0
    ):
        observed_tree_retention_status = "recommended"
    add_check(
        "slurm-tree-retention:overall-status",
        surface="slurm-tree-retention",
        condition=str(snapshot.slurm_tree_retention_summary["overall_policy_status"])
        == observed_tree_retention_status,
        expected=snapshot.slurm_tree_retention_summary["overall_policy_status"],
        observed=observed_tree_retention_status,
        detail="tree-retention summary overall_policy_status matches the written policy counts",
    )
    slurm_merge_variant_ids = review_context.slurm_merge_variant_ids
    add_check(
        "slurm-merge:job-coverage",
        surface="slurm-merge",
        condition=config_variant_ids == slurm_merge_variant_ids,
        expected=config_variant_ids,
        observed=slurm_merge_variant_ids,
        detail="merge-variant rows cover the configured variant ids",
    )
    add_check(
        "slurm-merge:check-count",
        surface="slurm-merge",
        condition=int(snapshot.slurm_merge_summary["check_count"])
        == len(snapshot.slurm_merge_check_rows),
        expected=snapshot.slurm_merge_summary["check_count"],
        observed=len(snapshot.slurm_merge_check_rows),
        detail="merge summary check_count matches the written merge checks",
    )
    merged_variant_count = sum(
        1
        for row in snapshot.slurm_merge_variant_rows
        if str(row["included_in_merge"]) == "true"
    )
    add_check(
        "slurm-merge:merged-variant-count",
        surface="slurm-merge",
        condition=int(snapshot.slurm_merge_summary["merged_variant_count"])
        == merged_variant_count,
        expected=snapshot.slurm_merge_summary["merged_variant_count"],
        observed=merged_variant_count,
        detail="merge summary merged_variant_count matches the written merge-variant rows",
    )
    add_check(
        "slurm-merge:merge-ready",
        surface="slurm-merge",
        condition=bool(snapshot.slurm_merge_summary["merge_ready"])
        == (
            int(snapshot.slurm_merge_summary["failed_check_count"]) == 0
            and merged_variant_count == len(config_variant_ids)
        ),
        expected=int(snapshot.slurm_merge_summary["failed_check_count"]) == 0
        and merged_variant_count == len(config_variant_ids),
        observed=snapshot.slurm_merge_summary["merge_ready"],
        detail="merge summary merge_ready matches the merge-check and merged-variant totals",
    )
    slurm_output_freshness_variant_ids = (
        review_context.slurm_output_freshness_variant_ids
    )
    add_check(
        "slurm-freshness:job-coverage",
        surface="slurm-freshness",
        condition=config_variant_ids == slurm_output_freshness_variant_ids,
        expected=config_variant_ids,
        observed=slurm_output_freshness_variant_ids,
        detail="output-freshness rows cover the configured variant ids",
    )
    freshness_failed_check_count = sum(
        1
        for row in snapshot.slurm_output_freshness_check_rows
        if str(row["status"]) == "failed"
    )
    freshness_stale_job_count = sum(
        1
        for row in snapshot.slurm_output_freshness_rows
        if str(row["freshness_status"]) == "stale"
    )
    add_check(
        "slurm-freshness:check-count",
        surface="slurm-freshness",
        condition=int(snapshot.slurm_output_freshness_summary["check_count"])
        == len(snapshot.slurm_output_freshness_check_rows),
        expected=snapshot.slurm_output_freshness_summary["check_count"],
        observed=len(snapshot.slurm_output_freshness_check_rows),
        detail="output-freshness summary check_count matches the written check rows",
    )
    add_check(
        "slurm-freshness:job-count",
        surface="slurm-freshness",
        condition=(
            int(snapshot.slurm_output_freshness_summary["fresh_job_count"])
            + int(snapshot.slurm_output_freshness_summary["stale_job_count"])
        )
        == len(snapshot.slurm_output_freshness_rows),
        expected=len(snapshot.slurm_output_freshness_rows),
        observed=(
            int(snapshot.slurm_output_freshness_summary["fresh_job_count"])
            + int(snapshot.slurm_output_freshness_summary["stale_job_count"])
        ),
        detail="output-freshness summary job counts match the written per-job ledger",
    )
    add_check(
        "slurm-freshness:all-outputs-fresh",
        surface="slurm-freshness",
        condition=bool(snapshot.slurm_output_freshness_summary["all_outputs_fresh"])
        == (freshness_failed_check_count == 0 and freshness_stale_job_count == 0),
        expected=freshness_failed_check_count == 0 and freshness_stale_job_count == 0,
        observed=snapshot.slurm_output_freshness_summary["all_outputs_fresh"],
        detail="output-freshness summary all_outputs_fresh matches the written checks and per-job statuses",
    )
    for script_path_text in sorted(slurm_script_paths):
        script_path = snapshot.bundle_root / script_path_text
        add_check(
            f"slurm-arrays:script:{Path(script_path_text).name}",
            surface="slurm-arrays",
            condition=script_path.is_file(),
            expected="script file exists",
            observed="present" if script_path.is_file() else "missing",
            detail="array partition table references an existing sbatch script",
        )
    slurm_job_status_variant_ids = review_context.slurm_job_status_variant_ids
    slurm_partition_status_ids = review_context.slurm_partition_status_ids
    slurm_output_freshness_rows_by_variant = (
        review_context.slurm_output_freshness_rows_by_variant
    )
    slurm_job_evidence_rows_by_variant = (
        review_context.slurm_job_evidence_rows_by_variant
    )
    slurm_merge_rows_by_variant = review_context.slurm_merge_rows_by_variant
    slurm_job_status_rows_by_variant = review_context.slurm_job_status_rows_by_variant
    add_check(
        "slurm-status:job-coverage",
        surface="slurm-status",
        condition=config_variant_ids == slurm_job_status_variant_ids,
        expected=config_variant_ids,
        observed=slurm_job_status_variant_ids,
        detail="job-status rows cover the configured variant ids",
    )
    add_check(
        "slurm-status:partition-coverage",
        surface="slurm-status",
        condition=sorted(slurm_partition_ids) == slurm_partition_status_ids,
        expected=sorted(slurm_partition_ids),
        observed=slurm_partition_status_ids,
        detail="partition-status rows cover the same partition ids as the array partition table",
    )
    add_check(
        "slurm-status:workflow-job-count",
        surface="slurm-status",
        condition=int(snapshot.slurm_workflow_status["job_count"])
        == len(snapshot.slurm_job_status_rows),
        expected=snapshot.slurm_workflow_status["job_count"],
        observed=len(snapshot.slurm_job_status_rows),
        detail="workflow status job_count matches the number of job-status rows",
    )
    add_check(
        "slurm-status:workflow-partition-count",
        surface="slurm-status",
        condition=int(snapshot.slurm_workflow_status["partition_count"])
        == len(snapshot.slurm_partition_status_rows),
        expected=snapshot.slurm_workflow_status["partition_count"],
        observed=len(snapshot.slurm_partition_status_rows),
        detail="workflow status partition_count matches the number of partition-status rows",
    )
    slurm_failure_recovery_variant_ids = (
        review_context.slurm_failure_recovery_variant_ids
    )
    add_check(
        "slurm-failure-recovery:variant-coverage",
        surface="slurm-failure-recovery",
        condition=config_variant_ids == slurm_failure_recovery_variant_ids,
        expected=config_variant_ids,
        observed=slurm_failure_recovery_variant_ids,
        detail="failure-recovery job rows cover the configured variant ids",
    )
    add_check(
        "slurm-failure-recovery:partition-coverage",
        surface="slurm-failure-recovery",
        condition=int(snapshot.slurm_failure_recovery_summary["partition_count"])
        == len(snapshot.slurm_failure_recovery_partition_rows),
        expected=snapshot.slurm_failure_recovery_summary["partition_count"],
        observed=len(snapshot.slurm_failure_recovery_partition_rows),
        detail="failure-recovery summary partition_count matches the written partition rows",
    )
    observed_rerunnable_job_count = sum(
        1
        for row in snapshot.slurm_failure_recovery_job_rows
        if str(row["rerunnable"]) == "true"
    )
    add_check(
        "slurm-failure-recovery:rerunnable-count",
        surface="slurm-failure-recovery",
        condition=int(snapshot.slurm_failure_recovery_summary["rerunnable_job_count"])
        == observed_rerunnable_job_count,
        expected=snapshot.slurm_failure_recovery_summary["rerunnable_job_count"],
        observed=observed_rerunnable_job_count,
        detail="failure-recovery summary rerunnable_job_count matches the written job rows",
    )
    observed_failure_recovery_status = "clean"
    if observed_rerunnable_job_count > 0:
        observed_failure_recovery_status = "recovery_needed"
    elif int(snapshot.slurm_failure_recovery_summary["blocked_job_count"]) > 0:
        observed_failure_recovery_status = "workflow_active"
    add_check(
        "slurm-failure-recovery:overall-status",
        surface="slurm-failure-recovery",
        condition=str(
            snapshot.slurm_failure_recovery_summary["overall_recovery_status"]
        )
        == observed_failure_recovery_status,
        expected=snapshot.slurm_failure_recovery_summary["overall_recovery_status"],
        observed=observed_failure_recovery_status,
        detail="failure-recovery summary overall_recovery_status matches the written rerun decisions",
    )
    for variant_id in config_variant_ids:
        freshness_row = slurm_output_freshness_rows_by_variant.get(variant_id)
        job_evidence_row = slurm_job_evidence_rows_by_variant.get(variant_id)
        merge_row = slurm_merge_rows_by_variant.get(variant_id)
        job_status_row = slurm_job_status_rows_by_variant.get(variant_id)
        add_check(
            f"slurm-freshness:status-link:{variant_id}",
            surface="slurm-freshness",
            condition=(
                freshness_row is not None
                and job_status_row is not None
                and str(freshness_row["freshness_status"])
                == str(job_status_row["output_freshness_status"])
            ),
            expected=None
            if freshness_row is None
            else freshness_row["freshness_status"],
            observed=None
            if job_status_row is None
            else job_status_row["output_freshness_status"],
            detail="job-status rows expose the same output-freshness status as the freshness ledger",
        )
        add_check(
            f"slurm-job-evidence:status-link:{variant_id}",
            surface="slurm-job-evidence",
            condition=(
                job_evidence_row is not None
                and job_status_row is not None
                and str(job_evidence_row["status"])
                == str(job_status_row["task_status"])
            ),
            expected=None if job_status_row is None else job_status_row["task_status"],
            observed=None if job_evidence_row is None else job_evidence_row["status"],
            detail="job-evidence rows expose the same terminal task status recorded in the job-status ledger",
        )
        add_check(
            f"slurm-merge:status-link:{variant_id}",
            surface="slurm-merge",
            condition=(
                merge_row is not None
                and job_status_row is not None
                and str(merge_row["job_status"]) == str(job_status_row["status"])
            ),
            expected=None if job_status_row is None else job_status_row["status"],
            observed=None if merge_row is None else merge_row["job_status"],
            detail="merge-variant rows expose the same job status as the job-status ledger",
        )
        add_check(
            f"slurm-merge:freshness-link:{variant_id}",
            surface="slurm-merge",
            condition=(
                merge_row is not None
                and freshness_row is not None
                and str(merge_row["output_freshness_status"])
                == str(freshness_row["freshness_status"])
            ),
            expected=None
            if freshness_row is None
            else freshness_row["freshness_status"],
            observed=None
            if merge_row is None
            else merge_row["output_freshness_status"],
            detail="merge-variant rows expose the same freshness status as the freshness ledger",
        )
        if job_evidence_row is None:
            continue
        evidence_json_path = snapshot.bundle_root / str(
            job_evidence_row["evidence_json_path"]
        )
        evidence_html_path = snapshot.bundle_root / str(
            job_evidence_row["evidence_html_path"]
        )
        add_check(
            f"slurm-job-evidence:json:{variant_id}",
            surface="slurm-job-evidence",
            condition=evidence_json_path.is_file(),
            expected="file exists",
            observed="present" if evidence_json_path.is_file() else "missing",
            detail="job-evidence index references an existing JSON evidence package",
        )
        add_check(
            f"slurm-job-evidence:html:{variant_id}",
            surface="slurm-job-evidence",
            condition=evidence_html_path.is_file(),
            expected="file exists",
            observed="present" if evidence_html_path.is_file() else "missing",
            detail="job-evidence index references an existing HTML evidence package",
        )
