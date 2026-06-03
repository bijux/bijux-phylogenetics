from __future__ import annotations

from .contracts import RabiesMethodSensitivitySlurmOutputExplosionReport
from .inputs import _add_check, load_output_explosion_inputs
from .shared import (
    _BOOTSTRAP_HIGH_REPLICATES,
    _BOOTSTRAP_WARNING_REPLICATES,
    _POSTERIOR_HIGH_BYTES,
    _POSTERIOR_WARNING_BYTES,
    _REPORT_HIGH_BYTES,
    _REPORT_WARNING_BYTES,
    _TOTAL_OUTPUT_HIGH_MIB,
    _TOTAL_OUTPUT_WARNING_MIB,
    _TOTAL_STORAGE_HIGH_MIB,
    _TOTAL_STORAGE_WARNING_MIB,
    _TREE_HIGH_BYTES,
    _TREE_WARNING_BYTES,
    _TREE_WARNING_FILES,
)
from .variant_risk import build_output_explosion_variant_rows, escalate_risk, raise_risk


def build_rabies_method_sensitivity_slurm_output_explosion_report(
    bundle_root,
) -> RabiesMethodSensitivitySlurmOutputExplosionReport:
    """Assess whether retained workflow outputs are trending toward an explosion risk."""
    loaded_inputs = load_output_explosion_inputs(bundle_root)
    bundle_root = loaded_inputs.bundle_root
    config = loaded_inputs.config
    configured_variant_ids = loaded_inputs.configured_variant_ids
    storage_summary = loaded_inputs.storage_summary
    checks = list(loaded_inputs.checks)
    storage_category_by_id = loaded_inputs.storage_category_by_id
    total_estimated_output_mib = loaded_inputs.total_estimated_output_mib

    variant_rows = build_output_explosion_variant_rows(
        configured_variant_ids=configured_variant_ids,
        job_plan_by_variant=loaded_inputs.job_plan_by_variant,
        storage_variant_by_variant=loaded_inputs.storage_variant_by_variant,
        total_estimated_output_mib=total_estimated_output_mib,
    )

    storage_trees_row = storage_category_by_id["trees"]
    storage_posterior_row = storage_category_by_id["posterior_samples"]
    storage_reports_row = storage_category_by_id["reports"]
    global_issues: list[str] = []
    global_severity = "low"
    global_severity = raise_risk(
        global_severity,
        global_issues,
        value=total_estimated_output_mib,
        warning_threshold=_TOTAL_OUTPUT_WARNING_MIB,
        high_threshold=_TOTAL_OUTPUT_HIGH_MIB,
        warning_detail="total retained output MiB suggests scaling pressure",
        high_detail="total retained output MiB suggests a severe retained-output burden",
    )
    global_severity = raise_risk(
        global_severity,
        global_issues,
        value=int(storage_summary["total_estimated_storage_mib"]),
        warning_threshold=_TOTAL_STORAGE_WARNING_MIB,
        high_threshold=_TOTAL_STORAGE_HIGH_MIB,
        warning_detail="total retained storage MiB suggests scaling pressure",
        high_detail="total retained storage MiB suggests severe retained-storage pressure",
    )
    global_severity = raise_risk(
        global_severity,
        global_issues,
        value=int(storage_trees_row["total_byte_count"]),
        warning_threshold=_TREE_WARNING_BYTES,
        high_threshold=_TREE_HIGH_BYTES,
        warning_detail="total tree artifact bytes are large enough to pressure retained storage",
        high_detail="total tree artifact bytes are dominating retained storage",
    )
    global_severity = raise_risk(
        global_severity,
        global_issues,
        value=int(storage_posterior_row["total_byte_count"]),
        warning_threshold=_POSTERIOR_WARNING_BYTES,
        high_threshold=_POSTERIOR_HIGH_BYTES,
        warning_detail="posterior sample bytes are becoming large enough to pressure retained storage",
        high_detail="posterior sample bytes are dominating retained storage",
    )
    global_severity = raise_risk(
        global_severity,
        global_issues,
        value=int(storage_reports_row["total_byte_count"]),
        warning_threshold=_REPORT_WARNING_BYTES,
        high_threshold=_REPORT_HIGH_BYTES,
        warning_detail="review artifact bytes are becoming large enough to pressure retained storage",
        high_detail="review artifact bytes are dominating retained storage",
    )
    bootstrap_replicates = loaded_inputs.bootstrap_replicates
    total_tree_file_count = int(storage_trees_row["total_file_count"])
    if (
        bootstrap_replicates >= _BOOTSTRAP_HIGH_REPLICATES
        and total_tree_file_count >= _TREE_WARNING_FILES
    ):
        global_severity = escalate_risk(global_severity, "high")
        global_issues.append(
            "bootstrap replicates and retained tree-file counts together suggest a severe tree-output explosion risk"
        )
    elif (
        bootstrap_replicates >= _BOOTSTRAP_WARNING_REPLICATES
        and total_tree_file_count >= _TREE_WARNING_FILES
    ):
        global_severity = escalate_risk(global_severity, "warning")
        global_issues.append(
            "bootstrap replicates and retained tree-file counts together suggest tree-output growth pressure"
        )

    low_risk_variant_count = sum(1 for row in variant_rows if row.risk_status == "low")
    warning_variant_count = sum(
        1 for row in variant_rows if row.risk_status == "warning"
    )
    high_risk_variant_count = sum(
        1 for row in variant_rows if row.risk_status == "high"
    )
    failed_check_count = sum(1 for row in checks if row.status == "failed")
    overall_risk_status = _derive_overall_risk_status(
        failed_check_count=failed_check_count,
        warning_variant_count=warning_variant_count,
        high_risk_variant_count=high_risk_variant_count,
        global_severity=global_severity,
    )

    largest_variant = max(
        variant_rows,
        key=lambda row: (row.output_share, row.variant_id),
    )
    _add_check(
        checks,
        "risk-summary:variant-counts",
        surface="risk-summary",
        condition=low_risk_variant_count
        + warning_variant_count
        + high_risk_variant_count
        == len(variant_rows),
        expected=len(variant_rows),
        observed=low_risk_variant_count
        + warning_variant_count
        + high_risk_variant_count,
        detail="risk-status counts cover every configured variant exactly once",
    )
    _add_check(
        checks,
        "risk-summary:overall-status",
        surface="risk-summary",
        condition=overall_risk_status
        == _derive_overall_risk_status(
            failed_check_count=failed_check_count,
            warning_variant_count=warning_variant_count,
            high_risk_variant_count=high_risk_variant_count,
            global_severity=global_severity,
        ),
        expected=_derive_overall_risk_status(
            failed_check_count=failed_check_count,
            warning_variant_count=warning_variant_count,
            high_risk_variant_count=high_risk_variant_count,
            global_severity=global_severity,
        ),
        observed=overall_risk_status,
        detail="overall explosion-risk status matches the global and per-variant risk counts",
    )

    return RabiesMethodSensitivitySlurmOutputExplosionReport(
        dataset_id=str(config["dataset_id"]),
        workflow_prefix=str(config["workflow_prefix"]),
        bundle_root=bundle_root,
        bootstrap_replicates=bootstrap_replicates,
        overall_risk_status=overall_risk_status,
        variant_count=len(configured_variant_ids),
        check_count=len(checks),
        failed_check_count=failed_check_count,
        low_risk_variant_count=low_risk_variant_count,
        warning_variant_count=warning_variant_count,
        high_risk_variant_count=high_risk_variant_count,
        global_issue_count=len(global_issues),
        total_estimated_output_mib=total_estimated_output_mib,
        total_estimated_storage_mib=int(storage_summary["total_estimated_storage_mib"]),
        total_tree_byte_count=int(storage_trees_row["total_byte_count"]),
        total_tree_file_count=int(storage_trees_row["total_file_count"]),
        total_posterior_sample_byte_count=int(
            storage_posterior_row["total_byte_count"]
        ),
        total_posterior_sample_file_count=int(
            storage_posterior_row["total_file_count"]
        ),
        total_report_byte_count=int(storage_reports_row["total_byte_count"]),
        largest_variant_id=largest_variant.variant_id,
        largest_variant_output_share=largest_variant.output_share,
        global_issues=tuple(global_issues),
        checks=tuple(checks),
        variants=tuple(variant_rows),
    )


def _derive_overall_risk_status(
    *,
    failed_check_count: int,
    warning_variant_count: int,
    high_risk_variant_count: int,
    global_severity: str,
) -> str:
    if (
        failed_check_count > 0
        or high_risk_variant_count > 0
        or global_severity == "high"
    ):
        return "high"
    if warning_variant_count > 0 or global_severity == "warning":
        return "warning"
    return "low"
