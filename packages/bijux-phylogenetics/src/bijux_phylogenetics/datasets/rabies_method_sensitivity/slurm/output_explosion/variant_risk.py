from __future__ import annotations

from .contracts import RabiesMethodSensitivitySlurmOutputExplosionVariantRow
from .shared import (
    _DOMINANT_VARIANT_OUTPUT_MIB,
    _DOMINANT_VARIANT_OUTPUT_SHARE,
    _OUTPUT_HIGH_MIB,
    _OUTPUT_WARNING_MIB,
    _POSTERIOR_HIGH_BYTES,
    _POSTERIOR_HIGH_FILES,
    _POSTERIOR_WARNING_BYTES,
    _POSTERIOR_WARNING_FILES,
    _REPORT_HIGH_BYTES,
    _REPORT_WARNING_BYTES,
    _STORAGE_HIGH_MIB,
    _STORAGE_WARNING_MIB,
    _TREE_HIGH_BYTES,
    _TREE_HIGH_FILES,
    _TREE_WARNING_BYTES,
    _TREE_WARNING_FILES,
)


def build_output_explosion_variant_rows(
    *,
    configured_variant_ids: list[str],
    job_plan_by_variant: dict[str, dict[str, str]],
    storage_variant_by_variant: dict[str, dict[str, str]],
    total_estimated_output_mib: int,
) -> list[RabiesMethodSensitivitySlurmOutputExplosionVariantRow]:
    variant_rows: list[RabiesMethodSensitivitySlurmOutputExplosionVariantRow] = []
    total_output_denom = max(1, total_estimated_output_mib)
    for variant_id in configured_variant_ids:
        plan_row = job_plan_by_variant.get(variant_id)
        storage_row = storage_variant_by_variant.get(variant_id)
        if plan_row is None or storage_row is None:
            issues = ["planning or storage row is missing"]
            variant_rows.append(
                RabiesMethodSensitivitySlurmOutputExplosionVariantRow(
                    variant_id=variant_id,
                    risk_status="high",
                    estimated_output_mib=0
                    if plan_row is None
                    else int(plan_row["estimated_output_mib"]),
                    estimated_storage_mib=0
                    if storage_row is None
                    else int(storage_row["estimated_storage_mib"]),
                    tree_file_count=0
                    if storage_row is None
                    else int(storage_row["tree_file_count"]),
                    tree_byte_count=0
                    if storage_row is None
                    else int(storage_row["tree_byte_count"]),
                    posterior_sample_file_count=0
                    if storage_row is None
                    else int(storage_row["posterior_sample_file_count"]),
                    posterior_sample_byte_count=0
                    if storage_row is None
                    else int(storage_row["posterior_sample_byte_count"]),
                    report_byte_count=0
                    if storage_row is None
                    else int(storage_row["report_byte_count"]),
                    output_share=0.0,
                    issue_count=1,
                    issues=tuple(issues),
                )
            )
            continue

        estimated_output_mib = int(plan_row["estimated_output_mib"])
        estimated_storage_mib = int(storage_row["estimated_storage_mib"])
        tree_file_count = int(storage_row["tree_file_count"])
        tree_byte_count = int(storage_row["tree_byte_count"])
        posterior_sample_file_count = int(storage_row["posterior_sample_file_count"])
        posterior_sample_byte_count = int(storage_row["posterior_sample_byte_count"])
        report_byte_count = int(storage_row["report_byte_count"])
        output_share = round(estimated_output_mib / total_output_denom, 4)
        issues: list[str] = []
        severity = "low"

        severity = raise_risk(
            severity,
            issues,
            value=estimated_output_mib,
            warning_threshold=_OUTPUT_WARNING_MIB,
            high_threshold=_OUTPUT_HIGH_MIB,
            warning_detail="estimated retained output MiB is large for one variant",
            high_detail="estimated retained output MiB is very large for one variant",
        )
        severity = raise_risk(
            severity,
            issues,
            value=estimated_storage_mib,
            warning_threshold=_STORAGE_WARNING_MIB,
            high_threshold=_STORAGE_HIGH_MIB,
            warning_detail="estimated retained storage MiB is large for one variant",
            high_detail="estimated retained storage MiB is very large for one variant",
        )
        severity = raise_risk(
            severity,
            issues,
            value=tree_byte_count,
            warning_threshold=_TREE_WARNING_BYTES,
            high_threshold=_TREE_HIGH_BYTES,
            warning_detail="tree artifact bytes are large enough to pressure retained storage",
            high_detail="tree artifact bytes are large enough to dominate retained storage",
        )
        severity = raise_risk(
            severity,
            issues,
            value=tree_file_count,
            warning_threshold=_TREE_WARNING_FILES,
            high_threshold=_TREE_HIGH_FILES,
            warning_detail="tree artifact file counts suggest tree-set growth pressure",
            high_detail="tree artifact file counts suggest a large retained tree-set burden",
        )
        severity = raise_risk(
            severity,
            issues,
            value=posterior_sample_byte_count,
            warning_threshold=_POSTERIOR_WARNING_BYTES,
            high_threshold=_POSTERIOR_HIGH_BYTES,
            warning_detail="posterior sample bytes are large enough to pressure retained storage",
            high_detail="posterior sample bytes are large enough to dominate retained storage",
        )
        severity = raise_risk(
            severity,
            issues,
            value=posterior_sample_file_count,
            warning_threshold=_POSTERIOR_WARNING_FILES,
            high_threshold=_POSTERIOR_HIGH_FILES,
            warning_detail="posterior sample file counts suggest chain or tree-set growth pressure",
            high_detail="posterior sample file counts suggest a very large retained chain or tree-set burden",
        )
        severity = raise_risk(
            severity,
            issues,
            value=report_byte_count,
            warning_threshold=_REPORT_WARNING_BYTES,
            high_threshold=_REPORT_HIGH_BYTES,
            warning_detail="review artifact bytes are becoming large for one variant",
            high_detail="review artifact bytes are dominating retained storage for one variant",
        )
        if (
            output_share >= _DOMINANT_VARIANT_OUTPUT_SHARE
            and estimated_output_mib >= _DOMINANT_VARIANT_OUTPUT_MIB
        ):
            severity = escalate_risk(severity, "warning")
            issues.append("one variant dominates most retained output MiB")
        variant_rows.append(
            RabiesMethodSensitivitySlurmOutputExplosionVariantRow(
                variant_id=variant_id,
                risk_status=severity,
                estimated_output_mib=estimated_output_mib,
                estimated_storage_mib=estimated_storage_mib,
                tree_file_count=tree_file_count,
                tree_byte_count=tree_byte_count,
                posterior_sample_file_count=posterior_sample_file_count,
                posterior_sample_byte_count=posterior_sample_byte_count,
                report_byte_count=report_byte_count,
                output_share=output_share,
                issue_count=len(issues),
                issues=tuple(issues),
            )
        )
    return variant_rows


def raise_risk(
    current_severity: str,
    issues: list[str],
    *,
    value: int,
    warning_threshold: int,
    high_threshold: int,
    warning_detail: str,
    high_detail: str,
) -> str:
    if value >= high_threshold:
        issues.append(high_detail)
        return escalate_risk(current_severity, "high")
    if value >= warning_threshold:
        issues.append(warning_detail)
        return escalate_risk(current_severity, "warning")
    return current_severity


def escalate_risk(current_severity: str, next_severity: str) -> str:
    order = {"low": 0, "warning": 1, "high": 2}
    return (
        next_severity
        if order[next_severity] > order[current_severity]
        else current_severity
    )
