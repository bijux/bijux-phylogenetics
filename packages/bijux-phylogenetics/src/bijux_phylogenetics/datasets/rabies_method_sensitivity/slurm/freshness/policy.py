from __future__ import annotations

import json

from .contracts import (
    RabiesMethodSensitivityOutputFreshnessCheckRow,
    RabiesMethodSensitivitySlurmOutputFreshnessRow,
)
from .interfaces import DatasetLike


def evaluate_freshness_checks(
    *,
    resolved_config: dict[str, object],
    dataset: DatasetLike,
    selected_variant_ids: tuple[str, ...],
    input_checksums: dict[str, str],
) -> tuple[
    list[RabiesMethodSensitivityOutputFreshnessCheckRow],
    bool,
    bool,
    dict[str, list[RabiesMethodSensitivityOutputFreshnessCheckRow]],
]:
    checks: list[RabiesMethodSensitivityOutputFreshnessCheckRow] = []

    def add_check(
        check_id: str,
        *,
        surface: str,
        scope: str,
        condition: bool,
        expected: object,
        observed: object,
        detail: str,
    ) -> None:
        checks.append(
            RabiesMethodSensitivityOutputFreshnessCheckRow(
                check_id=check_id,
                surface=surface,
                scope=scope,
                status="passed" if condition else "failed",
                expected="" if expected is None else str(expected),
                observed="" if observed is None else str(observed),
                detail=detail,
            )
        )

    recorded_input_checksums = {
        str(key): str(value)
        for key, value in dict(resolved_config.get("input_checksums", {})).items()
    }
    for filename, checksum in input_checksums.items():
        if filename not in recorded_input_checksums:
            continue
        add_check(
            f"input-checksum:{filename}",
            surface="input-checksum",
            scope="workflow",
            condition=recorded_input_checksums.get(filename) == checksum,
            expected=recorded_input_checksums.get(filename),
            observed=checksum,
            detail=(
                f"{filename} still matches the checksum recorded in "
                "workflow-config.resolved.json"
            ),
        )

    workflow_setting_pairs = (
        ("sequence_type", dataset.sequence_type),
        ("outgroup_taxa", list(dataset.outgroup_taxa)),
        ("iqtree_seed", dataset.iqtree_seed),
        ("iqtree_threads", dataset.iqtree_threads),
        ("bootstrap_replicates", dataset.bootstrap_replicates),
    )
    for setting_name, current_value in workflow_setting_pairs:
        recorded_value = resolved_config.get(setting_name)
        if recorded_value is None:
            continue
        add_check(
            f"workflow-setting:{setting_name}",
            surface="workflow-setting",
            scope="workflow",
            condition=recorded_value == current_value,
            expected=recorded_value,
            observed=current_value,
            detail=f"{setting_name} still matches the recorded workflow setting",
        )

    current_variants = {variant.variant_id: variant for variant in dataset.variants}
    if resolved_config.get("selected_variant_ids") is not None:
        add_check(
            "variant-selection:coverage",
            surface="variant-selection",
            scope="workflow",
            condition=all(
                variant_id in current_variants for variant_id in selected_variant_ids
            ),
            expected=sorted(selected_variant_ids),
            observed=sorted(
                variant_id
                for variant_id in selected_variant_ids
                if variant_id in current_variants
            ),
            detail="current packaged variants still cover every recorded selected variant id",
        )

    bundle_input_ok = all(
        row.status == "passed" and row.surface == "input-checksum"
        for row in checks
        if row.surface == "input-checksum"
    )
    bundle_workflow_settings_ok = all(
        row.status == "passed"
        for row in checks
        if row.surface in {"workflow-setting", "variant-selection"}
    )

    bundle_variant_rows = {
        str(row["variant_id"]): row for row in list(resolved_config.get("variants", []))
    }
    variant_check_rows_by_variant: dict[
        str, list[RabiesMethodSensitivityOutputFreshnessCheckRow]
    ] = {}
    for variant_id in selected_variant_ids:
        current_variant = current_variants.get(variant_id)
        recorded_variant = bundle_variant_rows.get(variant_id)
        normalized_recorded_variant = (
            None
            if recorded_variant is None
            else {
                "variant_id": str(recorded_variant["variant_id"]),
                "alignment_mode": str(recorded_variant["alignment_mode"]),
                "trimming_mode": str(recorded_variant["trimming_mode"]),
                "trim_gap_threshold": float(recorded_variant["trim_gap_threshold"]),
            }
        )
        if recorded_variant is None:
            continue
        if current_variant is None:
            row = RabiesMethodSensitivityOutputFreshnessCheckRow(
                check_id=f"variant-setting:{variant_id}",
                surface="variant-setting",
                scope=variant_id,
                status="failed",
                expected=json.dumps(normalized_recorded_variant, sort_keys=True),
                observed="missing",
                detail="current packaged workflow no longer exposes the recorded variant",
            )
            checks.append(row)
            variant_check_rows_by_variant.setdefault(variant_id, []).append(row)
            continue
        current_payload = {
            "variant_id": current_variant.variant_id,
            "alignment_mode": current_variant.alignment_mode,
            "trimming_mode": current_variant.trimming_mode,
            "trim_gap_threshold": current_variant.trim_gap_threshold,
        }
        row = RabiesMethodSensitivityOutputFreshnessCheckRow(
            check_id=f"variant-setting:{variant_id}",
            surface="variant-setting",
            scope=variant_id,
            status="passed"
            if normalized_recorded_variant == current_payload
            else "failed",
            expected=json.dumps(normalized_recorded_variant, sort_keys=True),
            observed=json.dumps(current_payload, sort_keys=True),
            detail="current packaged variant settings still match the recorded selection",
        )
        checks.append(row)
        variant_check_rows_by_variant.setdefault(variant_id, []).append(row)

    return (
        checks,
        bundle_input_ok,
        bundle_workflow_settings_ok,
        variant_check_rows_by_variant,
    )


def build_freshness_row(
    *,
    member_row: dict[str, str],
    bundle_input_ok: bool,
    bundle_workflow_settings_ok: bool,
    global_failed_reason_codes: tuple[str, ...],
    global_failed_reason_detail: str,
    variant_check_rows: list[RabiesMethodSensitivityOutputFreshnessCheckRow],
) -> RabiesMethodSensitivitySlurmOutputFreshnessRow:
    variant_settings_match = all(row.status == "passed" for row in variant_check_rows)
    stale_reason_codes = [
        *global_failed_reason_codes,
        *(row.check_id for row in variant_check_rows if row.status == "failed"),
    ]
    stale_reason_detail = "; ".join(
        part
        for part in [
            global_failed_reason_detail,
            "; ".join(
                row.detail for row in variant_check_rows if row.status == "failed"
            ),
        ]
        if part
    )
    return RabiesMethodSensitivitySlurmOutputFreshnessRow(
        partition_id=str(member_row["partition_id"]),
        array_index=int(member_row["array_index"]),
        variant_id=str(member_row["variant_id"]),
        freshness_status="fresh" if not stale_reason_codes else "stale",
        inputs_match=bundle_input_ok,
        workflow_settings_match=bundle_workflow_settings_ok,
        variant_settings_match=variant_settings_match,
        stale_reason_count=len(stale_reason_codes),
        stale_reason_codes=tuple(stale_reason_codes),
        stale_reason_detail=(
            "current inputs and output-affecting settings still match the bundle"
            if not stale_reason_codes
            else stale_reason_detail
        ),
    )
