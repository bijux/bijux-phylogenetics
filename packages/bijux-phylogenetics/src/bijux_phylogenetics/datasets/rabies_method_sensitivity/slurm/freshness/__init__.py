from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import hashlib
import json
from pathlib import Path
from typing import Protocol

__all__ = [
    "RabiesMethodSensitivityOutputFreshnessCheckRow",
    "RabiesMethodSensitivitySlurmOutputFreshnessRow",
    "RabiesMethodSensitivitySlurmOutputFreshnessReport",
    "build_rabies_method_sensitivity_slurm_output_freshness_report",
    "write_rabies_method_sensitivity_slurm_output_freshness_checks_table",
    "write_rabies_method_sensitivity_slurm_output_freshness_json",
    "write_rabies_method_sensitivity_slurm_output_freshness_table",
]

_CONFIG_FILENAME = "workflow-config.resolved.json"
_SLURM_ARRAY_MEMBERS_FILENAME = "slurm-array-members.tsv"


class _VariantLike(Protocol):
    variant_id: str
    alignment_mode: str
    trimming_mode: str
    trim_gap_threshold: float


class _DatasetLike(Protocol):
    dataset_id: str
    workflow_prefix: str
    sequence_type: str
    outgroup_taxa: tuple[str, ...]
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    sequences_path: Path
    metadata_path: Path
    variants: tuple[_VariantLike, ...]


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityOutputFreshnessCheckRow:
    """One freshness check over bundle inputs or output-affecting settings."""

    check_id: str
    surface: str
    scope: str
    status: str
    expected: str
    observed: str
    detail: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmOutputFreshnessRow:
    """One per-job output-freshness classification against current inputs/settings."""

    partition_id: str
    array_index: int
    variant_id: str
    freshness_status: str
    inputs_match: bool
    workflow_settings_match: bool
    variant_settings_match: bool
    stale_reason_count: int
    stale_reason_codes: tuple[str, ...]
    stale_reason_detail: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmOutputFreshnessReport:
    """One reviewer-facing freshness summary for rabies batch outputs."""

    dataset_id: str
    workflow_prefix: str
    bundle_root: Path
    all_outputs_fresh: bool
    selected_variant_ids: tuple[str, ...]
    check_count: int
    failed_check_count: int
    fresh_job_count: int
    stale_job_count: int
    checks: tuple[RabiesMethodSensitivityOutputFreshnessCheckRow, ...]
    jobs: tuple[RabiesMethodSensitivitySlurmOutputFreshnessRow, ...]


def build_rabies_method_sensitivity_slurm_output_freshness_report(
    bundle_root: Path,
    *,
    dataset: _DatasetLike | None = None,
) -> RabiesMethodSensitivitySlurmOutputFreshnessReport:
    """Detect whether bundle outputs still match the current packaged workflow state."""
    bundle_root = bundle_root.resolve()
    resolved_config = _load_json(bundle_root / _CONFIG_FILENAME)
    if dataset is None:
        from ...config import (
            load_rabies_method_sensitivity_panel_dataset,
        )

        dataset = load_rabies_method_sensitivity_panel_dataset()
    member_rows = _read_tsv_rows(bundle_root / _SLURM_ARRAY_MEMBERS_FILENAME)
    selected_variant_ids = tuple(
        str(value) for value in list(resolved_config.get("selected_variant_ids", []))
    )
    if not selected_variant_ids:
        selected_variant_ids = tuple(str(row["variant_id"]) for row in member_rows)
    bundle_variant_rows = {
        str(row["variant_id"]): row for row in list(resolved_config.get("variants", []))
    }
    current_variants = {variant.variant_id: variant for variant in dataset.variants}
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

    input_checksums = {
        "sequences.fasta": _sha256(dataset.sequences_path),
        "metadata.csv": _sha256(dataset.metadata_path),
    }
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

    global_failed_reason_codes = tuple(
        row.check_id
        for row in checks
        if row.scope == "workflow" and row.status == "failed"
    )
    global_failed_reason_detail = "; ".join(
        row.detail
        for row in checks
        if row.scope == "workflow" and row.status == "failed"
    )
    freshness_rows = tuple(
        _build_freshness_row(
            member_row=member_row,
            bundle_input_ok=bundle_input_ok,
            bundle_workflow_settings_ok=bundle_workflow_settings_ok,
            global_failed_reason_codes=global_failed_reason_codes,
            global_failed_reason_detail=global_failed_reason_detail,
            variant_check_rows=variant_check_rows_by_variant.get(
                str(member_row["variant_id"]),
                [],
            ),
        )
        for member_row in member_rows
    )
    failed_check_count = sum(1 for row in checks if row.status == "failed")
    fresh_job_count = sum(1 for row in freshness_rows if row.freshness_status == "fresh")
    stale_job_count = sum(1 for row in freshness_rows if row.freshness_status == "stale")

    return RabiesMethodSensitivitySlurmOutputFreshnessReport(
        dataset_id=dataset.dataset_id,
        workflow_prefix=dataset.workflow_prefix,
        bundle_root=bundle_root,
        all_outputs_fresh=failed_check_count == 0 and stale_job_count == 0,
        selected_variant_ids=selected_variant_ids,
        check_count=len(checks),
        failed_check_count=failed_check_count,
        fresh_job_count=fresh_job_count,
        stale_job_count=stale_job_count,
        checks=tuple(checks),
        jobs=freshness_rows,
    )


def write_rabies_method_sensitivity_slurm_output_freshness_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmOutputFreshnessReport,
) -> Path:
    """Write one per-job output-freshness ledger."""
    return _write_tsv(
        path,
        fieldnames=(
            "partition_id",
            "array_index",
            "variant_id",
            "freshness_status",
            "inputs_match",
            "workflow_settings_match",
            "variant_settings_match",
            "stale_reason_count",
            "stale_reason_codes",
            "stale_reason_detail",
        ),
        rows=[
            {
                "partition_id": row.partition_id,
                "array_index": row.array_index,
                "variant_id": row.variant_id,
                "freshness_status": row.freshness_status,
                "inputs_match": str(row.inputs_match).lower(),
                "workflow_settings_match": str(row.workflow_settings_match).lower(),
                "variant_settings_match": str(row.variant_settings_match).lower(),
                "stale_reason_count": row.stale_reason_count,
                "stale_reason_codes": ",".join(row.stale_reason_codes),
                "stale_reason_detail": row.stale_reason_detail,
            }
            for row in report.jobs
        ],
    )


def write_rabies_method_sensitivity_slurm_output_freshness_checks_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmOutputFreshnessReport,
) -> Path:
    """Write one check-level ledger for output freshness."""
    return _write_tsv(
        path,
        fieldnames=(
            "check_id",
            "surface",
            "scope",
            "status",
            "expected",
            "observed",
            "detail",
        ),
        rows=[
            {
                "check_id": row.check_id,
                "surface": row.surface,
                "scope": row.scope,
                "status": row.status,
                "expected": row.expected,
                "observed": row.observed,
                "detail": row.detail,
            }
            for row in report.checks
        ],
    )


def write_rabies_method_sensitivity_slurm_output_freshness_json(
    path: Path,
    report: RabiesMethodSensitivitySlurmOutputFreshnessReport,
) -> Path:
    """Write the structured output-freshness summary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(report)
    payload["bundle_root"] = "."
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def _build_freshness_row(
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
            "; ".join(row.detail for row in variant_check_rows if row.status == "failed"),
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


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_tsv(
    path: Path,
    *,
    fieldnames: tuple[str, ...],
    rows: list[dict[str, object]],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: "" if value is None else value
                    for key, value in row.items()
                }
            )
    return path
