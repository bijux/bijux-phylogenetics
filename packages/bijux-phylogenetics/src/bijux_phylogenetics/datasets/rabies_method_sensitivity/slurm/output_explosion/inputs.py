from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .contracts import RabiesMethodSensitivitySlurmOutputExplosionCheckRow
from .shared import (
    _CATEGORY_IDS,
    _CONFIG_FILENAME,
    _MEBIBYTE,
    _SLURM_JOB_PLAN_FILENAME,
    _SLURM_STORAGE_CATEGORIES_FILENAME,
    _SLURM_STORAGE_SUMMARY_FILENAME,
    _SLURM_STORAGE_VARIANTS_FILENAME,
    _load_json,
    _read_tsv_rows,
)


@dataclass(frozen=True, slots=True)
class OutputExplosionInputs:
    bundle_root: Path
    config: dict[str, object]
    configured_variant_ids: list[str]
    job_plan_rows: list[dict[str, str]]
    storage_category_rows: list[dict[str, str]]
    storage_variant_rows: list[dict[str, str]]
    storage_summary: dict[str, object]
    job_plan_by_variant: dict[str, dict[str, str]]
    storage_variant_by_variant: dict[str, dict[str, str]]
    storage_category_by_id: dict[str, dict[str, str]]
    total_estimated_output_mib: int
    bootstrap_replicates: int
    checks: tuple[RabiesMethodSensitivitySlurmOutputExplosionCheckRow, ...]


def load_output_explosion_inputs(bundle_root: Path) -> OutputExplosionInputs:
    bundle_root = bundle_root.resolve()
    config = _load_json(bundle_root / _CONFIG_FILENAME)
    job_plan_rows = _read_tsv_rows(bundle_root / _SLURM_JOB_PLAN_FILENAME)
    storage_category_rows = _read_tsv_rows(
        bundle_root / _SLURM_STORAGE_CATEGORIES_FILENAME
    )
    storage_variant_rows = _read_tsv_rows(
        bundle_root / _SLURM_STORAGE_VARIANTS_FILENAME
    )
    storage_summary = _load_json(bundle_root / _SLURM_STORAGE_SUMMARY_FILENAME)

    configured_variant_ids = sorted(
        str(row["variant_id"]) for row in list(config.get("variants", []))
    )
    job_plan_by_variant = {str(row["variant_id"]): row for row in job_plan_rows}
    storage_variant_by_variant = {
        str(row["variant_id"]): row for row in storage_variant_rows
    }
    storage_category_by_id = {
        str(row["category_id"]): row for row in storage_category_rows
    }
    total_estimated_output_mib = sum(
        int(row["estimated_output_mib"]) for row in job_plan_rows
    )

    checks: list[RabiesMethodSensitivitySlurmOutputExplosionCheckRow] = []
    _add_check(
        checks,
        "job-plan:variant-coverage",
        surface="job-plan",
        condition=sorted(job_plan_by_variant) == configured_variant_ids,
        expected=configured_variant_ids,
        observed=sorted(job_plan_by_variant),
        detail="job-plan rows cover the configured variant ids",
    )
    _add_check(
        checks,
        "storage:variant-coverage",
        surface="storage",
        condition=sorted(storage_variant_by_variant) == configured_variant_ids,
        expected=configured_variant_ids,
        observed=sorted(storage_variant_by_variant),
        detail="storage variant rows cover the configured variant ids",
    )
    _add_check(
        checks,
        "storage:category-coverage",
        surface="storage",
        condition=sorted(storage_category_by_id) == _CATEGORY_IDS,
        expected=_CATEGORY_IDS,
        observed=sorted(storage_category_by_id),
        detail="storage category rows cover the explicit retained-output categories",
    )
    _add_check(
        checks,
        "job-plan:total-output-mib",
        surface="job-plan",
        condition=int(storage_summary["variant_count"]) == len(configured_variant_ids),
        expected=len(configured_variant_ids),
        observed=storage_summary["variant_count"],
        detail="storage summary variant_count matches the configured variant count",
    )
    _add_check(
        checks,
        "storage:total-storage-mib",
        surface="storage",
        condition=int(storage_summary["total_estimated_storage_mib"])
        >= int(storage_summary["total_byte_count"]) // _MEBIBYTE,
        expected="storage summary rounds total retained bytes up into MiB",
        observed=storage_summary["total_estimated_storage_mib"],
        detail="storage summary total_estimated_storage_mib is a rounded-up retained-byte estimate",
    )
    _add_check(
        checks,
        "storage:largest-variant",
        surface="storage",
        condition=str(storage_summary["largest_variant_id"])
        == max(
            storage_variant_rows,
            key=lambda row: (int(row["total_byte_count"]), str(row["variant_id"])),
        )["variant_id"],
        expected=storage_summary["largest_variant_id"],
        observed=max(
            storage_variant_rows,
            key=lambda row: (int(row["total_byte_count"]), str(row["variant_id"])),
        )["variant_id"],
        detail="storage summary largest_variant_id matches the written storage-variant rows",
    )
    return OutputExplosionInputs(
        bundle_root=bundle_root,
        config=config,
        configured_variant_ids=configured_variant_ids,
        job_plan_rows=job_plan_rows,
        storage_category_rows=storage_category_rows,
        storage_variant_rows=storage_variant_rows,
        storage_summary=storage_summary,
        job_plan_by_variant=job_plan_by_variant,
        storage_variant_by_variant=storage_variant_by_variant,
        storage_category_by_id=storage_category_by_id,
        total_estimated_output_mib=total_estimated_output_mib,
        bootstrap_replicates=int(config["bootstrap_replicates"]),
        checks=tuple(checks),
    )


def _add_check(
    checks: list[RabiesMethodSensitivitySlurmOutputExplosionCheckRow],
    check_id: str,
    *,
    surface: str,
    condition: bool,
    expected: object,
    observed: object,
    detail: str,
) -> None:
    checks.append(
        RabiesMethodSensitivitySlurmOutputExplosionCheckRow(
            check_id=check_id,
            surface=surface,
            status="passed" if condition else "failed",
            expected="" if expected is None else str(expected),
            observed="" if observed is None else str(observed),
            detail=detail,
        )
    )
