from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .contracts import RabiesMethodSensitivitySlurmMergeCheckRow
from .shared import (
    _CHANGED_CLADES_FILENAME,
    _CONCLUSION_SUMMARY_FILENAME,
    _CONFIG_FILENAME,
    _PREPROCESSING_COMPARISONS_FILENAME,
    _SLURM_JOB_EVIDENCE_FILENAME,
    _SLURM_JOB_STATUS_FILENAME,
    _SLURM_OUTPUT_FRESHNESS_FILENAME,
    _STABLE_CLADES_FILENAME,
    _VARIANT_SUMMARY_FILENAME,
    _WORKFLOW_SUMMARY_FILENAME,
    _load_json,
    _read_tsv_rows,
)


@dataclass(frozen=True, slots=True)
class SlurmMergeInputs:
    bundle_root: Path
    config: dict[str, object]
    workflow_summary: dict[str, str]
    configured_variant_ids: list[str]
    variant_summary_rows: list[dict[str, str]]
    preprocessing_rows: list[dict[str, str]]
    stable_clade_rows: list[dict[str, str]]
    changed_clade_rows: list[dict[str, str]]
    conclusion_rows: list[dict[str, str]]
    variant_summary_by_variant: dict[str, dict[str, str]]
    job_status_by_variant: dict[str, dict[str, str]]
    freshness_by_variant: dict[str, dict[str, str]]
    job_evidence_by_variant: dict[str, dict[str, str]]
    checks: tuple[RabiesMethodSensitivitySlurmMergeCheckRow, ...]


def load_slurm_merge_inputs(bundle_root: Path) -> SlurmMergeInputs:
    bundle_root = bundle_root.resolve()
    config = _load_json(bundle_root / _CONFIG_FILENAME)
    workflow_summary_rows = _read_tsv_rows(bundle_root / _WORKFLOW_SUMMARY_FILENAME)
    variant_summary_rows = _read_tsv_rows(bundle_root / _VARIANT_SUMMARY_FILENAME)
    preprocessing_rows = _read_tsv_rows(
        bundle_root / _PREPROCESSING_COMPARISONS_FILENAME
    )
    stable_clade_rows = _read_tsv_rows(bundle_root / _STABLE_CLADES_FILENAME)
    changed_clade_rows = _read_tsv_rows(bundle_root / _CHANGED_CLADES_FILENAME)
    conclusion_rows = _read_tsv_rows(bundle_root / _CONCLUSION_SUMMARY_FILENAME)
    slurm_job_status_rows = _read_tsv_rows(bundle_root / _SLURM_JOB_STATUS_FILENAME)
    slurm_output_freshness_rows = _read_tsv_rows(
        bundle_root / _SLURM_OUTPUT_FRESHNESS_FILENAME
    )
    slurm_job_evidence_rows = _read_tsv_rows(bundle_root / _SLURM_JOB_EVIDENCE_FILENAME)

    configured_variant_ids = [
        str(row["variant_id"]) for row in list(config.get("variants", []))
    ]
    variant_summary_by_variant = {
        str(row["variant_id"]): row for row in variant_summary_rows
    }
    job_status_by_variant = {
        str(row["variant_id"]): row for row in slurm_job_status_rows
    }
    freshness_by_variant = {
        str(row["variant_id"]): row for row in slurm_output_freshness_rows
    }
    job_evidence_by_variant = {
        str(row["variant_id"]): row for row in slurm_job_evidence_rows
    }
    summary_variant_ids = sorted(variant_summary_by_variant)
    job_status_variant_ids = sorted(job_status_by_variant)
    freshness_variant_ids = sorted(freshness_by_variant)
    job_evidence_variant_ids = sorted(job_evidence_by_variant)
    expected_variant_ids = sorted(configured_variant_ids)

    checks = (
        _build_check(
            "variant-coverage:variant-summary",
            surface="variant-coverage",
            condition=expected_variant_ids == summary_variant_ids,
            expected=expected_variant_ids,
            observed=summary_variant_ids,
            detail="variant summary rows cover the configured variant ids",
        ),
        _build_check(
            "variant-coverage:job-status",
            surface="variant-coverage",
            condition=expected_variant_ids == job_status_variant_ids,
            expected=expected_variant_ids,
            observed=job_status_variant_ids,
            detail="job-status rows cover the configured variant ids",
        ),
        _build_check(
            "variant-coverage:freshness",
            surface="variant-coverage",
            condition=expected_variant_ids == freshness_variant_ids,
            expected=expected_variant_ids,
            observed=freshness_variant_ids,
            detail="freshness rows cover the configured variant ids",
        ),
        _build_check(
            "variant-coverage:job-evidence",
            surface="variant-coverage",
            condition=expected_variant_ids == job_evidence_variant_ids,
            expected=expected_variant_ids,
            observed=job_evidence_variant_ids,
            detail="job-evidence rows cover the configured variant ids",
        ),
    )
    return SlurmMergeInputs(
        bundle_root=bundle_root,
        config=config,
        workflow_summary=workflow_summary_rows[0],
        configured_variant_ids=configured_variant_ids,
        variant_summary_rows=variant_summary_rows,
        preprocessing_rows=preprocessing_rows,
        stable_clade_rows=stable_clade_rows,
        changed_clade_rows=changed_clade_rows,
        conclusion_rows=conclusion_rows,
        variant_summary_by_variant=variant_summary_by_variant,
        job_status_by_variant=job_status_by_variant,
        freshness_by_variant=freshness_by_variant,
        job_evidence_by_variant=job_evidence_by_variant,
        checks=checks,
    )


def _build_check(
    check_id: str,
    *,
    surface: str,
    condition: bool,
    expected: object,
    observed: object,
    detail: str,
) -> RabiesMethodSensitivitySlurmMergeCheckRow:
    return RabiesMethodSensitivitySlurmMergeCheckRow(
        check_id=check_id,
        surface=surface,
        status="passed" if condition else "failed",
        expected="" if expected is None else str(expected),
        observed="" if observed is None else str(observed),
        detail=detail,
    )
