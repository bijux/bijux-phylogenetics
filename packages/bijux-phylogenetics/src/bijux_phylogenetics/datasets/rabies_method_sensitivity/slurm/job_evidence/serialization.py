from __future__ import annotations

import json
import os
from pathlib import Path

from .contracts import (
    RabiesMethodSensitivitySlurmJobEvidenceReport,
    RabiesMethodSensitivitySlurmJobEvidenceRow,
)


def write_rabies_method_sensitivity_slurm_job_evidence_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmJobEvidenceReport,
) -> Path:
    """Write one workflow-wide index over per-job evidence packages."""
    rows = [
        [
            "partition_id",
            "array_index",
            "variant_id",
            "status",
            "selected_model",
            "total_runtime_seconds",
            "output_file_count",
            "output_byte_count",
            "artifact_file_count",
            "serious_conflict_count",
            "rooted_engine_rf_distance",
            "warning_count",
            "evidence_json_path",
            "evidence_html_path",
        ]
    ]
    for row in report.jobs:
        rows.append(
            [
                row.partition_id,
                str(row.array_index),
                row.variant_id,
                row.status,
                row.selected_model,
                _format_float(row.total_runtime_seconds),
                str(row.output_file_count),
                str(row.output_byte_count),
                str(row.artifact_file_count),
                str(row.serious_conflict_count),
                str(row.rooted_engine_rf_distance),
                str(row.warning_count),
                row.evidence_json_path,
                row.evidence_html_path,
            ]
        )
    return _write_tsv(path, rows)


def write_rabies_method_sensitivity_slurm_job_evidence_summary_json(
    path: Path,
    report: RabiesMethodSensitivitySlurmJobEvidenceReport,
) -> Path:
    """Write the structured workflow-wide job-evidence summary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "dataset_id": report.dataset_id,
        "workflow_prefix": report.workflow_prefix,
        "execution_mode": report.execution_mode,
        "parallel_workers": report.parallel_workers,
        "bundle_root": ".",
        "evidence_root": _relative_bundle_path(
            report.bundle_root, report.evidence_root
        ),
        "job_count": report.job_count,
        "completed_job_count": report.completed_job_count,
        "failed_job_count": report.failed_job_count,
        "total_runtime_seconds": report.total_runtime_seconds,
        "total_output_file_count": report.total_output_file_count,
        "total_output_byte_count": report.total_output_byte_count,
        "total_artifact_file_count": report.total_artifact_file_count,
        "jobs": [_serialize_row(row) for row in report.jobs],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def _serialize_row(
    row: RabiesMethodSensitivitySlurmJobEvidenceRow,
) -> dict[str, object]:
    return {
        "partition_id": row.partition_id,
        "array_index": row.array_index,
        "variant_id": row.variant_id,
        "label": row.label,
        "execution_mode": row.execution_mode,
        "status": row.status,
        "script_path": row.script_path,
        "output_root": row.output_root,
        "evidence_directory": row.evidence_directory,
        "evidence_json_path": row.evidence_json_path,
        "evidence_html_path": row.evidence_html_path,
        "task_log_copy_path": row.task_log_copy_path,
        "alignment_manifest_path": row.alignment_manifest_path,
        "trimming_manifest_path": row.trimming_manifest_path,
        "inference_manifest_path": row.inference_manifest_path,
        "model_selection_manifest_path": row.model_selection_manifest_path,
        "iqtree_support_manifest_path": row.iqtree_support_manifest_path,
        "fasttree_manifest_path": row.fasttree_manifest_path,
        "alignment_mode": row.alignment_mode,
        "trimming_mode": row.trimming_mode,
        "trim_gap_threshold": row.trim_gap_threshold,
        "selected_model": row.selected_model,
        "alignment_length": row.alignment_length,
        "trimmed_alignment_length": row.trimmed_alignment_length,
        "total_runtime_seconds": row.total_runtime_seconds,
        "alignment_runtime_seconds": row.alignment_runtime_seconds,
        "trimming_runtime_seconds": row.trimming_runtime_seconds,
        "inference_runtime_seconds": row.inference_runtime_seconds,
        "model_selection_runtime_seconds": row.model_selection_runtime_seconds,
        "iqtree_support_runtime_seconds": row.iqtree_support_runtime_seconds,
        "fasttree_runtime_seconds": row.fasttree_runtime_seconds,
        "rooted_engine_rf_distance": row.rooted_engine_rf_distance,
        "rooted_engine_normalized_rf": row.rooted_engine_normalized_rf,
        "rooted_same_taxa_different_rooting": row.rooted_same_taxa_different_rooting,
        "serious_conflict_count": row.serious_conflict_count,
        "stable_clade_count": row.stable_clade_count,
        "unstable_clade_count": row.unstable_clade_count,
        "engine_specific_clade_count": row.engine_specific_clade_count,
        "warning_count": row.warning_count,
        "output_file_count": row.output_file_count,
        "output_byte_count": row.output_byte_count,
        "artifact_file_count": row.artifact_file_count,
    }


def _relative_bundle_path(source: Path, target: Path) -> str:
    base = source if source.is_dir() else source.parent
    return os.path.relpath(target, start=base).replace(os.sep, "/")


def _write_tsv(path: Path, rows: list[list[str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join("\t".join(row) for row in rows) + "\n",
        encoding="utf-8",
    )
    return path


def _format_float(value: float) -> str:
    return format(value, ".12g")
