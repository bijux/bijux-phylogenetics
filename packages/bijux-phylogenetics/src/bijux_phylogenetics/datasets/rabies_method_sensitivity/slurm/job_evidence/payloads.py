from __future__ import annotations

import hashlib
import os
from pathlib import Path
import shutil

from .contracts import (
    RabiesMethodSensitivitySlurmJobEvidenceRow,
    TaskRecordLike,
    VariantRunLike,
)


def _build_job_payload(
    *,
    bundle_root: Path,
    row: RabiesMethodSensitivitySlurmJobEvidenceRow,
    task_record: TaskRecordLike,
    variant_run: VariantRunLike,
    member: object,
    dataset_id: str,
    workflow_prefix: str,
    execution_record_path: Path,
    workflow_manifest_path: Path,
    output_inventory: list[dict[str, object]],
    warnings: list[str],
) -> dict[str, object]:
    trimming_summary = variant_run.trimming_workflow.trimming_summary
    return {
        "dataset_id": dataset_id,
        "workflow_prefix": workflow_prefix,
        "partition_id": row.partition_id,
        "array_index": row.array_index,
        "variant_id": row.variant_id,
        "label": row.label,
        "execution_mode": row.execution_mode,
        "status": row.status,
        "script_path": row.script_path,
        "bundle_output_root": row.output_root,
        "task_log_path": row.task_log_copy_path,
        "workflow_manifest_path": _relative_bundle_path(
            bundle_root,
            workflow_manifest_path,
        ),
        "execution_record_path": _relative_bundle_path(
            bundle_root,
            execution_record_path,
        ),
        "error_code": task_record.error_code,
        "error_message": task_record.error_message,
        "variant": {
            "alignment_mode": row.alignment_mode,
            "trimming_mode": row.trimming_mode,
            "trim_gap_threshold": row.trim_gap_threshold,
            "selected_model": row.selected_model,
        },
        "alignment": {
            "alignment_length": row.alignment_length,
            "trimmed_alignment_length": row.trimmed_alignment_length,
            "retained_site_fraction": (
                None
                if trimming_summary is None
                else trimming_summary.retained_site_fraction
            ),
            "removed_site_fraction": (
                None
                if trimming_summary is None
                else trimming_summary.removed_site_fraction
            ),
        },
        "commands": {
            "alignment": variant_run.alignment_workflow.run.command,
            "trimming": variant_run.trimming_workflow.run.command,
            **variant_run.inference_comparison.commands,
        },
        "engine_versions": {
            "mafft": variant_run.alignment_workflow.run.version.text,
            "trimal": variant_run.trimming_workflow.run.version.text,
            **variant_run.inference_comparison.engine_versions,
        },
        "runtimes": {
            "alignment_seconds": row.alignment_runtime_seconds,
            "trimming_seconds": row.trimming_runtime_seconds,
            "inference_seconds": row.inference_runtime_seconds,
            "model_selection_seconds": row.model_selection_runtime_seconds,
            "iqtree_support_seconds": row.iqtree_support_runtime_seconds,
            "fasttree_seconds": row.fasttree_runtime_seconds,
            "total_seconds": row.total_runtime_seconds,
        },
        "rooting": {
            "fasttree": {
                "requested_taxa": list(variant_run.fasttree_rooting.requested_taxa),
                "matched_taxa": list(variant_run.fasttree_rooting.matched_taxa),
                "outgroup_monophyletic": variant_run.fasttree_rooting.outgroup_monophyletic,
                "rooted_outgroup_taxa": list(
                    variant_run.fasttree_rooting.rooted_outgroup_taxa
                ),
                "warning_count": len(variant_run.fasttree_rooting.warnings),
            },
            "iqtree": {
                "requested_taxa": list(variant_run.iqtree_rooting.requested_taxa),
                "matched_taxa": list(variant_run.iqtree_rooting.matched_taxa),
                "outgroup_monophyletic": variant_run.iqtree_rooting.outgroup_monophyletic,
                "rooted_outgroup_taxa": list(
                    variant_run.iqtree_rooting.rooted_outgroup_taxa
                ),
                "warning_count": len(variant_run.iqtree_rooting.warnings),
            },
        },
        "rooted_engine_comparison": {
            "robinson_foulds_distance": row.rooted_engine_rf_distance,
            "normalized_robinson_foulds": row.rooted_engine_normalized_rf,
            "same_taxa_different_rooting": row.rooted_same_taxa_different_rooting,
        },
        "unrooted_conclusion_summary": {
            "stable_clade_count": row.stable_clade_count,
            "unstable_clade_count": row.unstable_clade_count,
            "engine_specific_clade_count": row.engine_specific_clade_count,
            "serious_conflict_count": row.serious_conflict_count,
        },
        "warnings": warnings,
        "output_inventory": output_inventory,
        "manifests": {
            "alignment": row.alignment_manifest_path,
            "trimming": row.trimming_manifest_path,
            "inference_comparison": row.inference_manifest_path,
            "model_selection": row.model_selection_manifest_path,
            "iqtree_support": row.iqtree_support_manifest_path,
            "fasttree": row.fasttree_manifest_path,
        },
    }


def _collect_output_inventory(output_root: Path) -> list[dict[str, object]]:
    return [
        {
            "relative_path": path.relative_to(output_root).as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path),
        }
        for path in sorted(output_root.rglob("*"))
        if path.is_file()
    ]


def _collect_warnings(variant_run: VariantRunLike) -> list[str]:
    warnings: list[str] = []
    warning_groups = (
        variant_run.alignment_workflow.run.warning_lines,
        variant_run.trimming_workflow.run.warning_lines,
        variant_run.inference_comparison.warnings,
        variant_run.fasttree_rooting.warnings,
        variant_run.iqtree_rooting.warnings,
    )
    for group in warning_groups:
        for warning in group:
            if warning not in warnings:
                warnings.append(warning)
    return warnings


def _copy_output(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    return Path(shutil.copy2(source, destination))


def _relative_bundle_path(source: Path, target: Path) -> str:
    base = source if source.is_dir() else source.parent
    return os.path.relpath(target, start=base).replace(os.sep, "/")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()
