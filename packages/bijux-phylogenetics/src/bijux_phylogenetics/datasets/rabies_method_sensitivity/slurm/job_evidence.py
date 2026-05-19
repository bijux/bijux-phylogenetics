from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
from typing import Protocol

from bijux_phylogenetics.render.html import write_html_report

from .arrays import (
    RabiesMethodSensitivitySlurmArrayStrategyReport,
)

__all__ = [
    "RabiesMethodSensitivitySlurmJobEvidenceReport",
    "RabiesMethodSensitivitySlurmJobEvidenceRow",
    "write_rabies_method_sensitivity_slurm_job_evidence_bundle",
    "write_rabies_method_sensitivity_slurm_job_evidence_summary_json",
    "write_rabies_method_sensitivity_slurm_job_evidence_table",
]


class _TaskRecordLike(Protocol):
    variant_id: str
    label: str
    execution_mode: str
    status: str
    log_path: Path
    error_code: str | None
    error_message: str | None


class _VariantConfigLike(Protocol):
    variant_id: str
    label: str
    alignment_mode: str
    trimming_mode: str
    trim_gap_threshold: float


class _EngineVersionLike(Protocol):
    text: str


class _EngineRunLike(Protocol):
    version: _EngineVersionLike
    command: list[str]
    warning_lines: list[str]
    runtime_seconds: float
    started_at_utc: str
    ended_at_utc: str


class _TrimmingSummaryLike(Protocol):
    retained_site_fraction: float
    removed_site_fraction: float


class _EngineWorkflowLike(Protocol):
    manifest_path: Path
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    config: dict[str, object]
    run: _EngineRunLike
    trimming_summary: _TrimmingSummaryLike | None


class _ConclusionSummaryLike(Protocol):
    stable_clade_count: int
    unstable_clade_count: int
    engine_specific_clade_count: int
    serious_conflict_count: int


class _InferenceComparisonLike(Protocol):
    manifest_path: Path
    step_manifests: dict[str, Path]
    commands: dict[str, list[str]]
    engine_versions: dict[str, str]
    runtime_seconds: float
    selected_model: str
    conclusion_summary: _ConclusionSummaryLike
    warnings: list[str]


class _RootingLike(Protocol):
    requested_taxa: tuple[str, ...]
    matched_taxa: tuple[str, ...]
    outgroup_monophyletic: bool | None
    rooted_outgroup_taxa: tuple[str, ...]
    warnings: list[str]


class _RootedComparisonLike(Protocol):
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    same_taxa_different_rooting: bool


class _VariantRunLike(Protocol):
    config: _VariantConfigLike
    alignment_workflow: _EngineWorkflowLike
    trimming_workflow: _EngineWorkflowLike
    inference_comparison: _InferenceComparisonLike
    fasttree_rooting: _RootingLike
    iqtree_rooting: _RootingLike
    rooted_engine_comparison: _RootedComparisonLike
    alignment_length: int
    trimmed_alignment_length: int


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmJobEvidenceRow:
    """One independently reviewable provenance package for one Slurm job."""

    partition_id: str
    array_index: int
    variant_id: str
    label: str
    execution_mode: str
    status: str
    script_path: str
    output_root: str
    evidence_directory: str
    evidence_json_path: str
    evidence_html_path: str
    task_log_copy_path: str
    alignment_manifest_path: str
    trimming_manifest_path: str
    inference_manifest_path: str
    model_selection_manifest_path: str
    iqtree_support_manifest_path: str
    fasttree_manifest_path: str
    alignment_mode: str
    trimming_mode: str
    trim_gap_threshold: float
    selected_model: str
    alignment_length: int
    trimmed_alignment_length: int
    total_runtime_seconds: float
    alignment_runtime_seconds: float
    trimming_runtime_seconds: float
    inference_runtime_seconds: float
    model_selection_runtime_seconds: float
    iqtree_support_runtime_seconds: float
    fasttree_runtime_seconds: float
    rooted_engine_rf_distance: int
    rooted_engine_normalized_rf: float
    rooted_same_taxa_different_rooting: bool
    serious_conflict_count: int
    stable_clade_count: int
    unstable_clade_count: int
    engine_specific_clade_count: int
    warning_count: int
    output_file_count: int
    output_byte_count: int
    artifact_file_count: int


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmJobEvidenceReport:
    """One reviewer-facing per-job evidence surface over the full workflow."""

    dataset_id: str
    workflow_prefix: str
    execution_mode: str
    parallel_workers: int
    bundle_root: Path
    evidence_root: Path
    index_path: Path
    summary_path: Path
    job_count: int
    completed_job_count: int
    failed_job_count: int
    total_runtime_seconds: float
    total_output_file_count: int
    total_output_byte_count: int
    total_artifact_file_count: int
    jobs: tuple[RabiesMethodSensitivitySlurmJobEvidenceRow, ...]


def write_rabies_method_sensitivity_slurm_job_evidence_bundle(
    output_root: Path,
    *,
    bundle_root: Path,
    dataset_id: str,
    workflow_prefix: str,
    execution_mode: str,
    parallel_workers: int,
    task_records: tuple[_TaskRecordLike, ...],
    variant_runs: tuple[_VariantRunLike, ...],
    array_strategy_report: RabiesMethodSensitivitySlurmArrayStrategyReport,
    execution_record_path: Path,
    workflow_manifest_path: Path,
) -> RabiesMethodSensitivitySlurmJobEvidenceReport:
    """Write one independent provenance package per planned Slurm job."""
    output_root.mkdir(parents=True, exist_ok=True)
    task_records_by_variant = {
        task_record.variant_id: task_record for task_record in task_records
    }
    variant_runs_by_variant = {
        variant_run.config.variant_id: variant_run for variant_run in variant_runs
    }
    rows: list[RabiesMethodSensitivitySlurmJobEvidenceRow] = []
    for member in sorted(
        array_strategy_report.members,
        key=lambda row: (row.partition_id, row.array_index, row.variant_id),
    ):
        task_record = task_records_by_variant[member.variant_id]
        variant_run = variant_runs_by_variant[member.variant_id]
        row = _write_job_evidence_artifacts(
            output_root=output_root,
            bundle_root=bundle_root,
            task_record=task_record,
            variant_run=variant_run,
            member=member,
            dataset_id=dataset_id,
            workflow_prefix=workflow_prefix,
            execution_record_path=execution_record_path,
            workflow_manifest_path=workflow_manifest_path,
        )
        rows.append(row)
    report = RabiesMethodSensitivitySlurmJobEvidenceReport(
        dataset_id=dataset_id,
        workflow_prefix=workflow_prefix,
        execution_mode=execution_mode,
        parallel_workers=parallel_workers,
        bundle_root=bundle_root,
        evidence_root=output_root,
        index_path=bundle_root / "slurm-job-evidence.tsv",
        summary_path=bundle_root / "slurm-job-evidence-summary.json",
        job_count=len(rows),
        completed_job_count=sum(1 for row in rows if row.status == "succeeded"),
        failed_job_count=sum(1 for row in rows if row.status != "succeeded"),
        total_runtime_seconds=round(
            sum(row.total_runtime_seconds for row in rows),
            6,
        ),
        total_output_file_count=sum(row.output_file_count for row in rows),
        total_output_byte_count=sum(row.output_byte_count for row in rows),
        total_artifact_file_count=sum(row.artifact_file_count for row in rows),
        jobs=tuple(rows),
    )
    write_rabies_method_sensitivity_slurm_job_evidence_table(
        report.index_path,
        report,
    )
    write_rabies_method_sensitivity_slurm_job_evidence_summary_json(
        report.summary_path,
        report,
    )
    return report


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
        "evidence_root": _relative_bundle_path(report.bundle_root, report.evidence_root),
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


def _write_job_evidence_artifacts(
    *,
    output_root: Path,
    bundle_root: Path,
    task_record: _TaskRecordLike,
    variant_run: _VariantRunLike,
    member: object,
    dataset_id: str,
    workflow_prefix: str,
    execution_record_path: Path,
    workflow_manifest_path: Path,
) -> RabiesMethodSensitivitySlurmJobEvidenceRow:
    variant_id = variant_run.config.variant_id
    evidence_directory = output_root / variant_id
    evidence_directory.mkdir(parents=True, exist_ok=True)
    copied_task_log_path = _copy_output(
        bundle_root / member.task_log_path,
        evidence_directory / "task.log",
    )
    copied_manifest_paths = {
        "alignment_manifest_path": _copy_output(
            variant_run.alignment_workflow.manifest_path,
            evidence_directory / "alignment.manifest.json",
        ),
        "trimming_manifest_path": _copy_output(
            variant_run.trimming_workflow.manifest_path,
            evidence_directory / "trimming.manifest.json",
        ),
        "inference_manifest_path": _copy_output(
            variant_run.inference_comparison.manifest_path,
            evidence_directory / "inference-comparison.manifest.json",
        ),
        "model_selection_manifest_path": _copy_output(
            variant_run.inference_comparison.step_manifests["model_selection"],
            evidence_directory / "model-selection.manifest.json",
        ),
        "iqtree_support_manifest_path": _copy_output(
            variant_run.inference_comparison.step_manifests["iqtree_support"],
            evidence_directory / "iqtree-support.manifest.json",
        ),
        "fasttree_manifest_path": _copy_output(
            variant_run.inference_comparison.step_manifests["fasttree"],
            evidence_directory / "fasttree.manifest.json",
        ),
    }
    variant_output_root = bundle_root / member.bundle_output_directory
    output_inventory = _collect_output_inventory(variant_output_root)
    warnings = _collect_warnings(variant_run)
    total_runtime_seconds = round(
        variant_run.alignment_workflow.run.runtime_seconds
        + variant_run.trimming_workflow.run.runtime_seconds
        + variant_run.inference_comparison.runtime_seconds,
        6,
    )
    evidence_json_path = evidence_directory / "job-evidence.json"
    evidence_html_path = evidence_directory / "job-evidence.html"
    row = RabiesMethodSensitivitySlurmJobEvidenceRow(
        partition_id=member.partition_id,
        array_index=member.array_index,
        variant_id=variant_id,
        label=task_record.label,
        execution_mode=task_record.execution_mode,
        status=task_record.status,
        script_path=str(member.script_path),
        output_root=str(member.bundle_output_directory),
        evidence_directory=_relative_bundle_path(bundle_root, evidence_directory),
        evidence_json_path=_relative_bundle_path(bundle_root, evidence_json_path),
        evidence_html_path=_relative_bundle_path(bundle_root, evidence_html_path),
        task_log_copy_path=_relative_bundle_path(bundle_root, copied_task_log_path),
        alignment_manifest_path=_relative_bundle_path(
            bundle_root,
            copied_manifest_paths["alignment_manifest_path"],
        ),
        trimming_manifest_path=_relative_bundle_path(
            bundle_root,
            copied_manifest_paths["trimming_manifest_path"],
        ),
        inference_manifest_path=_relative_bundle_path(
            bundle_root,
            copied_manifest_paths["inference_manifest_path"],
        ),
        model_selection_manifest_path=_relative_bundle_path(
            bundle_root,
            copied_manifest_paths["model_selection_manifest_path"],
        ),
        iqtree_support_manifest_path=_relative_bundle_path(
            bundle_root,
            copied_manifest_paths["iqtree_support_manifest_path"],
        ),
        fasttree_manifest_path=_relative_bundle_path(
            bundle_root,
            copied_manifest_paths["fasttree_manifest_path"],
        ),
        alignment_mode=variant_run.config.alignment_mode,
        trimming_mode=variant_run.config.trimming_mode,
        trim_gap_threshold=variant_run.config.trim_gap_threshold,
        selected_model=variant_run.inference_comparison.selected_model,
        alignment_length=variant_run.alignment_length,
        trimmed_alignment_length=variant_run.trimmed_alignment_length,
        total_runtime_seconds=total_runtime_seconds,
        alignment_runtime_seconds=variant_run.alignment_workflow.run.runtime_seconds,
        trimming_runtime_seconds=variant_run.trimming_workflow.run.runtime_seconds,
        inference_runtime_seconds=variant_run.inference_comparison.runtime_seconds,
        model_selection_runtime_seconds=(
            variant_run.inference_comparison.model_selection_workflow.run.runtime_seconds
        ),
        iqtree_support_runtime_seconds=(
            variant_run.inference_comparison.iqtree_support_workflow.run.runtime_seconds
        ),
        fasttree_runtime_seconds=(
            variant_run.inference_comparison.fasttree_workflow.run.runtime_seconds
        ),
        rooted_engine_rf_distance=(
            variant_run.rooted_engine_comparison.robinson_foulds_distance
        ),
        rooted_engine_normalized_rf=(
            variant_run.rooted_engine_comparison.normalized_robinson_foulds
        ),
        rooted_same_taxa_different_rooting=(
            variant_run.rooted_engine_comparison.same_taxa_different_rooting
        ),
        serious_conflict_count=(
            variant_run.inference_comparison.conclusion_summary.serious_conflict_count
        ),
        stable_clade_count=(
            variant_run.inference_comparison.conclusion_summary.stable_clade_count
        ),
        unstable_clade_count=(
            variant_run.inference_comparison.conclusion_summary.unstable_clade_count
        ),
        engine_specific_clade_count=(
            variant_run.inference_comparison.conclusion_summary.engine_specific_clade_count
        ),
        warning_count=len(warnings),
        output_file_count=len(output_inventory),
        output_byte_count=sum(item["size_bytes"] for item in output_inventory),
        artifact_file_count=0,
    )
    payload = _build_job_payload(
        bundle_root=bundle_root,
        row=row,
        task_record=task_record,
        variant_run=variant_run,
        member=member,
        dataset_id=dataset_id,
        workflow_prefix=workflow_prefix,
        execution_record_path=execution_record_path,
        workflow_manifest_path=workflow_manifest_path,
        output_inventory=output_inventory,
        warnings=warnings,
    )
    evidence_json_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_job_html_report(
        path=evidence_html_path,
        row=row,
        payload=payload,
        warnings=warnings,
        output_inventory=output_inventory,
        execution_record_path=execution_record_path,
        workflow_manifest_path=workflow_manifest_path,
    )
    artifact_file_count = len([path for path in evidence_directory.iterdir() if path.is_file()])
    return replace(row, artifact_file_count=artifact_file_count)


def _build_job_payload(
    *,
    bundle_root: Path,
    row: RabiesMethodSensitivitySlurmJobEvidenceRow,
    task_record: _TaskRecordLike,
    variant_run: _VariantRunLike,
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


def _write_job_html_report(
    *,
    path: Path,
    row: RabiesMethodSensitivitySlurmJobEvidenceRow,
    payload: dict[str, object],
    warnings: list[str],
    output_inventory: list[dict[str, object]],
    execution_record_path: Path,
    workflow_manifest_path: Path,
) -> Path:
    commands = payload["commands"]
    sections = [
        (
            "job",
            "\n".join(
                [
                    f"variant_id: {row.variant_id}",
                    f"label: {row.label}",
                    f"status: {row.status}",
                    f"partition_id: {row.partition_id}",
                    f"array_index: {row.array_index}",
                    f"execution_mode: {row.execution_mode}",
                    f"script_path: {row.script_path}",
                ]
            ),
        ),
        (
            "configuration",
            "\n".join(
                [
                    f"alignment_mode: {row.alignment_mode}",
                    f"trimming_mode: {row.trimming_mode}",
                    f"trim_gap_threshold: {_format_float(row.trim_gap_threshold)}",
                    f"selected_model: {row.selected_model}",
                    f"alignment_length: {row.alignment_length}",
                    f"trimmed_alignment_length: {row.trimmed_alignment_length}",
                ]
            ),
        ),
        (
            "commands",
            "\n".join(
                [
                    f"{name}: {shlex.join(list(command))}"
                    for name, command in dict(commands).items()
                ]
            ),
        ),
        (
            "runtime-and-findings",
            "\n".join(
                [
                    f"total_runtime_seconds: {_format_float(row.total_runtime_seconds)}",
                    (
                        "model_selection_runtime_seconds: "
                        f"{_format_float(row.model_selection_runtime_seconds)}"
                    ),
                    (
                        "iqtree_support_runtime_seconds: "
                        f"{_format_float(row.iqtree_support_runtime_seconds)}"
                    ),
                    (
                        "fasttree_runtime_seconds: "
                        f"{_format_float(row.fasttree_runtime_seconds)}"
                    ),
                    f"serious_conflict_count: {row.serious_conflict_count}",
                    f"rooted_engine_rf_distance: {row.rooted_engine_rf_distance}",
                    (
                        "rooted_same_taxa_different_rooting: "
                        f"{str(row.rooted_same_taxa_different_rooting).lower()}"
                    ),
                ]
            ),
        ),
        (
            "warnings",
            "none" if not warnings else "\n".join(warnings),
        ),
        (
            "output-inventory",
            "\n".join(
                f"{item['relative_path']}: {item['size_bytes']} bytes {item['sha256']}"
                for item in output_inventory
            ),
        ),
    ]
    artifact_links = [
        ("task-log", Path("task.log").as_posix(), None),
        (
            "workflow-manifest",
            _relative_bundle_path(path, workflow_manifest_path),
            None,
        ),
        (
            "execution-record",
            _relative_bundle_path(path, execution_record_path),
            None,
        ),
        ("alignment-manifest", "alignment.manifest.json", None),
        ("trimming-manifest", "trimming.manifest.json", None),
        ("inference-comparison-manifest", "inference-comparison.manifest.json", None),
        ("model-selection-manifest", "model-selection.manifest.json", None),
        ("iqtree-support-manifest", "iqtree-support.manifest.json", None),
        ("fasttree-manifest", "fasttree.manifest.json", None),
    ]
    return write_html_report(
        title=f"Rabies Slurm Job Evidence: {row.variant_id}",
        sections=sections,
        out_path=path,
        embedded_json=payload,
        summary_metrics=[
            ("variant", row.variant_id),
            ("status", row.status),
            ("selected model", row.selected_model),
            ("runtime seconds", _format_float(row.total_runtime_seconds)),
            ("output files", row.output_file_count),
            ("warnings", row.warning_count),
            ("serious conflicts", row.serious_conflict_count),
            ("rooted RF", row.rooted_engine_rf_distance),
        ],
        artifact_links=artifact_links,
    )


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


def _collect_warnings(variant_run: _VariantRunLike) -> list[str]:
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


def _serialize_row(row: RabiesMethodSensitivitySlurmJobEvidenceRow) -> dict[str, object]:
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


def _copy_output(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    return Path(shutil.copy2(source, destination))


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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _format_float(value: float) -> str:
    return format(value, ".12g")
