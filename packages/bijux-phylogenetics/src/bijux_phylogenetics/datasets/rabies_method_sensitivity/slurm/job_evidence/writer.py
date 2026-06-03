from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path

from ..arrays import RabiesMethodSensitivitySlurmArrayStrategyReport
from .contracts import (
    RabiesMethodSensitivitySlurmJobEvidenceReport,
    RabiesMethodSensitivitySlurmJobEvidenceRow,
    TaskRecordLike,
    VariantRunLike,
)
from .payloads import (
    _build_job_payload,
    _collect_output_inventory,
    _collect_warnings,
    _copy_output,
    _relative_bundle_path,
)
from .presentation import _write_job_html_report
from .serialization import (
    write_rabies_method_sensitivity_slurm_job_evidence_summary_json,
    write_rabies_method_sensitivity_slurm_job_evidence_table,
)


def write_rabies_method_sensitivity_slurm_job_evidence_bundle(
    output_root: Path,
    *,
    bundle_root: Path,
    dataset_id: str,
    workflow_prefix: str,
    execution_mode: str,
    parallel_workers: int,
    task_records: tuple[TaskRecordLike, ...],
    variant_runs: tuple[VariantRunLike, ...],
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
        total_runtime_seconds=round(sum(row.total_runtime_seconds for row in rows), 6),
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


def _write_job_evidence_artifacts(
    *,
    output_root: Path,
    bundle_root: Path,
    task_record: TaskRecordLike,
    variant_run: VariantRunLike,
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
    artifact_file_count = len(
        [path for path in evidence_directory.iterdir() if path.is_file()]
    )
    return replace(row, artifact_file_count=artifact_file_count)
