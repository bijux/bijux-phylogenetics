from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment import AlignmentAlphabet

from ...common import (
    EngineVersionInfo,
    build_file_checksums,
    clear_incomplete_engine_run,
    load_engine_manifest,
    write_engine_manifest,
)
from ..models import CodonAwareAlignmentWorkflowReport, EngineWorkflowReport
from .report_restore import (
    _restore_codon_aware_alignment_report,
    _restore_workflow_report,
)


def _resume_existing_workflow(
    *,
    manifest_path: Path,
    input_paths: list[Path],
    expected_command: list[str],
    expected_version: EngineVersionInfo,
) -> EngineWorkflowReport | None:
    if not manifest_path.exists():
        return None
    payload = load_engine_manifest(manifest_path)
    report = _restore_workflow_report(payload)
    if report.run.command != expected_command:
        return None
    if report.run.version.text != expected_version.text:
        return None
    current_input_checksums = build_file_checksums(input_paths)
    if report.input_checksums != current_input_checksums:
        return None
    if any(not path.exists() for path in report.output_paths.values()):
        return None
    current_output_checksums = build_file_checksums(list(report.output_paths.values()))
    if report.output_checksums != current_output_checksums:
        return None
    clear_incomplete_engine_run(manifest_path)
    report.resumed = True
    return report


def _persist_workflow_report(report: EngineWorkflowReport) -> EngineWorkflowReport:
    report.output_checksums = build_file_checksums(list(report.output_paths.values()))
    clear_incomplete_engine_run(report.manifest_path)
    write_engine_manifest(report.manifest_path, report)
    return report


def _resume_existing_codon_aware_alignment(
    *,
    manifest_path: Path,
    input_path: Path,
    expected_command: list[str],
    expected_version: EngineVersionInfo,
    expected_sequence_type: AlignmentAlphabet,
    expected_genetic_code_id: int,
) -> CodonAwareAlignmentWorkflowReport | None:
    if not manifest_path.exists():
        return None
    payload = load_engine_manifest(manifest_path)
    report = _restore_codon_aware_alignment_report(payload)
    if report.run.command != expected_command:
        return None
    if report.run.version.text != expected_version.text:
        return None
    if report.sequence_type != expected_sequence_type:
        return None
    if report.genetic_code_id != expected_genetic_code_id:
        return None
    current_input_checksums = build_file_checksums([input_path])
    if report.input_checksums != current_input_checksums:
        return None
    if any(not path.exists() for path in report.output_paths.values()):
        return None
    current_output_checksums = build_file_checksums(list(report.output_paths.values()))
    if report.output_checksums != current_output_checksums:
        return None
    clear_incomplete_engine_run(manifest_path)
    report.resumed = True
    return report


def _resume_has_bootstrap_review_outputs(report: EngineWorkflowReport) -> bool:
    return (
        report.bootstrap_support_summary is not None
        and report.weak_backbone_report is not None
        and report.output_paths.get("support_table") is not None
        and report.output_paths.get("low_support_branches") is not None
        and report.output_paths.get("support_histogram") is not None
    )


def _resume_has_fasttree_review_outputs(report: EngineWorkflowReport) -> bool:
    return (
        report.fasttree_support_summary is not None
        and report.output_paths.get("support_table") is not None
        and report.output_paths.get("low_support_branches") is not None
        and report.output_paths.get("support_histogram") is not None
    )


def _resume_has_sh_alrt_review_outputs(report: EngineWorkflowReport) -> bool:
    return (
        report.sh_alrt_support_summary is not None
        and report.output_paths.get("support_table") is not None
        and report.output_paths.get("conflicting_support_branches") is not None
    )
