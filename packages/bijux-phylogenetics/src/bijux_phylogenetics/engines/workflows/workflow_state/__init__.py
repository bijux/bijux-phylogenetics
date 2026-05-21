# ruff: noqa: F401
from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.presentation import build_tree_comparison_report
from bijux_phylogenetics.io.fasta import (
    infer_alignment_alphabet,
    load_fasta_alignment,
    write_fasta_alignment,
)
from bijux_phylogenetics.io.fasta.coding import (
    back_translate_aligned_coding_sequences,
    prepare_coding_sequences_for_alignment,
    translate_prepared_coding_sequences,
)
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.phylo.alignment import (
    AlignmentAlphabet,
    AlignmentRecord,
    AlignmentSummary,
    CodingSequenceExclusion,
)
from bijux_phylogenetics.phylo.alignment.partitions import LocusPartition
from bijux_phylogenetics.runtime.errors import EngineWorkflowError, PhylogeneticsError

from ...artifacts.bootstrap import (
    build_bootstrap_support_histogram_rows,
    build_bootstrap_support_rows,
    build_low_support_bootstrap_rows,
    write_bootstrap_support_histogram,
    write_bootstrap_support_table,
)
from ...artifacts.fasttree import (
    build_fasttree_low_support_rows,
    build_fasttree_support_histogram_rows,
    build_fasttree_support_rows,
    write_fasttree_support_histogram,
    write_fasttree_support_table,
)
from ...artifacts.iqtree import (
    IqtreeModelCandidate,
    IqtreeModelSelectionSummary,
    parse_best_model_file,
    parse_iqtree_model_selection_summary,
    parse_log_likelihood_file,
    resolve_iqtree_model_sidecar,
    write_iqtree_model_candidates_table,
)
from ...artifacts.sh_alrt import (
    build_conflicting_sh_alrt_support_rows,
    build_sh_alrt_support_rows,
    write_sh_alrt_support_table,
)
from ...artifacts.support import (
    BootstrapSupportNode,
    BootstrapSupportSummaryReport,
    FastTreeSupportNode,
    FastTreeSupportSummaryReport,
    ShAlrtSupportNode,
    ShAlrtSupportSummaryReport,
    WeakBackboneReport,
)
from ...common import (
    EngineRunReport,
    EngineVersionInfo,
    active_engine_run_is_live,
    build_engine_output_error,
    build_file_checksums,
    cleanup_incomplete_engine_run,
    clear_incomplete_engine_run,
    engine_active_marker_path,
    engine_incomplete_marker_path,
    execute_engine_command,
    load_active_engine_run,
    load_engine_manifest,
    load_incomplete_engine_run,
    load_unaligned_fasta,
    observe_engine_outputs,
    read_engine_version,
    resolve_engine_executable,
    update_incomplete_engine_run,
    validate_timeout_seconds,
    write_engine_manifest,
)
from ...validation import (
    detect_weakly_supported_backbone,
    summarize_bootstrap_support_distribution,
    summarize_fasttree_support_distribution,
    summarize_sh_alrt_support_distribution,
)
from ..models import (
    AlignmentTrimmingSummary,
    CodonAwareAlignmentWorkflowReport,
    EngineWorkflowReport,
    ExternalTreeComparisonReport,
    IqtreeSupportValue,
    IqtreeWorkflowSummary,
)
from ..models import (
    PreparedIqtreePartitions as _PreparedIqtreePartitions,
)

_INCOMPLETE_RUN_POLICIES = {"reject", "clean"}

def _validate_incomplete_run_policy(policy: str) -> str:
    if policy not in _INCOMPLETE_RUN_POLICIES:
        available = ", ".join(sorted(_INCOMPLETE_RUN_POLICIES))
        raise ValueError(
            f"incomplete_run_policy must be one of: {available}; got {policy}"
        )
    return policy


def _resolve_incomplete_workflow_state(
    *,
    manifest_path: Path,
    incomplete_run_policy: str,
) -> list[str]:
    _validate_incomplete_run_policy(incomplete_run_policy)
    active_record = load_active_engine_run(manifest_path)
    if active_record is not None and active_engine_run_is_live(active_record):
        raise EngineWorkflowError(
            "engine workflow is already running for the requested output manifest",
            code="engine_workflow_already_running",
            details={
                "manifest_path": str(manifest_path),
                "marker_path": str(engine_active_marker_path(manifest_path)),
                "running_process_id": active_record.process_id,
                "running_workflow": active_record.workflow,
                "running_engine_name": active_record.engine_name,
            },
        )
    record = load_incomplete_engine_run(manifest_path)
    if record is None:
        return []
    if incomplete_run_policy == "clean":
        cleanup_incomplete_engine_run(manifest_path)
        return [
            "removed outputs from a previously incomplete engine run before restarting"
        ]
    marker_path = engine_incomplete_marker_path(manifest_path)
    observed_outputs = [
        {
            "output_name": observation.output_name,
            "path": str(observation.path),
            "exists": observation.exists,
            "path_kind": observation.path_kind,
            "size_bytes": observation.size_bytes,
            "sha256": observation.sha256,
        }
        for observation in record.observed_outputs
    ]
    raise EngineWorkflowError(
        "a previous engine run left incomplete outputs and resume could not safely "
        f"reuse them; marker: {marker_path}",
        code="engine_incomplete_outputs_present",
        details={
            "manifest_path": str(manifest_path),
            "marker_path": str(marker_path),
            "engine_name": record.engine_name,
            "workflow": record.workflow,
            "failure_reason": (
                record.failure_reason
                if record.failure_reason is not None
                else "engine_run_incomplete"
            ),
            "failure_message": record.failure_message,
            "timed_out": record.timed_out,
            "exit_code": record.exit_code,
            "timeout_seconds": record.timeout_seconds,
            "missing_output_names": list(record.missing_output_names),
            "observed_outputs": observed_outputs,
            "incomplete_run_policy": incomplete_run_policy,
            "available_actions": ["resume", "clean"],
        },
    )

def _record_output_validation_failure(
    manifest_path: Path,
    run: EngineRunReport,
    error: PhylogeneticsError,
) -> None:
    observations = observe_engine_outputs(run.output_paths)
    update_incomplete_engine_run(
        manifest_path,
        ended_at_utc=run.ended_at_utc,
        timed_out=run.timed_out,
        exit_code=run.exit_code,
        failure_reason=error.code,
        failure_message=(
            f"{run.engine_name} {run.workflow} produced outputs that failed "
            f"validation: {error.code}"
        ),
        missing_output_names=[
            observation.output_name
            for observation in observations
            if not observation.exists
        ],
        observed_outputs=observations,
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


from .paths import (
    _manifest_path_from_output,
    _partition_alignment_file_name,
    _partition_support_path,
    _prefix_path,
    _sidecar,
)
from .artifact_validation import (
    _ensure_inference_ready_alignment,
    _require_nonempty_text_output,
    _validate_alignment_output,
    _validate_iqtree_required_artifacts,
    _validate_support_value_count,
    _validate_tree_output,
    _validate_tree_set_output,
)
from .coding_alignment import (
    _build_alignment_trimming_summary,
    _write_coding_exclusion_table,
    _write_coding_summary_table,
)
from .report_restore import (
    _restore_codon_aware_alignment_report,
    _restore_workflow_report,
)
