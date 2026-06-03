# ruff: noqa: F401
from __future__ import annotations

import hashlib
from pathlib import Path
import re

from bijux_phylogenetics.compare.presentation import build_tree_comparison_report
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.fasta import (
    infer_alignment_alphabet,
    load_fasta_alignment,
    write_fasta_alignment,
)
from bijux_phylogenetics.io.fasta.coding import (
    back_translate_aligned_coding_sequences,
    classify_sequence_coding_behavior,
    prepare_coding_sequences_for_alignment,
    translate_prepared_coding_sequences,
)
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.alignment import (
    AlignmentAlphabet,
    AlignmentRecord,
    AlignmentSummary,
    CodingSequenceExclusion,
)
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    build_partition_summary_report,
    normalize_partition_data_type,
    parse_locus_partitions,
    slice_partition_sequence,
    write_locus_partitions,
    write_partition_summary_table,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError, PhylogeneticsError
from bijux_phylogenetics.trees import load_tree_set

from ..artifacts.bootstrap import (
    build_bootstrap_support_histogram_rows,
    build_bootstrap_support_rows,
    build_low_support_bootstrap_rows,
    write_bootstrap_support_histogram,
    write_bootstrap_support_table,
)
from ..artifacts.fasttree import (
    build_fasttree_low_support_rows,
    build_fasttree_support_histogram_rows,
    build_fasttree_support_rows,
    write_fasttree_support_histogram,
    write_fasttree_support_table,
)
from ..artifacts.iqtree import (
    IqtreeModelCandidate,
    IqtreeModelSelectionSummary,
    parse_best_model_file,
    parse_iqtree_model_selection_summary,
    parse_log_likelihood_file,
    resolve_iqtree_model_sidecar,
    write_iqtree_model_candidates_table,
)
from ..artifacts.sh_alrt import (
    build_conflicting_sh_alrt_support_rows,
    build_sh_alrt_support_rows,
    write_sh_alrt_support_table,
)
from ..artifacts.support import (
    BootstrapSupportNode,
    BootstrapSupportSummaryReport,
    FastTreeSupportNode,
    FastTreeSupportSummaryReport,
    ShAlrtSupportNode,
    ShAlrtSupportSummaryReport,
    WeakBackboneReport,
)
from ..common import (
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
from ..validation.audits import (
    detect_weakly_supported_backbone,
    summarize_bootstrap_support_distribution,
    summarize_fasttree_support_distribution,
    summarize_sh_alrt_support_distribution,
)
from ..validation.preflight import require_external_engine_surface
from .models import (
    AlignmentTrimmingSummary,
    CodonAwareAlignmentWorkflowReport,
    EngineWorkflowReport,
    ExternalTreeComparisonReport,
    IqtreeSupportValue,
    IqtreeWorkflowSummary,
)
from .models import (
    PreparedIqtreePartitions as _PreparedIqtreePartitions,
)
from .state import (
    _build_alignment_trimming_summary,
    _ensure_inference_ready_alignment,
    _manifest_path_from_output,
    _partition_alignment_file_name,
    _partition_support_path,
    _persist_workflow_report,
    _prefix_path,
    _record_output_validation_failure,
    _require_nonempty_text_output,
    _resolve_incomplete_workflow_state,
    _restore_codon_aware_alignment_report,
    _restore_workflow_report,
    _resume_existing_codon_aware_alignment,
    _resume_existing_workflow,
    _resume_has_bootstrap_review_outputs,
    _resume_has_fasttree_review_outputs,
    _resume_has_sh_alrt_review_outputs,
    _sidecar,
    _validate_alignment_output,
    _validate_incomplete_run_policy,
    _validate_iqtree_required_artifacts,
    _validate_support_value_count,
    _validate_tree_output,
    _validate_tree_set_output,
    _write_alignment_trimming_summary_table,
    _write_coding_exclusion_table,
    _write_coding_summary_table,
)

_MAFFT_ALIGNMENT_MODE_ARGUMENTS: dict[str, tuple[str, ...]] = {
    "auto": ("--auto",),
    "linsi": ("--localpair", "--maxiterate", "1000"),
    "ginsi": ("--globalpair", "--maxiterate", "1000"),
    "einsi": ("--ep", "0", "--genafpair", "--maxiterate", "1000"),
    "fast": ("--retree", "2", "--maxiterate", "0"),
}
_TRIMAL_TRIMMING_MODES: tuple[str, ...] = (
    "gap-threshold",
    "gappyout",
    "strict",
    "strictplus",
    "automated1",
)


def list_mafft_alignment_modes() -> tuple[str, ...]:
    """Return the supported named MAFFT alignment strategies."""
    return tuple(_MAFFT_ALIGNMENT_MODE_ARGUMENTS)


def resolve_mafft_alignment_mode(mode: str) -> tuple[str, ...]:
    """Resolve one named MAFFT alignment strategy into explicit engine arguments."""
    try:
        return _MAFFT_ALIGNMENT_MODE_ARGUMENTS[mode]
    except KeyError as error:
        available = ", ".join(sorted(_MAFFT_ALIGNMENT_MODE_ARGUMENTS))
        raise ValueError(
            f"unsupported mafft alignment mode '{mode}', expected one of: {available}"
        ) from error


def list_trimal_trimming_modes() -> tuple[str, ...]:
    """Return the supported named trimAl trimming strategies."""
    return _TRIMAL_TRIMMING_MODES


def resolve_trimal_trimming_mode(
    mode: str,
    *,
    gap_threshold: float,
) -> tuple[str, ...]:
    """Resolve one named trimAl strategy into explicit engine arguments."""
    if mode == "gap-threshold":
        return ("-gt", f"{gap_threshold:.6f}")
    if mode in {"gappyout", "strict", "strictplus", "automated1"}:
        return (f"-{mode}",)
    available = ", ".join(sorted(_TRIMAL_TRIMMING_MODES))
    raise ValueError(
        f"unsupported trimAl trimming mode '{mode}', expected one of: {available}"
    )


def run_multiple_sequence_alignment(
    input_path: Path,
    out_path: Path,
    *,
    executable: str | Path = "mafft",
    mode: str = "auto",
    extra_args: tuple[str, ...] = (),
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run a multiple-sequence alignment engine against an unaligned FASTA file."""
    load_unaligned_fasta(input_path)
    validate_timeout_seconds(timeout_seconds)
    require_external_engine_surface(
        workflow_id="multiple-sequence-alignment",
        summary="MAFFT multiple-sequence-alignment workflow.",
        required_engines=("mafft",),
        executables={"mafft": executable},
        preserve_missing_error=True,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path = _manifest_path_from_output(out_path)
    mode_args = resolve_mafft_alignment_mode(mode)
    version = read_engine_version(
        "mafft",
        executable,
        version_args=("--version",),
        timeout_seconds=timeout_seconds,
    )
    resolved = resolve_engine_executable(executable)
    command = [resolved, *mode_args, *extra_args, str(input_path.resolve())]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=[input_path],
            expected_command=command,
            expected_version=version,
        )
        if resumed is not None:
            return resumed
    incomplete_notes = _resolve_incomplete_workflow_state(
        manifest_path=manifest_path,
        incomplete_run_policy=incomplete_run_policy,
    )
    run = execute_engine_command(
        engine_name="mafft",
        workflow="multiple-sequence-alignment",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_path.parent,
        stdout_path=out_path,
        stderr_path=_sidecar(out_path, "stderr.log"),
        output_paths={"alignment": out_path},
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _validate_alignment_output(
            out_path,
            engine_name="mafft",
            workflow="multiple-sequence-alignment",
            output_name="alignment",
            artifact_kind="multiple-sequence-alignment",
        )
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    report = EngineWorkflowReport(
        workflow="multiple-sequence-alignment",
        engine_name="mafft",
        input_paths=[input_path],
        output_paths={"alignment": out_path},
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([input_path]),
        output_checksums={},
        config={
            "mode": mode,
            "extra_args": list(extra_args),
            "timeout_seconds": timeout_seconds,
        },
        notes=[
            f"mafft alignment mode: {mode}",
            "alignment output validated as deterministic equal-length FASTA",
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)


def run_codon_aware_multiple_sequence_alignment(
    input_path: Path,
    out_path: Path,
    *,
    executable: str | Path = "mafft",
    mode: str = "auto",
    sequence_type: AlignmentAlphabet | None = None,
    genetic_code: int | str | None = None,
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> CodonAwareAlignmentWorkflowReport:
    """Align coding nucleotide sequences through a translated amino-acid guide."""
    validate_timeout_seconds(timeout_seconds)
    require_external_engine_surface(
        workflow_id="codon-aware-multiple-sequence-alignment",
        summary="MAFFT codon-aware multiple-sequence-alignment workflow.",
        required_engines=("mafft",),
        executables={"mafft": executable},
        preserve_missing_error=True,
    )
    prepared_records, preparation = prepare_coding_sequences_for_alignment(
        input_path,
        sequence_type=sequence_type,
        genetic_code=genetic_code,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    guide_input_path = _sidecar(out_path, "guide-input.fasta")
    guide_alignment_path = _sidecar(out_path, "guide-alignment.fasta")
    exclusion_report_path = _sidecar(out_path, "excluded.tsv")
    coding_summary_path = _sidecar(out_path, "coding-summary.tsv")
    manifest_path = _manifest_path_from_output(out_path)
    guide_records = translate_prepared_coding_sequences(
        prepared_records,
        genetic_code=preparation.genetic_code_id,
    )
    write_fasta_alignment(guide_input_path, guide_records)
    mode_args = resolve_mafft_alignment_mode(mode)
    version = read_engine_version(
        "mafft",
        executable,
        version_args=("--version",),
        timeout_seconds=timeout_seconds,
    )
    resolved = resolve_engine_executable(executable)
    command = [resolved, *mode_args, str(guide_input_path.resolve())]
    if resume:
        resumed = _resume_existing_codon_aware_alignment(
            manifest_path=manifest_path,
            input_path=input_path,
            expected_command=command,
            expected_version=version,
            expected_sequence_type=preparation.sequence_type,
            expected_genetic_code_id=preparation.genetic_code_id,
        )
        if resumed is not None:
            return resumed
    incomplete_notes = _resolve_incomplete_workflow_state(
        manifest_path=manifest_path,
        incomplete_run_policy=incomplete_run_policy,
    )
    run = execute_engine_command(
        engine_name="mafft",
        workflow="codon-aware-multiple-sequence-alignment",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_path.parent,
        stdout_path=guide_alignment_path,
        stderr_path=_sidecar(out_path, "stderr.log"),
        output_paths={"guide_alignment": guide_alignment_path},
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _validate_alignment_output(
            guide_alignment_path,
            engine_name="mafft",
            workflow="codon-aware-multiple-sequence-alignment",
            output_name="guide_alignment",
            artifact_kind="mafft-guide-alignment",
        )
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    aligned_guide = load_fasta_alignment(guide_alignment_path)
    codon_records = back_translate_aligned_coding_sequences(
        aligned_guide,
        coding_records=prepared_records,
    )
    write_fasta_alignment(out_path, codon_records)
    try:
        codon_summary = _validate_alignment_output(
            out_path,
            engine_name="mafft",
            workflow="codon-aware-multiple-sequence-alignment",
            output_name="alignment",
            artifact_kind="codon-aware-alignment",
        )
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    if codon_summary.alignment_length % 3 != 0:
        raise EngineWorkflowError(
            "codon-aware alignment produced an alignment length that is not divisible by three"
        )
    _write_coding_exclusion_table(exclusion_report_path, preparation.excluded_sequences)
    _write_coding_summary_table(
        coding_summary_path,
        input_path=input_path,
        genetic_code=preparation.genetic_code_id,
        exclusions=preparation.excluded_sequences,
    )
    output_paths = {
        "alignment": out_path,
        "guide_input": guide_input_path,
        "guide_alignment": guide_alignment_path,
        "excluded_sequences": exclusion_report_path,
        "coding_summary": coding_summary_path,
    }
    report = CodonAwareAlignmentWorkflowReport(
        workflow="codon-aware-multiple-sequence-alignment",
        engine_name="mafft",
        input_path=input_path,
        output_paths=output_paths,
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([input_path]),
        output_checksums={},
        config={
            "mode": mode,
            "sequence_type": preparation.sequence_type,
            "genetic_code_id": preparation.genetic_code_id,
            "timeout_seconds": timeout_seconds,
        },
        sequence_type=preparation.sequence_type,
        genetic_code_id=preparation.genetic_code_id,
        genetic_code_name=preparation.genetic_code_name,
        input_sequence_count=preparation.input_sequence_count,
        accepted_sequence_count=preparation.accepted_sequence_count,
        invalid_codon_sequence_count=preparation.invalid_codon_sequence_count,
        excluded_sequences=preparation.excluded_sequences,
        terminal_stop_sequence_count=preparation.terminal_stop_sequence_count,
        notes=[
            "codon-aware alignment preserved nucleotide codon triplets through amino-acid guide alignment",
            f"mafft alignment mode: {mode}",
            f"genetic code: {preparation.genetic_code_name} ({preparation.genetic_code_id})",
            f"accepted coding sequences: {preparation.accepted_sequence_count} of {preparation.input_sequence_count}",
            f"retained nucleotide alignment length: {codon_summary.alignment_length}",
            *incomplete_notes,
        ],
        warnings=list(dict.fromkeys(run.warning_lines + preparation.warnings)),
        resumed=False,
    )
    report.output_checksums = build_file_checksums(list(output_paths.values()))
    write_engine_manifest(manifest_path, report)
    return report


def run_alignment_trimming(
    input_path: Path,
    out_path: Path,
    *,
    executable: str | Path = "trimal",
    mode: str = "gap-threshold",
    gap_threshold: float = 0.1,
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run an external alignment trimming engine against an aligned FASTA file."""
    load_fasta_alignment(input_path)
    validate_timeout_seconds(timeout_seconds)
    require_external_engine_surface(
        workflow_id="alignment-trimming",
        summary="trimAl alignment-trimming workflow.",
        required_engines=("trimal",),
        executables={"trimal": executable},
        preserve_missing_error=True,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    input_summary = summarise_fasta(input_path)
    trimming_summary_path = _sidecar(out_path, "retained-sites.tsv")
    manifest_path = _manifest_path_from_output(out_path)
    mode_args = resolve_trimal_trimming_mode(mode, gap_threshold=gap_threshold)
    version = read_engine_version(
        "trimal",
        executable,
        version_args=("--version",),
        timeout_seconds=timeout_seconds,
    )
    resolved = resolve_engine_executable(executable)
    command = [
        resolved,
        "-in",
        str(input_path.resolve()),
        "-out",
        str(out_path.resolve()),
        *mode_args,
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=[input_path],
            expected_command=command,
            expected_version=version,
        )
        if resumed is not None:
            return resumed
    incomplete_notes = _resolve_incomplete_workflow_state(
        manifest_path=manifest_path,
        incomplete_run_policy=incomplete_run_policy,
    )
    run = execute_engine_command(
        engine_name="trimal",
        workflow="alignment-trimming",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_path.parent,
        stdout_path=_sidecar(out_path, "stdout.log"),
        stderr_path=_sidecar(out_path, "stderr.log"),
        output_paths={"trimmed_alignment": out_path},
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        trimmed_summary = _validate_alignment_output(
            out_path,
            engine_name="trimal",
            workflow="alignment-trimming",
            output_name="trimmed_alignment",
            artifact_kind="trimmed-alignment",
        )
        trimming_summary = _build_alignment_trimming_summary(
            mode=mode,
            gap_threshold=gap_threshold,
            input_summary=input_summary,
            trimmed_summary=trimmed_summary,
        )
        _write_alignment_trimming_summary_table(
            trimming_summary_path,
            summary=trimming_summary,
        )
        _require_nonempty_text_output(
            trimming_summary_path,
            engine_name="trimal",
            workflow="alignment-trimming",
            output_name="trimming_summary",
            artifact_kind="trimmed-alignment-retained-sites",
        )
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    report = EngineWorkflowReport(
        workflow="alignment-trimming",
        engine_name="trimal",
        input_paths=[input_path],
        output_paths={
            "trimmed_alignment": out_path,
            "trimming_summary": trimming_summary_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([input_path]),
        output_checksums={},
        config={
            "mode": mode,
            "gap_threshold": gap_threshold,
            "timeout_seconds": timeout_seconds,
        },
        trimming_summary=trimming_summary,
        notes=[
            f"trimal trimming mode: {mode}",
            f"retained sites: {trimming_summary.retained_site_count} of {trimming_summary.input_alignment_length}",
            f"gap percentage: {trimming_summary.input_gap_percentage:.3f} -> {trimming_summary.trimmed_gap_percentage:.3f}",
            "trimmed alignment validated as nonempty equal-length FASTA",
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)
