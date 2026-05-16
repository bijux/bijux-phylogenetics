from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
from pathlib import Path
import re

from bijux_phylogenetics.compare.reports import (
    ComparisonReportBuildResult,
    build_tree_comparison_report,
)
from bijux_phylogenetics.core.alignment import (
    AlignmentAlphabet,
    AlignmentRecord,
    AlignmentSummary,
    CodingSequenceExclusion,
)
from bijux_phylogenetics.core.partitions import (
    LocusPartition,
    PartitionSummaryReport,
    build_partition_summary_report,
    normalize_partition_data_type,
    parse_locus_partitions,
    slice_partition_sequence,
    write_locus_partitions,
    write_partition_summary_table,
)
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.errors import EngineWorkflowError, PhylogeneticsError
from bijux_phylogenetics.io.fasta import (
    back_translate_aligned_coding_sequences,
    classify_sequence_coding_behavior,
    infer_alignment_alphabet,
    load_fasta_alignment,
    prepare_coding_sequences_for_alignment,
    summarise_fasta,
    translate_prepared_coding_sequences,
    write_fasta_alignment,
)
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.trees import load_tree_set

from .bootstrap_artifacts import (
    build_bootstrap_support_histogram_rows,
    build_bootstrap_support_rows,
    build_low_support_bootstrap_rows,
    write_bootstrap_support_histogram,
    write_bootstrap_support_table,
)
from .common import (
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
    read_engine_version,
    resolve_engine_executable,
    update_incomplete_engine_run,
    validate_timeout_seconds,
    write_engine_manifest,
)
from .fasttree_artifacts import (
    build_fasttree_low_support_rows,
    build_fasttree_support_histogram_rows,
    build_fasttree_support_rows,
    write_fasttree_support_histogram,
    write_fasttree_support_table,
)
from .iqtree_artifacts import (
    IqtreeModelCandidate,
    IqtreeModelSelectionSummary,
    parse_best_model_file,
    parse_iqtree_model_selection_summary,
    parse_log_likelihood_file,
    resolve_iqtree_model_sidecar,
    write_iqtree_model_candidates_table,
)
from .sh_alrt_artifacts import (
    build_conflicting_sh_alrt_support_rows,
    build_sh_alrt_support_rows,
    write_sh_alrt_support_table,
)
from .validation import (
    BootstrapSupportNode,
    BootstrapSupportSummaryReport,
    FastTreeSupportNode,
    FastTreeSupportSummaryReport,
    ShAlrtSupportNode,
    ShAlrtSupportSummaryReport,
    WeakBackboneReport,
    detect_weakly_supported_backbone,
    summarize_bootstrap_support_distribution,
    summarize_fasttree_support_distribution,
    summarize_sh_alrt_support_distribution,
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
_MINIMUM_UFBOOT_REPLICATES = 1000
_INCOMPLETE_RUN_POLICIES = {"reject", "clean"}


@dataclass(slots=True)
class EngineWorkflowReport:
    workflow: str
    engine_name: str
    input_paths: list[Path]
    output_paths: dict[str, Path]
    run: EngineRunReport
    manifest_path: Path
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    config: dict[str, object] = field(default_factory=dict)
    selected_model: str | None = None
    log_likelihood: float | None = None
    iqtree_summary: IqtreeWorkflowSummary | None = None
    model_selection_summary: IqtreeModelSelectionSummary | None = None
    bootstrap_support_summary: BootstrapSupportSummaryReport | None = None
    fasttree_support_summary: FastTreeSupportSummaryReport | None = None
    sh_alrt_support_summary: ShAlrtSupportSummaryReport | None = None
    weak_backbone_report: WeakBackboneReport | None = None
    trimming_summary: AlignmentTrimmingSummary | None = None
    resumed: bool = False
    notes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class IqtreeSupportValue:
    node: str
    descendant_taxa: list[str]
    support: float
    support_fraction: float
    is_backbone: bool


@dataclass(slots=True)
class IqtreeWorkflowSummary:
    selected_model: str | None
    log_likelihood: float | None
    support_value_count: int
    minimum_support: float | None
    maximum_support: float | None
    support_values: list[IqtreeSupportValue] = field(default_factory=list)


@dataclass(slots=True)
class AlignmentTrimmingSummary:
    mode: str
    gap_threshold: float | None
    input_alignment_length: int
    trimmed_alignment_length: int
    retained_site_count: int
    removed_site_count: int
    retained_site_fraction: float
    removed_site_fraction: float
    input_gap_fraction: float
    trimmed_gap_fraction: float
    input_gap_percentage: float
    trimmed_gap_percentage: float


@dataclass(slots=True)
class ExternalTreeComparisonReport:
    fast_tree_path: Path
    ml_tree_path: Path
    comparison_report: ComparisonReportBuildResult


@dataclass(slots=True)
class CodonAwareAlignmentWorkflowReport:
    workflow: str
    engine_name: str
    input_path: Path
    output_paths: dict[str, Path]
    run: EngineRunReport
    manifest_path: Path
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    config: dict[str, object]
    sequence_type: AlignmentAlphabet
    genetic_code_id: int
    genetic_code_name: str
    input_sequence_count: int
    accepted_sequence_count: int
    invalid_codon_sequence_count: int
    excluded_sequences: list[CodingSequenceExclusion]
    terminal_stop_sequence_count: int
    notes: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    resumed: bool = False


@dataclass(frozen=True, slots=True)
class _PreparedIqtreePartitions:
    command_args: list[str]
    summary: PartitionSummaryReport
    output_paths: dict[str, Path]
    notes: list[str]
    mixed_data_types: bool


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


def _sidecar(path: Path, label: str) -> Path:
    return path.parent / f"{path.name}.{label}"


def _prefix_path(out_dir: Path, prefix: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / prefix


def _manifest_path_from_output(path: Path) -> Path:
    return _sidecar(path, "manifest.json")


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
    raise EngineWorkflowError(
        "a previous engine run left incomplete outputs and resume could not safely "
        f"reuse them; marker: {marker_path}"
    )


def _partition_support_path(prefix_path: Path, suffix: str) -> Path:
    return prefix_path.parent / f"{prefix_path.name}.{suffix}"


def _record_output_validation_failure(
    manifest_path: Path,
    run: EngineRunReport,
    error: PhylogeneticsError,
) -> None:
    update_incomplete_engine_run(
        manifest_path,
        ended_at_utc=run.ended_at_utc,
        timed_out=run.timed_out,
        exit_code=run.exit_code,
        failure_message=(
            f"{run.engine_name} {run.workflow} produced outputs that failed "
            f"validation: {error.code}"
        ),
    )


def _require_nonempty_text_output(
    path: Path,
    *,
    engine_name: str,
    workflow: str,
    output_name: str,
    artifact_kind: str,
) -> str:
    if not path.exists():
        raise build_engine_output_error(
            f"{engine_name} {workflow} did not produce required output '{output_name}': {path}",
            code="engine_required_output_missing",
            engine_name=engine_name,
            workflow=workflow,
            path=path,
            output_name=output_name,
            artifact_kind=artifact_kind,
        )
    text = path.read_text(encoding="utf-8", errors="replace")
    if not text.strip():
        raise build_engine_output_error(
            f"{engine_name} {workflow} produced an empty required output '{output_name}': {path}",
            code="engine_output_empty",
            engine_name=engine_name,
            workflow=workflow,
            path=path,
            output_name=output_name,
            artifact_kind=artifact_kind,
        )
    return text


def _validate_alignment_output(
    path: Path,
    *,
    engine_name: str,
    workflow: str,
    output_name: str,
    artifact_kind: str,
) -> AlignmentSummary:
    _require_nonempty_text_output(
        path,
        engine_name=engine_name,
        workflow=workflow,
        output_name=output_name,
        artifact_kind=artifact_kind,
    )
    records = load_fasta_alignment(path)
    if not records or len(records[0].sequence) < 1:
        raise build_engine_output_error(
            f"{engine_name} {workflow} produced an empty alignment after validation: {path}",
            code="engine_output_empty",
            engine_name=engine_name,
            workflow=workflow,
            path=path,
            output_name=output_name,
            artifact_kind=artifact_kind,
        )
    return summarise_fasta(path)


def _write_coding_exclusion_table(
    path: Path, exclusions: list[CodingSequenceExclusion]
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "identifier",
                "comparable_length",
                "reason",
                "invalid_codon_count",
                "premature_stop_count",
                "terminal_stop_count",
                "trailing_bases",
                "note",
            ]
        )
    ]
    lines.extend(
        "\t".join(
            [
                row.identifier,
                str(row.comparable_length),
                row.reason,
                str(row.invalid_codon_count),
                str(row.premature_stop_count),
                str(row.terminal_stop_count),
                str(row.trailing_bases),
                row.note,
            ]
        )
        for row in exclusions
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_coding_summary_table(
    path: Path,
    *,
    input_path: Path,
    genetic_code: int,
    exclusions: list[CodingSequenceExclusion],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    exclusion_by_identifier = {row.identifier: row for row in exclusions}
    behaviors = classify_sequence_coding_behavior(
        input_path,
        genetic_code=genetic_code,
    )
    header = "\t".join(
        [
            "identifier",
            "status",
            "comparable_length",
            "divisible_by_three",
            "invalid_codon_count",
            "premature_stop_count",
            "terminal_stop_count",
            "exclusion_reason",
            "note",
        ]
    )
    rows = [header]
    for behavior in behaviors:
        exclusion = exclusion_by_identifier.get(behavior.identifier)
        rows.append(
            "\t".join(
                [
                    behavior.identifier,
                    "excluded" if exclusion is not None else "accepted",
                    str(behavior.comparable_length),
                    "yes" if behavior.divisible_by_three else "no",
                    str(behavior.invalid_codon_count),
                    str(behavior.premature_stop_count),
                    str(behavior.terminal_stop_count),
                    "" if exclusion is None else exclusion.reason,
                    behavior.note,
                ]
            )
        )
    ordered_rows = [rows[0], *sorted(rows[1:])]
    path.write_text("\n".join(ordered_rows) + "\n", encoding="utf-8")
    return path


def _build_alignment_trimming_summary(
    *,
    mode: str,
    gap_threshold: float,
    input_summary: AlignmentSummary,
    trimmed_summary: AlignmentSummary,
) -> AlignmentTrimmingSummary:
    if trimmed_summary.alignment_length > input_summary.alignment_length:
        raise EngineWorkflowError(
            "trimmed alignment is longer than the input alignment, which is not a valid trimAl result"
        )
    retained_site_count = trimmed_summary.alignment_length
    removed_site_count = (
        input_summary.alignment_length - trimmed_summary.alignment_length
    )
    retained_site_fraction = retained_site_count / input_summary.alignment_length
    removed_site_fraction = removed_site_count / input_summary.alignment_length
    return AlignmentTrimmingSummary(
        mode=mode,
        gap_threshold=gap_threshold if mode == "gap-threshold" else None,
        input_alignment_length=input_summary.alignment_length,
        trimmed_alignment_length=trimmed_summary.alignment_length,
        retained_site_count=retained_site_count,
        removed_site_count=removed_site_count,
        retained_site_fraction=retained_site_fraction,
        removed_site_fraction=removed_site_fraction,
        input_gap_fraction=input_summary.gap_fraction,
        trimmed_gap_fraction=trimmed_summary.gap_fraction,
        input_gap_percentage=input_summary.gap_fraction * 100.0,
        trimmed_gap_percentage=trimmed_summary.gap_fraction * 100.0,
    )


def _validate_tree_output(
    path: Path,
    *,
    engine_name: str,
    workflow: str,
    output_name: str,
    artifact_kind: str,
) -> None:
    tree_text = _require_nonempty_text_output(
        path,
        engine_name=engine_name,
        workflow=workflow,
        output_name=output_name,
        artifact_kind=artifact_kind,
    )
    loads_newick(tree_text)
    validate_tree_path(path)


def _validate_tree_set_output(
    path: Path,
    *,
    engine_name: str,
    workflow: str,
    output_name: str,
    artifact_kind: str,
) -> None:
    _require_nonempty_text_output(
        path,
        engine_name=engine_name,
        workflow=workflow,
        output_name=output_name,
        artifact_kind=artifact_kind,
    )
    load_tree_set(path)


def _validate_iqtree_required_artifacts(
    prefix_path: Path,
    *,
    workflow: str,
) -> None:
    _require_nonempty_text_output(
        prefix_path.with_suffix(".iqtree"),
        engine_name="iqtree",
        workflow=workflow,
        output_name="iqtree_report",
        artifact_kind="iqtree-report",
    )
    _require_nonempty_text_output(
        prefix_path.with_suffix(".log"),
        engine_name="iqtree",
        workflow=workflow,
        output_name="iqtree_log",
        artifact_kind="iqtree-log",
    )


def _validate_iqtree_model_result(
    prefix_path: Path,
    *,
    workflow: str,
    default_selected_model: str | None = None,
) -> str:
    selected_model = _parse_best_model_artifact(prefix_path)
    if selected_model is None and default_selected_model is not None:
        return default_selected_model
    if selected_model is None:
        raise build_engine_output_error(
            f"iqtree {workflow} did not expose a parsable best-fit model result",
            code="engine_model_result_missing",
            engine_name="iqtree",
            workflow=workflow,
            path=prefix_path.with_suffix(".iqtree"),
            output_name="iqtree_report",
            artifact_kind="iqtree-model-result",
            details={
                "model_sidecar_path": (
                    None
                    if resolve_iqtree_model_sidecar(prefix_path) is None
                    else str(resolve_iqtree_model_sidecar(prefix_path))
                )
            },
        )
    return selected_model


def _validate_support_value_count(
    *,
    engine_name: str,
    workflow: str,
    path: Path,
    output_name: str,
    artifact_kind: str,
    support_value_count: int,
    support_kind: str,
) -> None:
    if support_value_count > 0:
        return
    raise build_engine_output_error(
        f"{engine_name} {workflow} did not expose any parsable {support_kind} values in '{output_name}': {path}",
        code="engine_support_values_missing",
        engine_name=engine_name,
        workflow=workflow,
        path=path,
        output_name=output_name,
        artifact_kind=artifact_kind,
        details={"support_kind": support_kind},
    )


def _ensure_inference_ready_alignment(path: Path) -> None:
    records = load_fasta_alignment(path)
    if not records or len(records[0].sequence) < 1:
        raise EngineWorkflowError(
            f"inference alignment is empty after filtering: {path}"
        )
    summarise_fasta(path)


def _partition_alignment_file_name(partition: LocusPartition) -> str:
    normalized = re.sub(r"[^A-Za-z0-9._-]+", "-", partition.name.strip().lower())
    normalized = normalized.strip("-._") or "partition"
    digest = hashlib.sha1(partition.name.encode("utf-8")).hexdigest()[:8]  # nosec B324
    return f"{normalized}-{digest}.fasta"


def _iqtree_partition_supports_fixed_model(
    *,
    model: str,
    mixed_data_types: bool,
) -> bool:
    if not mixed_data_types:
        return True
    normalized = model.strip().upper()
    return normalized in {
        "TEST",
        "TESTONLY",
        "TESTNEWONLY",
        "MF",
        "MFP",
        "TESTMERGE",
        "TESTMERGEONLY",
        "MF+MERGE",
        "MFP+MERGE",
    }


def _prepare_iqtree_partitions(
    input_path: Path,
    partition_path: Path,
    *,
    prefix_path: Path,
) -> _PreparedIqtreePartitions:
    records = load_fasta_alignment(input_path)
    alignment_summary = summarise_fasta(input_path)
    partitions = parse_locus_partitions(partition_path)
    summary = build_partition_summary_report(
        partitions,
        alignment_length=alignment_summary.alignment_length,
    )
    summary_path = _partition_support_path(prefix_path, "partition-summary.tsv")
    write_partition_summary_table(summary_path, summary)
    notes = [
        f"validated {summary.partition_count} partitions across {summary.assigned_site_count} assigned sites",
    ]
    output_paths: dict[str, Path] = {
        "partition_summary": summary_path,
    }

    declared_types = {
        normalize_partition_data_type(partition.data_type)
        for partition in partitions
        if partition.data_type is not None
    }
    if len(declared_types) <= 1:
        normalized_partition_path = _partition_support_path(
            prefix_path, "partition-scheme.partitions"
        )
        write_locus_partitions(normalized_partition_path, partitions)
        output_paths["partition_scheme"] = normalized_partition_path
        notes.append(
            "prepared a normalized partition scheme for single-alignment IQ-TREE analysis"
        )
        return _PreparedIqtreePartitions(
            command_args=[
                "-s",
                str(input_path.resolve()),
                "-p",
                str(normalized_partition_path.resolve()),
            ],
            summary=summary,
            output_paths=output_paths,
            notes=notes,
            mixed_data_types=False,
        )

    if any(partition.data_type is None for partition in partitions):
        raise EngineWorkflowError(
            "mixed partition analyses require every partition to declare a data_type"
        )
    unsupported_types = sorted(
        {
            data_type
            for data_type in declared_types
            if data_type is not None and data_type not in {"DNA", "RNA", "PROTEIN"}
        }
    )
    if unsupported_types:
        raise EngineWorkflowError(
            "mixed partition analyses currently support only DNA, RNA, and PROTEIN datatypes; "
            f"got: {', '.join(unsupported_types)}"
        )

    partition_alignment_dir = _partition_support_path(
        prefix_path, "partition-alignments"
    )
    partition_alignment_dir.mkdir(parents=True, exist_ok=True)
    lines = ["#nexus", "begin sets;"]
    for partition in partitions:
        partition_alignment_path = (
            partition_alignment_dir / _partition_alignment_file_name(partition)
        )
        write_fasta_alignment(
            partition_alignment_path,
            [
                AlignmentRecord(
                    identifier=record.identifier,
                    sequence=slice_partition_sequence(record.sequence, partition),
                )
                for record in records
            ],
        )
        output_paths[f"partition_alignment_{partition.name}"] = partition_alignment_path
        lines.append(
            f"    charset {partition.name} = {partition_alignment_path.name}: *;"
        )
    lines.append("end;")
    mixed_partition_path = _partition_support_path(prefix_path, "partition-scheme.nex")
    mixed_partition_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    output_paths["partition_scheme"] = mixed_partition_path
    notes.append(
        "prepared a mixed-datatype NEXUS partition scheme with one extracted alignment per partition"
    )
    return _PreparedIqtreePartitions(
        command_args=["-p", str(mixed_partition_path.resolve())],
        summary=summary,
        output_paths=output_paths,
        notes=notes,
        mixed_data_types=True,
    )


def _restore_workflow_report(payload: dict[str, object]) -> EngineWorkflowReport:
    run_payload = payload["run"]
    if not isinstance(run_payload, dict):
        raise EngineWorkflowError("engine manifest contains an invalid run payload")
    version_payload = run_payload["version"]
    if not isinstance(version_payload, dict):
        raise EngineWorkflowError("engine manifest contains an invalid version payload")
    version = EngineVersionInfo(
        engine_name=str(version_payload["engine_name"]),
        executable=str(version_payload["executable"]),
        command=[str(item) for item in version_payload["command"]],
        text=str(version_payload["text"]),
    )
    run = EngineRunReport(
        engine_name=str(run_payload["engine_name"]),
        workflow=str(run_payload["workflow"]),
        executable=str(run_payload["executable"]),
        working_directory=Path(run_payload["working_directory"]),
        version=version,
        command=[str(item) for item in run_payload["command"]],
        exit_code=int(run_payload["exit_code"]),
        stdout_path=Path(run_payload["stdout_path"]),
        stderr_path=Path(run_payload["stderr_path"]),
        output_paths={
            str(key): Path(value)
            for key, value in dict(run_payload["output_paths"]).items()
        },
        warning_lines=[str(item) for item in run_payload["warning_lines"]],
        started_at_utc=str(run_payload.get("started_at_utc", "")),
        ended_at_utc=str(run_payload.get("ended_at_utc", "")),
        runtime_seconds=float(run_payload.get("runtime_seconds", 0.0)),
        timeout_seconds=(
            None
            if run_payload.get("timeout_seconds") is None
            else float(run_payload["timeout_seconds"])
        ),
        timed_out=bool(run_payload.get("timed_out", False)),
    )
    return EngineWorkflowReport(
        workflow=str(payload["workflow"]),
        engine_name=str(payload["engine_name"]),
        input_paths=[Path(item) for item in payload["input_paths"]],
        output_paths={
            str(key): Path(value)
            for key, value in dict(payload["output_paths"]).items()
        },
        run=run,
        manifest_path=Path(payload["manifest_path"]),
        input_checksums={
            str(key): str(value)
            for key, value in dict(payload.get("input_checksums", {})).items()
        },
        output_checksums={
            str(key): str(value)
            for key, value in dict(payload.get("output_checksums", {})).items()
        },
        config={
            str(key): value for key, value in dict(payload.get("config", {})).items()
        },
        selected_model=None
        if payload.get("selected_model") is None
        else str(payload["selected_model"]),
        log_likelihood=(
            None
            if payload.get("log_likelihood") is None
            else float(payload["log_likelihood"])
        ),
        iqtree_summary=(
            None
            if payload.get("iqtree_summary") is None
            else IqtreeWorkflowSummary(
                selected_model=(
                    None
                    if dict(payload["iqtree_summary"]).get("selected_model") is None
                    else str(dict(payload["iqtree_summary"])["selected_model"])
                ),
                log_likelihood=(
                    None
                    if dict(payload["iqtree_summary"]).get("log_likelihood") is None
                    else float(dict(payload["iqtree_summary"])["log_likelihood"])
                ),
                support_value_count=int(
                    dict(payload["iqtree_summary"])["support_value_count"]
                ),
                minimum_support=(
                    None
                    if dict(payload["iqtree_summary"]).get("minimum_support") is None
                    else float(dict(payload["iqtree_summary"])["minimum_support"])
                ),
                maximum_support=(
                    None
                    if dict(payload["iqtree_summary"]).get("maximum_support") is None
                    else float(dict(payload["iqtree_summary"])["maximum_support"])
                ),
                support_values=[
                    IqtreeSupportValue(
                        node=str(dict(item)["node"]),
                        descendant_taxa=[
                            str(taxon) for taxon in list(dict(item)["descendant_taxa"])
                        ],
                        support=float(dict(item)["support"]),
                        support_fraction=float(dict(item)["support_fraction"]),
                        is_backbone=bool(dict(item)["is_backbone"]),
                    )
                    for item in list(dict(payload["iqtree_summary"])["support_values"])
                ],
            )
        ),
        model_selection_summary=(
            None
            if payload.get("model_selection_summary") is None
            else IqtreeModelSelectionSummary(
                selected_model=(
                    None
                    if dict(payload["model_selection_summary"]).get("selected_model")
                    is None
                    else str(dict(payload["model_selection_summary"])["selected_model"])
                ),
                selected_criterion=(
                    None
                    if dict(payload["model_selection_summary"]).get(
                        "selected_criterion"
                    )
                    is None
                    else str(
                        dict(payload["model_selection_summary"])["selected_criterion"]
                    )
                ),
                best_model_aic=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_model_aic")
                    is None
                    else str(dict(payload["model_selection_summary"])["best_model_aic"])
                ),
                best_model_aicc=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_model_aicc")
                    is None
                    else str(
                        dict(payload["model_selection_summary"])["best_model_aicc"]
                    )
                ),
                best_model_bic=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_model_bic")
                    is None
                    else str(dict(payload["model_selection_summary"])["best_model_bic"])
                ),
                best_score_aic=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_score_aic")
                    is None
                    else float(
                        dict(payload["model_selection_summary"])["best_score_aic"]
                    )
                ),
                best_score_aicc=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_score_aicc")
                    is None
                    else float(
                        dict(payload["model_selection_summary"])["best_score_aicc"]
                    )
                ),
                best_score_bic=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_score_bic")
                    is None
                    else float(
                        dict(payload["model_selection_summary"])["best_score_bic"]
                    )
                ),
                candidate_count=int(
                    dict(payload["model_selection_summary"])["candidate_count"]
                ),
                candidates=[
                    IqtreeModelCandidate(
                        rank=int(dict(item)["rank"]),
                        model=str(dict(item)["model"]),
                        log_likelihood=float(dict(item)["log_likelihood"]),
                        parameter_count=(
                            None
                            if dict(item).get("parameter_count") is None
                            else int(dict(item)["parameter_count"])
                        ),
                        aic=float(dict(item)["aic"]),
                        aicc=float(dict(item)["aicc"]),
                        bic=float(dict(item)["bic"]),
                    )
                    for item in list(
                        dict(payload["model_selection_summary"])["candidates"]
                    )
                ],
                bic_near_best_models=[
                    str(item)
                    for item in list(
                        dict(payload["model_selection_summary"])["bic_near_best_models"]
                    )
                ],
            )
        ),
        bootstrap_support_summary=(
            None
            if payload.get("bootstrap_support_summary") is None
            else BootstrapSupportSummaryReport(
                tree_path=Path(dict(payload["bootstrap_support_summary"])["tree_path"]),
                internal_node_count=int(
                    dict(payload["bootstrap_support_summary"])["internal_node_count"]
                ),
                supported_node_count=int(
                    dict(payload["bootstrap_support_summary"])["supported_node_count"]
                ),
                minimum_support=(
                    None
                    if dict(payload["bootstrap_support_summary"]).get("minimum_support")
                    is None
                    else float(
                        dict(payload["bootstrap_support_summary"])["minimum_support"]
                    )
                ),
                maximum_support=(
                    None
                    if dict(payload["bootstrap_support_summary"]).get("maximum_support")
                    is None
                    else float(
                        dict(payload["bootstrap_support_summary"])["maximum_support"]
                    )
                ),
                median_support=(
                    None
                    if dict(payload["bootstrap_support_summary"]).get("median_support")
                    is None
                    else float(
                        dict(payload["bootstrap_support_summary"])["median_support"]
                    )
                ),
                weakly_supported_clade_count=int(
                    dict(payload["bootstrap_support_summary"])[
                        "weakly_supported_clade_count"
                    ]
                ),
                support_histogram={
                    str(key): int(value)
                    for key, value in dict(
                        dict(payload["bootstrap_support_summary"])["support_histogram"]
                    ).items()
                },
                nodes=[
                    BootstrapSupportNode(
                        node=str(dict(item)["node"]),
                        descendant_taxa=[
                            str(taxon) for taxon in list(dict(item)["descendant_taxa"])
                        ],
                        support=float(dict(item)["support"]),
                        support_fraction=float(dict(item)["support_fraction"]),
                        is_backbone=bool(dict(item)["is_backbone"]),
                    )
                    for item in list(
                        dict(payload["bootstrap_support_summary"])["nodes"]
                    )
                ],
                warnings=[
                    str(item)
                    for item in list(
                        dict(payload["bootstrap_support_summary"])["warnings"]
                    )
                ],
            )
        ),
        fasttree_support_summary=(
            None
            if payload.get("fasttree_support_summary") is None
            else FastTreeSupportSummaryReport(
                tree_path=Path(dict(payload["fasttree_support_summary"])["tree_path"]),
                internal_node_count=int(
                    dict(payload["fasttree_support_summary"])["internal_node_count"]
                ),
                annotated_node_count=int(
                    dict(payload["fasttree_support_summary"])["annotated_node_count"]
                ),
                minimum_local_support=(
                    None
                    if dict(payload["fasttree_support_summary"]).get(
                        "minimum_local_support"
                    )
                    is None
                    else float(
                        dict(payload["fasttree_support_summary"])[
                            "minimum_local_support"
                        ]
                    )
                ),
                maximum_local_support=(
                    None
                    if dict(payload["fasttree_support_summary"]).get(
                        "maximum_local_support"
                    )
                    is None
                    else float(
                        dict(payload["fasttree_support_summary"])[
                            "maximum_local_support"
                        ]
                    )
                ),
                median_local_support=(
                    None
                    if dict(payload["fasttree_support_summary"]).get(
                        "median_local_support"
                    )
                    is None
                    else float(
                        dict(payload["fasttree_support_summary"])[
                            "median_local_support"
                        ]
                    )
                ),
                weakly_supported_clade_count=int(
                    dict(payload["fasttree_support_summary"])[
                        "weakly_supported_clade_count"
                    ]
                ),
                support_histogram={
                    str(key): int(value)
                    for key, value in dict(
                        dict(payload["fasttree_support_summary"])["support_histogram"]
                    ).items()
                },
                approximate_method=bool(
                    dict(payload["fasttree_support_summary"])["approximate_method"]
                ),
                support_label_kind=str(
                    dict(payload["fasttree_support_summary"])["support_label_kind"]
                ),
                support_scale=str(
                    dict(payload["fasttree_support_summary"])["support_scale"]
                ),
                nodes=[
                    FastTreeSupportNode(
                        node=str(dict(item)["node"]),
                        descendant_taxa=[
                            str(taxon) for taxon in list(dict(item)["descendant_taxa"])
                        ],
                        local_support=float(dict(item)["local_support"]),
                        support_fraction=float(dict(item)["support_fraction"]),
                        is_backbone=bool(dict(item)["is_backbone"]),
                    )
                    for item in list(dict(payload["fasttree_support_summary"])["nodes"])
                ],
                warnings=[
                    str(item)
                    for item in list(
                        dict(payload["fasttree_support_summary"])["warnings"]
                    )
                ],
            )
        ),
        sh_alrt_support_summary=(
            None
            if payload.get("sh_alrt_support_summary") is None
            else ShAlrtSupportSummaryReport(
                tree_path=Path(dict(payload["sh_alrt_support_summary"])["tree_path"]),
                internal_node_count=int(
                    dict(payload["sh_alrt_support_summary"])["internal_node_count"]
                ),
                annotated_node_count=int(
                    dict(payload["sh_alrt_support_summary"])["annotated_node_count"]
                ),
                fully_scored_node_count=int(
                    dict(payload["sh_alrt_support_summary"])["fully_scored_node_count"]
                ),
                minimum_sh_alrt_support=(
                    None
                    if dict(payload["sh_alrt_support_summary"]).get(
                        "minimum_sh_alrt_support"
                    )
                    is None
                    else float(
                        dict(payload["sh_alrt_support_summary"])[
                            "minimum_sh_alrt_support"
                        ]
                    )
                ),
                maximum_sh_alrt_support=(
                    None
                    if dict(payload["sh_alrt_support_summary"]).get(
                        "maximum_sh_alrt_support"
                    )
                    is None
                    else float(
                        dict(payload["sh_alrt_support_summary"])[
                            "maximum_sh_alrt_support"
                        ]
                    )
                ),
                minimum_ufboot_support=(
                    None
                    if dict(payload["sh_alrt_support_summary"]).get(
                        "minimum_ufboot_support"
                    )
                    is None
                    else float(
                        dict(payload["sh_alrt_support_summary"])[
                            "minimum_ufboot_support"
                        ]
                    )
                ),
                maximum_ufboot_support=(
                    None
                    if dict(payload["sh_alrt_support_summary"]).get(
                        "maximum_ufboot_support"
                    )
                    is None
                    else float(
                        dict(payload["sh_alrt_support_summary"])[
                            "maximum_ufboot_support"
                        ]
                    )
                ),
                weak_sh_alrt_clade_count=int(
                    dict(payload["sh_alrt_support_summary"])["weak_sh_alrt_clade_count"]
                ),
                weak_ufboot_clade_count=int(
                    dict(payload["sh_alrt_support_summary"])["weak_ufboot_clade_count"]
                ),
                conflicting_support_signal_count=int(
                    dict(payload["sh_alrt_support_summary"])[
                        "conflicting_support_signal_count"
                    ]
                ),
                nodes=[
                    ShAlrtSupportNode(
                        node=str(dict(item)["node"]),
                        descendant_taxa=[
                            str(taxon) for taxon in list(dict(item)["descendant_taxa"])
                        ],
                        sh_alrt_support=(
                            None
                            if dict(item).get("sh_alrt_support") is None
                            else float(dict(item)["sh_alrt_support"])
                        ),
                        sh_alrt_support_fraction=(
                            None
                            if dict(item).get("sh_alrt_support_fraction") is None
                            else float(dict(item)["sh_alrt_support_fraction"])
                        ),
                        ufboot_support=(
                            None
                            if dict(item).get("ufboot_support") is None
                            else float(dict(item)["ufboot_support"])
                        ),
                        ufboot_support_fraction=(
                            None
                            if dict(item).get("ufboot_support_fraction") is None
                            else float(dict(item)["ufboot_support_fraction"])
                        ),
                        is_backbone=bool(dict(item)["is_backbone"]),
                        sh_alrt_strong=bool(dict(item)["sh_alrt_strong"]),
                        ufboot_strong=bool(dict(item)["ufboot_strong"]),
                        conflicting_support_signal=bool(
                            dict(item)["conflicting_support_signal"]
                        ),
                        support_agreement=str(dict(item)["support_agreement"]),
                    )
                    for item in list(dict(payload["sh_alrt_support_summary"])["nodes"])
                ],
                warnings=[
                    str(item)
                    for item in list(
                        dict(payload["sh_alrt_support_summary"])["warnings"]
                    )
                ],
            )
        ),
        weak_backbone_report=(
            None
            if payload.get("weak_backbone_report") is None
            else WeakBackboneReport(
                tree_path=Path(dict(payload["weak_backbone_report"])["tree_path"]),
                threshold=float(dict(payload["weak_backbone_report"])["threshold"]),
                evaluated_backbone_node_count=int(
                    dict(payload["weak_backbone_report"])[
                        "evaluated_backbone_node_count"
                    ]
                ),
                weak_backbone_node_count=int(
                    dict(payload["weak_backbone_report"])["weak_backbone_node_count"]
                ),
                weak_nodes=[
                    BootstrapSupportNode(
                        node=str(dict(item)["node"]),
                        descendant_taxa=[
                            str(taxon) for taxon in list(dict(item)["descendant_taxa"])
                        ],
                        support=float(dict(item)["support"]),
                        support_fraction=float(dict(item)["support_fraction"]),
                        is_backbone=bool(dict(item)["is_backbone"]),
                    )
                    for item in list(
                        dict(payload["weak_backbone_report"])["weak_nodes"]
                    )
                ],
                warnings=[
                    str(item)
                    for item in list(dict(payload["weak_backbone_report"])["warnings"])
                ],
            )
        ),
        trimming_summary=None
        if payload.get("trimming_summary") is None
        else AlignmentTrimmingSummary(
            mode=str(dict(payload["trimming_summary"])["mode"]),
            gap_threshold=(
                None
                if dict(payload["trimming_summary"]).get("gap_threshold") is None
                else float(dict(payload["trimming_summary"])["gap_threshold"])
            ),
            input_alignment_length=int(
                dict(payload["trimming_summary"])["input_alignment_length"]
            ),
            trimmed_alignment_length=int(
                dict(payload["trimming_summary"])["trimmed_alignment_length"]
            ),
            retained_site_count=int(
                dict(payload["trimming_summary"])["retained_site_count"]
            ),
            removed_site_count=int(
                dict(payload["trimming_summary"])["removed_site_count"]
            ),
            retained_site_fraction=float(
                dict(payload["trimming_summary"])["retained_site_fraction"]
            ),
            removed_site_fraction=float(
                dict(payload["trimming_summary"])["removed_site_fraction"]
            ),
            input_gap_fraction=float(
                dict(payload["trimming_summary"])["input_gap_fraction"]
            ),
            trimmed_gap_fraction=float(
                dict(payload["trimming_summary"])["trimmed_gap_fraction"]
            ),
            input_gap_percentage=float(
                dict(payload["trimming_summary"])["input_gap_percentage"]
            ),
            trimmed_gap_percentage=float(
                dict(payload["trimming_summary"])["trimmed_gap_percentage"]
            ),
        ),
        resumed=bool(payload.get("resumed", False)),
        notes=[str(item) for item in payload.get("notes", [])],
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


def _restore_codon_aware_alignment_report(
    payload: dict[str, object],
) -> CodonAwareAlignmentWorkflowReport:
    run_payload = payload["run"]
    if not isinstance(run_payload, dict):
        raise EngineWorkflowError("engine manifest contains an invalid run payload")
    version_payload = run_payload["version"]
    if not isinstance(version_payload, dict):
        raise EngineWorkflowError("engine manifest contains an invalid version payload")
    version = EngineVersionInfo(
        engine_name=str(version_payload["engine_name"]),
        executable=str(version_payload["executable"]),
        command=[str(item) for item in version_payload["command"]],
        text=str(version_payload["text"]),
    )
    run = EngineRunReport(
        engine_name=str(run_payload["engine_name"]),
        workflow=str(run_payload["workflow"]),
        executable=str(run_payload["executable"]),
        working_directory=Path(run_payload["working_directory"]),
        version=version,
        command=[str(item) for item in run_payload["command"]],
        exit_code=int(run_payload["exit_code"]),
        stdout_path=Path(run_payload["stdout_path"]),
        stderr_path=Path(run_payload["stderr_path"]),
        output_paths={
            str(key): Path(value)
            for key, value in dict(run_payload["output_paths"]).items()
        },
        warning_lines=[str(item) for item in run_payload["warning_lines"]],
        started_at_utc=str(run_payload.get("started_at_utc", "")),
        ended_at_utc=str(run_payload.get("ended_at_utc", "")),
        runtime_seconds=float(run_payload.get("runtime_seconds", 0.0)),
        timeout_seconds=(
            None
            if run_payload.get("timeout_seconds") is None
            else float(run_payload["timeout_seconds"])
        ),
        timed_out=bool(run_payload.get("timed_out", False)),
    )
    return CodonAwareAlignmentWorkflowReport(
        workflow=str(payload["workflow"]),
        engine_name=str(payload["engine_name"]),
        input_path=Path(payload["input_path"]),
        output_paths={
            str(key): Path(value)
            for key, value in dict(payload["output_paths"]).items()
        },
        run=run,
        manifest_path=Path(payload["manifest_path"]),
        input_checksums={
            str(key): str(value)
            for key, value in dict(payload.get("input_checksums", {})).items()
        },
        output_checksums={
            str(key): str(value)
            for key, value in dict(payload.get("output_checksums", {})).items()
        },
        config={
            str(key): value for key, value in dict(payload.get("config", {})).items()
        },
        sequence_type=str(payload["sequence_type"]),
        genetic_code_id=int(payload.get("genetic_code_id", 1)),
        genetic_code_name=str(payload.get("genetic_code_name", "Standard")),
        input_sequence_count=int(payload.get("input_sequence_count", 0)),
        accepted_sequence_count=int(payload["accepted_sequence_count"]),
        invalid_codon_sequence_count=int(
            payload.get("invalid_codon_sequence_count", 0)
        ),
        excluded_sequences=[
            CodingSequenceExclusion(
                identifier=str(item["identifier"]),
                comparable_length=int(item["comparable_length"]),
                reason=str(item["reason"]),
                invalid_codon_count=int(item.get("invalid_codon_count", 0)),
                premature_stop_count=int(item["premature_stop_count"]),
                terminal_stop_count=int(item["terminal_stop_count"]),
                trailing_bases=int(item["trailing_bases"]),
                note=str(item["note"]),
            )
            for item in payload.get("excluded_sequences", [])
        ],
        terminal_stop_sequence_count=int(payload["terminal_stop_sequence_count"]),
        notes=[str(item) for item in payload.get("notes", [])],
        warnings=[str(item) for item in payload.get("warnings", [])],
        resumed=bool(payload.get("resumed", False)),
    )


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


def _iqtree_sequence_type_flag(
    path: Path, sequence_type: AlignmentAlphabet | None
) -> list[str]:
    detected = sequence_type
    if detected is None:
        detected = infer_alignment_alphabet(load_fasta_alignment(path))
    if detected in {"dna", "rna"}:
        return ["-st", "DNA"]
    if detected == "protein":
        return ["-st", "AA"]
    return []


def _fasttree_args(path: Path, sequence_type: AlignmentAlphabet | None) -> list[str]:
    detected = sequence_type
    if detected is None:
        detected = infer_alignment_alphabet(load_fasta_alignment(path))
    if detected in {"dna", "rna"}:
        return ["-gtr", "-nt", str(path)]
    if detected == "protein":
        return ["-lg", str(path)]
    return [str(path)]


def _iqtree_execution_controls(*, seed: int, threads: int) -> list[str]:
    if seed < 1:
        raise ValueError(f"iqtree seed must be positive, got {seed}")
    if threads < 1:
        raise ValueError(f"iqtree threads must be positive, got {threads}")
    return ["-seed", str(seed), "-nt", str(threads)]


def _validate_ufboot_replicates(replicates: int) -> None:
    if replicates < _MINIMUM_UFBOOT_REPLICATES:
        raise EngineWorkflowError(
            "iqtree ultrafast bootstrap requires at least "
            f"{_MINIMUM_UFBOOT_REPLICATES} replicates, got {replicates}"
        )


def _validate_sh_alrt_replicates(replicates: int) -> None:
    if replicates < 1:
        raise ValueError(f"sh-alrt replicates must be positive, got {replicates}")


def _parse_best_model_artifact(prefix_path: Path) -> str | None:
    for candidate in (
        prefix_path.with_suffix(".iqtree"),
        resolve_iqtree_model_sidecar(prefix_path),
    ):
        if candidate is None:
            continue
        model = parse_best_model_file(candidate)
        if model is not None:
            return model
    return None


def _parse_log_likelihood_artifact(prefix_path: Path) -> float | None:
    for candidate in (
        prefix_path.with_suffix(".iqtree"),
        prefix_path.with_suffix(".log"),
    ):
        if not candidate.exists():
            continue
        log_likelihood = parse_log_likelihood_file(candidate)
        if log_likelihood is not None:
            return log_likelihood
    return None


def _existing_iqtree_outputs(
    prefix_path: Path,
    *,
    include_tree: bool = False,
    include_bootstrap: bool = False,
    include_consensus: bool = False,
) -> dict[str, Path]:
    outputs: dict[str, Path] = {}
    for key, candidate in (
        ("iqtree_report", prefix_path.with_suffix(".iqtree")),
        ("iqtree_log", prefix_path.with_suffix(".log")),
    ):
        if candidate.exists():
            outputs[key] = candidate
    model_sidecar = resolve_iqtree_model_sidecar(prefix_path)
    if model_sidecar is not None:
        outputs["model_selection_sidecar"] = model_sidecar
    tree_candidate = prefix_path.with_suffix(".treefile")
    if include_tree and tree_candidate.exists():
        outputs["tree"] = tree_candidate
    bootstrap_candidate = prefix_path.with_suffix(".ufboot")
    if include_bootstrap and bootstrap_candidate.exists():
        outputs["bootstrap_trees"] = bootstrap_candidate
    consensus_candidate = prefix_path.with_suffix(".contree")
    if include_consensus and consensus_candidate.exists():
        outputs["consensus_tree"] = consensus_candidate
    return outputs


def _build_iqtree_summary(
    prefix_path: Path,
    *,
    default_selected_model: str | None,
    support_tree_path: Path | None = None,
) -> IqtreeWorkflowSummary:
    selected_model = _parse_best_model_artifact(prefix_path) or default_selected_model
    log_likelihood = _parse_log_likelihood_artifact(prefix_path)
    support_values: list[IqtreeSupportValue] = []
    minimum_support: float | None = None
    maximum_support: float | None = None
    if support_tree_path is not None and support_tree_path.exists():
        support_summary = summarize_bootstrap_support_distribution(support_tree_path)
        support_values = [
            IqtreeSupportValue(
                node=node.node,
                descendant_taxa=list(node.descendant_taxa),
                support=node.support,
                support_fraction=node.support_fraction,
                is_backbone=node.is_backbone,
            )
            for node in support_summary.nodes
        ]
        minimum_support = support_summary.minimum_support
        maximum_support = support_summary.maximum_support
    return IqtreeWorkflowSummary(
        selected_model=selected_model,
        log_likelihood=log_likelihood,
        support_value_count=len(support_values),
        minimum_support=minimum_support,
        maximum_support=maximum_support,
        support_values=support_values,
    )


def _build_iqtree_model_selection_summary(
    prefix_path: Path,
) -> IqtreeModelSelectionSummary | None:
    return parse_iqtree_model_selection_summary(
        iqtree_report_path=prefix_path.with_suffix(".iqtree"),
        model_sidecar_path=resolve_iqtree_model_sidecar(prefix_path),
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
    out_path.parent.mkdir(parents=True, exist_ok=True)
    input_summary = summarise_fasta(input_path)
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
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    report = EngineWorkflowReport(
        workflow="alignment-trimming",
        engine_name="trimal",
        input_paths=[input_path],
        output_paths={"trimmed_alignment": out_path},
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


def run_model_selection(
    input_path: Path,
    *,
    out_dir: Path,
    prefix: str = "model-selection",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
    partition_path: Path | None = None,
    resume: bool = False,
    seed: int = 1,
    threads: int = 1,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run a model-selection workflow on an aligned FASTA file."""
    _ensure_inference_ready_alignment(input_path)
    validate_timeout_seconds(timeout_seconds)
    prefix_path = _prefix_path(out_dir, prefix)
    manifest_path = prefix_path.with_suffix(".manifest.json")
    version = read_engine_version(
        "iqtree",
        executable,
        version_args=("--version",),
        timeout_seconds=timeout_seconds,
    )
    resolved = resolve_engine_executable(executable)
    iqtree_report_path = prefix_path.with_suffix(".iqtree")
    iqtree_log_path = prefix_path.with_suffix(".log")
    prepared_partitions = (
        None
        if partition_path is None
        else _prepare_iqtree_partitions(
            input_path,
            partition_path,
            prefix_path=prefix_path,
        )
    )
    command = [
        resolved,
        *(
            prepared_partitions.command_args
            if prepared_partitions is not None
            else ["-s", str(input_path.resolve())]
        ),
        *(
            []
            if prepared_partitions is not None and prepared_partitions.mixed_data_types
            else _iqtree_sequence_type_flag(input_path, sequence_type)
        ),
        *_iqtree_execution_controls(seed=seed, threads=threads),
        "-m",
        "MF",
        "-pre",
        str(prefix_path.resolve()),
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=(
                [input_path] if partition_path is None else [input_path, partition_path]
            ),
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
        engine_name="iqtree",
        workflow="model-selection",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "iqtree_report": iqtree_report_path,
            "iqtree_log": iqtree_log_path,
        },
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _validate_iqtree_required_artifacts(prefix_path, workflow="model-selection")
        selected_model = _validate_iqtree_model_result(
            prefix_path,
            workflow="model-selection",
        )
        iqtree_summary = _build_iqtree_summary(
            prefix_path,
            default_selected_model=selected_model,
        )
        model_selection_summary = _build_iqtree_model_selection_summary(prefix_path)
        if (
            model_selection_summary is None
            or model_selection_summary.candidate_count < 1
        ):
            raise build_engine_output_error(
                "iqtree model-selection did not expose a parsable candidate-model table",
                code="iqtree_model_candidates_missing",
                engine_name="iqtree",
                workflow="model-selection",
                path=iqtree_report_path,
                output_name="iqtree_report",
                artifact_kind="iqtree-model-candidates",
            )
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    selected_model_path = prefix_path.with_suffix(".selected-model.txt")
    selected_model_path.write_text(selected_model + "\n", encoding="utf-8")
    model_candidates_path = prefix_path.with_suffix(".model-candidates.tsv")
    write_iqtree_model_candidates_table(model_candidates_path, model_selection_summary)
    report = EngineWorkflowReport(
        workflow="model-selection",
        engine_name="iqtree",
        input_paths=(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_paths={
            **({} if prepared_partitions is None else prepared_partitions.output_paths),
            **_existing_iqtree_outputs(prefix_path),
            "selected_model": selected_model_path,
            "model_candidates": model_candidates_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_checksums={},
        config={
            "sequence_type": sequence_type,
            "partition_path": None if partition_path is None else str(partition_path),
            "seed": seed,
            "threads": threads,
            "timeout_seconds": timeout_seconds,
        },
        selected_model=selected_model,
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        model_selection_summary=model_selection_summary,
        notes=[
            *([] if prepared_partitions is None else prepared_partitions.notes),
            f"iqtree random seed: {seed}",
            f"iqtree threads: {threads}",
            "best-fit substitution model parsed from engine output",
            f"parsed {model_selection_summary.candidate_count} candidate substitution models from iqtree output",
            *(
                []
                if model_selection_summary.selected_criterion is None
                else [
                    "model-selection workflow exposed the governing information criterion"
                ]
            ),
            *(
                []
                if iqtree_summary.log_likelihood is None
                else [
                    "model-selection workflow exposed a parsable log-likelihood score"
                ]
            ),
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)


def run_maximum_likelihood_tree_inference(
    input_path: Path,
    *,
    out_dir: Path,
    model: str,
    prefix: str = "maximum-likelihood",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
    partition_path: Path | None = None,
    resume: bool = False,
    seed: int = 1,
    threads: int = 1,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run an external maximum-likelihood tree inference workflow."""
    _ensure_inference_ready_alignment(input_path)
    validate_timeout_seconds(timeout_seconds)
    prefix_path = _prefix_path(out_dir, prefix)
    manifest_path = prefix_path.with_suffix(".manifest.json")
    version = read_engine_version(
        "iqtree",
        executable,
        version_args=("--version",),
        timeout_seconds=timeout_seconds,
    )
    resolved = resolve_engine_executable(executable)
    prepared_partitions = (
        None
        if partition_path is None
        else _prepare_iqtree_partitions(
            input_path,
            partition_path,
            prefix_path=prefix_path,
        )
    )
    if prepared_partitions is not None and not _iqtree_partition_supports_fixed_model(
        model=model,
        mixed_data_types=prepared_partitions.mixed_data_types,
    ):
        raise EngineWorkflowError(
            "mixed DNA/protein partition analyses require a model-selection keyword such as MF, MFP, TEST, or TESTMERGE"
        )
    tree_path = prefix_path.with_suffix(".treefile")
    report_path = prefix_path.with_suffix(".iqtree")
    log_path = prefix_path.with_suffix(".log")
    command = [
        resolved,
        *(
            prepared_partitions.command_args
            if prepared_partitions is not None
            else ["-s", str(input_path.resolve())]
        ),
        *(
            []
            if prepared_partitions is not None and prepared_partitions.mixed_data_types
            else _iqtree_sequence_type_flag(input_path, sequence_type)
        ),
        *_iqtree_execution_controls(seed=seed, threads=threads),
        "-m",
        model,
        "-pre",
        str(prefix_path.resolve()),
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=(
                [input_path] if partition_path is None else [input_path, partition_path]
            ),
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
        engine_name="iqtree",
        workflow="maximum-likelihood-tree",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "tree": tree_path,
            "iqtree_report": report_path,
            "iqtree_log": log_path,
        },
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _validate_iqtree_required_artifacts(
            prefix_path,
            workflow="maximum-likelihood-tree",
        )
        _validate_tree_output(
            tree_path,
            engine_name="iqtree",
            workflow="maximum-likelihood-tree",
            output_name="tree",
            artifact_kind="maximum-likelihood-tree",
        )
        selected_model = _validate_iqtree_model_result(
            prefix_path,
            workflow="maximum-likelihood-tree",
            default_selected_model=model,
        )
        iqtree_summary = _build_iqtree_summary(
            prefix_path,
            default_selected_model=selected_model,
            support_tree_path=tree_path,
        )
        model_selection_summary = _build_iqtree_model_selection_summary(prefix_path)
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    report = EngineWorkflowReport(
        workflow="maximum-likelihood-tree",
        engine_name="iqtree",
        input_paths=(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_paths={
            **({} if prepared_partitions is None else prepared_partitions.output_paths),
            **_existing_iqtree_outputs(prefix_path, include_tree=True),
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_checksums={},
        config={
            "model": model,
            "sequence_type": sequence_type,
            "partition_path": None if partition_path is None else str(partition_path),
            "seed": seed,
            "threads": threads,
            "timeout_seconds": timeout_seconds,
        },
        selected_model=iqtree_summary.selected_model,
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        model_selection_summary=model_selection_summary,
        notes=[
            *([] if prepared_partitions is None else prepared_partitions.notes),
            f"iqtree random seed: {seed}",
            f"iqtree threads: {threads}",
            "maximum-likelihood tree validated as parseable Newick output",
            *(
                []
                if iqtree_summary.log_likelihood is None
                else ["log-likelihood parsed from iqtree inference artifacts"]
            ),
            *(
                []
                if iqtree_summary.support_value_count == 0
                else ["support values parsed from the inferred maximum-likelihood tree"]
            ),
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)


def run_bootstrap_support_estimation(
    input_path: Path,
    *,
    out_dir: Path,
    model: str,
    replicates: int = 1000,
    prefix: str = "bootstrap-support",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
    partition_path: Path | None = None,
    resume: bool = False,
    seed: int = 1,
    threads: int = 1,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run external bootstrap support estimation and retain bootstrap trees."""
    if replicates < 1:
        raise ValueError(f"replicates must be positive, got {replicates}")
    _validate_ufboot_replicates(replicates)
    _ensure_inference_ready_alignment(input_path)
    validate_timeout_seconds(timeout_seconds)
    prefix_path = _prefix_path(out_dir, prefix)
    manifest_path = prefix_path.with_suffix(".manifest.json")
    version = read_engine_version(
        "iqtree",
        executable,
        version_args=("--version",),
        timeout_seconds=timeout_seconds,
    )
    resolved = resolve_engine_executable(executable)
    prepared_partitions = (
        None
        if partition_path is None
        else _prepare_iqtree_partitions(
            input_path,
            partition_path,
            prefix_path=prefix_path,
        )
    )
    if prepared_partitions is not None and not _iqtree_partition_supports_fixed_model(
        model=model,
        mixed_data_types=prepared_partitions.mixed_data_types,
    ):
        raise EngineWorkflowError(
            "mixed DNA/protein partition analyses require a model-selection keyword such as MF, MFP, TEST, or TESTMERGE"
        )
    support_tree_path = prefix_path.with_suffix(".treefile")
    bootstrap_tree_path = prefix_path.with_suffix(".ufboot")
    support_table_path = prefix_path.with_suffix(".support.tsv")
    low_support_branches_path = prefix_path.with_suffix(".low-support.tsv")
    support_histogram_path = prefix_path.with_suffix(".support-histogram.tsv")
    report_path = prefix_path.with_suffix(".iqtree")
    log_path = prefix_path.with_suffix(".log")
    command = [
        resolved,
        *(
            prepared_partitions.command_args
            if prepared_partitions is not None
            else ["-s", str(input_path.resolve())]
        ),
        *(
            []
            if prepared_partitions is not None and prepared_partitions.mixed_data_types
            else _iqtree_sequence_type_flag(input_path, sequence_type)
        ),
        *_iqtree_execution_controls(seed=seed, threads=threads),
        "-m",
        model,
        "-bb",
        str(replicates),
        "-wbt",
        "-pre",
        str(prefix_path.resolve()),
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=(
                [input_path] if partition_path is None else [input_path, partition_path]
            ),
            expected_command=command,
            expected_version=version,
        )
        if resumed is not None and _resume_has_bootstrap_review_outputs(resumed):
            return resumed
    incomplete_notes = _resolve_incomplete_workflow_state(
        manifest_path=manifest_path,
        incomplete_run_policy=incomplete_run_policy,
    )
    run = execute_engine_command(
        engine_name="iqtree",
        workflow="bootstrap-support",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "support_tree": support_tree_path,
            "bootstrap_trees": bootstrap_tree_path,
            "iqtree_report": report_path,
            "iqtree_log": log_path,
        },
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _validate_iqtree_required_artifacts(prefix_path, workflow="bootstrap-support")
        _validate_tree_output(
            support_tree_path,
            engine_name="iqtree",
            workflow="bootstrap-support",
            output_name="support_tree",
            artifact_kind="bootstrap-supported-tree",
        )
        _validate_tree_set_output(
            bootstrap_tree_path,
            engine_name="iqtree",
            workflow="bootstrap-support",
            output_name="bootstrap_trees",
            artifact_kind="bootstrap-tree-set",
        )
        selected_model = _validate_iqtree_model_result(
            prefix_path,
            workflow="bootstrap-support",
            default_selected_model=model,
        )
        iqtree_summary = _build_iqtree_summary(
            prefix_path,
            default_selected_model=selected_model,
            support_tree_path=support_tree_path,
        )
        bootstrap_support_summary = summarize_bootstrap_support_distribution(
            support_tree_path
        )
        _validate_support_value_count(
            engine_name="iqtree",
            workflow="bootstrap-support",
            path=support_tree_path,
            output_name="support_tree",
            artifact_kind="bootstrap-supported-tree",
            support_value_count=bootstrap_support_summary.supported_node_count,
            support_kind="bootstrap support",
        )
        weak_backbone_report = detect_weakly_supported_backbone(support_tree_path)
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    write_bootstrap_support_table(
        support_table_path,
        build_bootstrap_support_rows(bootstrap_support_summary),
    )
    write_bootstrap_support_table(
        low_support_branches_path,
        build_low_support_bootstrap_rows(bootstrap_support_summary),
    )
    write_bootstrap_support_histogram(
        support_histogram_path,
        build_bootstrap_support_histogram_rows(bootstrap_support_summary),
    )
    model_selection_summary = _build_iqtree_model_selection_summary(prefix_path)
    report = EngineWorkflowReport(
        workflow="bootstrap-support",
        engine_name="iqtree",
        input_paths=(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_paths={
            **({} if prepared_partitions is None else prepared_partitions.output_paths),
            **_existing_iqtree_outputs(
                prefix_path,
                include_tree=False,
                include_bootstrap=True,
                include_consensus=True,
            ),
            "support_tree": support_tree_path,
            "support_table": support_table_path,
            "low_support_branches": low_support_branches_path,
            "support_histogram": support_histogram_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_checksums={},
        config={
            "model": model,
            "replicates": replicates,
            "sequence_type": sequence_type,
            "partition_path": None if partition_path is None else str(partition_path),
            "seed": seed,
            "threads": threads,
            "timeout_seconds": timeout_seconds,
        },
        selected_model=iqtree_summary.selected_model,
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        model_selection_summary=model_selection_summary,
        bootstrap_support_summary=bootstrap_support_summary,
        weak_backbone_report=weak_backbone_report,
        notes=[
            *([] if prepared_partitions is None else prepared_partitions.notes),
            f"iqtree random seed: {seed}",
            f"iqtree threads: {threads}",
            "bootstrap tree set retained for downstream consensus construction",
            *(
                []
                if iqtree_summary.log_likelihood is None
                else ["log-likelihood parsed from iqtree bootstrap inference artifacts"]
            ),
            *(
                []
                if iqtree_summary.support_value_count == 0
                else [
                    "support values parsed from the bootstrap-supported tree artifact"
                ]
            ),
            "branch-level support table exported for bootstrap review",
            "low-support branch ledger exported for weak-clade review",
            "support histogram exported for reviewer-facing support distribution checks",
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)


def run_sh_alrt_support_estimation(
    input_path: Path,
    *,
    out_dir: Path,
    model: str,
    sh_alrt_replicates: int = 1000,
    bootstrap_replicates: int = 1000,
    prefix: str = "sh-alrt-support",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
    partition_path: Path | None = None,
    resume: bool = False,
    seed: int = 1,
    threads: int = 1,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run combined SH-aLRT and UFBoot branch-support estimation."""
    _validate_sh_alrt_replicates(sh_alrt_replicates)
    _validate_ufboot_replicates(bootstrap_replicates)
    _ensure_inference_ready_alignment(input_path)
    validate_timeout_seconds(timeout_seconds)
    prefix_path = _prefix_path(out_dir, prefix)
    manifest_path = prefix_path.with_suffix(".manifest.json")
    version = read_engine_version(
        "iqtree",
        executable,
        version_args=("--version",),
        timeout_seconds=timeout_seconds,
    )
    resolved = resolve_engine_executable(executable)
    prepared_partitions = (
        None
        if partition_path is None
        else _prepare_iqtree_partitions(
            input_path,
            partition_path,
            prefix_path=prefix_path,
        )
    )
    if prepared_partitions is not None and not _iqtree_partition_supports_fixed_model(
        model=model,
        mixed_data_types=prepared_partitions.mixed_data_types,
    ):
        raise EngineWorkflowError(
            "mixed DNA/protein partition analyses require a model-selection keyword such as MF, MFP, TEST, or TESTMERGE"
        )
    support_tree_path = prefix_path.with_suffix(".treefile")
    bootstrap_tree_path = prefix_path.with_suffix(".ufboot")
    support_table_path = prefix_path.with_suffix(".support.tsv")
    conflicting_support_branches_path = prefix_path.with_suffix(
        ".conflicting-support.tsv"
    )
    report_path = prefix_path.with_suffix(".iqtree")
    log_path = prefix_path.with_suffix(".log")
    command = [
        resolved,
        *(
            prepared_partitions.command_args
            if prepared_partitions is not None
            else ["-s", str(input_path.resolve())]
        ),
        *(
            []
            if prepared_partitions is not None and prepared_partitions.mixed_data_types
            else _iqtree_sequence_type_flag(input_path, sequence_type)
        ),
        *_iqtree_execution_controls(seed=seed, threads=threads),
        "-m",
        model,
        "-alrt",
        str(sh_alrt_replicates),
        "-bb",
        str(bootstrap_replicates),
        "-wbt",
        "-pre",
        str(prefix_path.resolve()),
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=(
                [input_path] if partition_path is None else [input_path, partition_path]
            ),
            expected_command=command,
            expected_version=version,
        )
        if resumed is not None and _resume_has_sh_alrt_review_outputs(resumed):
            return resumed
    incomplete_notes = _resolve_incomplete_workflow_state(
        manifest_path=manifest_path,
        incomplete_run_policy=incomplete_run_policy,
    )
    run = execute_engine_command(
        engine_name="iqtree",
        workflow="sh-alrt-support",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "support_tree": support_tree_path,
            "bootstrap_trees": bootstrap_tree_path,
            "iqtree_report": report_path,
            "iqtree_log": log_path,
        },
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _validate_iqtree_required_artifacts(prefix_path, workflow="sh-alrt-support")
        _validate_tree_output(
            support_tree_path,
            engine_name="iqtree",
            workflow="sh-alrt-support",
            output_name="support_tree",
            artifact_kind="sh-alrt-supported-tree",
        )
        _validate_tree_set_output(
            bootstrap_tree_path,
            engine_name="iqtree",
            workflow="sh-alrt-support",
            output_name="bootstrap_trees",
            artifact_kind="bootstrap-tree-set",
        )
        selected_model = _validate_iqtree_model_result(
            prefix_path,
            workflow="sh-alrt-support",
            default_selected_model=model,
        )
        iqtree_summary = _build_iqtree_summary(
            prefix_path,
            default_selected_model=selected_model,
            support_tree_path=support_tree_path,
        )
        bootstrap_support_summary = summarize_bootstrap_support_distribution(
            support_tree_path
        )
        _validate_support_value_count(
            engine_name="iqtree",
            workflow="sh-alrt-support",
            path=support_tree_path,
            output_name="support_tree",
            artifact_kind="sh-alrt-supported-tree",
            support_value_count=bootstrap_support_summary.supported_node_count,
            support_kind="ultrafast bootstrap support",
        )
        sh_alrt_support_summary = summarize_sh_alrt_support_distribution(
            support_tree_path
        )
        _validate_support_value_count(
            engine_name="iqtree",
            workflow="sh-alrt-support",
            path=support_tree_path,
            output_name="support_tree",
            artifact_kind="sh-alrt-supported-tree",
            support_value_count=sh_alrt_support_summary.annotated_node_count,
            support_kind="sh-alrt support",
        )
        _validate_support_value_count(
            engine_name="iqtree",
            workflow="sh-alrt-support",
            path=support_tree_path,
            output_name="support_tree",
            artifact_kind="sh-alrt-supported-tree",
            support_value_count=sh_alrt_support_summary.fully_scored_node_count,
            support_kind="joint sh-alrt and ultrafast bootstrap support",
        )
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    write_sh_alrt_support_table(
        support_table_path,
        build_sh_alrt_support_rows(sh_alrt_support_summary),
    )
    write_sh_alrt_support_table(
        conflicting_support_branches_path,
        build_conflicting_sh_alrt_support_rows(sh_alrt_support_summary),
    )
    model_selection_summary = _build_iqtree_model_selection_summary(prefix_path)
    report = EngineWorkflowReport(
        workflow="sh-alrt-support",
        engine_name="iqtree",
        input_paths=(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_paths={
            **({} if prepared_partitions is None else prepared_partitions.output_paths),
            **_existing_iqtree_outputs(
                prefix_path,
                include_tree=False,
                include_bootstrap=True,
                include_consensus=True,
            ),
            "support_tree": support_tree_path,
            "support_table": support_table_path,
            "conflicting_support_branches": conflicting_support_branches_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_checksums={},
        config={
            "model": model,
            "sh_alrt_replicates": sh_alrt_replicates,
            "bootstrap_replicates": bootstrap_replicates,
            "sequence_type": sequence_type,
            "partition_path": None if partition_path is None else str(partition_path),
            "seed": seed,
            "threads": threads,
            "timeout_seconds": timeout_seconds,
        },
        selected_model=iqtree_summary.selected_model,
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        model_selection_summary=model_selection_summary,
        bootstrap_support_summary=bootstrap_support_summary,
        sh_alrt_support_summary=sh_alrt_support_summary,
        notes=[
            *([] if prepared_partitions is None else prepared_partitions.notes),
            f"iqtree random seed: {seed}",
            f"iqtree threads: {threads}",
            f"sh-alrt replicates: {sh_alrt_replicates}",
            f"ultrafast bootstrap replicates: {bootstrap_replicates}",
            "combined sh-alrt and ufboot branch-support table exported for review",
            "conflicting sh-alrt versus ufboot support signals exported for review",
            *(
                []
                if iqtree_summary.log_likelihood is None
                else ["log-likelihood parsed from iqtree support inference artifacts"]
            ),
            *(
                []
                if sh_alrt_support_summary.annotated_node_count == 0
                else [
                    "combined sh-alrt and ufboot support values parsed from the supported tree artifact"
                ]
            ),
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)


def run_bootstrap_consensus_tree(
    bootstrap_trees_path: Path,
    *,
    out_dir: Path,
    prefix: str = "bootstrap-consensus",
    executable: str | Path = "iqtree2",
    minimum_support: float = 0.5,
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Construct a consensus tree from bootstrap trees."""
    if not 0.0 <= minimum_support <= 1.0:
        raise ValueError(
            f"minimum_support must be between 0 and 1 inclusive, got {minimum_support}"
        )
    if not bootstrap_trees_path.exists():
        raise FileNotFoundError(bootstrap_trees_path)
    validate_timeout_seconds(timeout_seconds)
    prefix_path = _prefix_path(out_dir, prefix)
    manifest_path = prefix_path.with_suffix(".manifest.json")
    version = read_engine_version(
        "iqtree",
        executable,
        version_args=("--version",),
        timeout_seconds=timeout_seconds,
    )
    resolved = resolve_engine_executable(executable)
    consensus_tree_path = prefix_path.with_suffix(".contree")
    log_path = prefix_path.with_suffix(".log")
    command = [
        resolved,
        "-t",
        str(bootstrap_trees_path.resolve()),
        "-con",
        "-minsup",
        str(minimum_support),
        "-pre",
        str(prefix_path.resolve()),
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=[bootstrap_trees_path],
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
        engine_name="iqtree",
        workflow="bootstrap-consensus",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "consensus_tree": consensus_tree_path,
            "iqtree_log": log_path,
        },
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _require_nonempty_text_output(
            log_path,
            engine_name="iqtree",
            workflow="bootstrap-consensus",
            output_name="iqtree_log",
            artifact_kind="iqtree-log",
        )
        _validate_tree_output(
            consensus_tree_path,
            engine_name="iqtree",
            workflow="bootstrap-consensus",
            output_name="consensus_tree",
            artifact_kind="bootstrap-consensus-tree",
        )
        iqtree_summary = _build_iqtree_summary(
            prefix_path,
            default_selected_model=None,
            support_tree_path=consensus_tree_path,
        )
        _validate_support_value_count(
            engine_name="iqtree",
            workflow="bootstrap-consensus",
            path=consensus_tree_path,
            output_name="consensus_tree",
            artifact_kind="bootstrap-consensus-tree",
            support_value_count=iqtree_summary.support_value_count,
            support_kind="bootstrap consensus support",
        )
        model_selection_summary = _build_iqtree_model_selection_summary(prefix_path)
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    report = EngineWorkflowReport(
        workflow="bootstrap-consensus",
        engine_name="iqtree",
        input_paths=[bootstrap_trees_path],
        output_paths=_existing_iqtree_outputs(prefix_path, include_consensus=True),
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([bootstrap_trees_path]),
        output_checksums={},
        config={
            "minimum_support": minimum_support,
            "timeout_seconds": timeout_seconds,
        },
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        model_selection_summary=model_selection_summary,
        notes=[
            "consensus tree validated as parseable Newick output",
            *(
                []
                if iqtree_summary.support_value_count == 0
                else ["support values parsed from the bootstrap consensus tree"]
            ),
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)


def run_fast_tree_inference(
    input_path: Path,
    out_path: Path,
    *,
    executable: str | Path = "FastTree",
    sequence_type: AlignmentAlphabet | None = None,
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run a fast approximate tree inference engine against an aligned FASTA file."""
    _ensure_inference_ready_alignment(input_path)
    validate_timeout_seconds(timeout_seconds)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    support_table_path = _sidecar(out_path, "support.tsv")
    low_support_branches_path = _sidecar(out_path, "low-support.tsv")
    support_histogram_path = _sidecar(out_path, "support-histogram.tsv")
    version = read_engine_version(
        "FastTree",
        executable,
        version_args=("-help",),
        timeout_seconds=timeout_seconds,
    )
    resolved = resolve_engine_executable(executable)
    manifest_path = _manifest_path_from_output(out_path)
    command = [resolved, *_fasttree_args(input_path.resolve(), sequence_type)]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=[input_path],
            expected_command=command,
            expected_version=version,
        )
        if resumed is not None and _resume_has_fasttree_review_outputs(resumed):
            return resumed
    incomplete_notes = _resolve_incomplete_workflow_state(
        manifest_path=manifest_path,
        incomplete_run_policy=incomplete_run_policy,
    )
    run = execute_engine_command(
        engine_name="FastTree",
        workflow="fast-approximate-tree",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_path.parent,
        stdout_path=out_path,
        stderr_path=_sidecar(out_path, "stderr.log"),
        output_paths={"tree": out_path},
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _validate_tree_output(
            out_path,
            engine_name="FastTree",
            workflow="fast-approximate-tree",
            output_name="tree",
            artifact_kind="fast-approximate-tree",
        )
        fasttree_support_summary = summarize_fasttree_support_distribution(out_path)
        _validate_support_value_count(
            engine_name="FastTree",
            workflow="fast-approximate-tree",
            path=out_path,
            output_name="tree",
            artifact_kind="fast-approximate-tree",
            support_value_count=fasttree_support_summary.annotated_node_count,
            support_kind="FastTree local support",
        )
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    write_fasttree_support_table(
        support_table_path,
        build_fasttree_support_rows(fasttree_support_summary),
    )
    write_fasttree_support_table(
        low_support_branches_path,
        build_fasttree_low_support_rows(fasttree_support_summary),
    )
    write_fasttree_support_histogram(
        support_histogram_path,
        build_fasttree_support_histogram_rows(fasttree_support_summary),
    )
    report = EngineWorkflowReport(
        workflow="fast-approximate-tree",
        engine_name="FastTree",
        input_paths=[input_path],
        output_paths={
            "tree": out_path,
            "support_table": support_table_path,
            "low_support_branches": low_support_branches_path,
            "support_histogram": support_histogram_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([input_path]),
        output_checksums={},
        config={
            "sequence_type": sequence_type,
            "timeout_seconds": timeout_seconds,
        },
        fasttree_support_summary=fasttree_support_summary,
        notes=[
            "fast approximate tree validated as parseable Newick output",
            "FastTree is an approximately maximum-likelihood engine and should be reviewed as an exploratory or large-alignment inference method",
            "FastTree local support values are SH-like support proportions on the documented 0-1 scale",
            "branch-level FastTree local support table exported for review",
            "low-support FastTree branch ledger exported for weak-clade review",
            "FastTree local support histogram exported for reviewer-facing distribution checks",
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)


def compare_fast_and_ml_trees(
    fast_tree_path: Path,
    ml_tree_path: Path,
    *,
    out_path: Path,
) -> ExternalTreeComparisonReport:
    """Compare a fast approximate tree against a maximum-likelihood tree."""
    comparison_report = build_tree_comparison_report(
        fast_tree_path, ml_tree_path, out_path=out_path
    )
    return ExternalTreeComparisonReport(
        fast_tree_path=fast_tree_path,
        ml_tree_path=ml_tree_path,
        comparison_report=comparison_report,
    )
