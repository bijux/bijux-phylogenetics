from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, datetime
import hashlib
import json
from pathlib import Path

from bijux_phylogenetics.core.alignment import (
    AlignmentAlphabet,
    AlignmentRecord,
    FastaInputValidationReport,
    FastaRepairReport,
)
from bijux_phylogenetics.core.manifest import build_run_manifest, write_run_manifest
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
from bijux_phylogenetics.io.fasta import (
    detect_fasta_sequence_type,
    write_fasta_alignment,
)
from bijux_phylogenetics.io.fasta.records import (
    repair_fasta_input,
    validate_fasta_input,
)
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    MethodTierAssessment,
    fasta_to_tree_method_tier,
)

from .common import build_file_checksums, write_engine_manifest
from .validation import (
    BootstrapSupportSummaryReport,
    ModelSelectionValidationReport,
    summarize_bootstrap_support_distribution,
    validate_model_selection_against_engine_outputs,
)
from .workflows.alignment import (
    run_alignment_trimming,
    run_multiple_sequence_alignment,
)
from .workflows.iqtree import (
    run_bootstrap_support_estimation,
    run_maximum_likelihood_tree_inference,
    run_model_selection,
)
from .workflows.models import (
    EngineWorkflowReport,
)

__all__ = [
    "FastaToTreeModelRow",
    "FastaToTreeStageFingerprint",
    "FastaToTreeSupportRow",
    "FastaToTreeWorkflowReport",
    "build_fasta_to_tree_model_rows",
    "build_fasta_to_tree_support_rows",
    "infer_unaligned_sequence_type",
    "run_fasta_to_tree_workflow",
    "write_fasta_to_tree_log",
    "write_fasta_to_tree_model_table",
    "write_fasta_to_tree_support_table",
]


@dataclass(frozen=True, slots=True)
class FastaToTreeModelRow:
    """One reviewer-facing record describing the selected substitution model."""

    workflow: str
    engine_name: str
    sequence_type: AlignmentAlphabet
    selected_model: str
    report_selected_model: str | None
    artifact_selected_model: str | None
    model_consistent: bool
    alignment_path: Path
    trimmed_alignment_path: Path
    manifest_path: Path


@dataclass(frozen=True, slots=True)
class FastaToTreeSupportRow:
    """One reviewer-facing branch-support record from the final tree."""

    node: str
    descendant_taxa: tuple[str, ...]
    support: float
    support_fraction: float
    is_backbone: bool


@dataclass(frozen=True, slots=True)
class FastaToTreeStageFingerprint:
    """One deterministic fingerprint record for one workflow stage."""

    stage: str
    fingerprint: str
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    config: dict[str, object]
    engine_versions: dict[str, str]
    upstream_fingerprints: dict[str, str]
    resumed: bool


@dataclass(slots=True)
class FastaToTreeWorkflowReport:
    """End-to-end result for one raw-FASTA-to-tree workflow run."""

    workflow: str
    input_path: Path
    prepared_input_path: Path
    out_dir: Path
    prefix: str
    sequence_type: AlignmentAlphabet
    selected_model: str
    alignment_mode: str
    trimming_mode: str
    trim_gap_threshold: float
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    started_at_utc: str
    ended_at_utc: str
    runtime_seconds: float
    engine_artifact_dir: Path
    manifest_path: Path
    run_manifest_path: Path
    output_paths: dict[str, Path]
    config: dict[str, object]
    commands: dict[str, list[str]]
    engine_versions: dict[str, str]
    step_manifests: dict[str, Path]
    stage_fingerprints: dict[str, FastaToTreeStageFingerprint]
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    input_validation: FastaInputValidationReport
    repaired_input_validation: FastaInputValidationReport | None
    input_repair: FastaRepairReport | None
    alignment_workflow: EngineWorkflowReport
    trimming_workflow: EngineWorkflowReport
    model_selection_workflow: EngineWorkflowReport
    maximum_likelihood_workflow: EngineWorkflowReport
    bootstrap_workflow: EngineWorkflowReport
    model_validation: ModelSelectionValidationReport
    support_summary: BootstrapSupportSummaryReport
    model_rows: list[FastaToTreeModelRow]
    support_rows: list[FastaToTreeSupportRow]
    method_tier: MethodTierAssessment
    warnings: list[str]
    notes: list[str]


def infer_unaligned_sequence_type(records: list[tuple[str, str]]) -> AlignmentAlphabet:
    """Infer a stable sequence type from raw FASTA records before alignment."""
    report = detect_fasta_sequence_type(
        Path("<memory>"),
        records=[
            AlignmentRecord(identifier=identifier, sequence=sequence)
            for identifier, sequence in records
        ],
    )
    return "unknown" if report.selected_type is None else report.selected_type


def build_fasta_to_tree_model_rows(
    workflow_report: EngineWorkflowReport,
    *,
    validation: ModelSelectionValidationReport,
    sequence_type: AlignmentAlphabet,
    alignment_path: Path,
    trimmed_alignment_path: Path,
) -> list[FastaToTreeModelRow]:
    """Convert one model-selection workflow report into a TSV-ready row set."""
    if workflow_report.selected_model is None:
        raise ValueError("model-selection workflow report must expose selected_model")
    return [
        FastaToTreeModelRow(
            workflow=workflow_report.workflow,
            engine_name=workflow_report.engine_name,
            sequence_type=sequence_type,
            selected_model=workflow_report.selected_model,
            report_selected_model=validation.report_selected_model,
            artifact_selected_model=validation.artifact_selected_model,
            model_consistent=validation.valid,
            alignment_path=alignment_path,
            trimmed_alignment_path=trimmed_alignment_path,
            manifest_path=workflow_report.manifest_path,
        )
    ]


def build_fasta_to_tree_support_rows(
    support_summary: BootstrapSupportSummaryReport,
) -> list[FastaToTreeSupportRow]:
    """Convert branch-support summaries into TSV-ready support rows."""
    return [
        FastaToTreeSupportRow(
            node=node.node,
            descendant_taxa=tuple(node.descendant_taxa),
            support=node.support,
            support_fraction=node.support_fraction,
            is_backbone=node.is_backbone,
        )
        for node in support_summary.nodes
    ]


def _serialize_support_taxa(taxa: tuple[str, ...]) -> str:
    return ",".join(taxa)


def _write_tsv(path: Path, *, header: list[str], rows: list[list[str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _copy_output(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(source.read_bytes())
    return destination


def _normalize_stage_fingerprint_payload(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {
            str(key): _normalize_stage_fingerprint_payload(item)
            for key, item in sorted(value.items())
        }
    if isinstance(value, (list, tuple)):
        return [_normalize_stage_fingerprint_payload(item) for item in value]
    return value


def _build_stage_fingerprint(
    *,
    stage: str,
    input_checksums: dict[str, str],
    output_checksums: dict[str, str],
    config: dict[str, object],
    engine_versions: dict[str, str],
    upstream_fingerprints: dict[str, str],
    resumed: bool,
) -> FastaToTreeStageFingerprint:
    payload = _normalize_stage_fingerprint_payload(
        {
            "stage": stage,
            "input_checksums": input_checksums,
            "output_checksums": output_checksums,
            "config": config,
            "engine_versions": engine_versions,
            "upstream_fingerprints": upstream_fingerprints,
        }
    )
    digest = hashlib.sha256(  # nosec B324
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return FastaToTreeStageFingerprint(
        stage=stage,
        fingerprint=digest,
        input_checksums=input_checksums,
        output_checksums=output_checksums,
        config=config,
        engine_versions=engine_versions,
        upstream_fingerprints=upstream_fingerprints,
        resumed=resumed,
    )


def _display_path(path: Path, *, root_dir: Path | None) -> str:
    if root_dir is None:
        return str(path)
    try:
        return str(path.relative_to(root_dir))
    except ValueError:
        return str(path)


def _display_command(
    command: list[str],
    *,
    root_dir: Path | None,
    aliases: dict[str, str] | None = None,
) -> str:
    rendered: list[str] = []
    token_aliases = {} if aliases is None else aliases
    for token in command:
        if token in token_aliases:
            rendered.append(token_aliases[token])
            continue
        if token.startswith("/"):
            rendered.append(_display_path(Path(token), root_dir=root_dir))
        else:
            rendered.append(token)
    return " ".join(rendered)


def _artifact_prefix(out_dir: Path, prefix: str, step_name: str) -> Path:
    return out_dir / "engine-artifacts" / prefix / step_name / step_name


def _final_output_paths(out_dir: Path, prefix: str) -> dict[str, Path]:
    root = out_dir / prefix
    return {
        "alignment": root.with_suffix(".aln"),
        "trimmed_alignment": root.with_suffix(".trimmed.aln"),
        "tree": root.with_suffix(".tree"),
        "log": root.with_suffix(".log"),
        "methods_summary": root.with_suffix(".methods-summary.md"),
        "model_table": root.with_suffix(".model.tsv"),
        "support_table": root.with_suffix(".support.tsv"),
        "manifest": root.with_suffix(".manifest.json"),
        "run_manifest": root.with_suffix(".run.json"),
    }


def write_fasta_to_tree_model_table(
    path: Path,
    rows: list[FastaToTreeModelRow],
    *,
    root_dir: Path | None = None,
) -> Path:
    """Write the selected-model table for one FASTA-to-tree workflow."""
    return _write_tsv(
        path,
        header=[
            "workflow",
            "engine_name",
            "sequence_type",
            "selected_model",
            "report_selected_model",
            "artifact_selected_model",
            "model_consistent",
            "alignment_path",
            "trimmed_alignment_path",
            "manifest_path",
        ],
        rows=[
            [
                row.workflow,
                row.engine_name,
                row.sequence_type,
                row.selected_model,
                "" if row.report_selected_model is None else row.report_selected_model,
                ""
                if row.artifact_selected_model is None
                else row.artifact_selected_model,
                "true" if row.model_consistent else "false",
                _display_path(row.alignment_path, root_dir=root_dir),
                _display_path(row.trimmed_alignment_path, root_dir=root_dir),
                _display_path(row.manifest_path, root_dir=root_dir),
            ]
            for row in rows
        ],
    )


def write_fasta_to_tree_support_table(
    path: Path,
    rows: list[FastaToTreeSupportRow],
) -> Path:
    """Write the branch-support table for one FASTA-to-tree workflow."""
    return _write_tsv(
        path,
        header=[
            "node",
            "descendant_taxa",
            "support",
            "support_fraction",
            "is_backbone",
        ],
        rows=[
            [
                row.node,
                _serialize_support_taxa(row.descendant_taxa),
                format(row.support, ".12g"),
                format(row.support_fraction, ".12g"),
                "true" if row.is_backbone else "false",
            ]
            for row in rows
        ],
    )


def write_fasta_to_tree_log(
    path: Path,
    report: FastaToTreeWorkflowReport,
    *,
    root_dir: Path | None = None,
) -> Path:
    """Write a reviewer-facing plain-text log for one workflow run."""
    lines = [
        "workflow: fasta-to-tree",
        f"input_path: {_display_path(report.input_path, root_dir=root_dir)}",
        "prepared_input_path: "
        f"{_display_path(report.prepared_input_path, root_dir=root_dir)}",
        f"sequence_type: {report.sequence_type}",
        f"selected_model: {report.selected_model}",
        f"started_at_utc: {report.started_at_utc}",
        f"ended_at_utc: {report.ended_at_utc}",
        f"runtime_seconds: {report.runtime_seconds:.3f}",
        "run_manifest_path: "
        f"{_display_path(report.run_manifest_path, root_dir=root_dir)}",
        "",
        "[input-validation]",
        f"sequence_count: {report.input_validation.summary.sequence_count}",
    ]
    if report.input_validation.warnings:
        lines.append("warnings:")
        lines.extend(f"- {warning}" for warning in report.input_validation.warnings)
    lines.append("")
    if report.input_repair is not None:
        lines.extend(
            [
                "[input-repair]",
                f"output_path: {report.input_repair.output_path}",
                f"normalized_identifier_count: {len(report.input_repair.normalized_identifiers)}",
                f"removed_record_count: {len(report.input_repair.removed_records)}",
                "",
            ]
        )
    steps = [
        ("alignment", report.alignment_workflow),
        ("trimming", report.trimming_workflow),
        ("model_selection", report.model_selection_workflow),
        ("maximum_likelihood", report.maximum_likelihood_workflow),
        ("bootstrap_support", report.bootstrap_workflow),
    ]
    command_aliases = {str(report.input_path.resolve()): str(report.input_path)}
    for label, workflow in steps:
        lines.extend(
            [
                f"[{label}]",
                f"engine: {workflow.engine_name}",
                "command: "
                + _display_command(
                    workflow.run.command,
                    root_dir=root_dir,
                    aliases=command_aliases,
                ),
                f"version: {workflow.run.version.text}",
                f"manifest: {_display_path(workflow.manifest_path, root_dir=root_dir)}",
            ]
        )
        for output_label, output_path in workflow.output_paths.items():
            lines.append(
                f"output.{output_label}: {_display_path(output_path, root_dir=root_dir)}"
            )
        if workflow.run.warning_lines:
            lines.append("warnings:")
            lines.extend(f"- {warning}" for warning in workflow.run.warning_lines)
        lines.append("")
    if report.support_summary.warnings:
        lines.append("[support-summary]")
        lines.extend(f"- {warning}" for warning in report.support_summary.warnings)
        lines.append("")
    if report.notes:
        lines.append("[notes]")
        lines.extend(f"- {note}" for note in report.notes)
        lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    return path


def run_fasta_to_tree_workflow(
    input_path: Path,
    *,
    out_dir: Path,
    prefix: str | None = None,
    sequence_type: AlignmentAlphabet | None = None,
    mafft_executable: str | Path = "mafft",
    alignment_mode: str = "auto",
    trimal_executable: str | Path = "trimal",
    trimming_mode: str = "gap-threshold",
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = 1,
    iqtree_threads: int = 1,
    trim_gap_threshold: float = 0.1,
    bootstrap_replicates: int = 1000,
    normalize_identifiers: bool = False,
    remove_invalid_records: bool = False,
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> FastaToTreeWorkflowReport:
    """Run alignment, trimming, model selection, ML inference, and bootstrap support in one workflow."""
    started_at = datetime.now(UTC)
    workflow_prefix = input_path.stem if prefix is None else prefix
    engine_artifact_dir = out_dir / "engine-artifacts" / workflow_prefix
    input_validation = validate_fasta_input(input_path, sequence_type=sequence_type)
    has_input_blockers = bool(
        input_validation.duplicate_identifiers
        or input_validation.illegal_characters
        or input_validation.empty_sequences
    )
    prepared_input_path = input_path
    repaired_input_validation: FastaInputValidationReport | None = None
    input_repair: FastaRepairReport | None = None
    if has_input_blockers and not (normalize_identifiers or remove_invalid_records):
        raise InvalidAlignmentError(
            "fasta-to-tree input contains duplicate identifiers, empty sequences, "
            "or illegal characters; use normalize_identifiers or "
            "remove_invalid_records to repair the input explicitly"
        )
    if normalize_identifiers or remove_invalid_records:
        repaired_records, input_repair = repair_fasta_input(
            input_path,
            sequence_type=sequence_type,
            normalize_identifiers=normalize_identifiers,
            remove_invalid_records=remove_invalid_records,
        )
        prepared_input_path = _artifact_prefix(
            out_dir, workflow_prefix, "input-curation"
        ).with_suffix(".fasta")
        write_fasta_alignment(prepared_input_path, repaired_records)
        input_repair.output_path = prepared_input_path
        repaired_input_validation = validate_fasta_input(
            prepared_input_path,
            sequence_type=sequence_type,
        )
        if (
            repaired_input_validation.duplicate_identifiers
            or repaired_input_validation.illegal_characters
            or repaired_input_validation.empty_sequences
        ):
            raise InvalidAlignmentError(
                "fasta-to-tree repair policy left unresolved duplicate identifiers, "
                "illegal characters, or empty sequences in the prepared input"
            )
    effective_input_validation = (
        input_validation
        if repaired_input_validation is None
        else repaired_input_validation
    )
    inferred_sequence_type = (
        effective_input_validation.sequence_type_report.selected_type
        if sequence_type is None
        else sequence_type
    )
    if inferred_sequence_type in {None, "unknown"}:
        raise InvalidAlignmentError(
            "fasta-to-tree workflow could not resolve one supported raw sequence type: "
            f"{effective_input_validation.sequence_type_report.note}"
        )
    final_outputs = _final_output_paths(out_dir, workflow_prefix)

    alignment_workflow = run_multiple_sequence_alignment(
        prepared_input_path,
        _artifact_prefix(out_dir, workflow_prefix, "alignment").with_suffix(".aln"),
        executable=mafft_executable,
        mode=alignment_mode,
        resume=resume,
        timeout_seconds=timeout_seconds,
        incomplete_run_policy=incomplete_run_policy,
    )
    trimming_workflow = run_alignment_trimming(
        alignment_workflow.output_paths["alignment"],
        _artifact_prefix(out_dir, workflow_prefix, "trimming").with_suffix(
            ".trimmed.aln"
        ),
        executable=trimal_executable,
        mode=trimming_mode,
        gap_threshold=trim_gap_threshold,
        resume=resume,
        timeout_seconds=timeout_seconds,
        incomplete_run_policy=incomplete_run_policy,
    )
    model_selection_workflow = run_model_selection(
        trimming_workflow.output_paths["trimmed_alignment"],
        out_dir=_artifact_prefix(out_dir, workflow_prefix, "model-selection").parent,
        prefix="model-selection",
        executable=iqtree_executable,
        sequence_type=inferred_sequence_type,
        resume=resume,
        seed=iqtree_seed,
        threads=iqtree_threads,
        timeout_seconds=timeout_seconds,
        incomplete_run_policy=incomplete_run_policy,
    )
    if model_selection_workflow.selected_model is None:
        raise ValueError("model-selection workflow did not expose a selected model")
    maximum_likelihood_workflow = run_maximum_likelihood_tree_inference(
        trimming_workflow.output_paths["trimmed_alignment"],
        out_dir=_artifact_prefix(out_dir, workflow_prefix, "maximum-likelihood").parent,
        model=model_selection_workflow.selected_model,
        prefix="maximum-likelihood",
        executable=iqtree_executable,
        sequence_type=inferred_sequence_type,
        resume=resume,
        seed=iqtree_seed,
        threads=iqtree_threads,
        timeout_seconds=timeout_seconds,
        incomplete_run_policy=incomplete_run_policy,
    )
    bootstrap_workflow = run_bootstrap_support_estimation(
        trimming_workflow.output_paths["trimmed_alignment"],
        out_dir=_artifact_prefix(out_dir, workflow_prefix, "bootstrap-support").parent,
        model=model_selection_workflow.selected_model,
        replicates=bootstrap_replicates,
        prefix="bootstrap-support",
        executable=iqtree_executable,
        sequence_type=inferred_sequence_type,
        resume=resume,
        seed=iqtree_seed,
        threads=iqtree_threads,
        timeout_seconds=timeout_seconds,
        incomplete_run_policy=incomplete_run_policy,
    )

    model_validation = validate_model_selection_against_engine_outputs(
        model_selection_workflow.manifest_path
    )
    _copy_output(
        alignment_workflow.output_paths["alignment"], final_outputs["alignment"]
    )
    _copy_output(
        trimming_workflow.output_paths["trimmed_alignment"],
        final_outputs["trimmed_alignment"],
    )
    _copy_output(bootstrap_workflow.output_paths["support_tree"], final_outputs["tree"])
    model_rows = build_fasta_to_tree_model_rows(
        model_selection_workflow,
        validation=model_validation,
        sequence_type=inferred_sequence_type,
        alignment_path=final_outputs["alignment"],
        trimmed_alignment_path=final_outputs["trimmed_alignment"],
    )
    write_fasta_to_tree_model_table(
        final_outputs["model_table"],
        model_rows,
        root_dir=out_dir,
    )
    support_summary = summarize_bootstrap_support_distribution(final_outputs["tree"])
    support_rows = build_fasta_to_tree_support_rows(support_summary)
    write_fasta_to_tree_support_table(final_outputs["support_table"], support_rows)

    warnings = list(
        dict.fromkeys(
            alignment_workflow.run.warning_lines
            + trimming_workflow.run.warning_lines
            + model_selection_workflow.run.warning_lines
            + maximum_likelihood_workflow.run.warning_lines
            + bootstrap_workflow.run.warning_lines
            + support_summary.warnings
            + input_validation.warnings
            + ([] if input_repair is None else input_repair.warnings)
            + (
                []
                if repaired_input_validation is None
                else repaired_input_validation.warnings
            )
        )
    )
    notes = [
        "final tree path contains the bootstrap-supported inference tree",
        "engine-specific intermediate artifacts remain under engine-artifacts/",
        f"mafft alignment mode: {alignment_mode}",
        f"trimal trimming mode: {trimming_mode}",
        f"iqtree random seed: {iqtree_seed}",
        f"iqtree threads: {iqtree_threads}",
        "raw sequence type detection: "
        f"{effective_input_validation.sequence_type_report.detected_type} "
        f"({effective_input_validation.sequence_type_report.confidence})",
        effective_input_validation.sequence_type_report.note,
    ]
    if input_repair is not None:
        notes.append(
            f"prepared input FASTA written to {prepared_input_path} before alignment"
        )
    run_manifest_arguments = [
        str(input_path),
        "--out-dir",
        str(out_dir),
        "--prefix",
        workflow_prefix,
        "--sequence-type",
        inferred_sequence_type,
        "--alignment-mode",
        alignment_mode,
        "--trimming-mode",
        trimming_mode,
        "--trim-gap-threshold",
        format(trim_gap_threshold, ".12g"),
        "--iqtree-seed",
        str(iqtree_seed),
        "--iqtree-threads",
        str(iqtree_threads),
        "--bootstrap-replicates",
        str(bootstrap_replicates),
        "--resume",
        "true" if resume else "false",
        "--incomplete-run-policy",
        incomplete_run_policy,
    ]
    if timeout_seconds is not None:
        run_manifest_arguments.extend(
            ["--timeout-seconds", format(timeout_seconds, ".12g")]
        )
    report = FastaToTreeWorkflowReport(
        workflow="fasta-to-tree",
        input_path=input_path,
        prepared_input_path=prepared_input_path,
        out_dir=out_dir,
        prefix=workflow_prefix,
        sequence_type=inferred_sequence_type,
        selected_model=model_selection_workflow.selected_model,
        alignment_mode=alignment_mode,
        trimming_mode=trimming_mode,
        trim_gap_threshold=trim_gap_threshold,
        iqtree_seed=iqtree_seed,
        iqtree_threads=iqtree_threads,
        bootstrap_replicates=bootstrap_replicates,
        started_at_utc=started_at.replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z"),
        ended_at_utc="",
        runtime_seconds=0.0,
        engine_artifact_dir=engine_artifact_dir,
        manifest_path=final_outputs["manifest"],
        run_manifest_path=final_outputs["run_manifest"],
        output_paths=final_outputs,
        config={
            "sequence_type": inferred_sequence_type,
            "alignment_mode": alignment_mode,
            "trimming_mode": trimming_mode,
            "trim_gap_threshold": trim_gap_threshold,
            "iqtree_seed": iqtree_seed,
            "iqtree_threads": iqtree_threads,
            "bootstrap_replicates": bootstrap_replicates,
            "timeout_seconds": timeout_seconds,
            "resume": resume,
            "incomplete_run_policy": incomplete_run_policy,
        },
        commands={
            "alignment": alignment_workflow.run.command,
            "trimming": trimming_workflow.run.command,
            "model_selection": model_selection_workflow.run.command,
            "maximum_likelihood": maximum_likelihood_workflow.run.command,
            "bootstrap_support": bootstrap_workflow.run.command,
        },
        engine_versions={
            "mafft": alignment_workflow.run.version.text,
            "trimal": trimming_workflow.run.version.text,
            "iqtree_model_selection": model_selection_workflow.run.version.text,
            "iqtree_maximum_likelihood": maximum_likelihood_workflow.run.version.text,
            "iqtree_bootstrap_support": bootstrap_workflow.run.version.text,
        },
        step_manifests={
            "alignment": alignment_workflow.manifest_path,
            "trimming": trimming_workflow.manifest_path,
            "model_selection": model_selection_workflow.manifest_path,
            "maximum_likelihood": maximum_likelihood_workflow.manifest_path,
            "bootstrap_support": bootstrap_workflow.manifest_path,
        },
        stage_fingerprints={},
        input_checksums=build_file_checksums(
            [input_path]
            if prepared_input_path == input_path
            else [input_path, prepared_input_path]
        ),
        output_checksums={},
        input_validation=input_validation,
        repaired_input_validation=repaired_input_validation,
        input_repair=input_repair,
        alignment_workflow=alignment_workflow,
        trimming_workflow=trimming_workflow,
        model_selection_workflow=model_selection_workflow,
        maximum_likelihood_workflow=maximum_likelihood_workflow,
        bootstrap_workflow=bootstrap_workflow,
        model_validation=model_validation,
        support_summary=support_summary,
        model_rows=model_rows,
        support_rows=support_rows,
        method_tier=fasta_to_tree_method_tier(),
        warnings=warnings,
        notes=notes,
    )
    ended_at = datetime.now(UTC)
    report.ended_at_utc = (
        ended_at.replace(microsecond=0).isoformat().replace("+00:00", "Z")
    )
    report.runtime_seconds = max(
        0.0,
        round((ended_at - started_at).total_seconds(), 6),
    )
    write_fasta_to_tree_log(final_outputs["log"], report, root_dir=out_dir)
    from bijux_phylogenetics.reports.methods import (
        write_tree_inference_methods_summary_text,
    )

    write_tree_inference_methods_summary_text(
        final_outputs["methods_summary"],
        workflow_report=report,
    )
    validation_output_checksums = (
        {}
        if prepared_input_path == input_path
        else build_file_checksums([prepared_input_path])
    )
    validation_stage = _build_stage_fingerprint(
        stage="fasta_validation",
        input_checksums=build_file_checksums([input_path]),
        output_checksums=validation_output_checksums,
        config={
            "declared_sequence_type": sequence_type,
            "normalize_identifiers": normalize_identifiers,
            "remove_invalid_records": remove_invalid_records,
        },
        engine_versions={},
        upstream_fingerprints={},
        resumed=False,
    )
    alignment_stage = _build_stage_fingerprint(
        stage="alignment",
        input_checksums=alignment_workflow.input_checksums,
        output_checksums=alignment_workflow.output_checksums,
        config={
            "alignment_mode": alignment_mode,
            "timeout_seconds": timeout_seconds,
        },
        engine_versions={"mafft": alignment_workflow.run.version.text},
        upstream_fingerprints={"fasta_validation": validation_stage.fingerprint},
        resumed=alignment_workflow.resumed,
    )
    trimming_stage = _build_stage_fingerprint(
        stage="trimming",
        input_checksums=trimming_workflow.input_checksums,
        output_checksums=trimming_workflow.output_checksums,
        config={
            "trimming_mode": trimming_mode,
            "trim_gap_threshold": trim_gap_threshold,
            "timeout_seconds": timeout_seconds,
        },
        engine_versions={"trimal": trimming_workflow.run.version.text},
        upstream_fingerprints={"alignment": alignment_stage.fingerprint},
        resumed=trimming_workflow.resumed,
    )
    model_selection_stage = _build_stage_fingerprint(
        stage="model_selection",
        input_checksums=model_selection_workflow.input_checksums,
        output_checksums=model_selection_workflow.output_checksums,
        config={
            "iqtree_seed": iqtree_seed,
            "iqtree_threads": iqtree_threads,
            "timeout_seconds": timeout_seconds,
            "sequence_type": inferred_sequence_type,
        },
        engine_versions={
            "iqtree_model_selection": model_selection_workflow.run.version.text
        },
        upstream_fingerprints={"trimming": trimming_stage.fingerprint},
        resumed=model_selection_workflow.resumed,
    )
    inference_stage = _build_stage_fingerprint(
        stage="inference",
        input_checksums=maximum_likelihood_workflow.input_checksums,
        output_checksums=maximum_likelihood_workflow.output_checksums,
        config={
            "selected_model": model_selection_workflow.selected_model,
            "iqtree_seed": iqtree_seed,
            "iqtree_threads": iqtree_threads,
            "timeout_seconds": timeout_seconds,
            "sequence_type": inferred_sequence_type,
        },
        engine_versions={
            "iqtree_maximum_likelihood": (maximum_likelihood_workflow.run.version.text)
        },
        upstream_fingerprints={
            "trimming": trimming_stage.fingerprint,
            "model_selection": model_selection_stage.fingerprint,
        },
        resumed=maximum_likelihood_workflow.resumed,
    )
    support_stage = _build_stage_fingerprint(
        stage="support",
        input_checksums=bootstrap_workflow.input_checksums,
        output_checksums=bootstrap_workflow.output_checksums,
        config={
            "selected_model": model_selection_workflow.selected_model,
            "bootstrap_replicates": bootstrap_replicates,
            "iqtree_seed": iqtree_seed,
            "iqtree_threads": iqtree_threads,
            "timeout_seconds": timeout_seconds,
            "sequence_type": inferred_sequence_type,
        },
        engine_versions={
            "iqtree_bootstrap_support": bootstrap_workflow.run.version.text
        },
        upstream_fingerprints={
            "trimming": trimming_stage.fingerprint,
            "model_selection": model_selection_stage.fingerprint,
        },
        resumed=bootstrap_workflow.resumed,
    )
    report_output_checksums = build_file_checksums(
        [
            final_outputs["log"],
            final_outputs["methods_summary"],
            final_outputs["model_table"],
            final_outputs["support_table"],
        ]
    )
    report_stage = _build_stage_fingerprint(
        stage="report",
        input_checksums=build_file_checksums(
            [
                final_outputs["alignment"],
                final_outputs["trimmed_alignment"],
                final_outputs["tree"],
            ]
        ),
        output_checksums=report_output_checksums,
        config={
            "selected_model": model_selection_workflow.selected_model,
            "warning_count": len(warnings),
            "note_count": len(notes),
        },
        engine_versions={},
        upstream_fingerprints={
            "fasta_validation": validation_stage.fingerprint,
            "alignment": alignment_stage.fingerprint,
            "trimming": trimming_stage.fingerprint,
            "model_selection": model_selection_stage.fingerprint,
            "inference": inference_stage.fingerprint,
            "support": support_stage.fingerprint,
        },
        resumed=False,
    )
    report.stage_fingerprints = {
        "fasta_validation": validation_stage,
        "alignment": alignment_stage,
        "trimming": trimming_stage,
        "model_selection": model_selection_stage,
        "inference": inference_stage,
        "support": support_stage,
        "report": report_stage,
    }
    report.output_checksums = build_file_checksums(
        [
            final_outputs["alignment"],
            final_outputs["trimmed_alignment"],
            final_outputs["tree"],
            final_outputs["log"],
            final_outputs["methods_summary"],
            final_outputs["model_table"],
            final_outputs["support_table"],
        ]
    )
    write_engine_manifest(report.manifest_path, report)
    write_run_manifest(
        report.run_manifest_path,
        build_run_manifest(
            command="run_fasta_to_tree_workflow",
            arguments=run_manifest_arguments,
            input_paths=(
                [input_path]
                if prepared_input_path == input_path
                else [input_path, prepared_input_path]
            ),
            output_paths=[
                final_outputs["alignment"],
                final_outputs["trimmed_alignment"],
                final_outputs["tree"],
                final_outputs["log"],
                final_outputs["methods_summary"],
                final_outputs["model_table"],
                final_outputs["support_table"],
                report.manifest_path,
            ],
        ),
    )
    return report
