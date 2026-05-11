from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.alignment import AlignmentAlphabet, AlignmentRecord
from bijux_phylogenetics.io.fasta import infer_alignment_alphabet

from .common import build_file_checksums, load_unaligned_fasta, write_engine_manifest
from .validation import (
    BootstrapSupportSummaryReport,
    ModelSelectionValidationReport,
    summarize_bootstrap_support_distribution,
    validate_model_selection_against_engine_outputs,
)
from .workflows import (
    EngineWorkflowReport,
    run_alignment_trimming,
    run_bootstrap_support_estimation,
    run_maximum_likelihood_tree_inference,
    run_model_selection,
    run_multiple_sequence_alignment,
)

__all__ = [
    "FastaToTreeModelRow",
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


@dataclass(slots=True)
class FastaToTreeWorkflowReport:
    """End-to-end result for one raw-FASTA-to-tree workflow run."""

    input_path: Path
    out_dir: Path
    prefix: str
    sequence_type: AlignmentAlphabet
    selected_model: str
    engine_artifact_dir: Path
    manifest_path: Path
    output_paths: dict[str, Path]
    step_manifests: dict[str, Path]
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    alignment_workflow: EngineWorkflowReport
    trimming_workflow: EngineWorkflowReport
    model_selection_workflow: EngineWorkflowReport
    maximum_likelihood_workflow: EngineWorkflowReport
    bootstrap_workflow: EngineWorkflowReport
    model_validation: ModelSelectionValidationReport
    support_summary: BootstrapSupportSummaryReport
    model_rows: list[FastaToTreeModelRow]
    support_rows: list[FastaToTreeSupportRow]
    warnings: list[str]
    notes: list[str]


def infer_unaligned_sequence_type(records: list[tuple[str, str]]) -> AlignmentAlphabet:
    """Infer a stable sequence type from raw FASTA records before alignment."""
    return infer_alignment_alphabet(
        [
            AlignmentRecord(identifier=identifier, sequence=sequence)
            for identifier, sequence in records
        ]
    )


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


def _artifact_prefix(out_dir: Path, prefix: str, step_name: str) -> Path:
    return out_dir / "engine-artifacts" / prefix / step_name / step_name


def _final_output_paths(out_dir: Path, prefix: str) -> dict[str, Path]:
    root = out_dir / prefix
    return {
        "alignment": root.with_suffix(".aln"),
        "trimmed_alignment": root.with_suffix(".trimmed.aln"),
        "tree": root.with_suffix(".tree"),
        "log": root.with_suffix(".log"),
        "model_table": root.with_suffix(".model.tsv"),
        "support_table": root.with_suffix(".support.tsv"),
        "manifest": root.with_suffix(".manifest.json"),
    }


def write_fasta_to_tree_model_table(
    path: Path,
    rows: list[FastaToTreeModelRow],
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
                str(row.alignment_path),
                str(row.trimmed_alignment_path),
                str(row.manifest_path),
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


def write_fasta_to_tree_log(path: Path, report: FastaToTreeWorkflowReport) -> Path:
    """Write a reviewer-facing plain-text log for one workflow run."""
    lines = [
        "workflow: fasta-to-tree",
        f"input_path: {report.input_path}",
        f"sequence_type: {report.sequence_type}",
        f"selected_model: {report.selected_model}",
        "",
    ]
    steps = [
        ("alignment", report.alignment_workflow),
        ("trimming", report.trimming_workflow),
        ("model_selection", report.model_selection_workflow),
        ("maximum_likelihood", report.maximum_likelihood_workflow),
        ("bootstrap_support", report.bootstrap_workflow),
    ]
    for label, workflow in steps:
        lines.extend(
            [
                f"[{label}]",
                f"engine: {workflow.engine_name}",
                f"command: {' '.join(workflow.run.command)}",
                f"version: {workflow.run.version.text}",
                f"manifest: {workflow.manifest_path}",
            ]
        )
        for output_label, output_path in workflow.output_paths.items():
            lines.append(f"output.{output_label}: {output_path}")
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
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    trim_gap_threshold: float = 0.1,
    bootstrap_replicates: int = 1000,
) -> FastaToTreeWorkflowReport:
    """Run alignment, trimming, model selection, ML inference, and bootstrap support in one workflow."""
    raw_records = load_unaligned_fasta(input_path)
    inferred_sequence_type = (
        infer_unaligned_sequence_type(raw_records)
        if sequence_type is None
        else sequence_type
    )
    if inferred_sequence_type == "unknown":
        raise ValueError(
            "fasta-to-tree workflow could not infer a supported sequence type from the input FASTA"
        )
    workflow_prefix = input_path.stem if prefix is None else prefix
    final_outputs = _final_output_paths(out_dir, workflow_prefix)

    alignment_workflow = run_multiple_sequence_alignment(
        input_path,
        _artifact_prefix(out_dir, workflow_prefix, "alignment").with_suffix(".aln"),
        executable=mafft_executable,
    )
    trimming_workflow = run_alignment_trimming(
        alignment_workflow.output_paths["alignment"],
        _artifact_prefix(out_dir, workflow_prefix, "trimming").with_suffix(
            ".trimmed.aln"
        ),
        executable=trimal_executable,
        gap_threshold=trim_gap_threshold,
    )
    model_selection_workflow = run_model_selection(
        trimming_workflow.output_paths["trimmed_alignment"],
        out_dir=_artifact_prefix(out_dir, workflow_prefix, "model-selection").parent,
        prefix="model-selection",
        executable=iqtree_executable,
        sequence_type=inferred_sequence_type,
    )
    if model_selection_workflow.selected_model is None:
        raise ValueError("model-selection workflow did not expose a selected model")
    maximum_likelihood_workflow = run_maximum_likelihood_tree_inference(
        trimming_workflow.output_paths["trimmed_alignment"],
        out_dir=_artifact_prefix(
            out_dir, workflow_prefix, "maximum-likelihood"
        ).parent,
        model=model_selection_workflow.selected_model,
        prefix="maximum-likelihood",
        executable=iqtree_executable,
        sequence_type=inferred_sequence_type,
    )
    bootstrap_workflow = run_bootstrap_support_estimation(
        trimming_workflow.output_paths["trimmed_alignment"],
        out_dir=_artifact_prefix(out_dir, workflow_prefix, "bootstrap-support").parent,
        model=model_selection_workflow.selected_model,
        replicates=bootstrap_replicates,
        prefix="bootstrap-support",
        executable=iqtree_executable,
        sequence_type=inferred_sequence_type,
    )

    model_validation = validate_model_selection_against_engine_outputs(
        model_selection_workflow.manifest_path
    )
    _copy_output(alignment_workflow.output_paths["alignment"], final_outputs["alignment"])
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
    write_fasta_to_tree_model_table(final_outputs["model_table"], model_rows)
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
        )
    )
    notes = [
        "final tree path contains the bootstrap-supported inference tree",
        "engine-specific intermediate artifacts remain under engine-artifacts/",
    ]
    report = FastaToTreeWorkflowReport(
        input_path=input_path,
        out_dir=out_dir,
        prefix=workflow_prefix,
        sequence_type=inferred_sequence_type,
        selected_model=model_selection_workflow.selected_model,
        engine_artifact_dir=out_dir / "engine-artifacts" / workflow_prefix,
        manifest_path=final_outputs["manifest"],
        output_paths=final_outputs,
        step_manifests={
            "alignment": alignment_workflow.manifest_path,
            "trimming": trimming_workflow.manifest_path,
            "model_selection": model_selection_workflow.manifest_path,
            "maximum_likelihood": maximum_likelihood_workflow.manifest_path,
            "bootstrap_support": bootstrap_workflow.manifest_path,
        },
        input_checksums=build_file_checksums([input_path]),
        output_checksums={},
        alignment_workflow=alignment_workflow,
        trimming_workflow=trimming_workflow,
        model_selection_workflow=model_selection_workflow,
        maximum_likelihood_workflow=maximum_likelihood_workflow,
        bootstrap_workflow=bootstrap_workflow,
        model_validation=model_validation,
        support_summary=support_summary,
        model_rows=model_rows,
        support_rows=support_rows,
        warnings=warnings,
        notes=notes,
    )
    write_fasta_to_tree_log(final_outputs["log"], report)
    report.output_checksums = build_file_checksums(
        [
            final_outputs["alignment"],
            final_outputs["trimmed_alignment"],
            final_outputs["tree"],
            final_outputs["log"],
            final_outputs["model_table"],
            final_outputs["support_table"],
        ]
    )
    write_engine_manifest(report.manifest_path, report)
    return report
