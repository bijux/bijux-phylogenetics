from __future__ import annotations

from pathlib import Path

from .contracts import (
    FastaToTreeModelRow,
    FastaToTreeSupportRow,
    FastaToTreeWorkflowReport,
)
from .workflow_layout import _display_command, _display_path, _write_tsv


def _serialize_support_taxa(taxa: tuple[str, ...]) -> str:
    return ",".join(taxa)


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
