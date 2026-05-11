from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.alignment import AlignmentAlphabet, AlignmentRecord
from bijux_phylogenetics.io.fasta import infer_alignment_alphabet

from .validation import BootstrapSupportSummaryReport, ModelSelectionValidationReport
from .workflows import EngineWorkflowReport

__all__ = [
    "FastaToTreeModelRow",
    "FastaToTreeSupportRow",
    "build_fasta_to_tree_model_rows",
    "build_fasta_to_tree_support_rows",
    "infer_unaligned_sequence_type",
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
