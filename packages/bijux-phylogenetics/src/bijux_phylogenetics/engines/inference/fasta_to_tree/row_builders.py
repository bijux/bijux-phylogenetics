from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta import detect_fasta_sequence_type
from bijux_phylogenetics.phylo.alignment import AlignmentAlphabet, AlignmentRecord

from ...validation import ModelSelectionValidationReport
from ...workflows.models import EngineWorkflowReport
from .contracts import FastaToTreeModelRow, FastaToTreeSupportRow


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
    support_summary,
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
