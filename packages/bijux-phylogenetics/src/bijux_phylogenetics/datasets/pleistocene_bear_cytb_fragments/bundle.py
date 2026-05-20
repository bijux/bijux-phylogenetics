from __future__ import annotations

from pathlib import Path
import shutil

from bijux_phylogenetics.io.fasta import write_fasta_alignment

from .models import (
    PleistoceneBearCytbFragmentWorkflowBundle,
    PleistoceneBearCytbFragmentWorkflowReport,
    PleistoceneBearMissingnessEffectRow,
)


def write_pleistocene_bear_cytb_fragment_workflow_bundle(
    output_root: Path,
    report: PleistoceneBearCytbFragmentWorkflowReport,
) -> PleistoceneBearCytbFragmentWorkflowBundle:
    """Write the governed degraded-sequence outputs for the packaged bear panel."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    workflow = report.workflow
    summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report,
    )
    missingness_effects_path = _write_missingness_effect_table(
        output_root / "missingness-effects.tsv",
        report.missingness_rows,
    )
    alignment_path = _copy_output(
        workflow.output_paths["alignment"],
        output_root / f"{report.dataset.workflow_prefix}.aln",
    )
    trimmed_alignment_path = _copy_output(
        workflow.output_paths["trimmed_alignment"],
        output_root / f"{report.dataset.workflow_prefix}.trimmed.aln",
    )
    cleaned_alignment_path = write_fasta_alignment(
        output_root / f"{report.dataset.workflow_prefix}.cleaned.aln",
        report.cleaned_records,
    )
    tree_path = _copy_output(
        workflow.output_paths["tree"],
        output_root / f"{report.dataset.workflow_prefix}.tree",
    )
    model_table_path = _copy_output(
        workflow.output_paths["model_table"],
        output_root / f"{report.dataset.workflow_prefix}.model.tsv",
    )
    support_table_path = _copy_output(
        workflow.output_paths["support_table"],
        output_root / f"{report.dataset.workflow_prefix}.support.tsv",
    )
    log_path = _copy_output(
        workflow.output_paths["log"],
        output_root / f"{report.dataset.workflow_prefix}.log",
    )
    manifest_path = _copy_output(
        workflow.manifest_path,
        output_root / f"{report.dataset.workflow_prefix}.manifest.json",
    )
    engine_artifact_root = (
        output_root / "engine-artifacts" / report.dataset.workflow_prefix
    )
    shutil.copytree(workflow.engine_artifact_dir, engine_artifact_root)
    return PleistoceneBearCytbFragmentWorkflowBundle(
        output_root=output_root,
        selected_model=workflow.selected_model,
        minimum_support=workflow.support_summary.minimum_support,
        maximum_support=workflow.support_summary.maximum_support,
        median_support=workflow.support_summary.median_support,
        weakly_supported_clade_count=workflow.support_summary.weakly_supported_clade_count,
        removed_column_count=len(report.missingness_cleanup.removed_columns),
        removed_sequence_count=len(report.missingness_cleanup.removed_sequences),
        cleaned_missing_data_fraction=report.cleaned_summary.missing_data_fraction,
        summary_path=summary_path,
        missingness_effects_path=missingness_effects_path,
        alignment_path=alignment_path,
        trimmed_alignment_path=trimmed_alignment_path,
        cleaned_alignment_path=cleaned_alignment_path,
        tree_path=tree_path,
        model_table_path=model_table_path,
        support_table_path=support_table_path,
        log_path=log_path,
        manifest_path=manifest_path,
        engine_artifact_root=engine_artifact_root,
    )


def _write_workflow_summary_table(
    path: Path,
    report: PleistoceneBearCytbFragmentWorkflowReport,
) -> Path:
    workflow = report.workflow
    support = workflow.support_summary
    rows = [
        "\t".join(
            [
                "dataset_id",
                "sequence_count",
                "degraded_sequence_count",
                "selected_model",
                "internal_node_count",
                "supported_node_count",
                "minimum_support",
                "maximum_support",
                "median_support",
                "aligned_missing_data_fraction",
                "cleaned_missing_data_fraction",
                "removed_column_count",
                "removed_sequence_count",
            ]
        ),
        "\t".join(
            [
                report.dataset.dataset_id,
                str(report.dataset.sequence_count),
                str(len(report.dataset.degraded_sequence_ids)),
                workflow.selected_model,
                str(support.internal_node_count),
                str(support.supported_node_count),
                _format_number(support.minimum_support),
                _format_number(support.maximum_support),
                _format_number(support.median_support),
                _format_number(report.aligned_summary.missing_data_fraction),
                _format_number(report.cleaned_summary.missing_data_fraction),
                str(len(report.missingness_cleanup.removed_columns)),
                str(len(report.missingness_cleanup.removed_sequences)),
            ]
        ),
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _write_missingness_effect_table(
    path: Path,
    rows: list[PleistoceneBearMissingnessEffectRow],
) -> Path:
    lines = [
        "\t".join(
            [
                "identifier",
                "raw_sequence_length",
                "degraded_sequence",
                "aligned_missing_fraction",
                "engine_trimmed_missing_fraction",
                "cleaned_missing_fraction",
                "removed_by_missingness_cleanup",
            ]
        )
    ]
    lines.extend(
        "\t".join(
            [
                row.identifier,
                str(row.raw_sequence_length),
                "true" if row.degraded_sequence else "false",
                _format_number(row.aligned_missing_fraction),
                _format_number(row.engine_trimmed_missing_fraction),
                _format_number(row.cleaned_missing_fraction),
                "true" if row.removed_by_missingness_cleanup else "false",
            ]
        )
        for row in rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _copy_output(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    return Path(shutil.copy2(source, destination))


def _format_number(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".12g")
