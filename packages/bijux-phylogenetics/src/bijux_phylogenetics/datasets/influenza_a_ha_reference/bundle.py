from __future__ import annotations

from pathlib import Path
import shutil

from .models import (
    InfluenzaAHAReferenceWorkflowBundle,
    InfluenzaAHAReferenceWorkflowReport,
)


def write_influenza_a_ha_reference_workflow_bundle(
    output_root: Path,
    report: InfluenzaAHAReferenceWorkflowReport,
) -> InfluenzaAHAReferenceWorkflowBundle:
    """Write the governed viral sequence-to-tree outputs for the packaged dataset."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    workflow = report.workflow
    summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report,
    )
    alignment_path = _copy_output(
        workflow.output_paths["alignment"],
        output_root / f"{report.dataset.workflow_prefix}.aln",
    )
    trimmed_alignment_path = _copy_output(
        workflow.output_paths["trimmed_alignment"],
        output_root / f"{report.dataset.workflow_prefix}.trimmed.aln",
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
    return InfluenzaAHAReferenceWorkflowBundle(
        output_root=output_root,
        selected_model=workflow.selected_model,
        minimum_support=workflow.support_summary.minimum_support,
        maximum_support=workflow.support_summary.maximum_support,
        median_support=workflow.support_summary.median_support,
        weakly_supported_clade_count=workflow.support_summary.weakly_supported_clade_count,
        summary_path=summary_path,
        alignment_path=alignment_path,
        trimmed_alignment_path=trimmed_alignment_path,
        tree_path=tree_path,
        model_table_path=model_table_path,
        support_table_path=support_table_path,
        log_path=log_path,
        manifest_path=manifest_path,
        engine_artifact_root=engine_artifact_root,
    )


def _copy_output(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    return Path(shutil.copy2(source, destination))


def _write_workflow_summary_table(
    path: Path,
    report: InfluenzaAHAReferenceWorkflowReport,
) -> Path:
    workflow = report.workflow
    support = workflow.support_summary
    rows = [
        "dataset_id\tsequence_count\tsequence_type\tselected_model\tinternal_node_count\tsupported_node_count\tminimum_support\tmaximum_support\tmedian_support\tweakly_supported_clade_count",
        "\t".join(
            [
                report.dataset.dataset_id,
                str(report.dataset.sequence_count),
                report.dataset.sequence_type,
                workflow.selected_model,
                str(support.internal_node_count),
                str(support.supported_node_count),
                _format_number(support.minimum_support),
                _format_number(support.maximum_support),
                _format_number(support.median_support),
                str(support.weakly_supported_clade_count),
            ]
        ),
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _format_number(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".12g")
