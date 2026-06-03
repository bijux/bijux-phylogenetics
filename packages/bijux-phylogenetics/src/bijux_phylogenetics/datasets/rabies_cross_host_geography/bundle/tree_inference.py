from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
import shutil

from bijux_phylogenetics.datasets.rabies_cross_host_geography.models import (
    RabiesCrossHostGeographyPanelWorkflowReport,
)
from bijux_phylogenetics.datasets.rabies_cross_host_geography.workflow.tree_transforms import (
    _stabilize_clade_report,
)
from bijux_phylogenetics.phylo.topology import write_tree_rooting_report
from bijux_phylogenetics.trees import (
    CladeTableReport,
    extract_tree_clades,
    write_clade_table,
)

from .input_artifacts import (
    _copy_output,
    _write_alignment_quality_table,
    _write_comparative_branch_repairs_table,
    _write_input_validation_table,
    _write_resolved_workflow_config,
    _write_sequence_ranking_table,
    _write_workflow_config_audit_table,
)


@dataclass(frozen=True)
class TreeInferenceArtifacts:
    config_audit_path: Path
    resolved_config_path: Path
    input_validation_path: Path
    alignment_quality_path: Path
    alignment_sequence_ranking_path: Path
    alignment_path: Path
    trimmed_alignment_path: Path
    tree_path: Path
    rooting_report_path: Path
    model_table_path: Path
    support_table_path: Path
    log_path: Path
    manifest_path: Path
    engine_artifact_root: Path
    clade_report: CladeTableReport
    clade_table_path: Path
    comparative_repairs_path: Path


def _write_tree_inference_artifacts(
    output_root: Path,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
) -> TreeInferenceArtifacts:
    workflow = report.fasta_to_tree
    config_audit_path = _write_workflow_config_audit_table(
        output_root / "workflow-config-audit.tsv",
        report.config_audit_rows,
    )
    resolved_config_path = _write_resolved_workflow_config(
        output_root / "workflow-config.resolved.json",
        report.config,
    )
    input_validation_path = _write_input_validation_table(
        output_root / "input-validation.tsv",
        workflow=workflow,
    )
    alignment_quality_path = _write_alignment_quality_table(
        output_root / "alignment-quality.tsv",
        aligned=report.aligned_quality,
        trimmed=report.trimmed_quality,
    )
    alignment_sequence_ranking_path = _write_sequence_ranking_table(
        output_root / "alignment-sequence-ranking.tsv",
        report.trimmed_sequence_ranking,
    )
    alignment_path = _copy_output(
        workflow.output_paths["alignment"],
        output_root / workflow.output_paths["alignment"].name,
    )
    trimmed_alignment_path = _copy_output(
        workflow.output_paths["trimmed_alignment"],
        output_root / workflow.output_paths["trimmed_alignment"].name,
    )
    tree_path = _copy_output(
        report.rooted_tree_path,
        output_root / report.rooted_tree_path.name,
    )
    stable_rooting_report = replace(
        report.rooting_report,
        tree_path=Path(tree_path.name),
    )
    rooting_report_path = write_tree_rooting_report(
        output_root / f"{report.dataset.workflow_prefix}.rooting.tsv",
        stable_rooting_report,
    )
    model_table_path = _copy_output(
        workflow.output_paths["model_table"],
        output_root / workflow.output_paths["model_table"].name,
    )
    support_table_path = _copy_output(
        workflow.output_paths["support_table"],
        output_root / workflow.output_paths["support_table"].name,
    )
    log_path = _copy_output(
        workflow.output_paths["log"],
        output_root / workflow.output_paths["log"].name,
    )
    manifest_path = _copy_output(
        workflow.manifest_path,
        output_root / workflow.manifest_path.name,
    )
    engine_artifact_root = (
        output_root / "engine-artifacts" / report.dataset.workflow_prefix
    )
    shutil.copytree(workflow.engine_artifact_dir, engine_artifact_root)

    clade_report = extract_tree_clades(
        report.rooted_tree_path,
        metadata_path=report.dataset.metadata_path,
        taxon_column="taxon",
        metadata_columns=list(report.dataset.clade_metadata_columns),
    )
    stable_clade_report = _stabilize_clade_report(
        clade_report,
        stable_source_path=Path(tree_path.name),
    )
    clade_table_path = write_clade_table(
        output_root / "clade-table.tsv",
        stable_clade_report,
    )
    comparative_repairs_path = _write_comparative_branch_repairs_table(
        output_root / "comparative-tree-adjustments.tsv",
        report.comparative_branch_repairs,
    )
    return TreeInferenceArtifacts(
        config_audit_path=config_audit_path,
        resolved_config_path=resolved_config_path,
        input_validation_path=input_validation_path,
        alignment_quality_path=alignment_quality_path,
        alignment_sequence_ranking_path=alignment_sequence_ranking_path,
        alignment_path=alignment_path,
        trimmed_alignment_path=trimmed_alignment_path,
        tree_path=tree_path,
        rooting_report_path=rooting_report_path,
        model_table_path=model_table_path,
        support_table_path=support_table_path,
        log_path=log_path,
        manifest_path=manifest_path,
        engine_artifact_root=engine_artifact_root,
        clade_report=stable_clade_report,
        clade_table_path=clade_table_path,
        comparative_repairs_path=comparative_repairs_path,
    )
