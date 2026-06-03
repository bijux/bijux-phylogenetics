from __future__ import annotations

from pathlib import Path
import shutil

from bijux_phylogenetics.io.fasta import load_fasta_alignment

from .models import (
    GnathostomeOrthologProteinBenchmarkWorkflowBundle,
    GnathostomeOrthologProteinBenchmarkWorkflowReport,
)


def write_gnathostome_ortholog_protein_benchmark_workflow_bundle(
    output_root: Path,
    report: GnathostomeOrthologProteinBenchmarkWorkflowReport,
) -> GnathostomeOrthologProteinBenchmarkWorkflowBundle:
    """Write the governed public benchmark bundle for the protein workflow."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    workflow = report.workflow
    aligned_records = load_fasta_alignment(workflow.output_paths["alignment"])
    trimmed_records = load_fasta_alignment(workflow.output_paths["trimmed_alignment"])
    summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report,
        alignment_length=len(aligned_records[0].sequence),
        trimmed_alignment_length=len(trimmed_records[0].sequence),
    )
    assumptions_path = _write_molecular_assumption_table(
        output_root / "molecular-assumptions.tsv"
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
    return GnathostomeOrthologProteinBenchmarkWorkflowBundle(
        output_root=output_root,
        selected_model=workflow.selected_model,
        sequence_count=report.dataset.sequence_count,
        alignment_length=len(aligned_records[0].sequence),
        trimmed_alignment_length=len(trimmed_records[0].sequence),
        minimum_support=workflow.support_summary.minimum_support,
        maximum_support=workflow.support_summary.maximum_support,
        median_support=workflow.support_summary.median_support,
        weakly_supported_clade_count=workflow.support_summary.weakly_supported_clade_count,
        summary_path=summary_path,
        assumptions_path=assumptions_path,
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
    report: GnathostomeOrthologProteinBenchmarkWorkflowReport,
    *,
    alignment_length: int,
    trimmed_alignment_length: int,
) -> Path:
    support = report.workflow.support_summary
    rows = [
        (
            "dataset_id\tsequence_count\tsequence_type\tselected_model\t"
            "alignment_length\ttrimmed_alignment_length\tinternal_node_count\t"
            "supported_node_count\tminimum_support\tmaximum_support\t"
            "median_support\tweakly_supported_clade_count\tstate_space\t"
            "model_selection_scope"
        ),
        "\t".join(
            [
                report.dataset.dataset_id,
                str(report.dataset.sequence_count),
                report.dataset.sequence_type,
                report.workflow.selected_model,
                str(alignment_length),
                str(trimmed_alignment_length),
                str(support.internal_node_count),
                str(support.supported_node_count),
                _format_number(support.minimum_support),
                _format_number(support.maximum_support),
                _format_number(support.median_support),
                str(support.weakly_supported_clade_count),
                "amino-acid",
                "protein-models-only",
            ]
        ),
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _write_molecular_assumption_table(path: Path) -> Path:
    rows = [
        "assumption_id\tvalue\texplanation",
        (
            "sequence_type\tprotein\tRaw FASTA residues are interpreted as "
            "amino acids rather than nucleotides."
        ),
        (
            "iqtree_sequence_type_keyword\tAA\tIQ-TREE runs with -st AA, so "
            "model search and inference stay in amino-acid state space."
        ),
        (
            "model_selection_scope\tprotein_exchangeability_models\tModelFinder "
            "searches protein substitution models instead of nucleotide models."
        ),
        (
            "translation_required\tfalse\tThe benchmark starts from protein FASTA "
            "directly and does not translate coding DNA."
        ),
        (
            "dna_specific_assumptions_applied\tfalse\tGC content, codon position, "
            "and nucleotide substitution assumptions are not part of this workflow."
        ),
        (
            "branch_support_interpretation\tamino_acid_bootstrap\tBootstrap support "
            "values are estimated from the trimmed amino-acid alignment."
        ),
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _format_number(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".12g")
