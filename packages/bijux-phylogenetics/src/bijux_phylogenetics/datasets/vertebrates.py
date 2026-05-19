from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from bijux_phylogenetics.engines.fasta_to_tree import (
    FastaToTreeWorkflowReport,
    run_fasta_to_tree_workflow,
)
from bijux_phylogenetics.io.fasta import (
    load_fasta_alignment,
    load_fasta_records,
)
from bijux_phylogenetics.io.fasta.records import validate_fasta_input

_DATASET_ID = "gnathostome_ortholog_protein_benchmark"
_DATASET_LABEL = "Gnathostome ortholog protein benchmark"
_SEQUENCE_TYPE = "protein"
_WORKFLOW_PREFIX = "gnathostome-ortholog-protein-benchmark"
_IQTREE_SEED = 1
_IQTREE_THREADS = 1
_BOOTSTRAP_REPLICATES = 1000
_SOURCE_REFERENCE = "trimAl governed reference corpus example.009.AA"
_SOURCE_TRANSFORMATION = (
    "Removed alignment gap characters and placeholder missing-data marks from "
    "the aligned reference FASTA to recover the raw protein inputs used for "
    "the packaged end-to-end benchmark."
)


@dataclass(slots=True)
class GnathostomeOrthologProteinBenchmarkDataset:
    """Packaged protein FASTA benchmark for end-to-end phylogeny review."""

    dataset_id: str
    label: str
    dataset_root: Path
    readme_path: Path
    sequences_path: Path
    reference_output_root: Path
    sequence_count: int
    sequence_type: str
    workflow_prefix: str
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    source_reference: str
    source_transformation: str
    minimum_sequence_length: int
    maximum_sequence_length: int
    source_summary: str


@dataclass(slots=True)
class GnathostomeOrthologProteinBenchmarkExportResult:
    """Materialized copy of the packaged protein benchmark dataset."""

    output_root: Path
    readme_path: Path
    sequences_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class GnathostomeOrthologProteinBenchmarkWorkflowReport:
    """One governed protein FASTA-to-tree workflow run over the packaged panel."""

    dataset: GnathostomeOrthologProteinBenchmarkDataset
    workflow: FastaToTreeWorkflowReport


@dataclass(slots=True)
class GnathostomeOrthologProteinBenchmarkWorkflowBundle:
    """Written alignment, inference, support, and assumption outputs."""

    output_root: Path
    selected_model: str
    sequence_count: int
    alignment_length: int
    trimmed_alignment_length: int
    minimum_support: float | None
    maximum_support: float | None
    median_support: float | None
    weakly_supported_clade_count: int
    summary_path: Path
    assumptions_path: Path
    alignment_path: Path
    trimmed_alignment_path: Path
    tree_path: Path
    model_table_path: Path
    support_table_path: Path
    log_path: Path
    manifest_path: Path
    engine_artifact_root: Path


@dataclass(slots=True)
class GnathostomeOrthologProteinBenchmarkDemoResult:
    """Dataset export plus rerun workflow outputs for the public protein demo."""

    output_root: Path
    dataset: GnathostomeOrthologProteinBenchmarkDataset
    dataset_export: GnathostomeOrthologProteinBenchmarkExportResult
    workflow_bundle: GnathostomeOrthologProteinBenchmarkWorkflowBundle
    overview_path: Path


def load_gnathostome_ortholog_protein_benchmark_dataset() -> (
    GnathostomeOrthologProteinBenchmarkDataset
):
    """Expose the packaged gnathostome protein panel as a first-class surface."""
    dataset_root = _resource_root()
    sequences_path = dataset_root / "sequences.fasta"
    validation = validate_fasta_input(sequences_path, sequence_type=_SEQUENCE_TYPE)
    sequence_lengths = [
        len(record.sequence) for record in load_fasta_records(sequences_path)
    ]
    return GnathostomeOrthologProteinBenchmarkDataset(
        dataset_id=_DATASET_ID,
        label=_DATASET_LABEL,
        dataset_root=dataset_root,
        readme_path=dataset_root / "README.md",
        sequences_path=sequences_path,
        reference_output_root=dataset_root / "expected",
        sequence_count=validation.summary.sequence_count,
        sequence_type=_SEQUENCE_TYPE,
        workflow_prefix=_WORKFLOW_PREFIX,
        iqtree_seed=_IQTREE_SEED,
        iqtree_threads=_IQTREE_THREADS,
        bootstrap_replicates=_BOOTSTRAP_REPLICATES,
        source_reference=_SOURCE_REFERENCE,
        source_transformation=_SOURCE_TRANSFORMATION,
        minimum_sequence_length=min(sequence_lengths),
        maximum_sequence_length=max(sequence_lengths),
        source_summary=(
            "Nine ungapped gnathostome ortholog proteins packaged as one public "
            "amino-acid benchmark for MAFFT alignment, trimAl trimming, IQ-TREE "
            "protein model selection, maximum-likelihood inference, and "
            "bootstrap support review."
        ),
    )


def export_gnathostome_ortholog_protein_benchmark_dataset(
    destination: Path,
) -> GnathostomeOrthologProteinBenchmarkExportResult:
    """Copy the packaged protein benchmark dataset and governed outputs."""
    dataset = load_gnathostome_ortholog_protein_benchmark_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = Path(shutil.copy2(dataset.readme_path, destination / "README.md"))
    sequences_path = Path(
        shutil.copy2(dataset.sequences_path, destination / "sequences.fasta")
    )
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return GnathostomeOrthologProteinBenchmarkExportResult(
        output_root=destination,
        readme_path=readme_path,
        sequences_path=sequences_path,
        expected_output_root=expected_output_root,
    )


def run_gnathostome_ortholog_protein_benchmark_workflow(
    out_dir: Path,
    *,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = _IQTREE_SEED,
    iqtree_threads: int = _IQTREE_THREADS,
    bootstrap_replicates: int = _BOOTSTRAP_REPLICATES,
) -> GnathostomeOrthologProteinBenchmarkWorkflowReport:
    """Run the owned protein FASTA-to-tree workflow over the packaged panel."""
    dataset = load_gnathostome_ortholog_protein_benchmark_dataset()
    workflow = run_fasta_to_tree_workflow(
        dataset.sequences_path,
        out_dir=out_dir,
        prefix=dataset.workflow_prefix,
        sequence_type=dataset.sequence_type,
        mafft_executable=mafft_executable,
        trimal_executable=trimal_executable,
        iqtree_executable=iqtree_executable,
        iqtree_seed=iqtree_seed,
        iqtree_threads=iqtree_threads,
        bootstrap_replicates=bootstrap_replicates,
    )
    return GnathostomeOrthologProteinBenchmarkWorkflowReport(
        dataset=dataset,
        workflow=workflow,
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


def run_gnathostome_ortholog_protein_benchmark_demo(
    output_root: Path,
    *,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = _IQTREE_SEED,
    iqtree_threads: int = _IQTREE_THREADS,
    bootstrap_replicates: int = _BOOTSTRAP_REPLICATES,
) -> GnathostomeOrthologProteinBenchmarkDemoResult:
    """Materialize the packaged protein benchmark and rerun its governed outputs."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    dataset = load_gnathostome_ortholog_protein_benchmark_dataset()
    dataset_export = export_gnathostome_ortholog_protein_benchmark_dataset(
        output_root / "dataset"
    )
    with TemporaryDirectory(prefix="gnathostome-ortholog-protein-benchmark-") as root:
        workflow_report = run_gnathostome_ortholog_protein_benchmark_workflow(
            Path(root),
            mafft_executable=mafft_executable,
            trimal_executable=trimal_executable,
            iqtree_executable=iqtree_executable,
            iqtree_seed=iqtree_seed,
            iqtree_threads=iqtree_threads,
            bootstrap_replicates=bootstrap_replicates,
        )
        workflow_bundle = write_gnathostome_ortholog_protein_benchmark_workflow_bundle(
            output_root / "workflow",
            workflow_report,
        )
    overview_path = _write_overview(
        output_root / "overview.md", dataset, workflow_bundle
    )
    return GnathostomeOrthologProteinBenchmarkDemoResult(
        output_root=output_root,
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
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


def _write_overview(
    path: Path,
    dataset: GnathostomeOrthologProteinBenchmarkDataset,
    workflow_bundle: GnathostomeOrthologProteinBenchmarkWorkflowBundle,
) -> Path:
    lines = [
        "# Gnathostome Ortholog Protein Benchmark",
        "",
        f"- dataset id: `{dataset.dataset_id}`",
        f"- sequence count: `{dataset.sequence_count}`",
        f"- sequence type: `{dataset.sequence_type}`",
        f"- selected protein model: `{workflow_bundle.selected_model}`",
        "",
        "Protein-specific assumptions:",
        "",
        "- the workflow starts from amino-acid FASTA and does not translate coding DNA",
        "- IQ-TREE model selection is restricted to protein models with `-st AA`",
        "- bootstrap support values summarize the trimmed amino-acid alignment",
        "- nucleotide-specific assumptions such as codon position and GC interpretation are not used",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{workflow_bundle.summary_path.name}`",
        f"- molecular assumptions: `{workflow_bundle.assumptions_path.name}`",
        f"- final supported tree: `{workflow_bundle.tree_path.name}`",
        f"- selected model table: `{workflow_bundle.model_table_path.name}`",
        f"- branch support table: `{workflow_bundle.support_table_path.name}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _format_number(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".12g")


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent
        / "resources"
        / "datasets"
        / "vertebrates"
        / _DATASET_ID
    )
