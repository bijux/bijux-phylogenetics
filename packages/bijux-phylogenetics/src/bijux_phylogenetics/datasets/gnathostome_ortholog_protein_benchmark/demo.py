from __future__ import annotations

from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from .bundle import write_gnathostome_ortholog_protein_benchmark_workflow_bundle
from .export import export_gnathostome_ortholog_protein_benchmark_dataset
from .models import (
    GnathostomeOrthologProteinBenchmarkDataset,
    GnathostomeOrthologProteinBenchmarkDemoResult,
    GnathostomeOrthologProteinBenchmarkWorkflowBundle,
)
from .panel import (
    BOOTSTRAP_REPLICATES,
    IQTREE_SEED,
    IQTREE_THREADS,
    load_gnathostome_ortholog_protein_benchmark_dataset,
)
from .workflow import run_gnathostome_ortholog_protein_benchmark_workflow


def run_gnathostome_ortholog_protein_benchmark_demo(
    output_root: Path,
    *,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = IQTREE_SEED,
    iqtree_threads: int = IQTREE_THREADS,
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
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
        output_root / "overview.md",
        dataset,
        workflow_bundle,
    )
    return GnathostomeOrthologProteinBenchmarkDemoResult(
        output_root=output_root,
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


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
