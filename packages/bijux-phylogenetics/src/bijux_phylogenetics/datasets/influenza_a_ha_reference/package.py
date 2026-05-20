from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from bijux_phylogenetics.engines.fasta_to_tree import (
    FastaToTreeWorkflowReport,
    run_fasta_to_tree_workflow,
)
from bijux_phylogenetics.io.fasta.records import validate_fasta_input

_DATASET_ID = "influenza_a_ha_reference_panel"
_DATASET_LABEL = "Influenza A hemagglutinin reference panel"
_SEQUENCE_TYPE = "dna"
_WORKFLOW_PREFIX = "influenza-a-ha-reference-panel"
_IQTREE_SEED = 1
_IQTREE_THREADS = 1
_BOOTSTRAP_REPLICATES = 1000
_SOURCE_ACCESSIONS = (
    "NC_002017.1",
    "CY033655.1",
    "CY046787.1",
    "NC_007366.1",
    "NC_007374.1",
    "AY653200.1",
)


@dataclass(slots=True)
class InfluenzaAHAReferenceDataset:
    """Packaged viral FASTA dataset for sequence-to-tree workflow review."""

    dataset_id: str
    label: str
    dataset_root: Path
    sequences_path: Path
    reference_output_root: Path
    sequence_count: int
    sequence_type: str
    workflow_prefix: str
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    source_accessions: tuple[str, ...]
    source_summary: str


@dataclass(slots=True)
class InfluenzaAHAReferenceDatasetExportResult:
    """Materialized copy of the packaged viral dataset."""

    output_root: Path
    readme_path: Path
    sequences_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class InfluenzaAHAReferenceWorkflowReport:
    """One governed FASTA-to-tree workflow run over the packaged viral dataset."""

    dataset: InfluenzaAHAReferenceDataset
    workflow: FastaToTreeWorkflowReport


@dataclass(slots=True)
class InfluenzaAHAReferenceWorkflowBundle:
    """Written sequence-to-tree outputs for the packaged viral dataset."""

    output_root: Path
    selected_model: str
    minimum_support: float | None
    maximum_support: float | None
    median_support: float | None
    weakly_supported_clade_count: int
    summary_path: Path
    alignment_path: Path
    trimmed_alignment_path: Path
    tree_path: Path
    model_table_path: Path
    support_table_path: Path
    log_path: Path
    manifest_path: Path
    engine_artifact_root: Path


@dataclass(slots=True)
class InfluenzaAHAReferenceDemoResult:
    """Dataset export plus workflow outputs for the public viral demo."""

    output_root: Path
    dataset: InfluenzaAHAReferenceDataset
    dataset_export: InfluenzaAHAReferenceDatasetExportResult
    workflow_bundle: InfluenzaAHAReferenceWorkflowBundle
    overview_path: Path


def load_influenza_a_ha_reference_dataset() -> InfluenzaAHAReferenceDataset:
    """Expose the packaged influenza A HA panel as a first-class runtime surface."""
    dataset_root = _resource_root()
    sequences_path = dataset_root / "sequences.fasta"
    validation = validate_fasta_input(sequences_path, sequence_type=_SEQUENCE_TYPE)
    return InfluenzaAHAReferenceDataset(
        dataset_id=_DATASET_ID,
        label=_DATASET_LABEL,
        dataset_root=dataset_root,
        sequences_path=sequences_path,
        reference_output_root=dataset_root / "expected",
        sequence_count=validation.summary.sequence_count,
        sequence_type=_SEQUENCE_TYPE,
        workflow_prefix=_WORKFLOW_PREFIX,
        iqtree_seed=_IQTREE_SEED,
        iqtree_threads=_IQTREE_THREADS,
        bootstrap_replicates=_BOOTSTRAP_REPLICATES,
        source_accessions=_SOURCE_ACCESSIONS,
        source_summary=(
            "Published influenza A hemagglutinin segment-4 panel assembled from "
            "stable NCBI GenBank and RefSeq accessions spanning H1N1, H2N2, "
            "H3N2, and H5N1 lineages."
        ),
    )


def export_influenza_a_ha_reference_dataset(
    destination: Path,
) -> InfluenzaAHAReferenceDatasetExportResult:
    """Copy the packaged influenza A HA dataset and reference outputs to one directory."""
    dataset = load_influenza_a_ha_reference_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = shutil.copy2(
        dataset.dataset_root / "README.md", destination / "README.md"
    )
    sequences_path = shutil.copy2(
        dataset.sequences_path, destination / "sequences.fasta"
    )
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return InfluenzaAHAReferenceDatasetExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        sequences_path=Path(sequences_path),
        expected_output_root=expected_output_root,
    )


def run_influenza_a_ha_reference_workflow(
    out_dir: Path,
    *,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = _IQTREE_SEED,
    iqtree_threads: int = _IQTREE_THREADS,
    bootstrap_replicates: int = _BOOTSTRAP_REPLICATES,
) -> InfluenzaAHAReferenceWorkflowReport:
    """Run the owned FASTA-to-tree workflow over the packaged influenza A HA panel."""
    dataset = load_influenza_a_ha_reference_dataset()
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
    return InfluenzaAHAReferenceWorkflowReport(dataset=dataset, workflow=workflow)


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
        output_root / "workflow-summary.tsv", report
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


def run_influenza_a_ha_reference_demo(
    output_root: Path,
    *,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = _IQTREE_SEED,
    iqtree_threads: int = _IQTREE_THREADS,
    bootstrap_replicates: int = _BOOTSTRAP_REPLICATES,
) -> InfluenzaAHAReferenceDemoResult:
    """Materialize the packaged viral dataset and rerun the governed workflow bundle."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    dataset = load_influenza_a_ha_reference_dataset()
    dataset_export = export_influenza_a_ha_reference_dataset(output_root / "dataset")
    with TemporaryDirectory(prefix="influenza-a-ha-reference-") as temporary_root:
        workflow_report = run_influenza_a_ha_reference_workflow(
            Path(temporary_root),
            mafft_executable=mafft_executable,
            trimal_executable=trimal_executable,
            iqtree_executable=iqtree_executable,
            iqtree_seed=iqtree_seed,
            iqtree_threads=iqtree_threads,
            bootstrap_replicates=bootstrap_replicates,
        )
        workflow_bundle = write_influenza_a_ha_reference_workflow_bundle(
            output_root / "workflow",
            workflow_report,
        )
    overview_path = _write_overview(
        output_root / "overview.md", dataset, workflow_bundle
    )
    return InfluenzaAHAReferenceDemoResult(
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


def _write_overview(
    path: Path,
    dataset: InfluenzaAHAReferenceDataset,
    workflow_bundle: InfluenzaAHAReferenceWorkflowBundle,
) -> Path:
    lines = [
        "# Influenza A HA Reference Demo",
        "",
        f"- dataset id: `{dataset.dataset_id}`",
        f"- sequence count: `{dataset.sequence_count}`",
        f"- sequence type: `{dataset.sequence_type}`",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{workflow_bundle.summary_path.name}`",
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
        Path(__file__).resolve().parent.parent.parent
        / "resources"
        / "datasets"
        / "viruses"
        / _DATASET_ID
    )
