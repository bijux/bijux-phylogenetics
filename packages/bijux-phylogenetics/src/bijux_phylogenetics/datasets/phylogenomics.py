from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.engines.workflows.iqtree import (
    run_bootstrap_support_estimation,
    run_model_selection,
)
from bijux_phylogenetics.engines.workflows.models import EngineWorkflowReport
from bijux_phylogenetics.io.fasta import (
    load_fasta_alignment,
    write_fasta_alignment,
)
from bijux_phylogenetics.io.fasta.records import validate_fasta_input
from bijux_phylogenetics.phylo.alignment import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.concatenation import (
    ConcatenatedAlignmentReport,
    concatenate_locus_alignments,
)
from bijux_phylogenetics.phylo.alignment.occupancy import (
    LocusOccupancyReport,
    write_locus_partitions,
)
from bijux_phylogenetics.phylo.alignment.partitions import LocusPartition

_DATASET_ID = "catarrhine_mitogenome_five_locus_panel"
_DATASET_LABEL = "Catarrhine mitogenome five-locus panel"
_SEQUENCE_TYPE = "dna"
_WORKFLOW_PREFIX = "catarrhine-mitogenome-five-locus-panel"
_IQTREE_SEED = 1
_IQTREE_THREADS = 1
_BOOTSTRAP_REPLICATES = 1000
_SOURCE_ACCESSIONS = (
    "NC_012920.1",
    "NC_001643.1",
    "NC_011120.1",
    "NC_001646.1",
    "NC_002082.1",
    "NC_005943.1",
)


@dataclass(slots=True)
class CatarrhineMitogenomeFiveLocusPanelDataset:
    """Packaged multi-locus dataset for concatenation and partitioned inference review."""

    dataset_id: str
    label: str
    dataset_root: Path
    taxa_path: Path
    locus_alignment_root: Path
    reference_output_root: Path
    taxon_count: int
    locus_count: int
    locus_names: tuple[str, ...]
    sequence_type: str
    workflow_prefix: str
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    source_accessions: tuple[str, ...]
    source_summary: str


@dataclass(slots=True)
class CatarrhineMitogenomeFiveLocusPanelExportResult:
    """Materialized copy of the packaged multi-locus dataset."""

    output_root: Path
    readme_path: Path
    taxa_path: Path
    locus_alignment_root: Path
    expected_output_root: Path


@dataclass(slots=True)
class CatarrhineMitogenomeFiveLocusPanelWorkflowReport:
    """One governed multi-locus workflow run over the packaged catarrhine panel."""

    dataset: CatarrhineMitogenomeFiveLocusPanelDataset
    supermatrix_records: list[AlignmentRecord]
    partitions: tuple[LocusPartition, ...]
    concatenation_report: ConcatenatedAlignmentReport
    supermatrix_path: Path
    partition_path: Path
    model_selection: EngineWorkflowReport
    partitioned_support: EngineWorkflowReport


@dataclass(slots=True)
class CatarrhineMitogenomeFiveLocusPanelWorkflowBundle:
    """Written concatenation and partitioned-inference outputs for the panel."""

    output_root: Path
    selected_model: str
    taxon_count: int
    locus_count: int
    alignment_length: int
    partition_count: int
    minimum_support: float | None
    maximum_support: float | None
    weakly_supported_clade_count: int
    workflow_summary_path: Path
    supermatrix_path: Path
    partitions_path: Path
    occupancy_taxa_path: Path
    occupancy_loci_path: Path
    occupancy_matrix_path: Path
    partition_summary_path: Path
    model_candidates_path: Path
    support_tree_path: Path
    support_table_path: Path


@dataclass(slots=True)
class CatarrhineMitogenomeFiveLocusPanelDemoResult:
    """Dataset export plus workflow outputs for the public multi-locus demo."""

    output_root: Path
    dataset: CatarrhineMitogenomeFiveLocusPanelDataset
    dataset_export: CatarrhineMitogenomeFiveLocusPanelExportResult
    workflow_bundle: CatarrhineMitogenomeFiveLocusPanelWorkflowBundle
    overview_path: Path


def load_catarrhine_mitogenome_five_locus_panel_dataset() -> (
    CatarrhineMitogenomeFiveLocusPanelDataset
):
    """Expose the packaged catarrhine multi-locus panel as a first-class surface."""
    dataset_root = _resource_root()
    taxa_path = dataset_root / "taxa.csv"
    locus_alignment_root = dataset_root / "loci"
    alignment_paths = _locus_alignment_paths(locus_alignment_root)
    _validate_locus_taxa(alignment_paths)
    taxon_count = validate_fasta_input(
        alignment_paths[0], sequence_type=_SEQUENCE_TYPE
    ).summary.sequence_count
    return CatarrhineMitogenomeFiveLocusPanelDataset(
        dataset_id=_DATASET_ID,
        label=_DATASET_LABEL,
        dataset_root=dataset_root,
        taxa_path=taxa_path,
        locus_alignment_root=locus_alignment_root,
        reference_output_root=dataset_root / "expected",
        taxon_count=taxon_count,
        locus_count=len(alignment_paths),
        locus_names=tuple(path.stem for path in alignment_paths),
        sequence_type=_SEQUENCE_TYPE,
        workflow_prefix=_WORKFLOW_PREFIX,
        iqtree_seed=_IQTREE_SEED,
        iqtree_threads=_IQTREE_THREADS,
        bootstrap_replicates=_BOOTSTRAP_REPLICATES,
        source_accessions=_SOURCE_ACCESSIONS,
        source_summary=(
            "Five mitochondrial coding-gene alignments extracted from stable "
            "catarrhine mitochondrial genome RefSeq accessions and packaged for "
            "concatenation, occupancy review, and partitioned inference."
        ),
    )


def export_catarrhine_mitogenome_five_locus_panel_dataset(
    destination: Path,
) -> CatarrhineMitogenomeFiveLocusPanelExportResult:
    """Copy the packaged multi-locus dataset and reference outputs."""
    dataset = load_catarrhine_mitogenome_five_locus_panel_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = shutil.copy2(
        dataset.dataset_root / "README.md", destination / "README.md"
    )
    taxa_path = shutil.copy2(dataset.taxa_path, destination / "taxa.csv")
    locus_alignment_root = destination / "loci"
    shutil.copytree(dataset.locus_alignment_root, locus_alignment_root)
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return CatarrhineMitogenomeFiveLocusPanelExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        taxa_path=Path(taxa_path),
        locus_alignment_root=locus_alignment_root,
        expected_output_root=expected_output_root,
    )


def run_catarrhine_mitogenome_five_locus_panel_workflow(
    out_dir: Path,
    *,
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = _IQTREE_SEED,
    iqtree_threads: int = _IQTREE_THREADS,
    bootstrap_replicates: int = _BOOTSTRAP_REPLICATES,
) -> CatarrhineMitogenomeFiveLocusPanelWorkflowReport:
    """Run the owned multi-locus workflow over the packaged catarrhine panel."""
    dataset = load_catarrhine_mitogenome_five_locus_panel_dataset()
    assembled_root = out_dir / "assembled"
    engine_root = out_dir / "engine"
    assembled_root.mkdir(parents=True, exist_ok=True)
    engine_root.mkdir(parents=True, exist_ok=True)

    alignment_paths = _locus_alignment_paths(dataset.locus_alignment_root)
    supermatrix_records, partitions, concatenation_report = (
        concatenate_locus_alignments(
            alignment_paths,
            concatenated_alignment_path=assembled_root
            / f"{dataset.workflow_prefix}.supermatrix.fasta",
            concatenated_partition_path=assembled_root
            / f"{dataset.workflow_prefix}.partitions.txt",
        )
    )
    supermatrix_path = write_fasta_alignment(
        assembled_root / f"{dataset.workflow_prefix}.supermatrix.fasta",
        supermatrix_records,
    )
    partition_path = write_locus_partitions(
        assembled_root / f"{dataset.workflow_prefix}.partitions.txt",
        partitions,
    )
    model_selection = run_model_selection(
        supermatrix_path,
        out_dir=engine_root / "model-selection",
        prefix=dataset.workflow_prefix,
        executable=iqtree_executable,
        sequence_type=dataset.sequence_type,
        partition_path=partition_path,
        seed=iqtree_seed,
        threads=iqtree_threads,
    )
    partitioned_support = run_bootstrap_support_estimation(
        supermatrix_path,
        out_dir=engine_root / "partitioned-support",
        model=model_selection.selected_model,
        replicates=bootstrap_replicates,
        prefix=dataset.workflow_prefix,
        executable=iqtree_executable,
        sequence_type=dataset.sequence_type,
        partition_path=partition_path,
        seed=iqtree_seed,
        threads=iqtree_threads,
    )
    return CatarrhineMitogenomeFiveLocusPanelWorkflowReport(
        dataset=dataset,
        supermatrix_records=supermatrix_records,
        partitions=partitions,
        concatenation_report=concatenation_report,
        supermatrix_path=supermatrix_path,
        partition_path=partition_path,
        model_selection=model_selection,
        partitioned_support=partitioned_support,
    )


def write_catarrhine_mitogenome_five_locus_panel_workflow_bundle(
    output_root: Path,
    report: CatarrhineMitogenomeFiveLocusPanelWorkflowReport,
) -> CatarrhineMitogenomeFiveLocusPanelWorkflowBundle:
    """Write the governed multi-locus outputs for the packaged catarrhine panel."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    occupancy = report.concatenation_report.occupancy_report
    support = report.partitioned_support.bootstrap_support_summary
    if support is None:
        raise ValueError(
            "partitioned support workflow did not expose bootstrap summary"
        )

    workflow_summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report,
    )
    supermatrix_path = _copy_output(
        report.supermatrix_path,
        output_root / f"{report.dataset.workflow_prefix}.supermatrix.fasta",
    )
    partitions_path = _copy_output(
        report.partition_path,
        output_root / f"{report.dataset.workflow_prefix}.partitions.txt",
    )
    occupancy_taxa_path = _write_occupancy_taxa_table(
        output_root / "occupancy-taxa.tsv",
        occupancy,
    )
    occupancy_loci_path = _write_occupancy_loci_table(
        output_root / "occupancy-loci.tsv",
        occupancy,
    )
    occupancy_matrix_path = _write_occupancy_matrix_table(
        output_root / "occupancy-matrix.tsv",
        occupancy,
    )
    partition_summary_path = _copy_output(
        report.partitioned_support.output_paths["partition_summary"],
        output_root / f"{report.dataset.workflow_prefix}.partition-summary.tsv",
    )
    model_candidates_path = _copy_output(
        report.model_selection.output_paths["model_candidates"],
        output_root / f"{report.dataset.workflow_prefix}.model-candidates.tsv",
    )
    support_tree_path = _copy_output(
        report.partitioned_support.output_paths["support_tree"],
        output_root / f"{report.dataset.workflow_prefix}.supported.tree",
    )
    support_table_path = _copy_output(
        report.partitioned_support.output_paths["support_table"],
        output_root / f"{report.dataset.workflow_prefix}.support.tsv",
    )
    return CatarrhineMitogenomeFiveLocusPanelWorkflowBundle(
        output_root=output_root,
        selected_model=report.model_selection.selected_model,
        taxon_count=report.concatenation_report.taxon_count,
        locus_count=report.concatenation_report.locus_count,
        alignment_length=report.concatenation_report.alignment_length,
        partition_count=len(report.partitions),
        minimum_support=support.minimum_support,
        maximum_support=support.maximum_support,
        weakly_supported_clade_count=support.weakly_supported_clade_count,
        workflow_summary_path=workflow_summary_path,
        supermatrix_path=supermatrix_path,
        partitions_path=partitions_path,
        occupancy_taxa_path=occupancy_taxa_path,
        occupancy_loci_path=occupancy_loci_path,
        occupancy_matrix_path=occupancy_matrix_path,
        partition_summary_path=partition_summary_path,
        model_candidates_path=model_candidates_path,
        support_tree_path=support_tree_path,
        support_table_path=support_table_path,
    )


def run_catarrhine_mitogenome_five_locus_panel_demo(
    output_root: Path,
    *,
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = _IQTREE_SEED,
    iqtree_threads: int = _IQTREE_THREADS,
    bootstrap_replicates: int = _BOOTSTRAP_REPLICATES,
) -> CatarrhineMitogenomeFiveLocusPanelDemoResult:
    """Materialize the packaged multi-locus dataset and rerun the workflow."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    dataset = load_catarrhine_mitogenome_five_locus_panel_dataset()
    dataset_export = export_catarrhine_mitogenome_five_locus_panel_dataset(
        output_root / "dataset"
    )
    with TemporaryDirectory(
        prefix="catarrhine-mitogenome-five-locus-"
    ) as temporary_root:
        workflow_report = run_catarrhine_mitogenome_five_locus_panel_workflow(
            Path(temporary_root),
            iqtree_executable=iqtree_executable,
            iqtree_seed=iqtree_seed,
            iqtree_threads=iqtree_threads,
            bootstrap_replicates=bootstrap_replicates,
        )
        workflow_bundle = write_catarrhine_mitogenome_five_locus_panel_workflow_bundle(
            output_root / "workflow",
            workflow_report,
        )
    overview_path = _write_overview(
        output_root / "overview.md", dataset, workflow_bundle
    )
    return CatarrhineMitogenomeFiveLocusPanelDemoResult(
        output_root=output_root,
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _locus_alignment_paths(locus_alignment_root: Path) -> list[Path]:
    paths = sorted(locus_alignment_root.glob("*.fasta"))
    if not paths:
        raise FileNotFoundError(f"no locus alignments found in {locus_alignment_root}")
    return paths


def _validate_locus_taxa(alignment_paths: list[Path]) -> None:
    expected: set[str] | None = None
    for path in alignment_paths:
        identifiers = {record.identifier for record in load_fasta_alignment(path)}
        if expected is None:
            expected = identifiers
            continue
        if identifiers != expected:
            raise ValueError(
                f"locus alignment taxa do not match across packaged inputs: {path}"
            )


def _copy_output(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    return Path(shutil.copy2(source, destination))


def _write_workflow_summary_table(
    path: Path,
    report: CatarrhineMitogenomeFiveLocusPanelWorkflowReport,
) -> Path:
    occupancy = report.concatenation_report.occupancy_report
    support = report.partitioned_support.bootstrap_support_summary
    rows = [
        "\t".join(
            [
                "dataset_id",
                "taxon_count",
                "locus_count",
                "alignment_length",
                "partition_count",
                "selected_model",
                "minimum_support",
                "maximum_support",
                "weakly_supported_clade_count",
                "low_coverage_taxon_count",
                "low_coverage_locus_count",
                "warning_count",
            ]
        ),
        "\t".join(
            [
                report.dataset.dataset_id,
                str(report.concatenation_report.taxon_count),
                str(report.concatenation_report.locus_count),
                str(report.concatenation_report.alignment_length),
                str(len(report.partitions)),
                report.model_selection.selected_model or "",
                _format_number(None if support is None else support.minimum_support),
                _format_number(None if support is None else support.maximum_support),
                str(0 if support is None else support.weakly_supported_clade_count),
                str(len(occupancy.low_coverage_taxa)),
                str(len(occupancy.low_coverage_loci)),
                str(
                    len(report.concatenation_report.warnings)
                    + len(report.model_selection.notes)
                    + len(report.partitioned_support.notes)
                ),
            ]
        ),
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _write_occupancy_taxa_table(path: Path, report: LocusOccupancyReport) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "covered_locus_count",
            "total_locus_count",
            "locus_coverage_fraction",
            "observed_site_count",
            "total_site_count",
            "site_coverage_fraction",
            "low_coverage",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "covered_locus_count": str(row.covered_locus_count),
                "total_locus_count": str(row.total_locus_count),
                "locus_coverage_fraction": str(row.locus_coverage_fraction),
                "observed_site_count": str(row.observed_site_count),
                "total_site_count": str(row.total_site_count),
                "site_coverage_fraction": str(row.site_coverage_fraction),
                "low_coverage": str(row.low_coverage).lower(),
            }
            for row in report.taxa
        ],
    )


def _write_occupancy_loci_table(path: Path, report: LocusOccupancyReport) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "locus_name",
            "covered_taxon_count",
            "total_taxa",
            "taxon_coverage_fraction",
            "observed_site_count",
            "total_site_count",
            "site_coverage_fraction",
            "low_coverage",
        ],
        rows=[
            {
                "locus_name": row.locus_name,
                "covered_taxon_count": str(row.covered_taxon_count),
                "total_taxa": str(row.total_taxa),
                "taxon_coverage_fraction": str(row.taxon_coverage_fraction),
                "observed_site_count": str(row.observed_site_count),
                "total_site_count": str(row.total_site_count),
                "site_coverage_fraction": str(row.site_coverage_fraction),
                "low_coverage": str(row.low_coverage).lower(),
            }
            for row in report.loci
        ],
    )


def _write_occupancy_matrix_table(path: Path, report: LocusOccupancyReport) -> Path:
    locus_names = [partition.name for partition in report.partitions]
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            *locus_names,
            "covered_locus_count",
            "total_locus_count",
            "locus_coverage_fraction",
            "observed_site_count",
            "total_site_count",
            "site_coverage_fraction",
            "low_coverage",
        ],
        rows=[
            {
                "taxon": row.taxon,
                **{
                    locus_name: str(row.occupancies[locus_name])
                    for locus_name in locus_names
                },
                "covered_locus_count": str(row.covered_locus_count),
                "total_locus_count": str(row.total_locus_count),
                "locus_coverage_fraction": str(row.locus_coverage_fraction),
                "observed_site_count": str(row.observed_site_count),
                "total_site_count": str(row.total_site_count),
                "site_coverage_fraction": str(row.site_coverage_fraction),
                "low_coverage": str(row.low_coverage).lower(),
            }
            for row in report.taxa
        ],
    )


def _write_overview(
    path: Path,
    dataset: CatarrhineMitogenomeFiveLocusPanelDataset,
    workflow_bundle: CatarrhineMitogenomeFiveLocusPanelWorkflowBundle,
) -> Path:
    lines = [
        "# Catarrhine Mitogenome Five-Locus Demo",
        "",
        f"- dataset id: `{dataset.dataset_id}`",
        f"- taxon count: `{dataset.taxon_count}`",
        f"- locus count: `{dataset.locus_count}`",
        f"- selected model: `{workflow_bundle.selected_model}`",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{workflow_bundle.workflow_summary_path.name}`",
        f"- concatenated supermatrix: `{workflow_bundle.supermatrix_path.name}`",
        f"- partition file: `{workflow_bundle.partitions_path.name}`",
        f"- occupancy matrix: `{workflow_bundle.occupancy_matrix_path.name}`",
        f"- model candidates: `{workflow_bundle.model_candidates_path.name}`",
        f"- supported tree: `{workflow_bundle.support_tree_path.name}`",
        f"- support table: `{workflow_bundle.support_table_path.name}`",
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
        / "phylogenomics"
        / _DATASET_ID
    )
