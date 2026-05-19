from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from bijux_phylogenetics.core.alignment import (
    AlignmentQualityReport,
    AlignmentSummary,
    AlignmentTrimReport,
)
from bijux_phylogenetics.engines.fasta_to_tree import (
    FastaToTreeWorkflowReport,
    run_fasta_to_tree_workflow,
)
from bijux_phylogenetics.io.fasta import (
    AlignmentRecord,
    load_fasta_records,
    write_fasta_alignment,
)
from bijux_phylogenetics.io.fasta.cleaning import trim_alignment
from bijux_phylogenetics.io.fasta.quality import build_alignment_quality_report
from bijux_phylogenetics.io.fasta.records import (
    summarise_fasta,
    summarise_records_as_alignment_summary,
    validate_fasta_input,
)

_DATASET_ID = "pleistocene_bear_cytb_fragments"
_DATASET_LABEL = "Pleistocene bear CYTB fragment panel"
_SEQUENCE_TYPE = "dna"
_WORKFLOW_PREFIX = "pleistocene-bear-cytb-fragments"
_IQTREE_SEED = 1
_IQTREE_THREADS = 1
_BOOTSTRAP_REPLICATES = 1000
_SITE_MISSINGNESS_THRESHOLD = 0.15
_SEQUENCE_MISSINGNESS_THRESHOLD = 0.15
_DEGRADED_SEQUENCE_IDS = (
    "cave_bear_ud1838_fragment",
    "cave_bear_wk01_fragment",
)
_SOURCE_ACCESSIONS = (
    "OQ318974.1",
    "NC_003428.1",
    "OQ318956.1",
    "KX641337.1",
    "KX641335.1",
)


@dataclass(slots=True)
class PleistoceneBearCytbFragmentDataset:
    """Packaged ancient-DNA-style dataset for degraded sequence workflow review."""

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
    site_missingness_threshold: float
    sequence_missingness_threshold: float
    degraded_sequence_ids: tuple[str, ...]
    source_accessions: tuple[str, ...]
    source_summary: str


@dataclass(slots=True)
class PleistoceneBearCytbFragmentDatasetExportResult:
    """Materialized copy of the packaged ancient-DNA-style dataset."""

    output_root: Path
    readme_path: Path
    sequences_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class PleistoceneBearMissingnessEffectRow:
    """Reviewer-facing missingness comparison for one sequence across workflow stages."""

    identifier: str
    raw_sequence_length: int
    degraded_sequence: bool
    aligned_missing_fraction: float
    engine_trimmed_missing_fraction: float
    cleaned_missing_fraction: float
    removed_by_missingness_cleanup: bool


@dataclass(slots=True)
class PleistoceneBearCytbFragmentWorkflowReport:
    """One governed degraded-sequence workflow run over the packaged bear panel."""

    dataset: PleistoceneBearCytbFragmentDataset
    workflow: FastaToTreeWorkflowReport
    aligned_summary: AlignmentSummary
    trimmed_summary: AlignmentSummary
    cleaned_summary: AlignmentSummary
    aligned_quality: AlignmentQualityReport
    trimmed_quality: AlignmentQualityReport
    cleaned_quality: AlignmentQualityReport
    missingness_cleanup: AlignmentTrimReport
    cleaned_records: list[AlignmentRecord]
    missingness_rows: list[PleistoceneBearMissingnessEffectRow]


@dataclass(slots=True)
class PleistoceneBearCytbFragmentWorkflowBundle:
    """Written degraded-sequence workflow outputs for the packaged bear panel."""

    output_root: Path
    selected_model: str
    minimum_support: float | None
    maximum_support: float | None
    median_support: float | None
    weakly_supported_clade_count: int
    removed_column_count: int
    removed_sequence_count: int
    cleaned_missing_data_fraction: float
    summary_path: Path
    missingness_effects_path: Path
    alignment_path: Path
    trimmed_alignment_path: Path
    cleaned_alignment_path: Path
    tree_path: Path
    model_table_path: Path
    support_table_path: Path
    log_path: Path
    manifest_path: Path
    engine_artifact_root: Path


@dataclass(slots=True)
class PleistoceneBearCytbFragmentDemoResult:
    """Dataset export plus workflow outputs for the public degraded-sequence demo."""

    output_root: Path
    dataset: PleistoceneBearCytbFragmentDataset
    dataset_export: PleistoceneBearCytbFragmentDatasetExportResult
    workflow_bundle: PleistoceneBearCytbFragmentWorkflowBundle
    overview_path: Path


def load_pleistocene_bear_cytb_fragment_dataset() -> PleistoceneBearCytbFragmentDataset:
    """Expose the packaged bear fragment panel as a first-class runtime surface."""
    dataset_root = _resource_root()
    sequences_path = dataset_root / "sequences.fasta"
    validation = validate_fasta_input(sequences_path, sequence_type=_SEQUENCE_TYPE)
    return PleistoceneBearCytbFragmentDataset(
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
        site_missingness_threshold=_SITE_MISSINGNESS_THRESHOLD,
        sequence_missingness_threshold=_SEQUENCE_MISSINGNESS_THRESHOLD,
        degraded_sequence_ids=_DEGRADED_SEQUENCE_IDS,
        source_accessions=_SOURCE_ACCESSIONS,
        source_summary=(
            "Modern bear CYTB references paired with real ancient cave-bear CYTB "
            "sequences reduced to short fragment-style inputs with explicit "
            "internal missing-data blocks."
        ),
    )


def export_pleistocene_bear_cytb_fragment_dataset(
    destination: Path,
) -> PleistoceneBearCytbFragmentDatasetExportResult:
    """Copy the packaged degraded bear dataset and reference outputs to one directory."""
    dataset = load_pleistocene_bear_cytb_fragment_dataset()
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
    return PleistoceneBearCytbFragmentDatasetExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        sequences_path=Path(sequences_path),
        expected_output_root=expected_output_root,
    )


def run_pleistocene_bear_cytb_fragment_workflow(
    out_dir: Path,
    *,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = _IQTREE_SEED,
    iqtree_threads: int = _IQTREE_THREADS,
    bootstrap_replicates: int = _BOOTSTRAP_REPLICATES,
) -> PleistoceneBearCytbFragmentWorkflowReport:
    """Run the owned degraded-sequence workflow over the packaged bear panel."""
    dataset = load_pleistocene_bear_cytb_fragment_dataset()
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
    aligned_summary = summarise_fasta(workflow.output_paths["alignment"])
    trimmed_summary = summarise_fasta(workflow.output_paths["trimmed_alignment"])
    aligned_quality = build_alignment_quality_report(workflow.output_paths["alignment"])
    trimmed_quality = build_alignment_quality_report(
        workflow.output_paths["trimmed_alignment"]
    )
    cleaned_records, missingness_cleanup = trim_alignment(
        workflow.output_paths["alignment"],
        site_missingness_threshold=dataset.site_missingness_threshold,
        sequence_missingness_threshold=dataset.sequence_missingness_threshold,
    )
    cleaned_summary = summarise_records_as_alignment_summary(
        path=workflow.output_paths["alignment"],
        records=cleaned_records,
    )
    with TemporaryDirectory(prefix="pleistocene-bear-cleaned-") as temporary_root:
        cleaned_alignment_path = Path(temporary_root) / "cleaned.aln"
        write_fasta_alignment(cleaned_alignment_path, cleaned_records)
        cleaned_quality = build_alignment_quality_report(cleaned_alignment_path)
    missingness_rows = _build_missingness_rows(
        dataset,
        aligned_summary=aligned_summary,
        trimmed_summary=trimmed_summary,
        cleaned_summary=cleaned_summary,
        removed_sequence_ids={
            row.identifier for row in missingness_cleanup.removed_sequences
        },
    )
    return PleistoceneBearCytbFragmentWorkflowReport(
        dataset=dataset,
        workflow=workflow,
        aligned_summary=aligned_summary,
        trimmed_summary=trimmed_summary,
        cleaned_summary=cleaned_summary,
        aligned_quality=aligned_quality,
        trimmed_quality=trimmed_quality,
        cleaned_quality=cleaned_quality,
        missingness_cleanup=missingness_cleanup,
        cleaned_records=cleaned_records,
        missingness_rows=missingness_rows,
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
        output_root / "workflow-summary.tsv", report
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


def run_pleistocene_bear_cytb_fragment_demo(
    output_root: Path,
    *,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = _IQTREE_SEED,
    iqtree_threads: int = _IQTREE_THREADS,
    bootstrap_replicates: int = _BOOTSTRAP_REPLICATES,
) -> PleistoceneBearCytbFragmentDemoResult:
    """Materialize the packaged degraded bear dataset and rerun the governed workflow."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    dataset = load_pleistocene_bear_cytb_fragment_dataset()
    dataset_export = export_pleistocene_bear_cytb_fragment_dataset(
        output_root / "dataset"
    )
    with TemporaryDirectory(prefix="pleistocene-bear-workflow-") as temporary_root:
        workflow_report = run_pleistocene_bear_cytb_fragment_workflow(
            Path(temporary_root),
            mafft_executable=mafft_executable,
            trimal_executable=trimal_executable,
            iqtree_executable=iqtree_executable,
            iqtree_seed=iqtree_seed,
            iqtree_threads=iqtree_threads,
            bootstrap_replicates=bootstrap_replicates,
        )
        workflow_bundle = write_pleistocene_bear_cytb_fragment_workflow_bundle(
            output_root / "workflow",
            workflow_report,
        )
    overview_path = _write_overview(
        output_root / "overview.md", dataset, workflow_bundle
    )
    return PleistoceneBearCytbFragmentDemoResult(
        output_root=output_root,
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _build_missingness_rows(
    dataset: PleistoceneBearCytbFragmentDataset,
    *,
    aligned_summary: AlignmentSummary,
    trimmed_summary: AlignmentSummary,
    cleaned_summary: AlignmentSummary,
    removed_sequence_ids: set[str],
) -> list[PleistoceneBearMissingnessEffectRow]:
    raw_lengths = {
        record.identifier: len(record.sequence)
        for record in load_fasta_records(dataset.sequences_path)
    }
    aligned_missing = {
        row.identifier: row.missing_fraction
        for row in aligned_summary.per_sequence_missingness
    }
    trimmed_missing = {
        row.identifier: row.missing_fraction
        for row in trimmed_summary.per_sequence_missingness
    }
    cleaned_missing = {
        row.identifier: row.missing_fraction
        for row in cleaned_summary.per_sequence_missingness
    }
    rows: list[PleistoceneBearMissingnessEffectRow] = []
    for identifier in aligned_summary.ids:
        rows.append(
            PleistoceneBearMissingnessEffectRow(
                identifier=identifier,
                raw_sequence_length=raw_lengths[identifier],
                degraded_sequence=identifier in dataset.degraded_sequence_ids,
                aligned_missing_fraction=aligned_missing[identifier],
                engine_trimmed_missing_fraction=trimmed_missing[identifier],
                cleaned_missing_fraction=cleaned_missing.get(identifier, 1.0),
                removed_by_missingness_cleanup=identifier in removed_sequence_ids,
            )
        )
    return rows


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


def _write_overview(
    path: Path,
    dataset: PleistoceneBearCytbFragmentDataset,
    workflow_bundle: PleistoceneBearCytbFragmentWorkflowBundle,
) -> Path:
    lines = [
        "# Pleistocene Bear CYTB Fragment Demo",
        "",
        f"- dataset id: `{dataset.dataset_id}`",
        f"- sequence count: `{dataset.sequence_count}`",
        f"- degraded sequence count: `{len(dataset.degraded_sequence_ids)}`",
        f"- selected model: `{workflow_bundle.selected_model}`",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{workflow_bundle.summary_path.name}`",
        f"- missingness effects: `{workflow_bundle.missingness_effects_path.name}`",
        f"- final supported tree: `{workflow_bundle.tree_path.name}`",
        f"- cleaned alignment: `{workflow_bundle.cleaned_alignment_path.name}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _copy_output(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    return Path(shutil.copy2(source, destination))


def _format_number(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".12g")


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent
        / "resources"
        / "datasets"
        / "ancient_dna"
        / _DATASET_ID
    )
