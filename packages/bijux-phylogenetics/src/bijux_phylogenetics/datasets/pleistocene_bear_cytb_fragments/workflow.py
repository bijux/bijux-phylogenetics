from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

from bijux_phylogenetics.engines.inference import run_fasta_to_tree_workflow
from bijux_phylogenetics.io.fasta import load_fasta_records, write_fasta_alignment
from bijux_phylogenetics.io.fasta.cleaning import trim_alignment
from bijux_phylogenetics.io.fasta.quality import build_alignment_quality_report
from bijux_phylogenetics.io.fasta.records import (
    summarise_fasta,
    summarise_records_as_alignment_summary,
)
from bijux_phylogenetics.phylo.alignment import AlignmentSummary

from .models import (
    PleistoceneBearCytbFragmentDataset,
    PleistoceneBearCytbFragmentWorkflowReport,
    PleistoceneBearMissingnessEffectRow,
)
from .panel import (
    BOOTSTRAP_REPLICATES,
    IQTREE_SEED,
    IQTREE_THREADS,
    load_pleistocene_bear_cytb_fragment_dataset,
)


def run_pleistocene_bear_cytb_fragment_workflow(
    out_dir: Path,
    *,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = IQTREE_SEED,
    iqtree_threads: int = IQTREE_THREADS,
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
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
