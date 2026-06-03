from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.io.fasta import write_fasta_alignment
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.marginal_ancestral_probabilities import (
    evaluate_nucleotide_marginal_ancestral_probabilities,
    evaluate_nucleotide_marginal_ancestral_probabilities_from_alignment,
)
from bijux_phylogenetics.phylo.likelihood.marginal_ancestral_sites import (
    summarize_marginal_ancestral_sites,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    MarginalAncestralSequenceExportRecord,
    MarginalAncestralSequenceFastaExportReport,
    MarginalAncestralSequenceUncertaintyRow,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


def reconstruct_nucleotide_marginal_ancestral_sequences(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model_name: str,
    posterior_probability_threshold: float = 0.5,
    low_confidence_state_symbol: str = "N",
    kappa: float | None = None,
    base_frequencies: dict[str, float] | None = None,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | list[float]
        | tuple[float, ...]
        | None
    ) = None,
) -> MarginalAncestralSequenceFastaExportReport:
    """Project marginal nucleotide posteriors into FASTA sequences plus uncertainty."""
    posterior_report = evaluate_nucleotide_marginal_ancestral_probabilities(
        tree,
        records,
        model_name=model_name,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
    )
    return reconstruct_nucleotide_marginal_ancestral_sequences_from_report(
        posterior_report,
        posterior_probability_threshold=posterior_probability_threshold,
        low_confidence_state_symbol=low_confidence_state_symbol,
    )


def reconstruct_nucleotide_marginal_ancestral_sequences_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    posterior_probability_threshold: float = 0.5,
    low_confidence_state_symbol: str = "N",
    kappa: float | None = None,
    base_frequencies: dict[str, float] | None = None,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | list[float]
        | tuple[float, ...]
        | None
    ) = None,
) -> MarginalAncestralSequenceFastaExportReport:
    """Build marginal ancestral FASTA export from file paths."""
    posterior_report = (
        evaluate_nucleotide_marginal_ancestral_probabilities_from_alignment(
            tree_path,
            alignment_path,
            model_name=model_name,
            kappa=kappa,
            base_frequencies=base_frequencies,
            exchangeabilities=exchangeabilities,
        )
    )
    return reconstruct_nucleotide_marginal_ancestral_sequences_from_report(
        posterior_report,
        posterior_probability_threshold=posterior_probability_threshold,
        low_confidence_state_symbol=low_confidence_state_symbol,
    )


def reconstruct_nucleotide_marginal_ancestral_sequences_from_report(
    posterior_report,
    *,
    posterior_probability_threshold: float = 0.5,
    low_confidence_state_symbol: str = "N",
) -> MarginalAncestralSequenceFastaExportReport:
    """Project one existing posterior report into FASTA and uncertainty artifacts."""
    _validate_posterior_probability_threshold(posterior_probability_threshold)
    validated_low_confidence_symbol = _validate_low_confidence_state_symbol(
        low_confidence_state_symbol
    )
    summary_rows = summarize_marginal_ancestral_sites(posterior_report)
    sequence_by_node_id: dict[str, list[str]] = {}
    record_metadata_by_node_id: dict[str, tuple[str | None, list[str]]] = {}
    uncertainty_rows: list[MarginalAncestralSequenceUncertaintyRow] = []
    for row in summary_rows:
        low_confidence = row.max_posterior_probability < posterior_probability_threshold
        exported_state = (
            validated_low_confidence_symbol if low_confidence else row.most_likely_state
        )
        sequence_by_node_id.setdefault(row.node_id, []).append(exported_state)
        record_metadata_by_node_id.setdefault(
            row.node_id,
            (row.node_name, row.descendant_taxa),
        )
        uncertainty_rows.append(
            MarginalAncestralSequenceUncertaintyRow(
                node_id=row.node_id,
                node_name=row.node_name,
                descendant_taxa=row.descendant_taxa,
                pattern_id=row.pattern_id,
                site_position=row.site_position,
                exported_state=exported_state,
                most_likely_state=row.most_likely_state,
                max_posterior_probability=row.max_posterior_probability,
                low_confidence=low_confidence,
                posterior_probability_a=row.posterior_probability_a,
                posterior_probability_c=row.posterior_probability_c,
                posterior_probability_g=row.posterior_probability_g,
                posterior_probability_t=row.posterior_probability_t,
            )
        )
    sequence_records = [
        MarginalAncestralSequenceExportRecord(
            node_id=node_id,
            node_name=record_metadata_by_node_id[node_id][0],
            descendant_taxa=record_metadata_by_node_id[node_id][1],
            sequence="".join(sequence_by_node_id[node_id]),
        )
        for node_id in sorted(sequence_by_node_id)
    ]
    return MarginalAncestralSequenceFastaExportReport(
        model_name=posterior_report.model_name,
        taxa=posterior_report.taxa,
        site_count=posterior_report.site_count,
        pattern_count=posterior_report.pattern_count,
        internal_node_count=posterior_report.internal_node_count,
        compression_used=posterior_report.compression_used,
        expansion_policy="thresholded-fasta-plus-uncertainty-table",
        tree_newick=posterior_report.tree_newick,
        parameter_values=posterior_report.parameter_values,
        posterior_probability_threshold=posterior_probability_threshold,
        low_confidence_state_symbol=validated_low_confidence_symbol,
        sequence_records=sequence_records,
        uncertainty_rows=uncertainty_rows,
    )


def write_marginal_ancestral_sequence_fasta(
    path: Path,
    report: MarginalAncestralSequenceFastaExportReport,
) -> Path:
    """Write thresholded marginal ancestral FASTA sequences for all internal nodes."""
    return write_fasta_alignment(
        path,
        [
            AlignmentRecord(
                identifier=record.node_id,
                sequence=record.sequence,
            )
            for record in report.sequence_records
        ],
    )


def write_marginal_ancestral_sequence_uncertainty_table(
    path: Path,
    report: MarginalAncestralSequenceFastaExportReport,
) -> Path:
    """Write one TSV uncertainty table aligned with the marginal ancestral FASTA export."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "model_name",
                "node_id",
                "node_name",
                "descendant_taxa",
                "pattern_id",
                "site_position",
                "exported_state",
                "most_likely_state",
                "max_posterior_probability",
                "low_confidence",
                "posterior_probability_threshold",
                "low_confidence_state_symbol",
                "posterior_probability_a",
                "posterior_probability_c",
                "posterior_probability_g",
                "posterior_probability_t",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.uncertainty_rows:
            writer.writerow(
                {
                    "model_name": report.model_name,
                    "node_id": row.node_id,
                    "node_name": row.node_name or "",
                    "descendant_taxa": "|".join(row.descendant_taxa),
                    "pattern_id": row.pattern_id,
                    "site_position": row.site_position,
                    "exported_state": row.exported_state,
                    "most_likely_state": row.most_likely_state,
                    "max_posterior_probability": repr(row.max_posterior_probability),
                    "low_confidence": str(row.low_confidence).lower(),
                    "posterior_probability_threshold": repr(
                        report.posterior_probability_threshold
                    ),
                    "low_confidence_state_symbol": report.low_confidence_state_symbol,
                    "posterior_probability_a": repr(row.posterior_probability_a),
                    "posterior_probability_c": repr(row.posterior_probability_c),
                    "posterior_probability_g": repr(row.posterior_probability_g),
                    "posterior_probability_t": repr(row.posterior_probability_t),
                }
            )
    return path


def _validate_posterior_probability_threshold(
    posterior_probability_threshold: float,
) -> None:
    if not 0.0 < posterior_probability_threshold <= 1.0:
        raise ValueError(
            "posterior_probability_threshold must be greater than 0 and at most 1"
        )


def _validate_low_confidence_state_symbol(low_confidence_state_symbol: str) -> str:
    symbol = low_confidence_state_symbol.strip().upper()
    if len(symbol) != 1:
        raise ValueError("low_confidence_state_symbol must be exactly one character")
    return symbol
