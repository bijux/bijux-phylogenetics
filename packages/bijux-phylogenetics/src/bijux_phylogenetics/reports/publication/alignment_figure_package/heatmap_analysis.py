from __future__ import annotations

from bijux_phylogenetics.phylo.alignment import (
    AlignmentRecord,
    AlignmentSummary,
    SequenceQualityRankingRow,
)

from .contracts import AlignmentHeatmapCell


def site_bins(alignment_length: int, *, maximum_bins: int) -> list[tuple[int, int]]:
    """Partition alignment sites into reviewer-facing heatmap bins."""
    if alignment_length <= 0:
        return []
    if maximum_bins <= 0:
        raise ValueError(f"maximum_bins must be positive, got {maximum_bins}")
    if alignment_length <= maximum_bins:
        return [(position, position) for position in range(1, alignment_length + 1)]
    bin_width = -(-alignment_length // maximum_bins)
    bins: list[tuple[int, int]] = []
    start = 1
    while start <= alignment_length:
        end = min(start + bin_width - 1, alignment_length)
        bins.append((start, end))
        start = end + 1
    return bins


def mean(values: list[float]) -> float:
    """Return a stable reviewer-facing mean for alignment burden metrics."""
    if not values:
        return 0.0
    return round(sum(values) / len(values), 15)


def classify_residue(
    residue: str,
    *,
    alphabet: str,
) -> tuple[float, float, float]:
    """Classify one residue into gap, missing, and ambiguity burden fractions."""
    if residue == "-":
        return 1.0, 0.0, 0.0
    if residue == "?":
        return 0.0, 1.0, 0.0
    canonical = (
        {"A", "C", "G", "T"}
        if alphabet == "dna"
        else {"A", "C", "G", "U"}
        if alphabet == "rna"
        else set("ACDEFGHIKLMNPQRSTVWY")
    )
    return (0.0, 0.0, 0.0) if residue.upper() in canonical else (0.0, 0.0, 1.0)


def build_heatmap_cells(
    summary: AlignmentSummary,
    records: list[AlignmentRecord],
    ranking_rows: list[SequenceQualityRankingRow],
    *,
    maximum_bins: int,
) -> tuple[list[AlignmentHeatmapCell], int, int]:
    """Build one heatmap cell per ranked sequence and site bin."""
    bins = site_bins(summary.alignment_length, maximum_bins=maximum_bins)
    records_by_id = {record.identifier: record for record in records}
    cells: list[AlignmentHeatmapCell] = []
    for row in ranking_rows:
        record = records_by_id[row.identifier]
        for start, end in bins:
            gap_values: list[float] = []
            missing_values: list[float] = []
            ambiguity_values: list[float] = []
            for position in range(start - 1, end):
                gap_fraction, missing_fraction, ambiguity_fraction = classify_residue(
                    record.sequence[position],
                    alphabet=summary.inferred_alphabet,
                )
                gap_values.append(gap_fraction)
                missing_values.append(missing_fraction)
                ambiguity_values.append(ambiguity_fraction)
            cells.append(
                AlignmentHeatmapCell(
                    identifier=row.identifier,
                    bin_start=start,
                    bin_end=end,
                    uncertainty_fraction=mean(
                        [
                            gap_fraction + missing_fraction + ambiguity_fraction
                            for gap_fraction, missing_fraction, ambiguity_fraction in zip(
                                gap_values,
                                missing_values,
                                ambiguity_values,
                                strict=True,
                            )
                        ]
                    ),
                    gap_fraction=mean(gap_values),
                    missing_fraction=mean(missing_values),
                    ambiguity_fraction=mean(ambiguity_values),
                )
            )
        if not bins:
            cells.append(
                AlignmentHeatmapCell(
                    identifier=row.identifier,
                    bin_start=1,
                    bin_end=0,
                    uncertainty_fraction=0.0,
                    gap_fraction=0.0,
                    missing_fraction=0.0,
                    ambiguity_fraction=0.0,
                )
            )
    return cells, len(ranking_rows), len(bins)
