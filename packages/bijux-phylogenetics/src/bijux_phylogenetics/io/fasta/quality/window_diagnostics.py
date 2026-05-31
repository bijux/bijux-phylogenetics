from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment import (
    AlignmentRecord,
    AlignmentSummary,
    AlignmentSuspiciousRegion,
    AlignmentWindowSummary,
)

from ..core import (
    _GAP_CHARACTERS,
    _is_ambiguity_character,
    _is_explicit_missing,
    load_fasta_alignment,
)
from ..records import (
    summarise_fasta,
)


def summarize_alignment_windows_from_records(
    summary: AlignmentSummary,
    records: list[AlignmentRecord],
    *,
    window_size: int,
    step_size: int,
) -> list[AlignmentWindowSummary]:
    windows: list[AlignmentWindowSummary] = []
    for start_index in range(0, summary.alignment_length, step_size):
        end_index = min(start_index + window_size, summary.alignment_length)
        if end_index <= start_index:
            continue
        columns = list(
            zip(
                *(record.sequence[start_index:end_index] for record in records),
                strict=True,
            )
        )
        if not columns:
            continue

        gap_count = 0
        missing_count = 0
        ambiguity_count = 0
        variable_sites = 0
        disagreement_fractions: list[float] = []
        comparable_sites = 0
        total_residues = len(records) * len(columns)
        alphabet = summary.inferred_alphabet
        for column in columns:
            gap_count += sum(1 for residue in column if residue in _GAP_CHARACTERS)
            missing_count += sum(
                1 for residue in column if _is_explicit_missing(residue)
            )
            ambiguity_count += sum(
                1
                for residue in column
                if _is_ambiguity_character(residue, alphabet=alphabet)
            )
            comparable = [
                residue.upper()
                for residue in column
                if residue not in _GAP_CHARACTERS
                and not _is_explicit_missing(residue)
                and not _is_ambiguity_character(residue, alphabet=alphabet)
            ]
            if comparable:
                comparable_sites += len(comparable)
                state_counts = {
                    state: comparable.count(state) for state in set(comparable)
                }
                if len(state_counts) > 1:
                    variable_sites += 1
                disagreement_fractions.append(
                    1.0 - (max(state_counts.values()) / len(comparable))
                )
            else:
                disagreement_fractions.append(0.0)

        windows.append(
            AlignmentWindowSummary(
                start=start_index + 1,
                end=end_index,
                site_count=len(columns),
                gap_fraction=round(gap_count / total_residues, 15),
                missing_fraction=round(missing_count / total_residues, 15),
                ambiguity_fraction=round(ambiguity_count / total_residues, 15),
                variable_fraction=round(variable_sites / len(columns), 15),
                disagreement_fraction=round(
                    sum(disagreement_fractions) / len(disagreement_fractions),
                    15,
                ),
                comparable_fraction=round(comparable_sites / total_residues, 15),
            )
        )
        if end_index == summary.alignment_length:
            break
    return windows


def detect_over_aligned_regions_from_windows(
    windows: list[AlignmentWindowSummary],
) -> list[AlignmentSuspiciousRegion]:
    regions: list[AlignmentSuspiciousRegion] = []
    for window in windows:
        uncertainty_fraction = (
            window.gap_fraction + window.missing_fraction + window.ambiguity_fraction
        )
        if uncertainty_fraction >= 0.4 and window.variable_fraction <= 0.2:
            regions.append(
                AlignmentSuspiciousRegion(
                    start=window.start,
                    end=window.end,
                    kind="over_aligned",
                    score=round(uncertainty_fraction - window.variable_fraction, 15),
                    note="gap- or ambiguity-heavy window with little residual variability; review for aggressive or artifactual alignment",
                )
            )
    return regions


def detect_under_aligned_regions_from_windows(
    windows: list[AlignmentWindowSummary],
) -> list[AlignmentSuspiciousRegion]:
    regions: list[AlignmentSuspiciousRegion] = []
    for window in windows:
        if window.variable_fraction >= 0.7 and window.disagreement_fraction >= 0.35:
            regions.append(
                AlignmentSuspiciousRegion(
                    start=window.start,
                    end=window.end,
                    kind="under_aligned",
                    score=round(
                        window.variable_fraction + window.disagreement_fraction,
                        15,
                    ),
                    note="high local mismatch and disagreement suggest the region may require realignment or masking",
                )
            )
    return regions


def summarize_alignment_windows(
    path: Path,
    *,
    window_size: int = 30,
    step_size: int = 10,
) -> list[AlignmentWindowSummary]:
    """Summarize an aligned FASTA file in sliding windows."""
    if window_size <= 0:
        raise ValueError(f"window_size must be positive, got {window_size}")
    if step_size <= 0:
        raise ValueError(f"step_size must be positive, got {step_size}")

    summary = summarise_fasta(path)
    records = load_fasta_alignment(path)
    return summarize_alignment_windows_from_records(
        summary,
        records,
        window_size=window_size,
        step_size=step_size,
    )


def detect_over_aligned_regions(
    path: Path,
    *,
    window_size: int = 30,
    step_size: int = 10,
) -> list[AlignmentSuspiciousRegion]:
    """Flag suspicious windows that look excessively gap-heavy or over-regularized."""
    return detect_over_aligned_regions_from_windows(
        summarize_alignment_windows(path, window_size=window_size, step_size=step_size)
    )


def detect_under_aligned_regions(
    path: Path,
    *,
    window_size: int = 30,
    step_size: int = 10,
) -> list[AlignmentSuspiciousRegion]:
    """Flag suspicious windows with strong local mismatch or gap disorder."""
    return detect_under_aligned_regions_from_windows(
        summarize_alignment_windows(path, window_size=window_size, step_size=step_size)
    )
