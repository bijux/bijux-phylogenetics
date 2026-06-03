from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta.coding import classify_sequence_coding_behavior
from bijux_phylogenetics.phylo.alignment import (
    AlignmentSummary,
    CodingSequenceExclusion,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError

from ..models import AlignmentTrimmingSummary


def _write_coding_exclusion_table(
    path: Path, exclusions: list[CodingSequenceExclusion]
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "identifier",
                "comparable_length",
                "reason",
                "invalid_codon_count",
                "premature_stop_count",
                "terminal_stop_count",
                "trailing_bases",
                "note",
            ]
        )
    ]
    lines.extend(
        "\t".join(
            [
                row.identifier,
                str(row.comparable_length),
                row.reason,
                str(row.invalid_codon_count),
                str(row.premature_stop_count),
                str(row.terminal_stop_count),
                str(row.trailing_bases),
                row.note,
            ]
        )
        for row in exclusions
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_coding_summary_table(
    path: Path,
    *,
    input_path: Path,
    genetic_code: int,
    exclusions: list[CodingSequenceExclusion],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    exclusion_by_identifier = {row.identifier: row for row in exclusions}
    behaviors = classify_sequence_coding_behavior(
        input_path,
        genetic_code=genetic_code,
    )
    header = "\t".join(
        [
            "identifier",
            "status",
            "comparable_length",
            "divisible_by_three",
            "invalid_codon_count",
            "premature_stop_count",
            "terminal_stop_count",
            "exclusion_reason",
            "note",
        ]
    )
    rows = [header]
    for behavior in behaviors:
        exclusion = exclusion_by_identifier.get(behavior.identifier)
        rows.append(
            "\t".join(
                [
                    behavior.identifier,
                    "excluded" if exclusion is not None else "accepted",
                    str(behavior.comparable_length),
                    "yes" if behavior.divisible_by_three else "no",
                    str(behavior.invalid_codon_count),
                    str(behavior.premature_stop_count),
                    str(behavior.terminal_stop_count),
                    "" if exclusion is None else exclusion.reason,
                    behavior.note,
                ]
            )
        )
    ordered_rows = [rows[0], *sorted(rows[1:])]
    path.write_text("\n".join(ordered_rows) + "\n", encoding="utf-8")
    return path


def _write_alignment_trimming_summary_table(
    path: Path,
    *,
    summary: AlignmentTrimmingSummary,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        "metric\tvalue",
        f"mode\t{summary.mode}",
        (
            "gap_threshold\t"
            if summary.gap_threshold is None
            else f"gap_threshold\t{summary.gap_threshold:.6f}"
        ),
        f"input_alignment_length\t{summary.input_alignment_length}",
        f"trimmed_alignment_length\t{summary.trimmed_alignment_length}",
        f"retained_site_count\t{summary.retained_site_count}",
        f"removed_site_count\t{summary.removed_site_count}",
        f"retained_site_fraction\t{summary.retained_site_fraction:.12g}",
        f"removed_site_fraction\t{summary.removed_site_fraction:.12g}",
        f"input_gap_fraction\t{summary.input_gap_fraction:.12g}",
        f"trimmed_gap_fraction\t{summary.trimmed_gap_fraction:.12g}",
        f"input_gap_percentage\t{summary.input_gap_percentage:.12g}",
        f"trimmed_gap_percentage\t{summary.trimmed_gap_percentage:.12g}",
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _build_alignment_trimming_summary(
    *,
    mode: str,
    gap_threshold: float,
    input_summary: AlignmentSummary,
    trimmed_summary: AlignmentSummary,
) -> AlignmentTrimmingSummary:
    if trimmed_summary.alignment_length > input_summary.alignment_length:
        raise EngineWorkflowError(
            "trimmed alignment is longer than the input alignment, which is not a valid trimAl result"
        )
    retained_site_count = trimmed_summary.alignment_length
    removed_site_count = (
        input_summary.alignment_length - trimmed_summary.alignment_length
    )
    retained_site_fraction = retained_site_count / input_summary.alignment_length
    removed_site_fraction = removed_site_count / input_summary.alignment_length
    return AlignmentTrimmingSummary(
        mode=mode,
        gap_threshold=gap_threshold if mode == "gap-threshold" else None,
        input_alignment_length=input_summary.alignment_length,
        trimmed_alignment_length=trimmed_summary.alignment_length,
        retained_site_count=retained_site_count,
        removed_site_count=removed_site_count,
        retained_site_fraction=retained_site_fraction,
        removed_site_fraction=removed_site_fraction,
        input_gap_fraction=input_summary.gap_fraction,
        trimmed_gap_fraction=trimmed_summary.gap_fraction,
        input_gap_percentage=input_summary.gap_fraction * 100.0,
        trimmed_gap_percentage=trimmed_summary.gap_fraction * 100.0,
    )
