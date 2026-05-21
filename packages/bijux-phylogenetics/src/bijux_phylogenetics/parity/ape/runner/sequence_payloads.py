from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.distance import (
    build_tree_from_imported_distance_matrix,
    compute_pairwise_genetic_distance_matrix,
)
from bijux_phylogenetics.io.fasta.coding import translate_coding_alignment
from bijux_phylogenetics.io.fasta.matrix import (
    compute_alignment_base_frequency_report,
    compute_alignment_segregating_site_report,
    load_dna_bin_alignment,
)

from .tree_payloads import _tree_structure_payload


def _build_bijux_dnabin_rows(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    matrix = load_dna_bin_alignment(input_fixture)
    return {
        "sequence_count": matrix.sequence_count,
        "alignment_length": matrix.alignment_length,
        "state_count": len(matrix.rows),
    }, [
        {
            "identifier": row.identifier,
            "position": row.position,
            "state": row.state,
        }
        for row in matrix.rows
    ]


def _ape_base_frequency_rows(input_fixture: Path) -> list[dict[str, object]]:
    report = compute_alignment_base_frequency_report(input_fixture)
    return [
        {"state": row.state, "frequency": row.frequency}
        for row in report.alignment_rows
    ]


def _build_bijux_base_frequency_summary(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_alignment_base_frequency_report(input_fixture)
    rows = [
        {"state": row.state, "frequency": row.frequency}
        for row in report.alignment_rows
    ]
    return {
        "sequence_count": report.sequence_count,
        "alignment_length": report.alignment_length,
        "state_count": len(rows),
    }, rows


def _build_bijux_segregating_site_rows(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_alignment_segregating_site_report(input_fixture)
    return {
        "sequence_count": report.sequence_count,
        "alignment_length": report.alignment_length,
        "segregating_site_count": len(report.segregating_site_positions),
    }, [{"position": row.position} for row in report.rows]


def _build_bijux_distance_rows(
    input_fixture: Path,
    *,
    pairwise_deletion: bool,
    distance_model: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_pairwise_genetic_distance_matrix(
        input_fixture,
        model=distance_model,
        gap_handling="pairwise-deletion" if pairwise_deletion else "complete-deletion",
        ambiguity_policy="ignore",
    )
    pair_lookup = {
        (row.left_identifier, row.right_identifier): row for row in report.pairs
    }
    rows: list[dict[str, object]] = []
    finite_distance_count = 0
    undefined_distance_count = 0
    infinite_distance_count = 0
    for left_identifier in report.identifiers:
        for right_identifier in report.identifiers:
            pair = pair_lookup.get((left_identifier, right_identifier))
            if pair is None:
                pair = pair_lookup.get((right_identifier, left_identifier))
            if pair is None:
                raise ValueError(
                    "distance parity rows require a complete symmetric pair lookup"
                )
            if pair.distance is None:
                if (
                    pair.saturation_reason
                    and "tends to infinity" in pair.saturation_reason
                ):
                    distance = ""
                    distance_status = "infinite"
                    infinite_distance_count += 1
                else:
                    distance = ""
                    distance_status = "undefined"
                    undefined_distance_count += 1
            else:
                distance = pair.distance
                distance_status = "finite"
                finite_distance_count += 1
            rows.append(
                {
                    "left_identifier": left_identifier,
                    "right_identifier": right_identifier,
                    "distance": distance,
                    "distance_status": distance_status,
                }
            )
    return {
        "sequence_count": len(report.identifiers),
        "alignment_length": report.alignment_length,
        "pairwise_deletion": pairwise_deletion,
        "distance_model": report.model,
        "finite_distance_count": finite_distance_count,
        "undefined_distance_count": undefined_distance_count,
        "infinite_distance_count": infinite_distance_count,
    }, rows


def _build_bijux_neighbor_joining_structure(
    input_fixture: Path,
) -> tuple[dict[str, object], None, str]:
    tree, _report = build_tree_from_imported_distance_matrix(
        input_fixture,
        method="neighbor-joining",
    )
    summary, _rows, normalized_text = _tree_structure_payload(tree, False, [])
    return summary, None, normalized_text


def _build_bijux_translation_rows(
    input_fixture: Path,
    *,
    genetic_code_id: int,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    translated, report = translate_coding_alignment(
        input_fixture, genetic_code=genetic_code_id
    )
    return {
        "sequence_count": report.translated_sequence_count,
        "translated_length": report.translated_alignment_length,
        "stop_codon_count": report.stop_codon_count,
        "dropped_trailing_nucleotide_count": report.dropped_trailing_nucleotide_count,
        "warning_count": len(report.warnings),
        "warnings": report.warnings,
    }, [
        {
            "identifier": row.identifier,
            "amino_acid_sequence": row.sequence,
        }
        for row in translated
    ]
