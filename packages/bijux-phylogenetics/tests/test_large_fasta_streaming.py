from __future__ import annotations

from pathlib import Path
import tracemalloc

import pytest

from bijux_phylogenetics.io.fasta import (
    build_alignment_quality_report,
    validate_fasta_input,
)


def _write_large_fasta(
    path: Path,
    *,
    sequence_count: int,
    sequence_length: int,
    alphabet: str = "ACGT",
) -> Path:
    with path.open("w", encoding="utf-8") as handle:
        for index in range(sequence_count):
            sequence = "".join(
                alphabet[(index + offset) % len(alphabet)]
                for offset in range(sequence_length)
            )
            handle.write(f">taxon_{index:04d}\n{sequence}\n")
    return path


@pytest.mark.slow
def test_validate_fasta_input_handles_thousand_sequence_raw_fasta_with_recorded_peak(
    tmp_path: Path,
) -> None:
    input_path = _write_large_fasta(
        tmp_path / "large-raw.fasta",
        sequence_count=1200,
        sequence_length=1024,
    )

    tracemalloc.start()
    report = validate_fasta_input(input_path, sequence_type="dna")
    _, peak_memory_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    assert report.summary.sequence_count == 1200
    assert report.summary.min_sequence_length == 1024
    assert report.summary.max_sequence_length == 1024
    assert report.summary.total_residue_count == 1200 * 1024
    assert report.sequence_type_report.detected_type == "dna"
    assert report.duplicate_identifiers == []
    assert report.illegal_characters == []
    assert report.empty_sequences == []
    assert report.length_outliers == []
    assert peak_memory_bytes > 0


@pytest.mark.slow
@pytest.mark.stress_small
def test_build_alignment_quality_report_handles_large_alignment_without_pairwise_blowup(
    tmp_path: Path,
) -> None:
    alignment_path = _write_large_fasta(
        tmp_path / "large-alignment.fasta",
        sequence_count=1024,
        sequence_length=1024,
    )

    tracemalloc.start()
    report = build_alignment_quality_report(alignment_path)
    _, peak_memory_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    assert report.sequence_count == 1024
    assert report.alignment_length == 1024
    assert report.near_duplicate_pairs == []
    assert report.near_duplicate_scan_performed is False
    assert any(
        "near-duplicate sequence scan was skipped" in warning
        for warning in report.warnings
    )
    assert peak_memory_bytes > 0
