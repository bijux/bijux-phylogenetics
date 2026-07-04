from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.reports import (
    write_supplementary_alignment_diagnostics_table,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_write_supplementary_alignment_diagnostics_table_tracks_filtering_status(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "supplementary-alignment.tsv"

    result = write_supplementary_alignment_diagnostics_table(
        output_path,
        alignment_path=fixture("example_taxon_workflow_alignment.fasta"),
        filtered_alignment_path=fixture(
            "example_taxon_workflow_filtered_alignment.fasta"
        ),
    )

    assert result.output_path == output_path
    assert result.row_count == 3
    assert result.retained_sequence_count == 2
    assert result.removed_sequence_count == 1
    assert result.filtered_only_sequence_count == 0

    rows = {row["sequence_id"]: row for row in read_tsv(output_path)}
    assert rows["A"]["filtering_status"] == "retained_after_filtering"
    assert rows["A"]["original_sequence_present"] == "True"
    assert rows["A"]["filtered_sequence_present"] == "True"
    assert rows["A"]["original_variable_site_count"] == "2"
    assert rows["A"]["original_parsimony_informative_site_count"] == "0"
    assert rows["A"]["filtered_alignment_length"] == "8"
    assert rows["B"]["filtering_status"] == "removed_during_filtering"
    assert rows["B"]["filtering_reason"] == "absent_from_filtered_alignment"
    assert rows["B"]["filtered_sequence_present"] == "False"
    assert rows["C"]["filtering_status"] == "retained_after_filtering"


def test_write_supplementary_alignment_diagnostics_table_marks_missing_filter_request(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "supplementary-alignment.tsv"

    result = write_supplementary_alignment_diagnostics_table(
        output_path,
        alignment_path=fixture("example_alignment.fasta"),
    )

    assert result.row_count == 4
    assert result.retained_sequence_count == 0
    assert result.removed_sequence_count == 0
    assert result.filtered_only_sequence_count == 0

    rows = {row["sequence_id"]: row for row in read_tsv(output_path)}
    assert rows["A"]["filtering_status"] == "not_requested"
    assert rows["A"]["original_low_information"] == "False"
    assert rows["A"]["filtered_alignment_length"] == ""
