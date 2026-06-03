from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.reports import write_supplementary_diversification_table

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


def test_write_supplementary_diversification_table_writes_clade_and_model_review_rows(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "supplementary-diversification.tsv"

    result = write_supplementary_diversification_table(
        output_path,
        tree_path=fixture("example_tree.nwk"),
        metadata_path=fixture("example_sampling_fractions.tsv"),
        clade_model="birth-death",
    )

    assert result.output_path == output_path
    assert result.row_count == 3
    assert result.clade_model == "birth-death"
    assert result.better_model in {"yule", "birth-death"}
    assert result.high_clade_count == 1
    assert result.low_clade_count == 1
    assert result.sampling_metadata_complete is True
    assert all(row.tree_source.endswith("example_tree.nwk") for row in result.rows)
    assert all(row.metadata_source is not None for row in result.rows)
    assert all(row.sampling_fraction == 0.75 for row in result.rows)
    assert all(row.sampling_metadata_complete is True for row in result.rows)
    assert {row.clade_classification for row in result.rows} == {
        "baseline",
        "high",
        "low",
    }
    assert any(row.clade_classification == "high" for row in result.rows)
    assert any(row.clade_classification == "low" for row in result.rows)
    assert all(row.yule_corrected_tip_count == 5.33333333333333 for row in result.rows)
    assert all(
        row.birth_death_corrected_tip_count == 5.33333333333333 for row in result.rows
    )
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 3
    assert rows[0]["sampling_metadata_complete"] == "true"
    assert rows[0]["sampling_fraction"] == "0.75"


def test_write_supplementary_diversification_table_carries_sampling_problems_into_rows(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "supplementary-diversification-incomplete.tsv"

    result = write_supplementary_diversification_table(
        output_path,
        tree_path=fixture("example_tree.nwk"),
        metadata_path=fixture("example_sampling_fractions_incomplete.tsv"),
        clade_model="yule",
    )

    assert result.row_count == 3
    assert result.clade_model == "yule"
    assert result.sampling_metadata_complete is False
    assert result.warning_count >= 2
    assert all(row.sampling_metadata_complete is False for row in result.rows)
    assert all(row.sampling_fraction == 0.75 for row in result.rows)
    assert any("D" in row.sampling_missing_taxa for row in result.rows)
    assert any(
        "B:missing-sampling-fraction:<missing>" in row.sampling_invalid_rows
        for row in result.rows
    )
    assert any(
        "C:out-of-range-sampling-fraction:1.2" in row.sampling_invalid_rows
        for row in result.rows
    )
    assert any(
        "sampling metadata does not cover every tree tip" in warning
        for row in result.rows
        for warning in row.warnings
    )
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 3
    assert "B:missing-sampling-fraction:<missing>" in rows[0]["sampling_invalid_rows"]
