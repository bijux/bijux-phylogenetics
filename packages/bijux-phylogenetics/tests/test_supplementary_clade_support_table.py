from __future__ import annotations

import csv
from pathlib import Path

import pytest

from bijux_phylogenetics.reports import write_supplementary_clade_support_table

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


def _rows_by_taxa(rows):
    return {tuple(row.descendant_taxa): row for row in rows}


def test_write_supplementary_clade_support_table_writes_tree_label_support_rows(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "supplementary-clade-support.tsv"

    result = write_supplementary_clade_support_table(
        output_path,
        tree_path=fixture("example_tree_support_left.nwk"),
    )

    assert result.output_path == output_path
    assert result.row_count == 3
    assert result.supported_clade_count == 3
    assert result.frequency_scored_clade_count == 0

    rows = _rows_by_taxa(result.rows)
    assert rows[("A", "B")].support == 95.0
    assert rows[("A", "B")].support_fraction == 0.95
    assert rows[("A", "B")].support_class == "strong"
    assert rows[("A", "B")].support_method == "tree-label"
    assert rows[("A", "B")].clade_frequency is None

    written_lookup = {
        tuple(row["descendant_taxa"].split("|")): row for row in read_tsv(output_path)
    }
    assert written_lookup[("A", "B")]["support_method"] == "tree-label"
    assert written_lookup[("A", "B")]["clade_frequency"] == ""


def test_write_supplementary_clade_support_table_joins_tree_set_frequency_rows(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "supplementary-clade-support.tsv"

    result = write_supplementary_clade_support_table(
        output_path,
        tree_path=fixture("example_tree_support_left.nwk"),
        comparison_tree_set_path=fixture("example_tree_set_left.nwk"),
    )

    assert result.row_count == 3
    assert result.frequency_scored_clade_count == 3
    assert result.frequency_partial_support_count == 2
    assert result.frequency_absent_clade_count == 0
    assert result.frequency_unscored_clade_count == 0

    rows = _rows_by_taxa(result.rows)
    assert rows[("A", "B")].support == 95.0
    assert rows[("A", "B")].supporting_tree_count == 2
    assert rows[("A", "B")].clade_frequency == pytest.approx(2.0 / 3.0)
    assert rows[("A", "B")].support_percent == pytest.approx(66.66666666666667)
    assert rows[("A", "B")].frequency_status == "partial-support"
    assert rows[("A", "B")].frequency_method == "reference-tree-clade-frequency"
    assert rows[("A", "B", "C", "D")].supporting_tree_count == 3
    assert rows[("A", "B", "C", "D")].frequency_status == "fixed"

    written_rows = read_tsv(output_path)
    written_lookup = {
        tuple(row["descendant_taxa"].split("|")): row for row in written_rows
    }
    assert written_lookup[("A", "B")]["supporting_tree_count"] == "2"
    assert written_lookup[("A", "B")]["frequency_status"] == "partial-support"
