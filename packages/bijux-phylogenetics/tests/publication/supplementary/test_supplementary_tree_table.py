from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.reports import write_supplementary_tree_diagnostics_table

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


def test_write_supplementary_tree_diagnostics_table_writes_tree_level_summary(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "supplementary-tree.tsv"

    result = write_supplementary_tree_diagnostics_table(
        output_path,
        tree_path=fixture("example_tree_support_left.nwk"),
    )

    assert result.output_path == output_path
    assert result.row_count == 1
    row = result.rows[0]
    assert row.tip_count == 4
    assert row.topology_shape == "binary"
    assert row.rooted is True
    assert row.ultrametric is True
    assert row.supported_branch_count == 3
    assert row.strong_support_branch_count == 2
    assert row.moderate_support_branch_count == 1
    assert row.missing_support_branch_count == 0
    assert row.warning_count == 0

    written = read_tsv(output_path)[0]
    assert written["topology_shape"] == "binary"
    assert written["supported_branch_count"] == "3"
    assert written["ultrametric"] == "True"
    assert written["warnings"] == ""


def test_write_supplementary_tree_diagnostics_table_preserves_warning_ledger(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "supplementary-tree.tsv"

    result = write_supplementary_tree_diagnostics_table(
        output_path,
        tree_path=fixture("example_tree_negative_length.nwk"),
    )

    row = result.rows[0]
    assert row.negative_branch_count == 1
    assert row.warning_count >= 1
    assert row.safe_for_publication is False
    written = read_tsv(output_path)[0]
    assert written["negative_branch_count"] == "1"
    assert "negative" in written["warnings"]
