from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.datasets.study_inputs import (
    check_tree_and_trait_taxon_names,
    write_tree_trait_name_mismatch_table,
)
from bijux_phylogenetics.runtime.errors import MetadataJoinError

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


def test_check_tree_and_trait_taxon_names_reports_both_mismatch_sides() -> None:
    report = check_tree_and_trait_taxon_names(
        fixture("example_tree.nwk"),
        fixture("example_traits.tsv"),
    )

    assert report.compatible is False
    assert report.reference_outcome == "mismatch"
    assert report.matching_policy == "case-sensitive-exact-label-matching"
    assert report.tree_not_data == ["D"]
    assert report.data_not_tree == ["E"]
    assert [(row.mismatch_side, row.taxon) for row in report.mismatch_rows] == [
        ("tree_not_data", "D"),
        ("data_not_tree", "E"),
    ]


def test_check_tree_and_trait_taxon_names_reports_clean_ok_surface() -> None:
    report = check_tree_and_trait_taxon_names(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative_reordered.tsv"),
    )

    assert report.compatible is True
    assert report.reference_outcome == "OK"
    assert report.tree_not_data == []
    assert report.data_not_tree == []
    assert report.mismatch_rows == []


def test_check_tree_and_trait_taxon_names_is_case_sensitive_by_policy(
    tmp_path: Path,
) -> None:
    table_path = tmp_path / "case-sensitive-traits.tsv"
    table_path.write_text(
        "taxon\tvalue\na\t1.0\nB\t2.0\nC\t3.0\nD\t4.0\n",
        encoding="utf-8",
    )

    report = check_tree_and_trait_taxon_names(
        fixture("example_tree.nwk"),
        table_path,
    )

    assert report.tree_not_data == ["A"]
    assert report.data_not_tree == ["a"]


def test_check_tree_and_trait_taxon_names_rejects_duplicate_taxa() -> None:
    with pytest.raises(MetadataJoinError, match="duplicate taxon 'A'"):
        check_tree_and_trait_taxon_names(
            fixture("example_tree.nwk"),
            fixture("example_metadata_duplicate.tsv"),
        )


def test_write_tree_trait_name_mismatch_table_writes_machine_readable_rows(
    tmp_path: Path,
) -> None:
    report = check_tree_and_trait_taxon_names(
        fixture("example_tree.nwk"),
        fixture("example_traits.tsv"),
    )

    output_path = tmp_path / "trait-name-mismatches.tsv"
    write_tree_trait_name_mismatch_table(output_path, report)

    assert output_path.read_text(encoding="utf-8") == (
        "mismatch_side\ttaxon\tpresent_in_tree\tpresent_in_traits\n"
        "tree_not_data\tD\ttrue\tfalse\n"
        "data_not_tree\tE\tfalse\ttrue\n"
    )
