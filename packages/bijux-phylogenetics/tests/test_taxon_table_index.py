from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import inspect_taxon_table_index

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


def test_inspect_taxon_table_index_reports_duplicate_taxa() -> None:
    audit = inspect_taxon_table_index(fixture("example_metadata_duplicate.tsv"))
    assert audit.row_count == 2
    assert audit.taxon_column == "taxon"
    assert audit.taxa == ["A"]
    assert audit.duplicate_taxa == ["A"]
    assert audit.empty_taxon_rows == []


def test_inspect_taxon_table_index_reports_unique_trait_rows() -> None:
    audit = inspect_taxon_table_index(fixture("example_traits_comparative.tsv"))
    assert audit.row_count == 4
    assert audit.taxon_column == "taxon"
    assert audit.taxa == ["A", "B", "C", "D"]
    assert audit.duplicate_taxa == []
    assert audit.empty_taxon_rows == []
