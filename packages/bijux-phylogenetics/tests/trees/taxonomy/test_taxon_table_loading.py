from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table

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


def test_load_taxon_table_normalizes_trailing_missing_cells_to_empty_strings() -> None:
    table = load_taxon_table(
        fixture("example_traits_phylogenetic_residuals_missing.tsv")
    )

    brain_mass_by_taxon = {
        row[table.taxon_column]: row["brain_mass"] for row in table.rows
    }

    assert brain_mass_by_taxon["E"] == ""
    assert brain_mass_by_taxon["F"] == "12.0"
