from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.comparative.pgls import (
    build_pgls_model_matrix,
    inspect_pgls_inputs,
    write_pgls_model_matrix_table,
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


def test_pgls_formula_supports_interceptless_numeric_predictor() -> None:
    report = inspect_pgls_inputs(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        formula="response ~ predictor_one - 1",
    )
    assert report.ready is True
    assert report.formula.include_intercept is False
    assert report.formula_audit.includes_intercept is False
    assert report.encoded_columns == ["predictor_one"]
    assert report.model_matrix.encoded_columns == ["predictor_one"]


def test_pgls_formula_supports_interceptless_categorical_encoding() -> None:
    report = inspect_pgls_inputs(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        formula="response ~ 0 + habitat",
    )
    assert report.ready is True
    assert report.formula.include_intercept is False
    assert report.encoded_columns == ["habitat[forest]", "habitat[tundra]"]
    habitat = next(row for row in report.predictors if row.name == "habitat")
    assert habitat.reference_level is None
    assert habitat.encoded_columns == ["habitat[forest]", "habitat[tundra]"]


def test_build_pgls_model_matrix_preserves_interaction_columns() -> None:
    report = build_pgls_model_matrix(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_interaction.tsv"),
        formula="response ~ 0 + predictor_one * habitat",
    )
    assert report.formula.include_intercept is False
    assert report.encoded_columns == [
        "predictor_one",
        "habitat[forest]",
        "habitat[tundra]",
        "predictor_one:habitat[forest]",
        "predictor_one:habitat[tundra]",
    ]
    first_row = report.rows[0]
    assert first_row.taxon == "A"
    assert first_row.encoded_values["predictor_one:habitat[forest]"] == 1.0
    assert first_row.encoded_values["predictor_one:habitat[tundra]"] == 0.0


def test_write_pgls_model_matrix_table_writes_encoded_rows(tmp_path: Path) -> None:
    report = build_pgls_model_matrix(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        formula="response ~ 0 + habitat",
    )
    out_path = tmp_path / "model-matrix.tsv"
    write_pgls_model_matrix_table(out_path, report)
    with out_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["taxon"] for row in rows] == ["A", "B", "C", "D"]
    assert rows[0]["response_value"] == "1.5"
    assert rows[0]["habitat[forest]"] == "1"
    assert rows[0]["habitat[tundra]"] == "0"
