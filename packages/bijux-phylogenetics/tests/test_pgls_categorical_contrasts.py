from __future__ import annotations

import json
import math
from pathlib import Path

from bijux_phylogenetics.comparative.pgls import build_pgls_model_matrix
from bijux_phylogenetics.comparative.pgls.categorical_contrasts import (
    summarize_pgls_categorical_contrasts,
    write_pgls_categorical_contrast_table,
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


def test_summarize_pgls_categorical_contrasts_reports_baseline_and_group_rows() -> None:
    report = summarize_pgls_categorical_contrasts(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one", "habitat"],
        lambda_value=0.0,
    )
    assert report.categorical_predictor_count == 1
    assert [row.level for row in report.rows] == ["forest", "tundra"]
    baseline_row = report.rows[0]
    assert baseline_row.is_reference_level is True
    assert baseline_row.baseline_level == "forest"
    assert baseline_row.coefficient_name is None
    contrast_row = report.rows[1]
    assert contrast_row.coefficient_name == "habitat[tundra]"
    assert math.isclose(contrast_row.coefficient_estimate or 0.0, -2.0, abs_tol=1e-12)
    assert contrast_row.observed_taxon_count == 2


def test_summarize_pgls_categorical_contrasts_reports_missing_category_taxa() -> None:
    report = summarize_pgls_categorical_contrasts(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative_missing_habitat.tsv"),
        response="response",
        predictors=["habitat"],
        lambda_value=0.0,
    )
    assert [row.missing_category_taxa for row in report.rows] == [["D"], ["D"]]


def test_pgls_categorical_model_matrix_matches_reference_fixture() -> None:
    reference = json.loads(
        fixture("pgls_categorical_model_matrix_reference.json").read_text(
            encoding="utf-8"
        )
    )
    for observation in reference["observations"]:
        report = build_pgls_model_matrix(
            fixture("example_tree.nwk"),
            fixture("example_traits_comparative.tsv"),
            formula=observation["formula"],
        )
        assert report.formula.include_intercept is observation["include_intercept"]
        assert report.encoded_columns == observation["encoded_columns"]
        rows_by_taxon = {row.taxon: row for row in report.rows}
        for expected_row in observation["rows"]:
            actual = rows_by_taxon[expected_row["taxon"]]
            for column in observation["encoded_columns"]:
                assert actual.encoded_values[column] == expected_row[column]


def test_write_pgls_categorical_contrast_table_writes_rows(tmp_path: Path) -> None:
    report = summarize_pgls_categorical_contrasts(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one", "habitat"],
        lambda_value=0.0,
    )
    out_path = tmp_path / "categorical-contrasts.tsv"
    write_pgls_categorical_contrast_table(out_path, report)
    contents = out_path.read_text(encoding="utf-8")
    assert "predictor\tsource_column\tencoding_scheme" in contents
    assert "habitat\thabitat\treference-level\tforest\tforest\ttrue" in contents
