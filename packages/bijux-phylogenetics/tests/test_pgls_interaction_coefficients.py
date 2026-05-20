from __future__ import annotations

import json
import math
from pathlib import Path

from bijux_phylogenetics.comparative.pgls import build_pgls_model_matrix, run_pgls
from bijux_phylogenetics.comparative.pgls.interaction_coefficients import (
    summarize_pgls_interaction_coefficients,
    write_pgls_interaction_coefficient_table,
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


def test_summarize_pgls_interaction_coefficients_reports_continuous_by_categorical() -> (
    None
):
    report = summarize_pgls_interaction_coefficients(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_interaction.tsv"),
        formula="response ~ predictor_one * habitat",
        lambda_value=0.0,
    )
    assert report.interaction_term_count == 1
    row = report.rows[0]
    assert row.interaction_term == "predictor_one:habitat"
    assert row.interaction_kind == "continuous-by-categorical"
    assert row.coefficient_name == "predictor_one:habitat[tundra]"
    assert row.component_terms == ["predictor_one", "habitat"]
    assert row.component_levels == [None, "tundra"]
    assert row.omitted_reference_levels == ["habitat=forest"]
    assert math.isclose(row.estimate, 0.5, abs_tol=1e-12)


def test_summarize_pgls_interaction_coefficients_supports_continuous_by_continuous() -> (
    None
):
    report = summarize_pgls_interaction_coefficients(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_continuous_interaction.tsv"),
        formula="response ~ predictor_one * predictor_two",
        lambda_value=0.0,
    )
    row = report.rows[0]
    assert row.interaction_kind == "continuous-by-continuous"
    assert row.component_columns == ["predictor_one", "predictor_two"]
    assert row.component_levels == [None, None]
    assert row.omitted_reference_levels == []


def test_summarize_pgls_interaction_coefficients_supports_categorical_by_categorical() -> (
    None
):
    report = summarize_pgls_interaction_coefficients(
        fixture("example_tree_eight_taxa.nwk"),
        fixture("example_traits_comparative_categorical_interaction.tsv"),
        formula="response ~ habitat * diet",
        lambda_value=0.0,
    )
    row = report.rows[0]
    assert row.interaction_kind == "categorical-by-categorical"
    assert row.component_columns == ["habitat[tundra]", "diet[herbivore]"]
    assert row.component_levels == ["tundra", "herbivore"]
    assert row.omitted_reference_levels == ["habitat=forest", "diet=carnivore"]
    assert math.isclose(row.estimate, 3.7, abs_tol=1e-12)


def test_pgls_interaction_outputs_match_reference_fixture() -> None:
    reference = json.loads(
        fixture("pgls_interaction_reference.json").read_text(encoding="utf-8")
    )
    cases = {
        "continuous-by-continuous": (
            fixture("example_tree_six_taxa.nwk"),
            fixture("example_traits_comparative_continuous_interaction.tsv"),
        ),
        "continuous-by-categorical": (
            fixture("example_tree_six_taxa.nwk"),
            fixture("example_traits_comparative_interaction.tsv"),
        ),
        "categorical-by-categorical": (
            fixture("example_tree_eight_taxa.nwk"),
            fixture("example_traits_comparative_categorical_interaction.tsv"),
        ),
    }
    for observation in reference["observations"]:
        tree_path, traits_path = cases[observation["case"]]
        matrix_report = build_pgls_model_matrix(
            tree_path,
            traits_path,
            formula=observation["formula"],
        )
        assert matrix_report.encoded_columns == observation["encoded_columns"]
        rows_by_taxon = {row.taxon: row for row in matrix_report.rows}
        for expected_row in observation["rows"]:
            actual_row = rows_by_taxon[expected_row["taxon"]]
            for column in observation["encoded_columns"]:
                assert actual_row.encoded_values[column] == expected_row[column]
        model_report = run_pgls(
            tree_path,
            traits_path,
            formula=observation["formula"],
            lambda_value=0.0,
        )
        coefficients = {
            coefficient.name: coefficient.estimate
            for coefficient in model_report.coefficients
        }
        for coefficient_name, expected_estimate in observation["coefficients"].items():
            assert math.isclose(
                coefficients[coefficient_name],
                expected_estimate,
                rel_tol=1e-9,
                abs_tol=1e-9,
            )


def test_write_pgls_interaction_coefficient_table_writes_rows(tmp_path: Path) -> None:
    report = summarize_pgls_interaction_coefficients(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_interaction.tsv"),
        formula="response ~ predictor_one * habitat",
        lambda_value=0.0,
    )
    out_path = tmp_path / "interaction-coefficients.tsv"
    write_pgls_interaction_coefficient_table(out_path, report)
    contents = out_path.read_text(encoding="utf-8")
    assert "interaction_term\tinteraction_kind\tcoefficient_name" in contents
    assert "predictor_one:habitat\tcontinuous-by-categorical" in contents
