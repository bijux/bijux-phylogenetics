from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative.pgls import inspect_pgls_inputs, run_pgls


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


def test_pgls_input_inspection_encodes_categorical_predictors() -> None:
    report = inspect_pgls_inputs(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one", "habitat"],
    )
    assert report.ready is True
    assert report.categorical_predictors == ["habitat"]
    habitat = next(row for row in report.predictors if row.name == "habitat")
    assert habitat.reference_level == "forest"
    assert habitat.encoded_columns == ["habitat[tundra]"]
    assert report.encoded_columns == ["intercept", "predictor_one", "habitat[tundra]"]


def test_pgls_input_inspection_rejects_missing_branch_lengths() -> None:
    report = inspect_pgls_inputs(
        fixture("example_tree_no_lengths.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
    )
    assert report.ready is False
    assert "PGLS requires complete tree branch lengths" in report.blockers


def test_run_pgls_supports_one_predictor() -> None:
    report = run_pgls(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        lambda_value=0.0,
    )
    assert report.predictors == ["predictor_one"]
    assert report.taxon_count == 4
    assert len(report.coefficients) == 2
    assert math.isclose(report.coefficients[0].estimate, 1.0, abs_tol=1e-6)
    assert math.isclose(report.coefficients[1].estimate, 0.7, abs_tol=1e-6)


def test_run_pgls_supports_multiple_predictors() -> None:
    report = run_pgls(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one", "predictor_two"],
        lambda_value=0.0,
    )
    coefficients = {coefficient.name: coefficient.estimate for coefficient in report.coefficients}
    assert math.isclose(coefficients["intercept"], 1.0, abs_tol=1e-6)
    assert math.isclose(coefficients["predictor_one"], 0.5, abs_tol=1e-6)
    assert math.isclose(coefficients["predictor_two"], 1.0, abs_tol=1e-6)


def test_run_pgls_supports_categorical_predictors_with_dummy_encoding() -> None:
    report = run_pgls(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one", "habitat"],
        lambda_value=0.0,
    )
    coefficients = {coefficient.name: coefficient.estimate for coefficient in report.coefficients}
    assert report.encoded_columns == ["intercept", "predictor_one", "habitat[tundra]"]
    assert math.isclose(coefficients["predictor_one"], 1.5, abs_tol=1e-6)
    assert math.isclose(coefficients["habitat[tundra]"], -2.0, abs_tol=1e-6)
    assert report.diagnostics.outlier_taxa == []
