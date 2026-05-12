from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from bijux_phylogenetics.comparative.phylogenetic_logistic import (
    summarize_phylogenetic_logistic,
    write_phylogenetic_logistic_coefficient_table,
    write_phylogenetic_logistic_excluded_taxa_table,
    write_phylogenetic_logistic_fitted_table,
)
from bijux_phylogenetics.errors import ComparativeMethodError

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


def test_summarize_phylogenetic_logistic_fits_binary_response() -> None:
    report = summarize_phylogenetic_logistic(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_phylogenetic_logistic.tsv"),
        response="presence",
        predictors=["body_size"],
    )
    coefficients = {row.name: row for row in report.coefficients}
    assert report.approximation_method == "phylogenetic-working-correlation-gee"
    assert report.taxon_count == 6
    assert report.success_count == 3
    assert report.failure_count == 3
    assert report.lambda_value == 1.0
    assert report.converged is True
    assert report.iteration_count >= 1
    assert report.separation_detected is False
    assert len(report.coefficients) == 2
    assert coefficients["body_size"].estimate > 0.0
    assert coefficients["body_size"].inference_distribution == "wald-normal"
    assert all(0.0 < row.fitted_probability < 1.0 for row in report.fitted_rows)
    assert math.isfinite(report.binomial_log_likelihood)


def test_summarize_phylogenetic_logistic_rejects_non_binary_response() -> None:
    with pytest.raises(ComparativeMethodError) as error:
        summarize_phylogenetic_logistic(
            fixture("example_tree.nwk"),
            fixture("example_traits_comparative.tsv"),
            response="response",
            predictors=["predictor_one"],
        )
    assert "requires a binary response encoded as 0 and 1" in str(error.value)


def test_summarize_phylogenetic_logistic_reports_separation_risk() -> None:
    report = summarize_phylogenetic_logistic(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_phylogenetic_logistic_separated.tsv"),
        formula="presence ~ habitat",
    )
    warning_codes = {row.code for row in report.warnings}
    assert report.separation_detected is True
    assert "large_coefficient_magnitude" in warning_codes or (
        "fitted_probability_boundary" in warning_codes
    )


def test_summarize_phylogenetic_logistic_approximates_reference_fixture() -> None:
    reference = json.loads(
        fixture("phylogenetic_logistic_reference.json").read_text(encoding="utf-8")
    )
    report = summarize_phylogenetic_logistic(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_phylogenetic_logistic.tsv"),
        response="presence",
        predictors=["body_size"],
    )
    coefficients = {row.name: row.estimate for row in report.coefficients}
    fitted_probabilities = {row.taxon: row.fitted_probability for row in report.fitted_rows}
    expected = reference["phyloglm"]
    assert abs(coefficients["intercept"] - expected["coefficients"]["intercept"]) < 0.5
    assert math.copysign(1.0, coefficients["intercept"]) == math.copysign(
        1.0, expected["coefficients"]["intercept"]
    )
    assert math.copysign(1.0, coefficients["body_size"]) == math.copysign(
        1.0, expected["coefficients"]["body_size"]
    )
    assert abs(coefficients["body_size"] - expected["coefficients"]["body_size"]) < 0.2
    assert all(
        abs(fitted_probabilities[taxon] - probability) < 0.1
        for taxon, probability in expected["fitted_probabilities"].items()
    )


def test_write_phylogenetic_logistic_tables_write_rows(tmp_path: Path) -> None:
    report = summarize_phylogenetic_logistic(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_phylogenetic_logistic.tsv"),
        response="presence",
        predictors=["body_size"],
    )
    coefficients_out = tmp_path / "phylogenetic-logistic-coefficients.tsv"
    fitted_out = tmp_path / "phylogenetic-logistic-fitted.tsv"
    excluded_out = tmp_path / "phylogenetic-logistic-excluded.tsv"
    write_phylogenetic_logistic_coefficient_table(coefficients_out, report)
    write_phylogenetic_logistic_fitted_table(fitted_out, report)
    write_phylogenetic_logistic_excluded_taxa_table(excluded_out, report)
    coefficient_rows = coefficients_out.read_text(encoding="utf-8").splitlines()
    fitted_rows = fitted_out.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    assert coefficient_rows[0].startswith("response\tterm\testimate")
    assert len(coefficient_rows) == 3
    assert fitted_rows[0].startswith(
        "taxon\tobserved_response\tfitted_probability\tlinear_predictor"
    )
    assert len(fitted_rows) == 7
    assert excluded_rows == ["taxon\treason\tdetails"]
