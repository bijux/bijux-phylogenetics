from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from bijux_phylogenetics.comparative.pgls import build_pgls_model_matrix, run_pgls

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


def _core_reference_case(case: str) -> dict[str, object]:
    payload = json.loads(
        fixture("reference_parity_core.json").read_text(encoding="utf-8")
    )
    return next(item for item in payload["observations"] if item["case"] == case)


def _extended_reference_case(case: str) -> dict[str, object]:
    payload = json.loads(
        fixture("reference_parity_extended_comparative.json").read_text(
            encoding="utf-8"
        )
    )
    return next(item for item in payload["observations"] if item["case"] == case)


def _assert_pgls_matches_reference(
    report,
    expected_output: dict[str, float],
    *,
    tolerance: float,
) -> None:
    assert math.isclose(
        report.log_likelihood,
        expected_output["log_likelihood"],
        rel_tol=tolerance,
        abs_tol=tolerance,
    )
    assert math.isclose(
        report.aic,
        expected_output["aic"],
        rel_tol=tolerance,
        abs_tol=tolerance,
    )
    assert math.isclose(
        report.lambda_value,
        expected_output["lambda_value"],
        rel_tol=tolerance,
        abs_tol=tolerance,
    )
    coefficients = {
        coefficient.name: coefficient for coefficient in report.coefficients
    }
    for coefficient in report.coefficients:
        prefix = f"coefficient.{coefficient.name}"
        assert math.isclose(
            coefficients[coefficient.name].estimate,
            expected_output[f"{prefix}.estimate"],
            rel_tol=tolerance,
            abs_tol=tolerance,
        )
        assert math.isclose(
            coefficients[coefficient.name].standard_error,
            expected_output[f"{prefix}.standard_error"],
            rel_tol=tolerance,
            abs_tol=tolerance,
        )
        assert math.isclose(
            coefficients[coefficient.name].p_value,
            expected_output[f"{prefix}.p_value"],
            rel_tol=tolerance,
            abs_tol=tolerance,
        )


def test_run_pgls_matches_simple_brownian_reference_case() -> None:
    reference = _core_reference_case("pgls-example-tree-brownian")
    report = run_pgls(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        lambda_value=1.0,
    )
    _assert_pgls_matches_reference(
        report,
        reference["expected_output"],
        tolerance=float(reference["tolerance"]),
    )


def test_run_pgls_matches_categorical_brownian_reference_case() -> None:
    reference = _core_reference_case("pgls-example-tree-brownian-categorical")
    report = run_pgls(
        fixture("example_tree_eight_taxa.nwk"),
        fixture("example_traits_comparative_categorical_interaction.tsv"),
        formula="response ~ habitat + diet",
        lambda_value=1.0,
    )
    matrix = build_pgls_model_matrix(
        fixture("example_tree_eight_taxa.nwk"),
        fixture("example_traits_comparative_categorical_interaction.tsv"),
        formula="response ~ habitat + diet",
    )
    _assert_pgls_matches_reference(
        report,
        reference["expected_output"],
        tolerance=float(reference["tolerance"]),
    )
    row_by_taxon = {row.taxon: row for row in matrix.rows}
    assert row_by_taxon["G"].encoded_values["habitat[tundra]"] == 1.0
    assert row_by_taxon["G"].encoded_values["diet[herbivore]"] == 1.0


def test_run_pgls_matches_interaction_brownian_reference_case() -> None:
    reference = _core_reference_case("pgls-example-tree-brownian-interaction")
    report = run_pgls(
        fixture("example_tree_eight_taxa.nwk"),
        fixture("example_traits_comparative_categorical_interaction.tsv"),
        formula="response ~ habitat * diet",
        lambda_value=1.0,
    )
    matrix = build_pgls_model_matrix(
        fixture("example_tree_eight_taxa.nwk"),
        fixture("example_traits_comparative_categorical_interaction.tsv"),
        formula="response ~ habitat * diet",
    )
    _assert_pgls_matches_reference(
        report,
        reference["expected_output"],
        tolerance=float(reference["tolerance"]),
    )
    row_by_taxon = {row.taxon: row for row in matrix.rows}
    assert row_by_taxon["H"].encoded_values["habitat[tundra]:diet[herbivore]"] == 1.0


@pytest.mark.slow
def test_run_pgls_matches_primate_fixed_reference_lambda_case() -> None:
    repository_root = Path(__file__).resolve().parents[3]
    reference = _extended_reference_case(
        "pgls-primate-longevity-fixed-reference-lambda"
    )
    report = run_pgls(
        repository_root
        / "evidence-book"
        / "studies"
        / "primate-longevity-signal"
        / "datasets"
        / "reference_trimmed_primatetree.nwk",
        repository_root
        / "evidence-book"
        / "studies"
        / "primate-longevity-signal"
        / "datasets"
        / "reference_primate.csv",
        response="longevity",
        predictors=["social_group_size"],
        taxon_column="species",
        lambda_value=float(reference["lambda_value"]),
    )
    _assert_pgls_matches_reference(
        report,
        reference["expected_output"],
        tolerance=float(reference["tolerance"]),
    )


@pytest.mark.slow
def test_run_pgls_matches_primate_estimated_lambda_reference_case() -> None:
    repository_root = Path(__file__).resolve().parents[3]
    reference = _extended_reference_case("pgls-primate-longevity-estimated-lambda")
    report = run_pgls(
        repository_root
        / "evidence-book"
        / "studies"
        / "primate-longevity-signal"
        / "datasets"
        / "reference_trimmed_primatetree.nwk",
        repository_root
        / "evidence-book"
        / "studies"
        / "primate-longevity-signal"
        / "datasets"
        / "reference_primate.csv",
        response="longevity",
        predictors=["social_group_size"],
        taxon_column="species",
        lambda_value="estimate",
    )
    _assert_pgls_matches_reference(
        report,
        reference["expected_output"],
        tolerance=float(reference["tolerance"]),
    )
