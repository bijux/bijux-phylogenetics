from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from bijux_phylogenetics.comparative.pgls.brownian_covariance import (
    summarize_brownian_covariance_pgls,
    write_brownian_covariance_table,
)
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

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


def test_summarize_brownian_covariance_pgls_matches_reference_fixture_case() -> None:
    fixture_path = fixture("comparative_reference_validation.json")
    reference = json.loads(fixture_path.read_text(encoding="utf-8"))
    expected = next(
        case
        for case in reference["observations"]
        if case["case"] == "pgls-example-tree-brownian"
    )
    report = summarize_brownian_covariance_pgls(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
    )
    coefficients = {row.name: row.estimate for row in report.model.coefficients}
    assert report.tree_is_ultrametric is True
    assert report.positive_definite_before_stabilization is True
    assert report.taxon_count == 4
    assert len(report.rows) == 16
    assert math.isclose(report.model.lambda_value, 1.0, abs_tol=1e-12)
    assert math.isclose(
        report.model.log_likelihood,
        expected["expected_parameters"]["log_likelihood"],
        rel_tol=1e-6,
        abs_tol=1e-6,
    )
    assert math.isclose(
        coefficients["intercept"],
        expected["expected_parameters"]["intercept"],
        rel_tol=1e-6,
        abs_tol=1e-6,
    )
    assert math.isclose(
        coefficients["predictor_one"],
        expected["expected_parameters"]["predictor_one"],
        rel_tol=1e-6,
        abs_tol=1e-6,
    )


def test_summarize_brownian_covariance_pgls_handles_rooted_non_ultrametric_tree() -> (
    None
):
    report = summarize_brownian_covariance_pgls(
        fixture("example_tree_internal_long_branch.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
    )
    diagonal_rows = [row for row in report.rows if row.is_diagonal]
    assert report.tree_is_ultrametric is False
    assert report.taxon_count == 4
    assert len(report.rows) == 16
    assert report.minimum_root_to_tip_depth == 0.2
    assert report.maximum_root_to_tip_depth == 1.1
    assert all(
        math.isclose(row.shared_path_length, row.left_root_depth, abs_tol=1e-12)
        for row in diagonal_rows
    )


def test_summarize_brownian_covariance_pgls_detects_invalid_covariance() -> None:
    with pytest.raises(ComparativeMethodError) as error:
        summarize_brownian_covariance_pgls(
            fixture("example_tree_zero_lengths.nwk"),
            fixture("example_traits_comparative.tsv"),
            response="response",
            predictors=["predictor_one"],
        )
    assert "Brownian covariance is invalid" in str(error.value)
    assert "non-positive root-to-tip path length" in str(error.value)


def test_summarize_brownian_covariance_pgls_reports_negative_branch_length_details() -> (
    None
):
    with pytest.raises(ComparativeMethodError) as error:
        summarize_brownian_covariance_pgls(
            fixture("example_tree_negative_length.nwk"),
            fixture("example_traits_comparative.tsv"),
            response="response",
            predictors=["predictor_one"],
        )

    assert (
        error.value.details["failure_reason"]
        == "brownian_covariance_negative_branch_lengths"
    )
    assert error.value.details["evidence"]["tree_path"].endswith(
        "example_tree_negative_length.nwk"
    )


def test_write_brownian_covariance_table_writes_pairwise_rows(tmp_path: Path) -> None:
    out_path = tmp_path / "brownian-covariance.tsv"
    report = summarize_brownian_covariance_pgls(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
    )
    write_brownian_covariance_table(out_path, report)
    written_rows = out_path.read_text(encoding="utf-8").splitlines()
    assert written_rows[0].startswith("left_taxon\tright_taxon\tis_diagonal")
    assert len(written_rows) == 17
    assert any(
        row.startswith("A\tA\ttrue\t0.3\t0.3\t0.3\ttrue") for row in written_rows[1:]
    )
