from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from bijux_phylogenetics.comparative.pgls.ou_covariance import (
    summarize_ou_covariance_pgls,
    write_ou_alpha_profile_table,
    write_ou_covariance_table,
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


def test_summarize_ou_covariance_pgls_matches_reference_fixture_case() -> None:
    expected = json.loads(
        fixture("ou_covariance_pgls_reference.json").read_text(encoding="utf-8")
    )
    report = summarize_ou_covariance_pgls(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        alpha=expected["alpha"],
    )
    coefficients = {row.name: row.estimate for row in report.model.coefficients}
    assert report.alpha_estimation_mode == "fixed"
    assert report.taxon_count == 4
    assert len(report.rows) == 16
    assert math.isclose(report.alpha, expected["alpha"], abs_tol=1e-12)
    assert math.isclose(
        report.model.log_likelihood,
        expected["log_likelihood"],
        rel_tol=1e-12,
        abs_tol=1e-12,
    )
    assert math.isclose(report.model.aic, expected["aic"], rel_tol=1e-12, abs_tol=1e-12)
    assert math.isclose(
        coefficients["intercept"], expected["coefficients"]["intercept"], abs_tol=1e-12
    )
    assert math.isclose(
        coefficients["predictor_one"],
        expected["coefficients"]["predictor_one"],
        abs_tol=1e-12,
    )


def test_summarize_ou_covariance_pgls_estimates_alpha_and_profile() -> None:
    report = summarize_ou_covariance_pgls(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        alpha="estimate",
    )
    assert report.alpha_estimation_mode == "estimated"
    assert len(report.alpha_profile_rows) == 8
    assert report.lower_95_confidence_interval is not None
    assert report.upper_95_confidence_interval is not None
    assert (
        report.lower_95_confidence_interval
        <= report.alpha
        <= report.upper_95_confidence_interval
    )
    assert any(row.within_95_confidence_interval for row in report.alpha_profile_rows)


def test_summarize_ou_covariance_pgls_rejects_non_positive_alpha() -> None:
    with pytest.raises(ComparativeMethodError, match="OU alpha must be positive"):
        summarize_ou_covariance_pgls(
            fixture("example_tree.nwk"),
            fixture("example_traits_comparative.tsv"),
            response="response",
            predictors=["predictor_one"],
            alpha=0.0,
        )


def test_summarize_ou_covariance_pgls_detects_invalid_covariance() -> None:
    with pytest.raises(ComparativeMethodError) as error:
        summarize_ou_covariance_pgls(
            fixture("example_tree_zero_lengths.nwk"),
            fixture("example_traits_comparative.tsv"),
            response="response",
            predictors=["predictor_one"],
            alpha=1.0,
        )
    assert "OU covariance is invalid" in str(error.value)


def test_summarize_ou_covariance_pgls_reports_negative_branch_length_details() -> None:
    with pytest.raises(ComparativeMethodError) as error:
        summarize_ou_covariance_pgls(
            fixture("example_tree_negative_length.nwk"),
            fixture("example_traits_comparative.tsv"),
            response="response",
            predictors=["predictor_one"],
            alpha=1.0,
        )

    assert (
        error.value.details["failure_reason"] == "ou_covariance_negative_branch_lengths"
    )
    assert error.value.details["evidence"]["tree_path"].endswith(
        "example_tree_negative_length.nwk"
    )


def test_write_ou_covariance_tables_write_rows(tmp_path: Path) -> None:
    covariance_out = tmp_path / "ou-covariance.tsv"
    profile_out = tmp_path / "ou-alpha-profile.tsv"
    report = summarize_ou_covariance_pgls(
        fixture("example_tree_internal_long_branch.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        alpha="estimate",
    )
    write_ou_covariance_table(covariance_out, report)
    write_ou_alpha_profile_table(profile_out, report)
    covariance_rows = covariance_out.read_text(encoding="utf-8").splitlines()
    profile_rows = profile_out.read_text(encoding="utf-8").splitlines()
    assert covariance_rows[0].startswith("left_taxon\tright_taxon\tis_diagonal")
    assert len(covariance_rows) == 17
    assert profile_rows[0].startswith("alpha_estimation_mode\talpha\tlog_likelihood")
    assert len(profile_rows) == 9
