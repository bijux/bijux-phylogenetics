from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from bijux_phylogenetics.comparative.pgls import run_pgls
from bijux_phylogenetics.comparative.pgls.lambda_fit import (
    summarize_pgls_lambda_fit,
    write_pgls_lambda_profile_table,
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


def test_run_pgls_reports_fixed_lambda_fit_surface() -> None:
    report = run_pgls(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        lambda_value=0.0,
    )
    assert report.lambda_fit.mode == "fixed"
    assert report.lambda_fit.lambda_value == 0.0
    assert report.lambda_fit.lower_95_confidence_interval is None
    assert report.lambda_fit.upper_95_confidence_interval is None
    assert len(report.lambda_fit.profile_rows) == 1
    assert report.lambda_fit.profile_rows[0].lambda_value == 0.0


@pytest.mark.slow
def test_run_pgls_estimated_lambda_matches_primate_reference_bundle() -> None:
    repository_root = Path(__file__).resolve().parents[3]
    parity_fixture = (
        repository_root
        / "evidence-book"
        / "studies"
        / "primate-pgls-and-signal"
        / "evidence-003"
        / "results"
        / "pagel-lambda-regression-parity.json"
    )
    reference = json.loads(parity_fixture.read_text(encoding="utf-8"))
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
    assert report.lambda_fit.mode == "estimated"
    assert math.isclose(
        report.lambda_value,
        reference["r_estimated_lambda"]["lambda_value"],
        abs_tol=0.05,
    )
    assert math.isclose(
        report.log_likelihood,
        reference["r_estimated_lambda"]["log_likelihood"],
        rel_tol=5e-4,
        abs_tol=0.15,
    )
    assert report.lambda_fit.lower_95_confidence_interval is not None
    assert report.lambda_fit.upper_95_confidence_interval is not None
    assert (
        report.lambda_fit.lower_95_confidence_interval
        <= reference["r_estimated_lambda"]["lambda_value"]
        <= report.lambda_fit.upper_95_confidence_interval
    )
    assert len(report.lambda_fit.profile_rows) == 101


def test_write_pgls_lambda_profile_table_writes_profile_rows(tmp_path: Path) -> None:
    report = summarize_pgls_lambda_fit(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        lambda_value="estimate",
    )
    out_path = tmp_path / "lambda-profile.tsv"
    write_pgls_lambda_profile_table(out_path, report)
    lines = out_path.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("mode\tlambda_value\tlog_likelihood")
    assert len(lines) == 102
    assert lines[1].startswith("estimated\t")
    assert any("\ttrue\t" in line for line in lines[1:])
