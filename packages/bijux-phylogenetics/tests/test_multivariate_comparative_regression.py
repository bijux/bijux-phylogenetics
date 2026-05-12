from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative.multivariate_regression import (
    run_multivariate_comparative_regression,
    write_multivariate_excluded_taxa_table,
    write_multivariate_residual_association_table,
    write_multivariate_residual_covariance_table,
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


def test_run_multivariate_comparative_regression_reports_covariance_and_association() -> (
    None
):
    report = run_multivariate_comparative_regression(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multiple.tsv"),
        responses=["response_growth", "response_range"],
        predictors=["predictor_one", "predictor_two"],
        lambda_value=0.0,
    )
    assert report.responses == ["response_growth", "response_range"]
    assert report.predictors == ["predictor_one", "predictor_two"]
    assert report.analysis_taxa == ["A", "B", "C", "D", "E", "F"]
    assert report.excluded_taxa == []
    assert len(report.response_models) == 2
    assert len(report.covariance_rows) == 4
    assert len(report.association_rows) == 1
    association = report.association_rows[0]
    assert association.left_response == "response_growth"
    assert association.right_response == "response_range"
    assert association.pair_count == 6
    assert math.isclose(association.correlation, 0.0, abs_tol=1e-12)
    assert association.p_value == 1.0
    diagonal = next(
        row
        for row in report.covariance_rows
        if row.left_response == "response_growth"
        and row.right_response == "response_growth"
    )
    assert diagonal.is_diagonal is True
    assert math.isclose(diagonal.correlation, 1.0)


def test_run_multivariate_comparative_regression_excludes_incomplete_taxa() -> None:
    report = run_multivariate_comparative_regression(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multivariate_missing.tsv"),
        responses=["response_growth", "response_range"],
        predictors=["predictor_one", "predictor_two"],
        lambda_value=0.0,
    )
    assert report.analysis_taxa == ["A", "C", "D", "E", "F"]
    assert len(report.excluded_taxa) == 1
    excluded = report.excluded_taxa[0]
    assert excluded.taxon == "B"
    assert excluded.reason == "missing_required_values"
    assert excluded.missing_columns == ["response_range"]


def test_write_multivariate_regression_tables_write_review_ledgers(
    tmp_path: Path,
) -> None:
    report = run_multivariate_comparative_regression(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multivariate_missing.tsv"),
        responses=["response_growth", "response_range"],
        predictors=["predictor_one", "predictor_two"],
        lambda_value=0.0,
    )
    covariance_path = tmp_path / "multivariate-residual-covariance.tsv"
    association_path = tmp_path / "multivariate-residual-associations.tsv"
    excluded_path = tmp_path / "multivariate-excluded-taxa.tsv"
    write_multivariate_residual_covariance_table(covariance_path, report)
    write_multivariate_residual_association_table(association_path, report)
    write_multivariate_excluded_taxa_table(excluded_path, report)
    covariance_rows = covariance_path.read_text(encoding="utf-8").splitlines()
    association_rows = association_path.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_path.read_text(encoding="utf-8").splitlines()
    assert covariance_rows[0].startswith(
        "left_response\tright_response\tpair_count\tis_diagonal"
    )
    assert association_rows[0].startswith(
        "left_response\tright_response\tpair_count\tcovariance\tcorrelation"
    )
    assert excluded_rows[0] == "taxon\treason\tmissing_columns"
    assert len(covariance_rows) == 5
    assert len(association_rows) == 2
    assert len(excluded_rows) == 2
