from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative.covariance import (
    summarize_comparative_covariance_audit,
    write_comparative_covariance_audit_candidate_table,
    write_comparative_covariance_audit_excluded_taxa_table,
    write_comparative_covariance_audit_summary_table,
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


def test_summarize_comparative_covariance_audit_reports_pgls_profile_surface() -> None:
    report = summarize_comparative_covariance_audit(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        analysis="pgls",
        response="response",
        predictors=["predictor_one"],
        lambda_value="estimate",
    )
    assert report.analysis == "pgls"
    assert report.covariance_model == "pagel-lambda"
    assert report.analysis_label == "response ~ predictor_one"
    assert report.tree_taxon_count == 4
    assert report.trait_taxon_count == 4
    assert report.matched_taxa == ["A", "B", "C", "D"]
    assert report.analysis_taxa == ["A", "B", "C", "D"]
    assert report.matrix_dimension == 4
    assert report.matrix_rank == 4
    assert report.fit_strategy == "exact"
    assert report.singular is False
    assert report.near_singular is False
    assert len(report.candidate_rows) == 101
    assert report.candidate_rows[0].candidate_label == "lambda=0.00"
    assert report.candidate_rows[-1].candidate_label == "lambda=1.00"
    assert all(row.fit_strategy == "exact" for row in report.candidate_rows)
    assert report.blockers == []


def test_summarize_comparative_covariance_audit_detects_duplicate_trait_taxa() -> None:
    report = summarize_comparative_covariance_audit(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative_duplicate.tsv"),
        analysis="pgls",
        response="response",
        predictors=["predictor_one"],
        lambda_value=1.0,
    )
    assert report.duplicate_trait_taxa == ["A"]
    assert report.fit_strategy == "failure"
    assert report.candidate_rows == []
    assert "trait table contains duplicate taxon keys" in report.blockers


def test_summarize_comparative_covariance_audit_detects_duplicate_tree_taxa() -> None:
    report = summarize_comparative_covariance_audit(
        fixture("example_tree_duplicate.nwk"),
        fixture("example_traits_comparative.tsv"),
        analysis="pgls",
        response="response",
        predictors=["predictor_one"],
        lambda_value=1.0,
    )
    assert report.duplicate_tree_taxa == ["A"]
    assert report.fit_strategy == "failure"
    assert report.candidate_rows == []
    assert "tree contains duplicate tip labels" in report.blockers


def test_summarize_comparative_covariance_audit_detects_invalid_branch_lengths() -> (
    None
):
    report = summarize_comparative_covariance_audit(
        fixture("example_tree_negative_length.nwk"),
        fixture("example_traits_comparative.tsv"),
        analysis="brownian-trait",
        trait="response",
    )
    assert report.analysis == "brownian-trait"
    assert report.fit_strategy == "failure"
    assert report.negative_branch_length_count == 1
    assert report.candidate_rows == []
    assert any("negative branch lengths" in blocker for blocker in report.blockers)


def test_summarize_comparative_covariance_audit_reports_taxon_overlap_mismatches() -> (
    None
):
    report = summarize_comparative_covariance_audit(
        fixture("example_tree.nwk"),
        fixture("example_traits.tsv"),
        analysis="brownian-trait",
        trait="value",
    )
    assert report.tree_taxon_count == 4
    assert report.trait_taxon_count == 4
    assert report.matched_taxa == ["A", "B", "C"]
    assert report.missing_from_traits == ["D"]
    assert report.extra_trait_taxa == ["E"]
    assert report.analysis_taxa == ["A", "B", "C"]
    assert report.fit_strategy == "exact"
    assert report.blockers == []


def test_summarize_comparative_covariance_audit_reports_zero_length_regularization() -> (
    None
):
    report = summarize_comparative_covariance_audit(
        fixture("example_tree_zero_lengths.nwk"),
        fixture("example_traits_comparative.tsv"),
        analysis="pgls",
        response="response",
        predictors=["predictor_one"],
        lambda_value=1.0,
    )
    assert report.fit_strategy == "regularization"
    assert report.zero_length_branch_count == 3
    assert report.matrix_dimension == 4
    assert report.matrix_rank == 3
    assert report.singular is True
    assert report.near_singular is True
    assert any("zero-length branches" in warning for warning in report.warnings)
    assert report.candidate_rows[0].fit_strategy == "regularization"
    assert report.candidate_rows[0].positive_definite_before_fit is False
    assert math.isfinite(report.candidate_rows[0].fit_condition_number or math.inf)


def test_summarize_comparative_covariance_audit_reports_ou_alpha_candidates() -> None:
    report = summarize_comparative_covariance_audit(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        analysis="ou-trait",
        trait="response",
        alpha="estimate",
    )
    assert report.analysis == "ou-trait"
    assert report.covariance_model == "ou"
    assert report.analysis_label == "response"
    assert report.matrix_dimension == 4
    assert report.matrix_rank == 4
    assert report.fit_strategy == "exact"
    assert len(report.candidate_rows) == 8
    assert all(row.parameter_name == "alpha" for row in report.candidate_rows)
    assert all(row.parameter_value is not None for row in report.candidate_rows)
    assert all(row.fit_strategy == "exact" for row in report.candidate_rows)


def test_write_comparative_covariance_audit_tables_write_rows(tmp_path: Path) -> None:
    report = summarize_comparative_covariance_audit(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative_missing_predictor.tsv"),
        analysis="pgls",
        response="response",
        predictors=["predictor_one"],
        lambda_value=1.0,
    )
    summary_out = tmp_path / "covariance-audit-summary.tsv"
    candidates_out = tmp_path / "covariance-audit-candidates.tsv"
    excluded_out = tmp_path / "covariance-audit-excluded.tsv"
    write_comparative_covariance_audit_summary_table(summary_out, report)
    write_comparative_covariance_audit_candidate_table(candidates_out, report)
    write_comparative_covariance_audit_excluded_taxa_table(excluded_out, report)
    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    candidate_rows = candidates_out.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("analysis\tcovariance_model\tanalysis_label")
    assert candidate_rows[0].startswith("candidate_label\tparameter_name")
    assert excluded_rows[0] == "taxon\treason\tdetails"
    assert len(summary_rows) == 2
    assert len(candidate_rows) == 2
    assert len(excluded_rows) == 2
