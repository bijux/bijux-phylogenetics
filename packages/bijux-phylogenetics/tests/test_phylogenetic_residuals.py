from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative.regression import (
    summarize_phylogenetic_residuals,
    write_phylogenetic_residual_coefficient_table,
    write_phylogenetic_residual_exclusion_table,
    write_phylogenetic_residual_summary_table,
    write_phylogenetic_residual_taxon_table,
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


def test_summarize_phylogenetic_residuals_flags_body_size_outlier() -> None:
    report = summarize_phylogenetic_residuals(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_phylogenetic_residuals.tsv"),
        response="brain_mass",
        predictor="body_mass",
        method="lambda",
    )

    assert report.tree_taxon_count == 6
    assert report.analyzed_taxa == ["A", "B", "C", "D", "E", "F"]
    assert report.analyzed_taxon_count == 6
    assert report.method == "lambda"
    assert report.lambda_estimation_mode == "estimated"
    assert report.coefficient_rows[0].name == "intercept"
    assert report.coefficient_rows[0].estimate > 0.0
    assert report.coefficient_rows[1].name == "body_mass"
    assert "F" in report.outlier_taxa
    assert report.max_abs_standardized_residual is not None
    assert report.max_abs_standardized_residual >= 2.0
    top_row = max(report.taxon_rows, key=lambda row: row.abs_standardized_residual)
    assert top_row.taxon == "F"
    assert top_row.outlier is True
    assert top_row.input_order == 6
    assert top_row.tree_tip_label == "F"
    assert top_row.fitted_value != top_row.observed_value


def test_summarize_phylogenetic_residuals_reports_missing_value_and_extra_taxon() -> (
    None
):
    report = summarize_phylogenetic_residuals(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_phylogenetic_residuals_missing.tsv"),
        response="brain_mass",
        predictor="body_mass",
        method="brownian",
    )

    assert report.method == "brownian"
    assert report.lambda_estimation_mode == "fixed"
    assert math.isclose(report.lambda_value, 1.0)
    assert report.analyzed_taxa == ["A", "B", "C", "D", "F"]
    assert {row.taxon: row.reason for row in report.excluded_taxa} == {
        "E": "missing_value",
        "G": "absent_from_tree",
    }
    assert "F" in report.outlier_taxa


def test_phylogenetic_residual_writers_emit_review_ledgers(tmp_path: Path) -> None:
    report = summarize_phylogenetic_residuals(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_phylogenetic_residuals.tsv"),
        response="brain_mass",
        predictor="body_mass",
        method="lambda",
    )
    summary_out = tmp_path / "phylogenetic-residual-summary.tsv"
    taxon_out = tmp_path / "phylogenetic-residuals.tsv"
    coefficient_out = tmp_path / "phylogenetic-residual-coefficients.tsv"
    excluded_out = tmp_path / "phylogenetic-residual-excluded.tsv"

    write_phylogenetic_residual_summary_table(summary_out, report)
    write_phylogenetic_residual_taxon_table(taxon_out, report)
    write_phylogenetic_residual_coefficient_table(coefficient_out, report)
    write_phylogenetic_residual_exclusion_table(excluded_out, report)

    assert (
        summary_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("response\tpredictor\tmethod")
    )
    assert (
        taxon_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith(
            "taxon\tinput_order\ttree_tip_label\tobserved_value\tfitted_value\tresidual"
        )
    )
    assert (
        coefficient_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("name\testimate\tstandard_error\tp_value")
    )
    assert excluded_out.read_text(encoding="utf-8").splitlines() == [
        "taxon\treason\tdetails"
    ]
