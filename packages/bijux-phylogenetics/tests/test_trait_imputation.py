from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.traits.imputation import (
    summarize_trait_imputation,
    write_trait_imputation_exclusion_table,
    write_trait_imputation_holdout_table,
    write_trait_imputation_summary_table,
    write_trait_imputation_table,
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


def test_summarize_trait_imputation_predicts_missing_taxa_and_validates_holdout() -> (
    None
):
    report = summarize_trait_imputation(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_trait_imputation.tsv"),
        trait="body_mass",
    )

    assert report.model == "brownian"
    assert report.tree_taxon_count == 6
    assert report.observed_taxa == ["A", "C", "D", "E", "F"]
    assert report.observed_taxon_count == 5
    assert len(report.imputation_rows) == 1
    assert report.imputation_rows[0].taxon == "B"
    assert report.imputation_rows[0].missing_reason == "missing_trait_value"
    assert report.imputation_rows[0].predicted_value > 0.0
    assert (
        report.imputation_rows[0].lower_95_confidence_interval
        < report.imputation_rows[0].predicted_value
        < report.imputation_rows[0].upper_95_confidence_interval
    )
    assert report.holdout_validation_status == "performed"
    assert len(report.holdout_rows) == 5
    assert report.holdout_mean_absolute_error is not None
    assert report.holdout_root_mean_squared_error is not None
    assert report.holdout_interval_coverage is not None
    assert 0.0 <= report.holdout_interval_coverage <= 1.0

    top_holdout = report.holdout_rows[0]
    assert top_holdout.rank == 1
    assert top_holdout.absolute_error >= 0.0
    assert (
        top_holdout.lower_95_confidence_interval
        < top_holdout.upper_95_confidence_interval
    )


def test_summarize_trait_imputation_tracks_excluded_taxa() -> None:
    report = summarize_trait_imputation(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_continuous_evolution_missing.tsv"),
        trait="response_growth",
    )

    assert report.observed_taxa == ["A", "D", "E", "F"]
    assert [row.taxon for row in report.imputation_rows] == ["B"]
    assert {row.taxon: row.reason for row in report.excluded_taxa} == {
        "C": "non_numeric_trait_value",
        "G": "absent_from_tree",
    }


def test_trait_imputation_writers_emit_review_ledgers(tmp_path: Path) -> None:
    report = summarize_trait_imputation(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_trait_imputation.tsv"),
        trait="body_mass",
    )
    summary_out = tmp_path / "trait-imputation-summary.tsv"
    imputations_out = tmp_path / "trait-imputations.tsv"
    holdout_out = tmp_path / "trait-imputation-holdout.tsv"
    excluded_out = tmp_path / "trait-imputation-excluded.tsv"

    write_trait_imputation_summary_table(summary_out, report)
    write_trait_imputation_table(imputations_out, report)
    write_trait_imputation_holdout_table(holdout_out, report)
    write_trait_imputation_exclusion_table(excluded_out, report)

    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    imputation_rows = imputations_out.read_text(encoding="utf-8").splitlines()
    holdout_rows = holdout_out.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("trait\ttaxon_column\tmodel")
    assert imputation_rows[0].startswith(
        "taxon\tmissing_reason\tobserved_support_taxon_count"
    )
    assert holdout_rows[0].startswith(
        "taxon\tobserved_value\tpredicted_value\tresidual"
    )
    assert excluded_rows == ["taxon\treason"]
    assert imputation_rows[1].startswith("B\tmissing_trait_value\t5\t")
