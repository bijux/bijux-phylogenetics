from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative.traits.outliers import (
    summarize_trait_outliers,
    write_trait_outlier_exclusion_table,
    write_trait_outlier_summary_table,
    write_trait_outlier_taxon_table,
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


def test_summarize_trait_outliers_ranks_focal_phylogenetic_exception() -> None:
    report = summarize_trait_outliers(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_trait_outliers.tsv"),
        trait="body_mass",
    )

    assert report.tree_taxon_count == 6
    assert report.analyzed_taxa == ["A", "B", "C", "D", "E", "F"]
    assert report.analyzed_taxon_count == 6
    assert report.top_outlier_taxon == "F"
    assert report.outlier_taxa == ["F"]
    assert report.top_abs_standardized_residual is not None
    assert report.top_abs_standardized_residual >= 2.0

    top_row = report.taxon_rows[0]
    assert top_row.taxon == "F"
    assert top_row.rank == 1
    assert top_row.outlier is True
    assert top_row.context_clade_id == "E|F"
    assert top_row.sibling_context_id == "E"
    assert top_row.context_taxa == ["E", "F"]
    assert top_row.sibling_taxa == ["E"]
    assert math.isclose(top_row.context_mean or 0.0, 14.0)
    assert math.isclose(top_row.sibling_mean or 0.0, 3.0)
    assert math.isclose(top_row.context_mean_shift or 0.0, 11.0)
    assert top_row.conditional_expected_value < top_row.observed_value


def test_summarize_trait_outliers_tracks_excluded_taxa() -> None:
    report = summarize_trait_outliers(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_continuous_evolution_missing.tsv"),
        trait="response_growth",
    )

    assert report.analyzed_taxa == ["A", "D", "E", "F"]
    assert {row.taxon: row.reason for row in report.excluded_taxa} == {
        "B": "missing_trait_value",
        "C": "non_numeric_trait_value",
        "G": "absent_from_tree",
    }


def test_trait_outlier_writers_emit_review_ledgers(tmp_path: Path) -> None:
    report = summarize_trait_outliers(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_trait_outliers.tsv"),
        trait="body_mass",
    )
    summary_out = tmp_path / "trait-outlier-summary.tsv"
    outliers_out = tmp_path / "trait-outliers.tsv"
    excluded_out = tmp_path / "trait-outlier-excluded.tsv"

    write_trait_outlier_summary_table(summary_out, report)
    write_trait_outlier_taxon_table(outliers_out, report)
    write_trait_outlier_exclusion_table(excluded_out, report)

    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    outlier_rows = outliers_out.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("trait\ttaxon_column\ttree_taxon_count")
    assert outlier_rows[0].startswith(
        "taxon\tobserved_value\tconditional_expected_value\tresidual"
    )
    assert outlier_rows[1].startswith("F\t25\t")
    assert excluded_rows == ["taxon\treason"]
