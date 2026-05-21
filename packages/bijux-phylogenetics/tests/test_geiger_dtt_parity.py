from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative.disparity import (
    summarize_continuous_clade_disparity,
    summarize_disparity_through_time,
)
from tests.support.geiger_dtt_reference import GEIGER_DTT_REFERENCE_PAYLOADS

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def test_summarize_continuous_clade_disparity_matches_governed_geiger_univariate_reference() -> (
    None
):
    report = summarize_continuous_clade_disparity(
        FIXTURES / "trees" / "example_tree_phytools_ultrametric_twenty_four_taxa.nwk",
        FIXTURES
        / "metadata"
        / "example_traits_geiger_continuous_model_panel_twenty_four_taxa.tsv",
        trait_columns=["ou_truth"],
    )
    reference = GEIGER_DTT_REFERENCE_PAYLOADS["ou_truth_univariate_twenty_four_taxa"]

    assert report.trait_columns == reference["trait_columns"]
    for observed, expected in zip(
        [row.disparity for row in report.clade_rows],
        reference["clade_disparity"],
        strict=True,
    ):
        assert math.isclose(observed, expected, rel_tol=1e-12, abs_tol=1e-12)


def test_summarize_disparity_through_time_matches_governed_geiger_univariate_reference() -> (
    None
):
    report = summarize_disparity_through_time(
        FIXTURES / "trees" / "example_tree_phytools_ultrametric_twenty_four_taxa.nwk",
        FIXTURES
        / "metadata"
        / "example_traits_geiger_continuous_model_panel_twenty_four_taxa.tsv",
        trait_columns=["ou_truth"],
    )
    reference = GEIGER_DTT_REFERENCE_PAYLOADS["ou_truth_univariate_twenty_four_taxa"]

    for row, expected_time, expected_disparity in zip(
        report.curve_rows,
        reference["times"],
        reference["relative_disparity"],
        strict=True,
    ):
        assert math.isclose(
            row.relative_time, expected_time, rel_tol=1e-12, abs_tol=1e-12
        )
        assert math.isclose(
            row.relative_disparity,
            expected_disparity,
            rel_tol=1e-12,
            abs_tol=1e-12,
        )


def test_summarize_disparity_through_time_matches_governed_geiger_multivariate_reference() -> (
    None
):
    report = summarize_disparity_through_time(
        FIXTURES / "trees" / "example_tree_phytools_ultrametric_twenty_four_taxa.nwk",
        FIXTURES
        / "metadata"
        / "example_traits_geiger_continuous_model_panel_twenty_four_taxa.tsv",
        trait_columns=["ou_truth", "early_burst_truth"],
    )
    reference = GEIGER_DTT_REFERENCE_PAYLOADS[
        "ou_truth_and_early_burst_truth_multivariate_twenty_four_taxa"
    ]

    for observed, expected in zip(
        [row.disparity for row in report.clade_rows],
        reference["clade_disparity"],
        strict=True,
    ):
        assert math.isclose(observed, expected, rel_tol=1e-12, abs_tol=1e-12)
    for row, expected_time, expected_disparity in zip(
        report.curve_rows,
        reference["times"],
        reference["relative_disparity"],
        strict=True,
    ):
        assert math.isclose(
            row.relative_time, expected_time, rel_tol=1e-12, abs_tol=1e-12
        )
        assert math.isclose(
            row.relative_disparity,
            expected_disparity,
            rel_tol=1e-12,
            abs_tol=1e-12,
        )
