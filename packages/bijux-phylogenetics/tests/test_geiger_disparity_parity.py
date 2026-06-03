from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative.disparity import (
    summarize_continuous_clade_disparity,
)
from tests.support.geiger_disparity_reference import (
    GEIGER_DISPARITY_REFERENCE_PAYLOADS,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def test_summarize_continuous_clade_disparity_matches_governed_geiger_univariate_reference() -> (
    None
):
    report = summarize_continuous_clade_disparity(
        FIXTURES / "trees" / "example_tree.nwk",
        FIXTURES / "metadata" / "example_traits_comparative.tsv",
        trait_columns=["response"],
    )
    reference = GEIGER_DISPARITY_REFERENCE_PAYLOADS["example_tree_response_univariate"]

    assert [row.ape_node_id for row in report.clade_rows] == reference[
        "ape_node_id_order"
    ]
    assert [row.descendant_taxa for row in report.clade_rows] == reference[
        "descendant_taxa"
    ]
    for observed, expected in zip(
        [row.disparity for row in report.clade_rows],
        reference["clade_disparity"],
        strict=True,
    ):
        assert math.isclose(observed, expected, rel_tol=1e-12, abs_tol=1e-12)


def test_summarize_continuous_clade_disparity_matches_governed_geiger_multivariate_reference() -> (
    None
):
    report = summarize_continuous_clade_disparity(
        FIXTURES / "trees" / "example_tree.nwk",
        FIXTURES / "metadata" / "example_traits_comparative.tsv",
        trait_columns=["response", "predictor_one"],
    )
    reference = GEIGER_DISPARITY_REFERENCE_PAYLOADS[
        "example_tree_response_predictor_one_multivariate"
    ]

    assert [row.ape_node_id for row in report.clade_rows] == reference[
        "ape_node_id_order"
    ]
    assert [row.descendant_taxa for row in report.clade_rows] == reference[
        "descendant_taxa"
    ]
    for observed, expected in zip(
        [row.disparity for row in report.clade_rows],
        reference["clade_disparity"],
        strict=True,
    ):
        assert math.isclose(observed, expected, rel_tol=1e-12, abs_tol=1e-12)
