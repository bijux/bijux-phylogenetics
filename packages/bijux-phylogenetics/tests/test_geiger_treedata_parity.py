from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import align_tree_and_trait_table
from tests.support.geiger_treedata_reference import (
    GEIGER_TREEDATA_REFERENCE_PAYLOADS,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures"
EXAMPLE_TREE = FIXTURES / "trees" / "example_tree.nwk"


def test_align_tree_and_trait_table_matches_governed_geiger_treedata_mismatch_reference() -> (
    None
):
    alignment = align_tree_and_trait_table(
        EXAMPLE_TREE,
        FIXTURES / "metadata" / "example_traits.tsv",
    )
    reference = GEIGER_TREEDATA_REFERENCE_PAYLOADS["example_tree_example_traits"]

    assert alignment.report.aligned_taxa == reference["aligned_taxa"]
    assert alignment.report.dropped_tree_taxa == reference["dropped_tree_taxa"]
    assert alignment.report.dropped_trait_taxa == reference["dropped_trait_taxa"]
    assert [row["taxon"] for row in alignment.rows] == reference["aligned_taxa"]


def test_align_tree_and_trait_table_matches_governed_geiger_treedata_reorder_reference() -> (
    None
):
    alignment = align_tree_and_trait_table(
        EXAMPLE_TREE,
        FIXTURES / "metadata" / "example_traits_comparative_reordered.tsv",
        required_trait_columns=("response",),
    )
    reference = GEIGER_TREEDATA_REFERENCE_PAYLOADS["example_tree_reordered_traits"]

    assert alignment.report.aligned_taxa == reference["aligned_taxa"]
    assert alignment.report.dropped_tree_taxa == reference["dropped_tree_taxa"]
    assert alignment.report.dropped_trait_taxa == reference["dropped_trait_taxa"]
    assert [row["taxon"] for row in alignment.rows] == reference["aligned_taxa"]


def test_align_tree_and_trait_table_matches_governed_geiger_treedata_missing_value_reference() -> (
    None
):
    alignment = align_tree_and_trait_table(
        EXAMPLE_TREE,
        FIXTURES / "metadata" / "example_traits_validate.tsv",
        required_trait_columns=("status",),
    )
    reference = GEIGER_TREEDATA_REFERENCE_PAYLOADS["example_tree_traits_validate"]

    assert alignment.report.aligned_taxa == reference["aligned_taxa"]
    assert alignment.report.dropped_tree_taxa == reference["dropped_tree_taxa"]
    assert alignment.report.dropped_trait_taxa == reference["dropped_trait_taxa"]
    assert [row["taxon"] for row in alignment.rows] == reference["aligned_taxa"]
    assert alignment.report.dropped_missing_value_taxa == []
