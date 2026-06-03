from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import check_tree_and_trait_taxon_names
from tests.support.geiger_name_check_reference import (
    GEIGER_NAME_CHECK_REFERENCE_PAYLOADS,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def test_check_tree_and_trait_taxon_names_matches_governed_geiger_mismatch_reference() -> (
    None
):
    report = check_tree_and_trait_taxon_names(
        FIXTURES / "trees" / "example_tree.nwk",
        FIXTURES / "metadata" / "example_traits.tsv",
    )
    reference = GEIGER_NAME_CHECK_REFERENCE_PAYLOADS["example_tree_example_traits"]

    assert report.reference_outcome == reference["reference_outcome"]
    assert report.tree_not_data == reference["tree_not_data"]
    assert report.data_not_tree == reference["data_not_tree"]


def test_check_tree_and_trait_taxon_names_matches_governed_geiger_ok_reference() -> (
    None
):
    report = check_tree_and_trait_taxon_names(
        FIXTURES / "trees" / "example_tree.nwk",
        FIXTURES / "metadata" / "example_traits_comparative_reordered.tsv",
    )
    reference = GEIGER_NAME_CHECK_REFERENCE_PAYLOADS["example_tree_reordered_traits"]

    assert report.reference_outcome == reference["reference_outcome"]
    assert report.tree_not_data == reference["tree_not_data"]
    assert report.data_not_tree == reference["data_not_tree"]


def test_check_tree_and_trait_taxon_names_matches_governed_geiger_case_policy(
    tmp_path: Path,
) -> None:
    table_path = tmp_path / "case-sensitive-traits.tsv"
    table_path.write_text(
        "taxon\tvalue\na\t1.0\nB\t2.0\nC\t3.0\nD\t4.0\n",
        encoding="utf-8",
    )

    report = check_tree_and_trait_taxon_names(
        FIXTURES / "trees" / "example_tree.nwk",
        table_path,
    )
    reference = GEIGER_NAME_CHECK_REFERENCE_PAYLOADS["case_sensitive_probe"]

    assert report.reference_outcome == reference["reference_outcome"]
    assert report.tree_not_data == reference["tree_not_data"]
    assert report.data_not_tree == reference["data_not_tree"]
