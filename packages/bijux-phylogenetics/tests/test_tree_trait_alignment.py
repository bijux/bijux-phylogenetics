from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.datasets.study_inputs import align_tree_and_trait_table
from bijux_phylogenetics.runtime.errors import MetadataJoinError

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


def test_align_tree_and_trait_table_drops_nonoverlap_in_tree_order() -> None:
    alignment = align_tree_and_trait_table(
        fixture("example_tree.nwk"),
        fixture("example_traits.tsv"),
    )

    assert alignment.tree.tip_names == ["A", "B", "C"]
    assert [row["taxon"] for row in alignment.rows] == ["A", "B", "C"]
    assert alignment.report.aligned_taxa == ["A", "B", "C"]
    assert alignment.report.dropped_tree_taxa == ["D"]
    assert alignment.report.dropped_trait_taxa == ["E"]
    assert alignment.report.dropped_missing_value_taxa == []
    assert alignment.report.missing_value_policy == "retain-overlapping-missing-values"


def test_align_tree_and_trait_table_reorders_rows_to_tree_tip_order() -> None:
    alignment = align_tree_and_trait_table(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative_reordered.tsv"),
        required_trait_columns=("response",),
    )

    assert alignment.tree.tip_names == ["A", "B", "C", "D"]
    assert [row["taxon"] for row in alignment.rows] == ["A", "B", "C", "D"]
    assert alignment.report.dropped_tree_taxa == []
    assert alignment.report.dropped_trait_taxa == []


def test_align_tree_and_trait_table_can_prune_requested_missing_values() -> None:
    alignment = align_tree_and_trait_table(
        fixture("example_tree.nwk"),
        fixture("example_traits_validate.tsv"),
        required_trait_columns=("status",),
        drop_missing_for_columns=("status",),
    )

    assert alignment.tree.tip_names == ["A", "B", "D"]
    assert [row["taxon"] for row in alignment.rows] == ["A", "B", "D"]
    assert alignment.report.dropped_missing_value_taxa == ["C"]
    assert [
        (item.taxon, item.trait) for item in alignment.report.missing_value_calls
    ] == [("C", "status")]
    assert (
        alignment.report.missing_value_policy
        == "drop-overlapping-missing-values-for-requested-traits"
    )


def test_align_tree_and_trait_table_rejects_duplicate_taxa() -> None:
    with pytest.raises(MetadataJoinError, match="duplicate taxon 'A'"):
        align_tree_and_trait_table(
            fixture("example_tree.nwk"),
            fixture("example_metadata_duplicate.tsv"),
        )
