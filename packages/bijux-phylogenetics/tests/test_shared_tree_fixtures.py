from __future__ import annotations

import pytest

from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_tree_fixture,
    list_shared_tree_fixtures,
)
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError, TreeParseError


def test_shared_tree_fixture_catalog_covers_required_tree_shapes() -> None:
    fixtures = list_shared_tree_fixtures()
    feature_tags = {tag for fixture in fixtures for tag in fixture.feature_tags}

    assert {
        "balanced-binary",
        "pectinate",
        "star",
        "polytomy",
        "rootable",
        "rooted",
        "unrooted",
        "ultrametric",
        "non-ultrametric",
        "zero-branch",
        "long-branch-outlier",
        "internal-node-labels",
        "branch-support-labels",
        "quoted-taxon-labels",
        "near-ultrametric",
        "malformed-newick",
    } <= feature_tags
    assert (
        max(
            fixture.tip_count or 0
            for fixture in fixtures
            if fixture.parse_expectation == "parseable"
        )
        > 100
    )


def test_shared_tree_fixture_catalog_preserves_durable_fixture_lookup() -> None:
    fixture = get_shared_tree_fixture("quoted_taxon_labels")

    assert fixture.relative_path == "trees/example_tree_labels.nwk"
    assert "quoted-taxon-labels" in fixture.feature_tags
    assert fixture.path.is_file()


@pytest.mark.parametrize(
    ("fixture_id", "expected_tip_count"),
    [
        ("two_tip_rooted_ultrametric", 2),
        ("balanced_rooted_ultrametric", 4),
        ("balanced_rooted_six_taxon", 6),
        ("near_ultrametric_branch_jitter", 4),
        ("cross_pairing_rooted_four_taxon", 4),
        ("pectinate_rooted_non_ultrametric", 4),
        ("star_unrooted_polytomy", 5),
        ("larger_binary_tree", 8),
        ("phytools_ultrametric_twenty_four_taxa", 24),
        ("phytools_ultrametric_one_hundred_twenty_eight_taxa", 128),
        ("phytools_non_ultrametric_twenty_four_taxa", 24),
        ("phytools_branch_edge_twenty_four_taxa", 24),
        ("quoted_taxon_labels", 3),
        ("outgroup_rootable_unrooted", 4),
        ("outgroup_rooted_on_d", 4),
    ],
)
def test_shared_tree_fixture_catalog_parseable_cases_load_as_expected(
    fixture_id: str, expected_tip_count: int
) -> None:
    fixture = get_shared_tree_fixture(fixture_id)

    tree = load_tree(fixture.path)

    assert len(tree.tip_names) == expected_tip_count


@pytest.mark.parametrize(
    "fixture_id",
    [
        "malformed_unbalanced_parentheses",
        "malformed_extra_closing_parenthesis",
    ],
)
def test_shared_tree_fixture_catalog_malformed_cases_raise_tree_parse_error(
    fixture_id: str,
) -> None:
    fixture = get_shared_tree_fixture(fixture_id)

    with pytest.raises(TreeParseError):
        load_tree(fixture.path)
    with pytest.raises(TreeParseError):
        validate_tree_path(fixture.path)


def test_shared_tree_fixture_catalog_negative_branch_fixture_fails_validation() -> None:
    fixture = get_shared_tree_fixture("negative_branch_length")

    tree = load_tree(fixture.path)

    assert len(tree.tip_names) == 3
    with pytest.raises(InvalidBranchLengthError):
        validate_tree_path(fixture.path)


def test_shared_tree_fixture_catalog_warning_rich_cases_report_warnings() -> None:
    fixture_ids = [
        "star_unrooted_polytomy",
        "rooted_polytomy",
        "unrooted_branch_length_tree",
        "zero_branch_lengths",
        "quoted_taxon_labels",
    ]

    for fixture_id in fixture_ids:
        fixture = get_shared_tree_fixture(fixture_id)
        report = validate_tree_path(fixture.path)
        assert report.validity_decision == "valid_with_warnings"
