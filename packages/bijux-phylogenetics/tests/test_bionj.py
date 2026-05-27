from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.distance import (
    build_bionj_tree,
    build_tree_from_imported_distance_matrix,
    load_imported_distance_matrix,
)
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.topology.neighbor_joining import (
    build_neighbor_joining_tree,
)
from bijux_phylogenetics.runtime.errors import InvalidDistanceMatrixError

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


def _lookup_from_imported_matrix(
    path: Path,
) -> tuple[list[str], dict[tuple[str, str], float]]:
    entries = load_imported_distance_matrix(path)
    identifiers = sorted(
        {entry.left_identifier for entry in entries}
        | {entry.right_identifier for entry in entries}
    )
    lookup = {(identifier, identifier): 0.0 for identifier in identifiers}
    for entry in entries:
        lookup[(entry.left_identifier, entry.right_identifier)] = entry.distance
    return identifiers, lookup


def test_build_bionj_tree_matches_noisy_fixture_and_differs_from_nj() -> None:
    identifiers, lookup = _lookup_from_imported_matrix(
        fixture("example_distance_matrix_bionj_noisy.tsv")
    )
    nj_tree = build_neighbor_joining_tree(identifiers, lookup)
    bionj_tree = build_bionj_tree(identifiers, lookup)
    assert dumps_newick(nj_tree) == (
        "((A:1.125,(B:1,C:2)Inner1:5.375)Inner2:5.375,D:4.625,E:-2.625)Inner3;"
    )
    assert dumps_newick(bionj_tree) == (
        "((A:2,(B:1,C:2)Inner1:5.66666666666667)Inner2:4.5,D:6.02,E:-4.02)Inner3;"
    )
    assert dumps_newick(bionj_tree) != dumps_newick(nj_tree)


def test_build_tree_from_imported_distance_matrix_supports_bionj() -> None:
    tree, report = build_tree_from_imported_distance_matrix(
        fixture("example_distance_matrix_bionj_noisy.tsv"),
        method="bionj",
    )
    assert dumps_newick(tree) == (
        "((A:2,(B:1,C:2)Inner1:5.66666666666667)Inner2:4.5,D:6.02,E:-4.02)Inner3;"
    )
    assert report.method == "bionj"
    assert report.method_policy.method == "bionj"


def test_build_tree_from_imported_distance_matrix_rejects_asymmetric_bionj_input() -> (
    None
):
    with pytest.raises(InvalidDistanceMatrixError) as error:
        build_tree_from_imported_distance_matrix(
            fixture("example_distance_matrix_asymmetric.tsv"),
            method="bionj",
        )
    assert error.value.code == "invalid_distance_matrix_error"
    assert "asymmetric directional entries" in str(error.value)
