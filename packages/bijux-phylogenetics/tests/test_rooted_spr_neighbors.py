from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.io.newick import loads_newick
import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    RootedSprMoveCandidate,
    RootedSprNeighborRow,
    RootedSprNeighborhoodReport,
    enumerate_rooted_spr_neighbors,
    validate_rooted_spr_tree,
)

FIXTURES = Path(__file__).parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_topology_gateway_exports_rooted_spr_neighbor_contracts() -> None:
    assert topology_api.RootedSprMoveCandidate is RootedSprMoveCandidate
    assert topology_api.RootedSprNeighborRow is RootedSprNeighborRow
    assert topology_api.RootedSprNeighborhoodReport is RootedSprNeighborhoodReport
    assert topology_api.enumerate_rooted_spr_neighbors is enumerate_rooted_spr_neighbors
    assert topology_api.validate_rooted_spr_tree is validate_rooted_spr_tree


def test_rooted_spr_neighbor_report_collapses_identity_and_duplicate_moves() -> None:
    report = enumerate_rooted_spr_neighbors(fixture("example_tree.nwk"))

    assert report.algorithm == "rooted-spr-neighbor-enumeration"
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert report.generated_move_candidate_count == 32
    assert report.identity_move_candidate_count == 8
    assert report.self_regraft_candidate_count == 0
    assert report.generated_neighbor_count == 12
    assert report.unique_neighbor_topology_count == 12
    assert len(report.duplicate_move_neighbor_topologies) == 4
    assert report.missing_tip_taxa == []
    assert report.unexpected_tip_taxa == []
    assert report.input_validation_errors == []


def test_rooted_spr_neighbor_report_matches_known_five_taxon_case() -> None:
    report = enumerate_rooted_spr_neighbors(loads_newick("((((A,D),B),C),E);"))

    assert report.generated_move_candidate_count == 52
    assert report.identity_move_candidate_count == 10
    assert report.self_regraft_candidate_count == 0
    assert report.generated_neighbor_count == 24
    assert report.unique_neighbor_topology_count == 24
    assert len(report.duplicate_move_neighbor_topologies) == 10


def test_rooted_spr_neighbors_preserve_taxa_and_exclude_input_topology() -> None:
    report = enumerate_rooted_spr_neighbors(loads_newick("(((A,C),B),D);"))

    assert all(sorted(row.tip_order) == ["A", "B", "C", "D"] for row in report.neighbor_rows)
    assert all(row.validation_errors == [] for row in report.neighbor_rows)
    assert all(row.supporting_move_count >= 1 for row in report.neighbor_rows)
    assert all(row.neighbor_tree_newick != "(((A,C),B),D);" for row in report.neighbor_rows)


def test_rooted_spr_validation_accepts_binary_root_representation_without_rooted_flag() -> None:
    validate_rooted_spr_tree(loads_newick("(((A,C),B),D);"))


def test_rooted_spr_validation_rejects_nonbinary_rooted_representation() -> None:
    with pytest.raises(ValueError, match="rooted SPR enumeration requires a binary root"):
        validate_rooted_spr_tree(loads_newick("(A,B,C,D);"))
