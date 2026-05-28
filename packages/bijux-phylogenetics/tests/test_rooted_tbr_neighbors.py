from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.topology import (
    enumerate_rooted_spr_neighbors,
    enumerate_rooted_tbr_neighbors,
    RootedTbrNeighborRow,
    RootedTbrNeighborhoodReport,
    validate_rooted_tbr_tree,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_topology_gateway_exports_rooted_tbr_neighbor_contracts() -> None:
    assert topology_api.RootedTbrNeighborRow is RootedTbrNeighborRow
    assert topology_api.RootedTbrNeighborhoodReport is RootedTbrNeighborhoodReport
    assert topology_api.enumerate_rooted_tbr_neighbors is enumerate_rooted_tbr_neighbors
    assert topology_api.validate_rooted_tbr_tree is validate_rooted_tbr_tree


def test_rooted_tbr_neighbor_report_matches_known_five_taxon_case() -> None:
    report = enumerate_rooted_tbr_neighbors(
        fixture("parsimony", "spr_search_start_tree_5_taxa.nwk")
    )

    assert report.algorithm == "rooted-tbr-neighbor-enumeration"
    assert report.tip_count == 5
    assert report.internal_node_count == 4
    assert report.generated_cut_edge_count == 3
    assert report.generated_reconnection_count == 52
    assert report.identity_reconnection_count == 6
    assert report.generated_neighbor_count == 10
    assert report.unique_neighbor_topology_count == 10
    assert len(report.duplicate_reconnection_neighbor_topologies) == 10
    assert report.missing_tip_taxa == []
    assert report.unexpected_tip_taxa == []
    assert report.input_validation_errors == []


def test_rooted_tbr_neighbor_surface_is_not_just_spr_renamed() -> None:
    tbr_report = enumerate_rooted_tbr_neighbors(loads_newick("(((A,C),B),D);"))
    spr_report = enumerate_rooted_spr_neighbors(loads_newick("(((A,C),B),D);"))

    assert tbr_report.generated_neighbor_count == 3
    assert spr_report.generated_neighbor_count == 10
    assert {
        row.neighbor_topology_fingerprint for row in tbr_report.neighbor_rows
    } != {
        row.neighbor_topology_fingerprint for row in spr_report.neighbor_rows
    }


def test_rooted_tbr_neighbors_preserve_taxa_and_generate_valid_topologies() -> None:
    report = enumerate_rooted_tbr_neighbors(
        fixture("parsimony", "spr_search_start_tree_5_taxa.nwk")
    )

    assert all(sorted(row.tip_order) == ["A", "B", "C", "D", "E"] for row in report.neighbor_rows)
    assert all(row.validation_errors == [] for row in report.neighbor_rows)
    assert all(row.supporting_reconnection_count >= 1 for row in report.neighbor_rows)


def test_rooted_tbr_validation_accepts_binary_root_representation_without_rooted_flag() -> None:
    validate_rooted_tbr_tree(loads_newick("(((A,C),B),D);"))


def test_rooted_tbr_validation_rejects_nonbinary_rooted_representation() -> None:
    with pytest.raises(ValueError, match="rooted TBR enumeration requires a binary root"):
        validate_rooted_tbr_tree(loads_newick("(A,B,C,D);"))
