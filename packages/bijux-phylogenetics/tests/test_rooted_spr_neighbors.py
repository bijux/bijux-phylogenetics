from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.io.newick import loads_newick
import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    RootedSprMoveCandidate,
    RootedSprNeighborRow,
    RootedSprNeighborhoodReport,
    enumerate_rooted_spr_neighbors,
    rooted_topology_fingerprint,
    validate_rooted_spr_tree,
    write_rooted_spr_artifacts,
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
    assert topology_api.write_rooted_spr_artifacts is write_rooted_spr_artifacts


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
    tree = loads_newick("(((A,C),B),D);")
    report = enumerate_rooted_spr_neighbors(tree)
    input_topology_fingerprint = rooted_topology_fingerprint(tree)

    assert all(sorted(row.tip_order) == ["A", "B", "C", "D"] for row in report.neighbor_rows)
    assert all(row.validation_errors == [] for row in report.neighbor_rows)
    assert all(row.supporting_move_count >= 1 for row in report.neighbor_rows)
    assert all(
        row.neighbor_topology_fingerprint != input_topology_fingerprint
        for row in report.neighbor_rows
    )


def test_rooted_spr_validation_accepts_binary_root_representation_without_rooted_flag() -> None:
    validate_rooted_spr_tree(loads_newick("(((A,C),B),D);"))


def test_rooted_spr_validation_rejects_nonbinary_rooted_representation() -> None:
    with pytest.raises(ValueError, match="rooted SPR enumeration requires a binary root"):
        validate_rooted_spr_tree(loads_newick("(A,B,C,D);"))


def test_rooted_spr_report_preserves_input_tree_path_for_file_inputs() -> None:
    report = enumerate_rooted_spr_neighbors(fixture("example_tree.nwk"))

    assert report.input_tree_path == fixture("example_tree.nwk")
    assert report.input_tree_newick == "((A:0.1,B:0.1):0.2,(C:0.2,D:0.2):0.1);"
    assert report.rooted is False


def test_rooted_spr_neighbor_enumeration_rejects_internal_polytomies() -> None:
    with pytest.raises(
        ValueError,
        match="rooted SPR enumeration requires a strictly bifurcating tree",
    ):
        enumerate_rooted_spr_neighbors(loads_newick("((A,B,C),D);"))


def test_write_rooted_spr_artifacts_materializes_governed_outputs(tmp_path: Path) -> None:
    report = enumerate_rooted_spr_neighbors(fixture("example_tree.nwk"))

    outputs = write_rooted_spr_artifacts(tmp_path / "rooted-spr-neighbors", report)

    assert set(outputs) == {
        "input_tree_path",
        "neighbors_path",
        "summary_path",
        "run_json_path",
    }
    assert outputs["neighbors_path"].read_text(encoding="utf-8").startswith(
        "neighbor_index\trepresentative_pruned_node_id\trepresentative_pruned_clade_id\trepresentative_pruned_descendant_taxa\trepresentative_regraft_target_branch_id\trepresentative_regraft_target_descendant_taxa\tsupporting_move_count\tneighbor_topology_fingerprint\ttip_order\tvalidation_errors\tneighbor_tree_newick\n"
    )
    assert outputs["summary_path"].read_text(encoding="utf-8").startswith(
        "neighborhood_family\talgorithm\tcandidate_count\tvalid_count\tduplicate_count\tskipped_count\tskipped_reason\tbudget_reason\n"
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "rooted-spr-neighbor-enumeration"
    assert payload["input_tree_path"].endswith("example_tree.nwk")
    assert payload["generated_move_candidate_count"] == 32
    assert payload["identity_move_candidate_count"] == 8
    assert payload["generated_neighbor_count"] == 12
