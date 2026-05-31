from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.io.newick import loads_newick
import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    RootedNniMoveCandidate,
    enumerate_rooted_nni_neighbors,
    expected_rooted_nni_neighbor_count,
    iter_rooted_nni_move_candidates,
    validate_rooted_nni_tree,
    write_rooted_nni_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_topology_gateway_exports_rooted_nni_neighbor_surface() -> None:
    assert topology_api.RootedNniMoveCandidate is RootedNniMoveCandidate
    assert topology_api.iter_rooted_nni_move_candidates is iter_rooted_nni_move_candidates
    assert topology_api.enumerate_rooted_nni_neighbors is enumerate_rooted_nni_neighbors
    assert (
        topology_api.expected_rooted_nni_neighbor_count
        is expected_rooted_nni_neighbor_count
    )
    assert topology_api.validate_rooted_nni_tree is validate_rooted_nni_tree
    assert topology_api.write_rooted_nni_artifacts is write_rooted_nni_artifacts


def test_rooted_nni_neighbor_count_matches_formula_on_four_taxon_tree() -> None:
    report = enumerate_rooted_nni_neighbors(loads_newick("(((A,C),B),D);"))

    assert report.algorithm == "rooted-nni-neighbor-enumeration"
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert report.expected_neighbor_count == 4
    assert report.generated_neighbor_count == 4
    assert report.unique_neighbor_topology_count == 4
    assert report.duplicate_neighbor_topologies == []
    assert report.missing_tip_taxa == []
    assert report.unexpected_tip_taxa == []
    assert report.input_validation_errors == []
    assert [row.neighbor_index for row in report.neighbor_rows] == [1, 2, 3, 4]


def test_rooted_nni_neighbor_count_matches_formula_on_five_taxon_tree() -> None:
    report = enumerate_rooted_nni_neighbors(loads_newick("((((A,B),C),D),E);"))

    assert report.tip_count == 5
    assert report.internal_node_count == 4
    assert expected_rooted_nni_neighbor_count(loads_newick("((((A,B),C),D),E);")) == 6
    assert report.expected_neighbor_count == 6
    assert report.generated_neighbor_count == 6


def test_rooted_nni_neighbors_preserve_taxa_and_generate_unique_valid_topologies() -> None:
    report = enumerate_rooted_nni_neighbors(loads_newick("(((A,C),B),D);"))

    assert len({row.neighbor_topology_fingerprint for row in report.neighbor_rows}) == 4
    assert all(sorted(row.tip_order) == ["A", "B", "C", "D"] for row in report.neighbor_rows)
    assert all(row.validation_errors == [] for row in report.neighbor_rows)


def test_rooted_nni_report_preserves_input_tree_path_for_file_inputs() -> None:
    report = enumerate_rooted_nni_neighbors(fixture("example_tree.nwk"))

    assert report.input_tree_path == fixture("example_tree.nwk")
    assert report.input_tree_newick == "((A:0.1,B:0.1):0.2,(C:0.2,D:0.2):0.1);"
    assert report.rooted is False


def test_rooted_nni_validation_accepts_binary_root_representation_without_rooted_flag() -> None:
    tree = loads_newick("(((A,C),B),D);")

    validate_rooted_nni_tree(tree)


def test_rooted_nni_validation_rejects_nonbinary_rooted_representation() -> None:
    with pytest.raises(ValueError, match="rooted NNI enumeration requires a binary root"):
        validate_rooted_nni_tree(loads_newick("(A,B,C,D);"))


def test_rooted_nni_neighbor_count_rejects_internal_polytomies() -> None:
    with pytest.raises(
        ValueError,
        match="rooted NNI enumeration requires a strictly bifurcating tree",
    ):
        expected_rooted_nni_neighbor_count(loads_newick("((A,B,C),D);"))


def test_write_rooted_nni_artifacts_materializes_governed_outputs(tmp_path: Path) -> None:
    report = enumerate_rooted_nni_neighbors(fixture("example_tree.nwk"))

    outputs = write_rooted_nni_artifacts(tmp_path / "rooted-nni-neighbors", report)

    assert set(outputs) == {
        "input_tree_path",
        "neighbors_path",
        "summary_path",
        "run_json_path",
    }
    assert outputs["neighbors_path"].read_text(encoding="utf-8").startswith(
        "neighbor_index\tparent_node_id\tchild_node_id\tsibling_node_id\texchanged_child_node_id\tpivot_branch_id\tsibling_clade_id\texchanged_clade_id\tneighbor_topology_fingerprint\ttip_order\tvalidation_errors\tneighbor_tree_newick\n"
    )
    assert outputs["summary_path"].read_text(encoding="utf-8").startswith(
        "neighborhood_family\talgorithm\tcandidate_count\tvalid_count\tduplicate_count\tskipped_count\tskipped_reason\tbudget_reason\n"
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "rooted-nni-neighbor-enumeration"
    assert payload["input_tree_path"].endswith("example_tree.nwk")
    assert payload["expected_neighbor_count"] == 4
    assert payload["generated_neighbor_count"] == 4
    assert len(payload["neighbor_rows"]) == 4
