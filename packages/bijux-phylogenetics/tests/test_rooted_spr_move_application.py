from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.io.newick import loads_newick
import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    RootedSprMoveApplicationReport,
    RootedSprMoveCandidate,
    apply_rooted_spr_move,
    resolve_rooted_spr_move_candidate,
    rooted_topology_fingerprint,
    summarize_rooted_spr_move_application,
)

FIXTURES = Path(__file__).parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_topology_gateway_exports_rooted_spr_move_application_surface() -> None:
    assert topology_api.RootedSprMoveApplicationReport is RootedSprMoveApplicationReport
    assert topology_api.RootedSprMoveCandidate is RootedSprMoveCandidate
    assert topology_api.resolve_rooted_spr_move_candidate is resolve_rooted_spr_move_candidate
    assert (
        topology_api.summarize_rooted_spr_move_application
        is summarize_rooted_spr_move_application
    )


def test_rooted_spr_move_application_preserves_taxa_and_builds_valid_tree() -> None:
    tree = loads_newick("(((A,C),B),D);")
    candidate, available_move_count = resolve_rooted_spr_move_candidate(tree, 1)

    moved_tree = apply_rooted_spr_move(tree, candidate)

    assert available_move_count == 30
    assert sorted(moved_tree.tip_names) == ["A", "B", "C", "D"]
    assert moved_tree.validation_errors() == []
    assert rooted_topology_fingerprint(moved_tree) != rooted_topology_fingerprint(tree)


def test_rooted_spr_move_application_report_emits_pruned_regraft_and_affected_clades() -> (
    None
):
    report = summarize_rooted_spr_move_application(loads_newick("(((A,C),B),D);"), 1)

    assert report.algorithm == "rooted-spr-move-application"
    assert report.selected_move_index == 1
    assert report.available_move_count == 30
    assert report.pruned_edge_id == report.selected_pruned_clade_id
    assert report.regraft_edge_id == report.selected_regraft_target_branch_id
    assert report.selected_pruned_descendant_taxa == ["A"]
    assert report.selected_regraft_target_descendant_taxa is None
    assert report.moved_topology_changed is True
    assert report.missing_tip_taxa == []
    assert report.unexpected_tip_taxa == []
    assert report.moved_validation_errors == []
    assert report.affected_subtree_report.affected_branch_clade_ids == [
        "A|C",
        "B|C",
        "A|B|C",
        "B|C|D",
    ]
    assert report.affected_subtree_report.unaffected_branch_clade_ids == [
        "A",
        "B",
        "C",
        "D",
    ]
    assert report.affected_clade_ids == ["A", "A|B|C", "A|C", "B|C", "B|C|D"]


def test_rooted_spr_move_application_report_preserves_input_tree_path() -> None:
    report = summarize_rooted_spr_move_application(fixture("example_tree.nwk"), 1)

    assert report.input_tree_path == fixture("example_tree.nwk")
    assert report.input_tree_newick == "((A:0.1,B:0.1):0.2,(C:0.2,D:0.2):0.1);"
    assert report.tip_count == 4


def test_rooted_spr_move_resolution_rejects_out_of_range_indices() -> None:
    with pytest.raises(
        ValueError,
        match="rooted SPR move index must be between 1 and 30",
    ):
        resolve_rooted_spr_move_candidate(loads_newick("(((A,C),B),D);"), 31)
