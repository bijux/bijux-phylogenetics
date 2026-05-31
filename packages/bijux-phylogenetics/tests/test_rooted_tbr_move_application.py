from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.topology.rooted_tbr import (
    RootedTbrMoveCandidate,
    apply_rooted_tbr_move,
    resolve_rooted_tbr_move_candidate,
    summarize_rooted_tbr_move_application,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_rooted_tbr_move_application_preserves_taxa_and_builds_valid_tree() -> None:
    tree = loads_newick("(((A,C),B),D);")
    candidate, available_move_count = resolve_rooted_tbr_move_candidate(tree, 1)

    moved_tree = apply_rooted_tbr_move(tree, candidate)

    assert isinstance(candidate, RootedTbrMoveCandidate)
    assert available_move_count == 16
    assert sorted(moved_tree.tip_names) == ["A", "B", "C", "D"]
    assert moved_tree.validation_errors() == []
    assert moved_tree.to_newick() == "((A,C),(B,D));"


def test_rooted_tbr_move_application_report_emits_cut_sizes_and_reconnection_edges() -> (
    None
):
    report = summarize_rooted_tbr_move_application(loads_newick("(((A,C),B),D);"), 1)

    assert report.algorithm == "rooted-tbr-move-application"
    assert report.selected_move_index == 1
    assert report.available_move_count == 16
    assert report.selected_cut_edge_id == "A|C"
    assert report.selected_cut_descendant_taxa == ["A", "C"]
    assert report.left_component_tip_count == 2
    assert report.right_component_tip_count == 2
    assert report.selected_left_attachment_branch_id == "interface"
    assert report.selected_left_attachment_descendant_taxa == ["A", "C"]
    assert report.selected_right_attachment_branch_id == "interface"
    assert report.selected_right_attachment_descendant_taxa == ["B", "D"]
    assert report.moved_tree_newick == "((A,C),(B,D));"
    assert report.moved_topology_changed is True
    assert report.reverse_move_available is False
    assert report.reverse_available_move_count == 0
    assert report.missing_tip_taxa == []
    assert report.unexpected_tip_taxa == []
    assert report.moved_validation_errors == []
    assert report.affected_subtree_report.affected_branch_clade_ids == [
        "B|D",
        "A|B|C",
    ]
    assert report.affected_subtree_report.unaffected_branch_clade_ids == [
        "A",
        "B",
        "C",
        "D",
        "A|C",
    ]


def test_rooted_tbr_move_application_reports_reverse_move_when_available() -> None:
    report = summarize_rooted_tbr_move_application(
        fixture("parsimony", "spr_search_start_tree_5_taxa.nwk"),
        39,
    )

    assert report.selected_cut_edge_id == "A|B|C|D"
    assert report.left_component_tip_count == 4
    assert report.right_component_tip_count == 1
    assert report.selected_left_attachment_branch_id == "A"
    assert report.selected_right_attachment_branch_id == "interface"
    assert report.moved_tree_newick == "((A,((B,C),D)),E);"
    assert report.reverse_move_available is True
    assert report.reverse_available_move_count == 2
    assert report.reverse_cut_edge_id == "A|B|C|D"
    assert report.reverse_left_attachment_branch_id == "C"
    assert report.reverse_right_attachment_branch_id == "interface"


def test_rooted_tbr_move_application_report_preserves_input_tree_path() -> None:
    report = summarize_rooted_tbr_move_application(
        fixture("parsimony", "spr_search_start_tree_5_taxa.nwk"),
        1,
    )

    assert report.input_tree_path == fixture("parsimony", "spr_search_start_tree_5_taxa.nwk")
    assert report.input_tree_newick == "((((A,D),B),C),E);"
    assert report.tip_count == 5


def test_rooted_tbr_move_resolution_rejects_out_of_range_indices() -> None:
    with pytest.raises(
        ValueError,
        match="rooted TBR move index must be between 1 and 16",
    ):
        resolve_rooted_tbr_move_candidate(loads_newick("(((A,C),B),D);"), 17)
