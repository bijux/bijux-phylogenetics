from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.newick import loads_newick
import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    RootedNniMoveCandidate,
    RootedSprMoveCandidate,
    RootedTbrMoveCandidate,
    TopologyMoveValidityReport,
    resolve_rooted_nni_move_candidate,
    resolve_rooted_spr_move_candidate,
    resolve_rooted_tbr_move_candidate,
    summarize_rooted_nni_move_validity,
    summarize_rooted_spr_move_validity,
    summarize_rooted_tbr_move_validity,
)

FIXTURES = Path(__file__).parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_topology_gateway_exports_move_validity_surface() -> None:
    assert topology_api.TopologyMoveValidityReport is TopologyMoveValidityReport
    assert topology_api.summarize_rooted_nni_move_validity is summarize_rooted_nni_move_validity
    assert topology_api.summarize_rooted_spr_move_validity is summarize_rooted_spr_move_validity
    assert topology_api.summarize_rooted_tbr_move_validity is summarize_rooted_tbr_move_validity


def test_rooted_nni_move_validity_accepts_available_candidate() -> None:
    tree = loads_newick("(((A,C),B),D);")
    candidate, available_move_count = resolve_rooted_nni_move_candidate(tree, 1)

    report = summarize_rooted_nni_move_validity(tree, candidate)

    assert report.algorithm == "topology-move-validity"
    assert report.move_family == "rooted-nni"
    assert report.available_move_count == available_move_count
    assert report.validity_decision == "accepted"
    assert report.rejection_code is None
    assert report.rejection_reason is None
    assert report.input_validation_errors == []
    assert report.evidence == {"candidate_available": True}


def test_rooted_spr_move_validity_rejects_self_regraft_candidates() -> None:
    tree = loads_newick("(((A,C),B),D);")
    resolved_candidate, available_move_count = resolve_rooted_spr_move_candidate(tree, 1)
    candidate = RootedSprMoveCandidate(
        pruned_node_id=resolved_candidate.pruned_node_id,
        pruned_clade_id=resolved_candidate.pruned_clade_id,
        pruned_descendant_taxa=resolved_candidate.pruned_descendant_taxa,
        regraft_target_branch_id=resolved_candidate.pruned_clade_id,
        regraft_target_descendant_taxa=resolved_candidate.pruned_descendant_taxa,
    )

    report = summarize_rooted_spr_move_validity(tree, candidate)

    assert report.move_family == "rooted-spr"
    assert report.available_move_count == available_move_count
    assert report.validity_decision == "rejected"
    assert report.rejection_code == "topology_move_self_regraft"
    assert "self-regraft" in report.rejection_reason
    assert report.evidence == {
        "pruned_clade_id": resolved_candidate.pruned_clade_id,
        "regraft_target_branch_id": resolved_candidate.pruned_clade_id,
    }


def test_rooted_nni_move_validity_rejects_too_small_tree() -> None:
    candidate = RootedNniMoveCandidate(
        parent_node_id="missing-parent",
        child_node_id="missing-child",
        sibling_node_id="missing-sibling",
        exchanged_child_node_id="missing-exchanged-child",
        pivot_branch_id="A",
        sibling_clade_id="B",
        exchanged_clade_id="A",
    )

    report = summarize_rooted_nni_move_validity(
        fixture("felsenstein_two_tip_tree.nwk"),
        candidate,
    )

    assert report.move_family == "rooted-nni"
    assert report.validity_decision == "rejected"
    assert report.available_move_count == 0
    assert report.rejection_code == "topology_move_tree_too_small"
    assert "at least one legal" in report.rejection_reason
    assert report.evidence["tip_count"] == 2


def test_rooted_tbr_move_validity_rejects_missing_branch_requests() -> None:
    tree = loads_newick("(((A,C),B),D);")
    candidate, available_move_count = resolve_rooted_tbr_move_candidate(tree, 1)
    missing_branch_candidate = RootedTbrMoveCandidate(
        cut_parent_node_id=candidate.cut_parent_node_id,
        cut_child_node_id=candidate.cut_child_node_id,
        cut_edge_id="missing-cut-edge",
        cut_descendant_taxa=candidate.cut_descendant_taxa,
        left_attachment_left_node_id=candidate.left_attachment_left_node_id,
        left_attachment_right_node_id=candidate.left_attachment_right_node_id,
        left_attachment_branch_id="missing-left-branch",
        left_attachment_descendant_taxa=candidate.left_attachment_descendant_taxa,
        right_attachment_left_node_id=candidate.right_attachment_left_node_id,
        right_attachment_right_node_id=candidate.right_attachment_right_node_id,
        right_attachment_branch_id=candidate.right_attachment_branch_id,
        right_attachment_descendant_taxa=candidate.right_attachment_descendant_taxa,
    )

    report = summarize_rooted_tbr_move_validity(tree, missing_branch_candidate)

    assert report.move_family == "rooted-tbr"
    assert report.available_move_count == available_move_count
    assert report.validity_decision == "rejected"
    assert report.rejection_code == "topology_move_missing_branch"
    assert report.evidence["missing_branch_fields"] == [
        "cut_edge_id",
        "left_attachment_branch_id",
    ]


def test_rooted_nni_move_validity_rejects_incompatible_rootedness_fixtures() -> None:
    candidate = RootedNniMoveCandidate(
        parent_node_id="missing-parent",
        child_node_id="missing-child",
        sibling_node_id="missing-sibling",
        exchanged_child_node_id="missing-exchanged-child",
        pivot_branch_id="A",
        sibling_clade_id="B",
        exchanged_clade_id="A",
    )

    report = summarize_rooted_nni_move_validity(
        fixture("example_tree_unrooted.nwk"),
        candidate,
    )

    assert report.move_family == "rooted-nni"
    assert report.validity_decision == "rejected"
    assert report.rejection_code == "topology_move_incompatible_rootedness"
    assert report.evidence == {"rooted": False, "root_child_count": 4}

