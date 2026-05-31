from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.io.newick import loads_newick
import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    PhyloTree,
    RootedNniMoveApplicationReport,
    TreeNode,
    apply_rooted_nni_move,
    derive_rooted_nni_reverse_move_candidate,
    resolve_rooted_nni_move_candidate,
    rooted_topology_fingerprint,
    summarize_rooted_nni_move_application,
)

FIXTURES = Path(__file__).parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def _build_metadata_rich_nni_tree() -> PhyloTree:
    return PhyloTree(
        root=TreeNode(
            name="root-label",
            metadata={"study": "rooted-nni"},
            edge_metadata={"edge_label": "root-edge"},
            children=[
                TreeNode(
                    name="inner-label",
                    branch_length=0.4,
                    metadata={"confidence": 0.88},
                    edge_metadata={"edge_label": "inner-edge"},
                    children=[
                        TreeNode(
                            name="ac-label",
                            branch_length=0.3,
                            metadata={"confidence": 0.91},
                            edge_metadata={"edge_label": "ac-edge"},
                            children=[
                                TreeNode(
                                    name="A",
                                    branch_length=0.1,
                                    metadata={"taxon_role": "left-tip"},
                                    edge_metadata={"edge_label": "A-edge"},
                                ),
                                TreeNode(
                                    name="C",
                                    branch_length=0.2,
                                    metadata={"taxon_role": "right-tip"},
                                    edge_metadata={"edge_label": "C-edge"},
                                ),
                            ],
                        ),
                        TreeNode(
                            name="B",
                            branch_length=0.15,
                            metadata={"taxon_role": "middle-tip"},
                            edge_metadata={"edge_label": "B-edge"},
                        ),
                    ],
                ),
                TreeNode(
                    name="D",
                    branch_length=0.25,
                    metadata={"taxon_role": "outer-tip"},
                    edge_metadata={"edge_label": "D-edge"},
                ),
            ],
        ),
        rooted=False,
    )


def test_topology_gateway_exports_rooted_nni_move_application_surface() -> None:
    assert (
        topology_api.RootedNniMoveApplicationReport
        is RootedNniMoveApplicationReport
    )
    assert (
        topology_api.resolve_rooted_nni_move_candidate
        is resolve_rooted_nni_move_candidate
    )
    assert (
        topology_api.derive_rooted_nni_reverse_move_candidate
        is derive_rooted_nni_reverse_move_candidate
    )
    assert (
        topology_api.summarize_rooted_nni_move_application
        is summarize_rooted_nni_move_application
    )


def test_rooted_nni_move_application_restores_original_topology_hash() -> None:
    tree = loads_newick("(((A,C),B),D);")
    candidate, available_move_count = resolve_rooted_nni_move_candidate(tree, 1)

    moved_tree = apply_rooted_nni_move(tree, candidate)
    reverse_candidate = derive_rooted_nni_reverse_move_candidate(
        tree,
        moved_tree,
        candidate,
    )
    restored_tree = apply_rooted_nni_move(moved_tree, reverse_candidate)

    assert available_move_count == 4
    assert rooted_topology_fingerprint(moved_tree) != rooted_topology_fingerprint(tree)
    assert (
        rooted_topology_fingerprint(restored_tree)
        == rooted_topology_fingerprint(tree)
    )


def test_rooted_nni_move_application_report_preserves_labels_metadata_and_branch_lengths(
) -> None:
    tree = _build_metadata_rich_nni_tree()
    original_topology_fingerprint = rooted_topology_fingerprint(tree)
    original_newick = tree.to_newick()

    report = summarize_rooted_nni_move_application(tree, 1)

    assert report.algorithm == "rooted-nni-move-application"
    assert report.available_move_count == 4
    assert report.selected_move_index == 1
    assert report.moved_topology_changed is True
    assert report.reverse_restores_original_topology is True
    assert report.node_names_preserved is True
    assert report.node_metadata_preserved is True
    assert report.edge_metadata_preserved is True
    assert report.branch_lengths_preserved is True
    assert report.total_branch_length_preserved is True
    assert report.missing_tip_taxa == []
    assert report.unexpected_tip_taxa == []
    assert report.moved_validation_errors == []
    assert report.reversed_validation_errors == []
    assert report.input_topology_fingerprint == original_topology_fingerprint
    assert report.reversed_topology_fingerprint == original_topology_fingerprint
    assert tree.to_newick() == original_newick


def test_rooted_nni_move_application_report_preserves_input_tree_path() -> None:
    report = summarize_rooted_nni_move_application(fixture("example_tree.nwk"), 1)

    assert report.input_tree_path == fixture("example_tree.nwk")
    assert report.input_tree_newick == "((A:0.1,B:0.1):0.2,(C:0.2,D:0.2):0.1);"
    assert report.tip_count == 4


def test_rooted_nni_move_resolution_rejects_out_of_range_indices() -> None:
    with pytest.raises(
        ValueError,
        match="rooted NNI move index must be between 1 and 4",
    ):
        resolve_rooted_nni_move_candidate(loads_newick("(((A,C),B),D);"), 5)
