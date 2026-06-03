from __future__ import annotations

from bijux_phylogenetics.io.biopython import tree_from_biophylo, tree_to_biophylo
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


def _example_tree() -> PhyloTree:
    return PhyloTree(
        root=TreeNode(
            name="Root",
            children=[
                TreeNode(
                    name="A",
                    branch_length=0.1,
                    metadata={"habitat": "forest"},
                    edge_metadata={"support_class": "tip"},
                ),
                TreeNode(
                    name="Mammals",
                    branch_length=0.2,
                    metadata={"confidence": 95.0},
                    edge_metadata={"support_class": "internal"},
                    children=[
                        TreeNode(name="B", branch_length=0.1),
                        TreeNode(name="C", branch_length=0.1),
                    ],
                ),
            ],
        ),
        rooted=True,
    )


def test_native_tree_model_assigns_stable_node_ids_and_parent_links() -> None:
    tree = _example_tree()

    node_ids = [node.node_id for node in tree.iter_nodes()]

    assert tree.root.parent is None
    assert all(node_id is not None for node_id in node_ids)
    assert len(set(node_ids)) == len(node_ids)
    for parent, child in tree.iter_edges():
        assert child.parent is parent
        assert tree.node_by_id(child.node_id or "") is child


def test_native_tree_model_copies_without_shared_mutation() -> None:
    tree = _example_tree()
    cloned = tree.copy()

    assert cloned is not tree
    assert cloned.root is not tree.root
    assert cloned.root.children[0] is not tree.root.children[0]
    assert cloned.root.children[1].node_id == tree.root.children[1].node_id

    cloned.root.metadata["analysis"] = "copy-only"
    cloned.root.children[0].metadata["habitat"] = "desert"
    cloned.root.children[1].children[0].name = "B_copy"
    cloned.refresh()

    assert "analysis" not in tree.root.metadata
    assert tree.root.children[0].metadata["habitat"] == "forest"
    assert tree.root.children[1].children[0].name == "B"
    assert (
        cloned.root.children[1].children[0].node_id
        != tree.root.children[1].children[0].node_id
    )


def test_native_tree_model_validation_detects_duplicate_parentage() -> None:
    shared_leaf = TreeNode(name="A", branch_length=0.1)
    tree = PhyloTree(
        root=TreeNode(
            children=[
                TreeNode(children=[shared_leaf]),
                TreeNode(children=[shared_leaf]),
            ]
        ),
        rooted=True,
    )

    errors = tree.validation_errors()

    assert any("duplicate parentage" in error for error in errors)


def test_native_tree_model_newick_helpers_roundtrip_deterministically() -> None:
    tree = PhyloTree.from_newick("((B:0.1,C:0.1)Mammals:0.2,A:0.1)Root;")

    assert tree.to_newick() == "(A:0.1,(B:0.1,C:0.1)Mammals:0.2)Root;"
    assert tree.root.name == "Root"
    assert [node.name for node in tree.iter_internal_nodes()] == ["Root", "Mammals"]


def test_native_tree_model_biopython_bridge_preserves_support_metadata() -> None:
    source = _example_tree()

    restored = tree_from_biophylo(tree_to_biophylo(source), source_format="newick")
    internal = next(
        node
        for node in restored.iter_internal_nodes()
        if node is not restored.root and set(node.descendant_taxa) == {"B", "C"}
    )

    assert internal.name == "Mammals"
    assert internal.metadata["confidence"] == 95.0
