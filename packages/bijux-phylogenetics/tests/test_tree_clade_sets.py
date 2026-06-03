from __future__ import annotations

from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_bipartition,
    informative_rooted_clades,
    informative_unrooted_splits,
    node_support_value,
    robinson_foulds_metrics,
    rooted_topology_fingerprint,
    rooted_topology_signature_ids,
    split_sort_key,
    tree_has_polytomy,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


def _balanced_subtree(start: int, stop: int) -> TreeNode:
    if stop - start == 1:
        return TreeNode(name=f"T{start:04d}")
    middle = start + ((stop - start) // 2)
    return TreeNode(
        children=[
            _balanced_subtree(start, middle),
            _balanced_subtree(middle, stop),
        ]
    )


def _balanced_tree(tip_count: int) -> PhyloTree:
    return PhyloTree(root=_balanced_subtree(0, tip_count), rooted=True)


def test_canonical_bipartition_normalizes_side_choice_and_child_order() -> None:
    taxa = {"A", "B", "C", "D", "E", "F"}

    assert canonical_bipartition({"A", "B"}, taxa) == frozenset({"A", "B"})
    assert canonical_bipartition({"C", "D", "E", "F"}, taxa) == frozenset({"A", "B"})
    assert canonical_bipartition({"E", "F"}, taxa) == frozenset({"E", "F"})


def test_rooted_clades_respect_shared_taxon_scope() -> None:
    tree = loads_newick("(((A:1,B:1):1,C:1):1,D:1);")

    assert informative_rooted_clades(tree) == {
        frozenset({"A", "B"}),
        frozenset({"A", "B", "C"}),
    }
    assert informative_rooted_clades(tree, {"A", "B", "C"}) == {frozenset({"A", "B"})}


def test_unrooted_splits_are_canonical_and_ignore_root_orientation() -> None:
    left = loads_newick("((A:1,B:1):1,(C:1,D:1):1);")
    right = loads_newick("((C:1,D:1):1,(B:1,A:1):1);")

    assert informative_unrooted_splits(left) == {frozenset({"A", "B"})}
    assert informative_unrooted_splits(right) == {frozenset({"A", "B"})}


def test_rooted_topology_fingerprint_ignores_branch_lengths_and_child_order() -> None:
    left = loads_newick("(((A:0.1,B:0.2):0.3,C:0.4):0.5,D:0.6);")
    right = loads_newick("((C:9.0,(B:8.0,A:7.0):6.0):5.0,D:4.0);")

    assert rooted_topology_signature_ids(left) == ("A|B", "A|B|C")
    assert rooted_topology_signature_ids(right) == ("A|B", "A|B|C")
    assert rooted_topology_fingerprint(left) == rooted_topology_fingerprint(right)
    assert len(rooted_topology_fingerprint(left)) == 64


def test_rooted_topology_fingerprint_changes_when_rooted_clades_change() -> None:
    left = loads_newick("(((A:1,B:1):1,C:1):1,D:1);")
    right = loads_newick("(((A:1,C:1):1,B:1):1,D:1);")

    assert rooted_topology_fingerprint(left) != rooted_topology_fingerprint(right)


def test_node_support_value_prefers_native_confidence_and_iqtree_composite() -> None:
    confident = loads_newick("((A:0.1,B:0.1)95:0.2,C:0.3);")
    confident_node = next(
        node
        for node in confident.iter_internal_nodes()
        if node is not confident.root and set(node.descendant_taxa) == {"A", "B"}
    )

    iqtree = loads_newick("((A:0.1,B:0.1)82/97:0.2,C:0.3);")
    iqtree_node = next(
        node
        for node in iqtree.iter_internal_nodes()
        if node is not iqtree.root and set(node.descendant_taxa) == {"A", "B"}
    )

    assert node_support_value(confident_node) == 95.0
    assert node_support_value(iqtree_node) == 97.0


def test_robinson_foulds_metrics_report_shared_and_unique_signatures() -> None:
    left = loads_newick("((A:1,B:1):1,(C:1,D:1):1);")
    right = loads_newick("(((A:1,C:1):1,B:1):1,D:1);")
    shared_taxa = {"A", "B", "C", "D"}

    rooted = robinson_foulds_metrics(left, right, shared_taxa, rf_mode="rooted")
    unrooted = robinson_foulds_metrics(left, right, shared_taxa, rf_mode="unrooted")

    assert rooted.shared_signatures == frozenset()
    assert rooted.left_only_signatures == frozenset(
        {frozenset({"A", "B"}), frozenset({"C", "D"})}
    )
    assert rooted.right_only_signatures == frozenset(
        {frozenset({"A", "C"}), frozenset({"A", "B", "C"})}
    )
    assert rooted.distance == 4
    assert rooted.normalized_distance == 1.0

    assert unrooted.shared_signatures == frozenset()
    assert unrooted.left_only_signatures == frozenset({frozenset({"A", "B"})})
    assert unrooted.right_only_signatures == frozenset({frozenset({"A", "C"})})
    assert unrooted.distance == 2
    assert unrooted.normalized_distance == 1.0


def test_tree_has_polytomy_reports_multifurcations() -> None:
    binary = loads_newick("((A:1,B:1):1,(C:1,D:1):1);")
    multifurcating = loads_newick("(A:1,B:1,C:1,D:1);")

    assert tree_has_polytomy(binary) is False
    assert tree_has_polytomy(multifurcating) is True


def test_native_clade_set_core_scales_to_thousand_tip_trees() -> None:
    tree = _balanced_tree(1024)

    rooted_clades = informative_rooted_clades(tree)
    unrooted_splits = informative_unrooted_splits(tree)

    assert len(tree.tip_names) == 1024
    assert len(rooted_clades) == 1022
    assert len(unrooted_splits) == 1021
    assert sorted(rooted_clades, key=split_sort_key)[0] == frozenset({"T0000", "T0001"})
