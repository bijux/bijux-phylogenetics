from __future__ import annotations

import pytest

import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.io.newick import dumps_newick, loads_newick
from bijux_phylogenetics.phylo.topology import (
    root_unrooted_tree_deterministically,
)


def test_topology_gateway_exports_deterministic_rooting_surface() -> None:
    assert (
        topology_api.root_unrooted_tree_deterministically
        is root_unrooted_tree_deterministically
    )


def test_root_unrooted_tree_deterministically_selects_largest_root_clade() -> None:
    unrooted_tree = loads_newick("((A:0.0625,B:0.0625):0.4375,C:0.0625,D:0.0625);")

    rooted_tree = root_unrooted_tree_deterministically(unrooted_tree)

    assert rooted_tree.rooted is True
    assert len(rooted_tree.root.children) == 2
    assert (
        dumps_newick(rooted_tree)
        == "((A:0.0625,B:0.0625):0,(C:0.0625,D:0.0625):0.4375);"
    )


def test_root_unrooted_tree_deterministically_preserves_binary_rooted_inputs() -> None:
    rooted_tree = loads_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
    rooted_tree.rooted = True

    rerooted_tree = root_unrooted_tree_deterministically(rooted_tree)

    assert dumps_newick(rerooted_tree) == dumps_newick(rooted_tree)


def test_root_unrooted_tree_deterministically_rejects_nontrifurcating_roots() -> None:
    with pytest.raises(
        ValueError,
        match="deterministic rooting requires an explicit unrooted trifurcation at the root",
    ):
        root_unrooted_tree_deterministically(loads_newick("(A:0.1,B:0.2,C:0.3,D:0.4);"))
