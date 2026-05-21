from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.io.newick import (
    dumps_newick_tree_set,
    load_newick_tree_set,
    write_newick_tree_set,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    TreeParseError,
    UnnamedTipError,
)


def test_write_newick_tree_set_roundtrips_deterministically(tmp_path: Path) -> None:
    trees = load_newick_tree_set(
        Path(__file__).parent / "fixtures" / "trees" / "example_tree_set_left.nwk"
    )

    first_text = dumps_newick_tree_set(trees)
    second_text = dumps_newick_tree_set(trees)
    output_path = write_newick_tree_set(tmp_path / "tree-set.nwk", trees)

    assert first_text == second_text
    assert output_path.read_text(encoding="utf-8") == first_text
    assert output_path.read_text(encoding="utf-8").count("\n") == len(trees)


def test_dumps_newick_tree_set_rejects_empty_tree_sets() -> None:
    with pytest.raises(TreeParseError):
        dumps_newick_tree_set([])


def test_write_newick_tree_set_rejects_unnamed_tips() -> None:
    malformed_tree = PhyloTree(
        root=TreeNode(
            children=[
                TreeNode(name="A", branch_length=0.1),
                TreeNode(branch_length=0.2),
            ]
        )
    )

    with pytest.raises(UnnamedTipError):
        dumps_newick_tree_set([malformed_tree])


def test_write_newick_tree_set_rejects_nonfinite_branch_lengths() -> None:
    malformed_tree = PhyloTree(
        root=TreeNode(
            children=[
                TreeNode(name="A", branch_length=0.1),
                TreeNode(name="B", branch_length=math.inf),
            ]
        )
    )

    with pytest.raises(InvalidBranchLengthError):
        dumps_newick_tree_set([malformed_tree])
