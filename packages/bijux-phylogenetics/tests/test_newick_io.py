from __future__ import annotations

import pytest

from bijux_phylogenetics.io.newick import (
    dumps_newick,
    iter_newick_tree_records,
    loads_newick,
    loads_newick_tree_set,
)
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError, TreeParseError


def test_loads_newick_preserves_numeric_support_metadata() -> None:
    tree = loads_newick("((A:0.1,B:0.1)95:0.2,C:0.3)Root;")

    internal = next(
        node
        for node in tree.iter_internal_nodes()
        if node is not tree.root and set(node.descendant_taxa) == {"A", "B"}
    )

    assert internal.name == "95"
    assert internal.metadata["confidence"] == 95.0
    assert dumps_newick(tree) == "((A:0.1,B:0.1)95:0.2,C:0.3)Root;"


def test_loads_newick_unescapes_embedded_quotes_in_labels() -> None:
    tree = loads_newick("('O''Brien taxon':0.1,B:0.2);")

    assert tree.tip_names == ["O'Brien taxon", "B"]
    assert dumps_newick(tree) == "(B:0.2,'O''Brien taxon':0.1);"


def test_loads_newick_allows_unnamed_tips_with_branch_lengths() -> None:
    tree = loads_newick("((:0.1,B:0.2):0.3,C:0.4);")

    unnamed_tip = next(node for node in tree.iter_leaves() if node.name is None)

    assert unnamed_tip.branch_length == 0.1
    assert tree.tip_count == 3


def test_loads_newick_tree_set_reads_multiple_records_from_text() -> None:
    trees = loads_newick_tree_set("((A:0.1,B:0.1):0.2,C:0.3);\n(A:0.1,B:0.2,C:0.3);\n")

    assert len(trees) == 2
    assert dumps_newick(trees[0]) == "((A:0.1,B:0.1):0.2,C:0.3);"
    assert dumps_newick(trees[1]) == "(A:0.1,B:0.2,C:0.3);"


def test_iter_newick_tree_records_normalizes_multiline_and_trailing_records() -> None:
    records = list(
        iter_newick_tree_records(
            "((A:0.1,B:0.1):0.2,\nC:0.3);\n\n(A:0.1,B:0.2,C:0.3)\n"
        )
    )

    assert records == [
        (1, "((A:0.1,B:0.1):0.2, C:0.3);"),
        (2, "(A:0.1,B:0.2,C:0.3)"),
    ]


def test_loads_newick_reports_location_for_malformed_structure() -> None:
    with pytest.raises(TreeParseError) as error_info:
        loads_newick("((A:0.1,B:0.2):0.3,C:0.4")

    error = error_info.value
    assert error.details == {"position": 24, "line": 1, "column": 25}
    assert "line 1, column 25" in error.message


def test_loads_newick_reports_location_for_invalid_branch_length() -> None:
    with pytest.raises(InvalidBranchLengthError) as error_info:
        loads_newick("((A:abc,B:0.2):0.3,C:0.4);")

    error = error_info.value
    assert error.details == {"position": 4, "line": 1, "column": 5}
    assert "line 1, column 5" in error.message
