from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.trees import (
    summarize_tree_set_shapes,
    summarize_tree_shape,
    write_tree_shape_table,
)


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_write_tree_shape_table_writes_single_tree_shape_row(tmp_path: Path) -> None:
    output = tmp_path / "tree-shape.tsv"
    report = summarize_tree_shape(tree_fixture("example_tree_ladderized.nwk"))

    write_tree_shape_table(output, report)

    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[0] == (
        "source_path\ttree_index\trooted\ttip_count\tinternal_node_count\tis_binary\t"
        "cherry_count\tsackin_imbalance_index\tcolless_imbalance_index\t"
        "normalized_colless_imbalance\ttree_height_edges\t"
        "tree_height_branch_length\tmean_tip_depth_edges\t"
        "mean_root_to_tip_branch_length\timbalance_summary\tladderized\t"
        "star_like\tcomb_like\tunusually_imbalanced"
    )
    assert (
        f"{tree_fixture('example_tree_ladderized.nwk')}\t\tFalse\t4\t3\tTrue\t1\t9\t3.0\t1.0\t3\t0.3\t2.25\t0.225\tladderized\tTrue\tFalse\tTrue\tTrue"
        in lines[1:]
    )


def test_write_tree_shape_table_writes_tree_set_rows_with_indices(
    tmp_path: Path,
) -> None:
    output = tmp_path / "tree-set-shape.tsv"
    report = summarize_tree_set_shapes(tree_fixture("example_tree_set_left.nwk"))

    write_tree_shape_table(output, report)

    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("source_path\ttree_index\trooted\t")
    assert any(
        line.startswith(
            f"{tree_fixture('example_tree_set_left.nwk')}\t3\tFalse\t4\t3\tTrue\t2\t8\t0.0\t0.0\t2\t0.3\t2.0\t0.3\tbalanced\tFalse\tFalse\tFalse\tFalse"
        )
        for line in lines[1:]
    )
