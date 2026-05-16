from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.trees import (
    analyze_branch_length_distribution,
    analyze_tree_set_branch_lengths,
    write_branch_length_table,
)


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_write_branch_length_table_writes_single_tree_rows(tmp_path: Path) -> None:
    output = tmp_path / "branch-lengths.tsv"
    report = analyze_branch_length_distribution(
        tree_fixture("example_tree_long_branch.nwk")
    )

    write_branch_length_table(output, report)

    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[0] == (
        "source_path\ttree_index\tnode\tbranch_type\ttip_taxon\t"
        "descendant_tip_count\tbranch_length\troot_depth\t"
        "tree_positive_branch_median\tzero_length\tnegative_length\t"
        "missing_length\tlong_outlier\tshort_outlier\toutlier_class"
    )
    assert any(
        line.startswith(
            f"{tree_fixture('example_tree_long_branch.nwk')}\t\tA\tterminal\tA\t1\t1.0\t1.1\t0.1\tFalse\tFalse\tFalse\tTrue\tFalse\tlong"
        )
        for line in lines[1:]
    )


def test_write_branch_length_table_writes_tree_set_indices(tmp_path: Path) -> None:
    output = tmp_path / "tree-set-branch-lengths.tsv"
    report = analyze_tree_set_branch_lengths(
        tree_fixture("example_tree_set_branch_lengths.nwk")
    )

    write_branch_length_table(output, report)

    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("source_path\ttree_index\tnode\tbranch_type\t")
    assert any(
        line.startswith(
            f"{tree_fixture('example_tree_set_branch_lengths.nwk')}\t2\tA\tterminal\tA\t1\t1.0\t1.1\t0.1\tFalse\tFalse\tFalse\tTrue\tFalse\tlong"
        )
        for line in lines[1:]
    )
