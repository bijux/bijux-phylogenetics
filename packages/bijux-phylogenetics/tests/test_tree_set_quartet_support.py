from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    TreeSetQuartetSupportReport,
    TreeSetQuartetSupportRow,
    compute_reference_tree_quartet_support,
    write_reference_tree_quartet_support_table,
)


def test_package_tree_gateway_exports_quartet_support_surface() -> None:
    assert trees_api.TreeSetQuartetSupportRow is TreeSetQuartetSupportRow
    assert trees_api.TreeSetQuartetSupportReport is TreeSetQuartetSupportReport
    assert (
        trees_api.compute_reference_tree_quartet_support
        is compute_reference_tree_quartet_support
    )
    assert (
        trees_api.write_reference_tree_quartet_support_table
        is write_reference_tree_quartet_support_table
    )


def test_compute_reference_tree_quartet_support_matches_hand_counted_fixture(
    tmp_path: Path,
) -> None:
    reference_tree = tmp_path / "quartet-reference-tree.nwk"
    comparison_tree_set = tmp_path / "quartet-support-tree-set.nwk"
    reference_tree.write_text("((A,B),(C,D));\n", encoding="utf-8")
    comparison_tree_set.write_text(
        "\n".join(
            [
                "((A,B),(C,D));",
                "((A,C),(B,D));",
                "((A,D),(B,C));",
                "(A,B,C,D);",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    report = compute_reference_tree_quartet_support(
        reference_tree,
        comparison_tree_set,
    )

    assert report.tree_count == 4
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.branch_count == 1
    assert report.total_quartet_count == 4
    assert report.concordant_quartet_count == 1
    assert report.discordant_quartet_count == 2
    assert report.uninformative_quartet_count == 1
    assert report.rows == [
        TreeSetQuartetSupportRow(
            branch_id="A|B::C|D",
            left_taxa=["A", "B"],
            right_taxa=["C", "D"],
            quartet_count_per_tree=1,
            concordant_quartet_count=1,
            discordant_quartet_count=2,
            uninformative_quartet_count=1,
            concordant_frequency=0.25,
            discordant_frequency=0.5,
            uninformative_frequency=0.25,
        )
    ]


def test_compute_reference_tree_quartet_support_requires_exact_taxon_match(
    tmp_path: Path,
) -> None:
    reference_tree = tmp_path / "quartet-reference-tree.nwk"
    comparison_tree_set = tmp_path / "quartet-support-tree-set-mismatch.nwk"
    reference_tree.write_text("((A,B),(C,D));\n", encoding="utf-8")
    comparison_tree_set.write_text("((A,B),(C,E));\n", encoding="utf-8")

    with pytest.raises(InvalidAlignmentError, match="exact same taxon set"):
        compute_reference_tree_quartet_support(
            reference_tree,
            comparison_tree_set,
        )


def test_write_reference_tree_quartet_support_table_writes_expected_columns(
    tmp_path: Path,
) -> None:
    reference_tree = tmp_path / "quartet-reference-tree.nwk"
    comparison_tree_set = tmp_path / "quartet-support-tree-set.nwk"
    output_path = tmp_path / "quartet-support.tsv"
    reference_tree.write_text("((A,B),(C,D));\n", encoding="utf-8")
    comparison_tree_set.write_text(
        "\n".join(
            [
                "((A,B),(C,D));",
                "((A,C),(B,D));",
                "((A,D),(B,C));",
                "(A,B,C,D);",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    report = compute_reference_tree_quartet_support(
        reference_tree,
        comparison_tree_set,
    )
    write_reference_tree_quartet_support_table(output_path, report)

    lines = output_path.read_text(encoding="utf-8").splitlines()

    assert lines[0] == (
        "branch_id\tleft_taxa\tright_taxa\tquartet_count_per_tree\t"
        "concordant_quartet_count\tdiscordant_quartet_count\t"
        "uninformative_quartet_count\tconcordant_frequency\t"
        "discordant_frequency\tuninformative_frequency"
    )
    assert lines[1] == "A|B::C|D\tA|B\tC|D\t1\t1\t2\t1\t0.25\t0.5\t0.25"


def test_compute_reference_tree_quartet_support_handles_multiple_quartets_per_branch(
    tmp_path: Path,
) -> None:
    reference_tree = tmp_path / "quartet-reference-tree-five-taxa.nwk"
    comparison_tree_set = tmp_path / "quartet-support-tree-set-five-taxa.nwk"
    reference_tree.write_text("(((A,B),C),(D,E));\n", encoding="utf-8")
    comparison_tree_set.write_text(
        "\n".join(
            [
                "(((A,B),C),(D,E));",
                "((E,D),(C,(B,A)));",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    report = compute_reference_tree_quartet_support(
        reference_tree,
        comparison_tree_set,
    )

    assert report.branch_count == 2
    assert report.total_quartet_count == 12
    assert report.concordant_quartet_count == 12
    assert report.discordant_quartet_count == 0
    assert report.uninformative_quartet_count == 0
    assert [
        (
            row.branch_id,
            row.left_taxa,
            row.right_taxa,
            row.quartet_count_per_tree,
            row.concordant_quartet_count,
        )
        for row in report.rows
    ] == [
        ("A|B::C|D|E", ["A", "B"], ["C", "D", "E"], 3, 6),
        ("D|E::A|B|C", ["D", "E"], ["A", "B", "C"], 3, 6),
    ]
