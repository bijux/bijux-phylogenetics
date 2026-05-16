from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
from bijux_phylogenetics.trees import (
    compute_reference_tree_clade_support,
    write_reference_tree_clade_support_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def _row_lookup(report):
    return {tuple(row.descendant_taxa): row for row in report.rows}


def test_compute_reference_tree_clade_support_maps_duplicate_and_conflicting_clades() -> (
    None
):
    report = compute_reference_tree_clade_support(
        fixture("example_tree.nwk"),
        fixture("example_tree_set_left.nwk"),
    )

    rows = _row_lookup(report)

    assert report.tree_count == 3
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.supported_clade_count == 2
    assert report.absent_clade_count == 0
    assert report.unscored_clade_count == 0
    assert rows[("A", "B", "C", "D")].supporting_tree_count == 3
    assert rows[("A", "B", "C", "D")].clade_frequency == pytest.approx(1.0)
    assert rows[("A", "B", "C", "D")].support_status == "fixed"
    assert rows[("A", "B")].supporting_tree_count == 2
    assert rows[("A", "B")].clade_frequency == pytest.approx(2.0 / 3.0)
    assert rows[("A", "B")].support_status == "partial-support"
    assert rows[("C", "D")].supporting_tree_count == 2
    assert rows[("C", "D")].clade_frequency == pytest.approx(2.0 / 3.0)


def test_compute_reference_tree_clade_support_is_child_order_insensitive() -> None:
    report = compute_reference_tree_clade_support(
        fixture("example_tree_support_iqtree_composite.nwk"),
        fixture("example_tree_set_topology_distance_rooted_child_order.nwk"),
    )

    rows = _row_lookup(report)

    assert report.tree_count == 2
    assert rows[("A", "B")].supporting_tree_count == 2
    assert rows[("A", "B")].support_status == "fixed"
    assert rows[("C", "D")].supporting_tree_count == 2
    assert rows[("C", "D")].support_percent == pytest.approx(100.0)


def test_compute_reference_tree_clade_support_matches_root_adjacent_split_support() -> (
    None
):
    report = compute_reference_tree_clade_support(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_posterior_tree_set_six_taxa.nwk"),
    )

    rows = _row_lookup(report)

    assert report.tree_count == 5
    assert rows[("A", "B", "C", "D")].supporting_tree_count == 3
    assert rows[("A", "B", "C", "D")].clade_frequency == pytest.approx(0.6)
    assert rows[("E", "F")].supporting_tree_count == 3
    assert rows[("A", "B")].supporting_tree_count == 2
    assert rows[("C", "D")].supporting_tree_count == 3


def test_compute_reference_tree_clade_support_leaves_absent_root_splits_unscored() -> (
    None
):
    report = compute_reference_tree_clade_support(
        fixture("example_tree_topology_diff.nwk"),
        fixture("example_tree_set_right.nwk"),
    )

    rows = _row_lookup(report)

    assert report.tree_count == 3
    assert report.supported_clade_count == 0
    assert report.absent_clade_count == 0
    assert report.unscored_clade_count == 2
    assert rows[("A", "C")].supporting_tree_count is None
    assert rows[("A", "C")].support_status == "not-counted"
    assert "never realizes the matching bipartition" in rows[("A", "C")].explanation
    assert rows[("B", "D")].supporting_tree_count is None
    assert rows[("B", "D")].clade_frequency is None


def test_compute_reference_tree_clade_support_reports_absent_non_root_clades(
    tmp_path: Path,
) -> None:
    comparison_tree = "(((A:1,C:1):1,(B:1,D:1):1):1,(E:1,F:1):2);"
    comparison_tree_set_path = tmp_path / "support-set-internal-absent.nwk"
    comparison_tree_set_path.write_text(
        "\n".join([comparison_tree] * 3) + "\n",
        encoding="utf-8",
    )

    report = compute_reference_tree_clade_support(
        fixture("example_tree_six_taxa.nwk"),
        comparison_tree_set_path,
    )

    rows = _row_lookup(report)

    assert report.tree_count == 3
    assert report.supported_clade_count == 2
    assert report.absent_clade_count == 2
    assert report.unscored_clade_count == 0
    assert rows[("A", "B")].supporting_tree_count == 0
    assert rows[("A", "B")].support_status == "absent"
    assert rows[("C", "D")].supporting_tree_count == 0
    assert rows[("C", "D")].clade_frequency == pytest.approx(0.0)
    assert rows[("A", "B", "C", "D")].supporting_tree_count == 3
    assert rows[("E", "F")].supporting_tree_count == 3


def test_compute_reference_tree_clade_support_marks_singleton_complements_unscored() -> (
    None
):
    report = compute_reference_tree_clade_support(
        fixture("example_tree_rooted_on_d.nwk"),
        fixture("example_tree_set_left.nwk"),
    )

    rows = _row_lookup(report)

    assert report.unscored_clade_count == 1
    assert rows[("A", "B", "C")].supporting_tree_count is None
    assert rows[("A", "B", "C")].clade_frequency is None
    assert rows[("A", "B", "C")].support_status == "not-counted"
    assert "singleton tip" in rows[("A", "B", "C")].explanation
    assert rows[("A", "B")].supporting_tree_count == 2


def test_compute_reference_tree_clade_support_requires_exact_taxon_match() -> None:
    with pytest.raises(
        InvalidAlignmentError,
        match="share the exact same taxon set",
    ):
        compute_reference_tree_clade_support(
            fixture("example_tree.nwk"),
            fixture("example_tree_set_mismatched.nwk"),
        )


def test_write_reference_tree_clade_support_table_writes_expected_columns(
    tmp_path: Path,
) -> None:
    report = compute_reference_tree_clade_support(
        fixture("example_tree.nwk"),
        fixture("example_tree_set_left.nwk"),
    )

    output_path = tmp_path / "reference-tree-support.tsv"
    write_reference_tree_clade_support_table(output_path, report)

    lines = output_path.read_text(encoding="utf-8").splitlines()

    assert lines[0].startswith(
        "node_id\tnode_kind\tnode_label\tdescendant_taxa\tsupporting_tree_count"
    )
    assert any(
        "\tA|B\t2\t0.666666666666667\t66.6666666666667\tpartial-support\t" in line
        for line in lines[1:]
    )


def test_compute_reference_tree_clade_support_scales_to_thousand_tree_stress_case(
    tmp_path: Path,
) -> None:
    balanced = fixture("example_tree.nwk").read_text(encoding="utf-8").strip()
    cross_pairing = (
        fixture("example_tree_topology_diff.nwk").read_text(encoding="utf-8").strip()
    )
    tree_set_path = tmp_path / "thousand-tree-support-set.nwk"
    tree_set_path.write_text(
        "\n".join([balanced] * 600 + [cross_pairing] * 400) + "\n",
        encoding="utf-8",
    )

    report = compute_reference_tree_clade_support(
        fixture("example_tree.nwk"),
        tree_set_path,
    )

    rows = _row_lookup(report)

    assert report.tree_count == 1000
    assert rows[("A", "B")].supporting_tree_count == 600
    assert rows[("A", "B")].clade_frequency == pytest.approx(0.6)
    assert rows[("C", "D")].supporting_tree_count == 600
