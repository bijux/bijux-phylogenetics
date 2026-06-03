from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
from bijux_phylogenetics.trees import (
    compute_posterior_clade_correlation_matrix,
    write_posterior_clade_correlation_artifacts,
    write_posterior_clade_correlation_matrix_table,
    write_posterior_clade_correlation_pair_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_compute_posterior_clade_correlation_matrix_reports_joint_structure() -> None:
    report = compute_posterior_clade_correlation_matrix(
        fixture("posterior_clade_correlation_tree_set.nwk")
    )
    row_map = {(row.left_clade, row.right_clade): row for row in report.rows}

    assert report.tree_count == 4
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.clade_count == 6
    assert report.pair_count == 21
    assert report.clade_order == ["A|B", "C|D", "A|C", "A|D", "B|C", "B|D"]

    a_b__c_d = row_map[("A|B", "C|D")]
    assert a_b__c_d.compatibility_relation == "compatible"
    assert a_b__c_d.compatibility_reason == "disjoint"
    assert a_b__c_d.left_tree_count == 2
    assert a_b__c_d.right_tree_count == 2
    assert a_b__c_d.cooccurrence_tree_count == 2
    assert a_b__c_d.cooccurrence_frequency == 0.5
    assert a_b__c_d.expected_cooccurrence_frequency == 0.25
    assert a_b__c_d.binary_correlation == pytest.approx(1.0)

    a_b__a_c = row_map[("A|B", "A|C")]
    assert a_b__a_c.compatibility_relation == "conflict"
    assert a_b__a_c.compatibility_reason == "overlap-without-containment"
    assert a_b__a_c.cooccurrence_tree_count == 0
    assert a_b__a_c.cooccurrence_frequency == 0.0
    assert a_b__a_c.expected_cooccurrence_frequency == 0.125
    assert a_b__a_c.binary_correlation == pytest.approx(-0.577350269189626)

    a_c__b_d = row_map[("A|C", "B|D")]
    assert a_c__b_d.cooccurrence_tree_count == 1
    assert a_c__b_d.cooccurrence_frequency == 0.25
    assert a_c__b_d.binary_correlation == pytest.approx(1.0)


def test_write_posterior_clade_correlation_tables_write_expected_rows(
    tmp_path: Path,
) -> None:
    report = compute_posterior_clade_correlation_matrix(
        fixture("posterior_clade_correlation_tree_set.nwk")
    )

    matrix_path = tmp_path / "posterior-clade-correlation-matrix.tsv"
    pair_path = tmp_path / "posterior-clade-correlation-pairs.tsv"
    write_posterior_clade_correlation_matrix_table(matrix_path, report)
    write_posterior_clade_correlation_pair_table(pair_path, report)

    assert matrix_path.read_text(encoding="utf-8").splitlines()[:4] == [
        "clade\tA|B\tC|D\tA|C\tA|D\tB|C\tB|D",
        "A|B\t1\t1\t-0.577350269189626\t-0.577350269189626\t-0.577350269189626\t-0.577350269189626",
        "C|D\t1\t1\t-0.577350269189626\t-0.577350269189626\t-0.577350269189626\t-0.577350269189626",
        "A|C\t-0.577350269189626\t-0.577350269189626\t1\t-0.333333333333333\t-0.333333333333333\t1",
    ]
    assert pair_path.read_text(encoding="utf-8").splitlines()[:4] == [
        "left_clade\tright_clade\tcompatibility_relation\tcompatibility_reason\tleft_tree_count\tright_tree_count\tleft_frequency\tright_frequency\tcooccurrence_tree_count\tcooccurrence_frequency\texpected_cooccurrence_frequency\tbinary_correlation",
        "A|B\tA|B\tidentical\tsame-clade\t2\t2\t0.5\t0.5\t2\t0.5\t0.25\t1",
        "A|B\tC|D\tcompatible\tdisjoint\t2\t2\t0.5\t0.5\t2\t0.5\t0.25\t1",
        "A|B\tA|C\tconflict\toverlap-without-containment\t2\t1\t0.5\t0.25\t0\t0\t0.125\t-0.577350269189626",
    ]


def test_write_posterior_clade_correlation_artifacts_writes_matrix_and_pairs(
    tmp_path: Path,
) -> None:
    report = compute_posterior_clade_correlation_matrix(
        fixture("posterior_clade_correlation_tree_set.nwk")
    )

    paths = write_posterior_clade_correlation_artifacts(tmp_path, report)

    assert sorted(paths) == [
        "posterior_clade_correlation_matrix_path",
        "posterior_clade_correlation_pair_path",
    ]
    assert (
        paths["posterior_clade_correlation_matrix_path"]
        .read_text(encoding="utf-8")
        .startswith("clade\tA|B\tC|D")
    )
    assert (
        paths["posterior_clade_correlation_pair_path"]
        .read_text(encoding="utf-8")
        .startswith("left_clade\tright_clade")
    )


def test_compute_posterior_clade_correlation_matrix_requires_exact_taxa() -> None:
    with pytest.raises(InvalidAlignmentError, match="exact same taxon set"):
        compute_posterior_clade_correlation_matrix(
            fixture("example_tree_set_mismatched.nwk")
        )
