from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.topology import (
    compare_support_values,
    write_support_comparison_table,
    write_tree_comparison_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_compare_support_values_classifies_conflicting_clades_by_support_strength(
) -> None:
    report = compare_support_values(
        fixture("example_tree_support_conflict_left.nwk"),
        fixture("example_tree_support_conflict_right.nwk"),
    )

    assert [
        (
            row.split_id,
            row.comparison_status,
            row.conflict_classification,
            row.support_strength,
        )
        for row in report.conflicting_clades
    ] == [
        ("C|D", "left_only", "low_support_disagreement", "low"),
        ("A|B|C", "right_only", "high_support_conflict", "strong"),
    ]


def test_write_support_comparison_table_writes_shared_and_conflicting_rows(
    tmp_path: Path,
) -> None:
    output = tmp_path / "support.tsv"

    write_support_comparison_table(
        output,
        fixture("example_tree_support_conflict_left.nwk"),
        fixture("example_tree_support_conflict_right.nwk"),
    )

    assert output.read_text(encoding="utf-8") == (
        "split_id\trow_kind\tcomparison_status\tleft_present\tright_present\tleft_support\tright_support\tleft_support_fraction\tright_support_fraction\tsupport_fraction_delta\tsupport_disagreement\tstrongest_support_fraction\tsupport_strength\tconflict_classification\tdetail\n"
        "A|B\tshared_clade\tshared\ttrue\ttrue\t95.0\t89.0\t0.95\t0.89\t0.05999999999999994\tfalse\t\t\t\tshared clade support is aligned across trees\n"
        "C|D\tconflicting_clade\tleft_only\ttrue\tfalse\t40.0\t\t0.4\t\t\tfalse\t0.4\tlow\tlow_support_disagreement\tconflicting clade was only weakly supported in the tree where it was present\n"
        "A|B|C\tconflicting_clade\tright_only\tfalse\ttrue\t\t92.0\t\t0.92\t\tfalse\t0.92\tstrong\thigh_support_conflict\tconflicting clade carried strong branch support in the tree where it was present\n"
    )


def test_write_tree_comparison_table_carries_support_conflict_columns(
    tmp_path: Path,
) -> None:
    output = tmp_path / "comparison.tsv"

    write_tree_comparison_table(
        output,
        fixture("example_tree_support_conflict_left.nwk"),
        fixture("example_tree_support_conflict_right.nwk"),
    )

    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[0] == (
        "split_id\tcomparison_status\tshared_clade\tleft_support\tright_support\t"
        "left_support_fraction\tright_support_fraction\tsupport_fraction_delta\t"
        "support_disagreement\tsupport_conflict_classification\t"
        "support_conflict_strength\tleft_length\tright_length\tlength_delta\t"
        "length_ratio\tbranch_score_status\tbranch_score_difference\t"
        "branch_score_squared_difference"
    )
    assert (
        "C|D\tleft_only\tfalse\t40.0\t\t0.4\t\t\tfalse\t"
        "low_support_disagreement\tlow\t\t\t\t\t\t\t"
        in lines
    )
    assert (
        "A|B|C\tright_only\tfalse\t\t92.0\t\t0.92\t\tfalse\t"
        "high_support_conflict\tstrong\t\t\t\t\t\t\t"
        in lines
    )
