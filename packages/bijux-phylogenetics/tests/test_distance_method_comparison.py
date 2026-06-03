from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.distance import (
    compare_distance_tree_methods_from_imported_distance_matrix,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "metadata" / name


def test_distance_method_comparison_reports_all_owned_trees_and_scores() -> None:
    report = compare_distance_tree_methods_from_imported_distance_matrix(
        fixture("example_distance_matrix_bionj_noisy.tsv")
    )

    assert report.source_kind == "imported-distance-matrix"
    assert report.missing_distance_policy == "reject"
    assert report.taxa == ["A", "B", "C", "D", "E"]
    assert report.compared_methods == [
        "neighbor-joining",
        "bionj",
        "upgma",
        "wpgma",
    ]
    assert [
        (row.method, row.patristic_residual_sum_squares) for row in report.rows
    ] == [
        ("neighbor-joining", 168.09375),
        ("bionj", 194.315733333333),
        ("upgma", 258.0),
        ("wpgma", 262.59375),
    ]
    assert [
        (row.method, row.balanced_minimum_evolution_score) for row in report.rows
    ] == [
        ("neighbor-joining", 28.25),
        ("bionj", 28.25),
        ("upgma", 16.875),
        ("wpgma", 16.875),
    ]
    assert [
        (
            row.method,
            row.ordinary_least_squares_residual_sum_squares,
            row.ordinary_least_squares_negative_branch_count,
        )
        for row in report.rows
    ] == [
        ("neighbor-joining", 165.75, 1),
        ("bionj", 165.75, 1),
        ("upgma", 165.75, 1),
        ("wpgma", 165.75, 1),
    ]
    assert report.rows[0].tree_newick == (
        "((A:1.125,(B:1,C:2)Inner1:5.375)Inner2:5.375,D:4.625,E:-2.625)Inner3;"
    )
    assert report.rows[1].tree_newick == (
        "((A:2,(B:1,C:2)Inner1:5.66666666666667)Inner2:4.5,D:6.02,E:-4.02)Inner3;"
    )


def test_distance_method_comparison_reports_rf_matrix_and_assumption_warnings() -> None:
    report = compare_distance_tree_methods_from_imported_distance_matrix(
        fixture("example_distance_matrix_bionj_noisy.tsv")
    )

    assert [
        (
            row.left_method,
            row.right_method,
            row.rooted_robinson_foulds_distance,
            row.rooted_normalized_robinson_foulds,
        )
        for row in report.rf_rows
    ] == [
        ("neighbor-joining", "bionj", 0, 0.0),
        ("neighbor-joining", "upgma", 3, 0.6),
        ("neighbor-joining", "wpgma", 3, 0.6),
        ("bionj", "upgma", 3, 0.6),
        ("bionj", "wpgma", 3, 0.6),
        ("upgma", "wpgma", 0, 0.0),
    ]
    assert [(row.scope, row.method, row.warning) for row in report.warning_rows] == [
        (
            "matrix",
            None,
            "distance matrix violates triangle inequality for one or more taxon triples",
        ),
        (
            "matrix",
            None,
            "pairwise distances are not ultrametric, so UPGMA's strict clock-like assumption is violated",
        ),
        (
            "method",
            "bionj",
            "bionj remains a distance-summary method rather than a full likelihood inference",
        ),
        (
            "method",
            "neighbor-joining",
            "neighbor-joining remains a distance-summary method rather than a full likelihood inference",
        ),
        (
            "method",
            "upgma",
            "source matrix is not ultrametric, so this clock-like clustering method is assumption-violating on this input",
        ),
        (
            "method",
            "upgma",
            "upgma assumes an ultrametric clock-like process and can misplace taxa when rates vary among lineages",
        ),
        (
            "method",
            "wpgma",
            "source matrix is not ultrametric, so this clock-like clustering method is assumption-violating on this input",
        ),
        (
            "method",
            "wpgma",
            "wpgma still assumes an ultrametric clock-like process and can overweight small clusters relative to upgma's taxon-count weighting",
        ),
    ]
