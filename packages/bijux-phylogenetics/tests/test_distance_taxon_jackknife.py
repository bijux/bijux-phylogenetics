from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.distance import (
    analyze_distance_taxon_jackknife_from_imported_distance_matrix,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "metadata" / name


def test_distance_taxon_jackknife_rebuilds_each_reduced_tree_and_reports_deltas() -> (
    None
):
    report = analyze_distance_taxon_jackknife_from_imported_distance_matrix(
        fixture("example_distance_matrix_taxon_influence_missing_noisy_five_taxon.tsv"),
        method="neighbor-joining",
        missing_distance_policy="mean-impute",
    )

    assert report.source_kind == "imported-distance-matrix"
    assert report.method == "neighbor-joining"
    assert report.missing_distance_policy == "mean-impute"
    assert report.taxa == ["A", "B", "C", "D", "E"]
    assert report.baseline_residual_sum_squares == 20.634259259263
    assert report.baseline_tree_newick == (
        "(((A:0.722222222222223,C:7.27777777777778)Inner1:5.58333333333333,"
        "E:5.58333333333333)Inner2:1.91666666666667,B:-3.75,D:4.75)Inner3;"
    )
    assert [row.removed_taxon for row in report.rows] == ["A", "B", "C", "D", "E"]
    assert all(row.rebuilt_tree_newick.endswith(";") for row in report.rows)

    first = report.rows[0]
    assert first.retained_taxa == ["B", "C", "D", "E"]
    assert first.rooted_robinson_foulds_distance == 2
    assert first.rooted_normalized_robinson_foulds == 1.0
    assert first.pruned_baseline_residual_sum_squares == 6.879629629632
    assert first.rebuilt_residual_sum_squares == 2.25
    assert first.residual_sum_squares_change == 4.629629629632
    assert first.reference_only_clades == ["C|E"]
    assert first.rebuilt_only_clades == ["B|C"]
    assert first.affected_clades == ["B|C", "C|E"]
    assert first.topology_changed is True

    second = report.rows[1]
    assert second.removed_taxon == "B"
    assert second.rooted_robinson_foulds_distance == 1
    assert second.rooted_normalized_robinson_foulds == 0.333333333333
    assert second.residual_sum_squares_change == 20.692808641978
    assert second.affected_clades == ["A|C|E"]
    assert second.topology_changed is True

    third = report.rows[2]
    assert third.removed_taxon == "C"
    assert third.rooted_robinson_foulds_distance == 0
    assert third.affected_clades == []
    assert third.topology_changed is False
    assert third.residual_sum_squares_change == 24.195555555557


def test_distance_taxon_jackknife_compares_against_pruned_baseline_clades() -> None:
    report = analyze_distance_taxon_jackknife_from_imported_distance_matrix(
        fixture("example_distance_matrix_taxon_influence_missing_noisy_five_taxon.tsv"),
        method="neighbor-joining",
        missing_distance_policy="mean-impute",
    )

    row = report.rows[3]
    assert row.removed_taxon == "D"
    assert row.reference_only_clades == ["A|C|E"]
    assert row.rebuilt_only_clades == []
    assert row.affected_clades == ["A|C|E"]
    assert row.pruned_baseline_tree_newick == (
        "(((A:0.722222222222223,C:7.27777777777778)Inner1:5.58333333333333,"
        "E:5.58333333333333)Inner2:1.91666666666667,B:-3.75)Inner3;"
    )
    assert row.rebuilt_tree_newick == (
        "((A:-0.45,C:8.45)Inner1:4.05,B:-1.55,E:6.55)Inner2;"
    )
