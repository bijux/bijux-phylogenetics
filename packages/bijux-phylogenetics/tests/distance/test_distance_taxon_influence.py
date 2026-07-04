from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.distance import (
    analyze_distance_taxon_influence_from_imported_distance_matrix,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "metadata" / name


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_distance_taxon_influence_ranks_leave_one_out_improvement_not_missingness() -> (
    None
):
    report = analyze_distance_taxon_influence_from_imported_distance_matrix(
        fixture("example_distance_matrix_taxon_influence_missing_noisy_five_taxon.tsv"),
        tree_fixture("example_tree_minimum_evolution_five_taxon.nwk"),
        method="neighbor-joining",
        missing_distance_policy="mean-impute",
    )

    assert report.taxa == ["A", "B", "C", "D", "E"]
    assert report.baseline_residual_sum_squares == 20.634259259263
    assert report.baseline_rooted_robinson_foulds_distance == 4
    assert [row.influence_rank for row in report.rows] == [1, 2, 3, 4, 5]
    assert [row.taxon for row in report.rows] == ["B", "A", "D", "E", "C"]
    assert report.rows[0].raw_missing_pair_count == 0
    assert report.rows[0].residual_sum_squares_improvement == 19.424259259263
    assert report.rows[0].rooted_robinson_foulds_improvement == 2
    assert report.rows[0].leave_one_out_residual_sum_squares == 1.21
    assert report.rows[1].raw_missing_pair_count == 1
    assert report.rows[3].raw_missing_pair_count == 1
    assert report.rows[0].raw_missing_pair_count < max(
        row.raw_missing_pair_count for row in report.rows
    )


def test_distance_taxon_influence_preserves_baseline_reference_surface() -> None:
    report = analyze_distance_taxon_influence_from_imported_distance_matrix(
        fixture("example_distance_matrix_taxon_influence_missing_noisy_five_taxon.tsv"),
        tree_fixture("example_tree_minimum_evolution_five_taxon.nwk"),
        method="neighbor-joining",
        missing_distance_policy="mean-impute",
    )

    assert report.source_kind == "imported-distance-matrix"
    assert report.method == "neighbor-joining"
    assert report.missing_distance_policy == "mean-impute"
    assert report.baseline_rooted_normalized_robinson_foulds == 1.0
    assert report.rows[-1].taxon == "C"
    assert report.rows[-1].rooted_robinson_foulds_improvement == 1
