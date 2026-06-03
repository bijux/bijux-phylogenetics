from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.distance import (
    compute_patristic_residual_diagnostics_from_imported_distance_matrix,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    for group in ("metadata", "trees"):
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_compute_patristic_residual_diagnostics_reports_ranked_unique_pairs() -> None:
    report = compute_patristic_residual_diagnostics_from_imported_distance_matrix(
        fixture("example_distance_matrix_minimum_evolution_five_taxon.tsv"),
        fixture("example_tree_minimum_evolution_five_taxon.nwk"),
    )

    assert report.taxa == ["A", "B", "C", "D", "E"]
    assert report.pair_count == 10
    assert report.residual_sum_squares == 2626.0
    assert report.max_absolute_residual == 20.0
    assert [
        (
            row.left_identifier,
            row.right_identifier,
            row.observed_distance,
            row.fitted_distance,
            row.residual,
            row.absolute_residual,
            row.rank,
        )
        for row in report.rows
    ] == [
        ("A", "D", 16.0, 36.0, -20.0, 20.0, 1),
        ("A", "C", 8.0, 27.0, -19.0, 19.0, 2),
        ("A", "E", 17.0, 36.0, -19.0, 19.0, 3),
        ("B", "D", 17.0, 36.0, -19.0, 19.0, 4),
        ("B", "C", 9.0, 27.0, -18.0, 18.0, 5),
        ("B", "E", 18.0, 36.0, -18.0, 18.0, 6),
        ("A", "B", 3.0, 18.0, -15.0, 15.0, 7),
        ("C", "D", 16.0, 27.0, -11.0, 11.0, 8),
        ("C", "E", 17.0, 27.0, -10.0, 10.0, 9),
        ("D", "E", 11.0, 18.0, -7.0, 7.0, 10),
    ]
