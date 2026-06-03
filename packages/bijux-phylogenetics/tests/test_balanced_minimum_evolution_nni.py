from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.distance import (
    load_imported_distance_matrix,
    search_balanced_minimum_evolution_nni,
    search_balanced_minimum_evolution_nni_from_imported_distance_matrix,
    validate_imported_distance_matrix,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    for group in ("metadata", "trees"):
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def distance_lookup(name: str) -> tuple[list[str], dict[tuple[str, str], float]]:
    matrix_path = fixture(name)
    validation = validate_imported_distance_matrix(matrix_path)
    entries = load_imported_distance_matrix(matrix_path)
    lookup = {(identifier, identifier): 0.0 for identifier in validation.identifiers}
    for entry in entries:
        lookup[(entry.left_identifier, entry.right_identifier)] = entry.distance
    return validation.identifiers, lookup


def test_search_balanced_minimum_evolution_nni_improves_from_nj_family_starts() -> None:
    identifiers, lookup = distance_lookup(
        "example_distance_matrix_balanced_minimum_evolution_nni_five_taxon.tsv"
    )

    for start_method in ("neighbor-joining", "bionj"):
        report = search_balanced_minimum_evolution_nni(
            identifiers,
            lookup,
            start_method=start_method,
        )

        assert report.algorithm == "balanced-minimum-evolution-nni-search"
        assert report.start_method == start_method
        assert report.taxon_count == 5
        assert report.pair_count == 10
        assert report.start_tree_newick == "(((A,E)Inner1,D)Inner2,B,C)Inner3;"
        assert report.final_tree_newick == "(((A,D)Inner1,E)Inner2,B,C)Inner3;"
        assert report.start_score == 34.0
        assert report.final_score == 33.75
        assert report.accepted_move_count == 1
        assert report.evaluated_neighbor_count == 4
        assert report.stopping_reason == "no-improving-neighbor"
        assert [row.event_kind for row in report.trace_rows] == [
            "start",
            "accepted-move",
            "final",
        ]
        assert report.trace_rows[1].score_before == 34.0
        assert report.trace_rows[1].score_after == 33.75
        assert report.trace_rows[1].score_delta == -0.25
        assert report.trace_rows[1].pivot_branch_id == "A|E"
        assert report.trace_rows[1].sibling_clade_id == "D"
        assert report.trace_rows[1].exchanged_clade_id == "E"
        assert report.trace_rows[-1].stopping_reason == "no-improving-neighbor"


def test_search_balanced_minimum_evolution_nni_from_imported_distance_matrix_carries_matrix_path() -> (
    None
):
    matrix_path = fixture(
        "example_distance_matrix_balanced_minimum_evolution_nni_five_taxon.tsv"
    )

    report = search_balanced_minimum_evolution_nni_from_imported_distance_matrix(
        matrix_path,
        start_method="bionj",
    )

    assert report.matrix_path == matrix_path
    assert report.accepted_move_count == 1
    assert report.start_score == 34.0
    assert report.final_score == 33.75
