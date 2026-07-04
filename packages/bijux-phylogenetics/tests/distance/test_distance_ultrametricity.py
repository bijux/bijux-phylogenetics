from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.distance import (
    build_tree_from_imported_distance_matrix,
    diagnose_distance_ultrametricity,
    diagnose_imported_distance_matrix_ultrametricity,
    load_imported_distance_matrix,
)
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    assess_tree_ultrametricity,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_diagnose_imported_distance_matrix_ultrametricity_accepts_ultrametric_fixture() -> (
    None
):
    report = diagnose_imported_distance_matrix_ultrametricity(
        fixture("example_distance_matrix_ultrametric.tsv")
    )

    assert report.ultrametric is True
    assert report.tested_triple_count == 4
    assert report.skipped_triple_count == 0
    assert report.max_violation == 0.0
    assert report.violating_triples == []
    assert report.tolerance == 1e-6


def test_diagnose_imported_distance_matrix_ultrametricity_reports_equal_row_sum_failures() -> (
    None
):
    entries = load_imported_distance_matrix(
        fixture("example_distance_matrix_equal_row_sum_nonultrametric.tsv")
    )
    row_sums: dict[str, float] = {}
    for entry in entries:
        if entry.left_identifier == entry.right_identifier:
            continue
        row_sums.setdefault(entry.left_identifier, 0.0)
        row_sums[entry.left_identifier] += entry.distance

    report = diagnose_imported_distance_matrix_ultrametricity(
        fixture("example_distance_matrix_equal_row_sum_nonultrametric.tsv")
    )

    assert sorted(row_sums.values()) == [12.0, 12.0, 12.0, 12.0]
    assert report.ultrametric is False
    assert report.tested_triple_count == 4
    assert report.skipped_triple_count == 0
    assert report.max_violation == 3.0
    assert [
        (
            row.left_identifier,
            row.middle_identifier,
            row.right_identifier,
            row.violation,
        )
        for row in report.violating_triples
    ] == [
        ("A", "B", "C", 3.0),
        ("A", "B", "D", 3.0),
        ("A", "C", "D", 3.0),
        ("B", "C", "D", 3.0),
    ]


def test_diagnose_distance_ultrametricity_respects_tolerance() -> None:
    strict = diagnose_distance_ultrametricity(
        fixture("example_alignment_distance.fasta"),
        model="p-distance",
        tolerance=1e-6,
    )
    loose = diagnose_distance_ultrametricity(
        fixture("example_alignment_distance.fasta"),
        model="p-distance",
        tolerance=0.2,
    )

    assert strict.ultrametric is False
    assert len(strict.violating_triples) == 4
    assert strict.max_violation == 0.125
    assert loose.ultrametric is True
    assert loose.violating_triples == []
    assert loose.max_violation == 0.125
    assert loose.tolerance == 0.2


def test_matrix_ultrametricity_is_not_tree_ultrametricity(tmp_path: Path) -> None:
    matrix_report = diagnose_imported_distance_matrix_ultrametricity(
        fixture("example_distance_matrix_equal_row_sum_nonultrametric.tsv")
    )
    tree, _ = build_tree_from_imported_distance_matrix(
        fixture("example_distance_matrix_equal_row_sum_nonultrametric.tsv"),
        method="upgma",
    )
    tree_path = tmp_path / "upgma-tree.nwk"
    write_newick(tree_path, tree)
    tree_report = assess_tree_ultrametricity(tree_path)

    assert matrix_report.ultrametric is False
    assert tree_report.ultrametric is True
