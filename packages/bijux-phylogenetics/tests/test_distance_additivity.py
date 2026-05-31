from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.distance import (
    diagnose_distance_additivity,
    diagnose_imported_distance_matrix_additivity,
    write_distance_additivity_artifacts,
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


def test_diagnose_imported_distance_matrix_additivity_accepts_additive_fixture() -> (
    None
):
    report = diagnose_imported_distance_matrix_additivity(
        fixture("example_distance_matrix_ultrametric.tsv")
    )

    assert report.additive is True
    assert report.tested_quartet_count == 1
    assert report.skipped_quartet_count == 0
    assert report.max_violation == 0.0
    assert report.violating_quartets == []


def test_diagnose_imported_distance_matrix_additivity_reports_four_point_failure() -> (
    None
):
    report = diagnose_imported_distance_matrix_additivity(
        fixture("example_distance_matrix_equal_row_sum_nonultrametric.tsv")
    )

    assert report.additive is False
    assert report.tested_quartet_count == 1
    assert report.skipped_quartet_count == 0
    assert report.max_violation == 6.0
    assert [
        (
            row.first_identifier,
            row.second_identifier,
            row.third_identifier,
            row.fourth_identifier,
            row.split_ab_cd_sum,
            row.split_ac_bd_sum,
            row.split_ad_bc_sum,
            row.best_split,
            row.violation_magnitude,
        )
        for row in report.violating_quartets
    ] == [
        ("A", "B", "C", "D", 2.0, 8.0, 14.0, "A,B|C,D", 6.0),
    ]


def test_diagnose_distance_additivity_reports_alignment_violation() -> None:
    report = diagnose_distance_additivity(
        fixture("example_alignment_distance.fasta"),
        model="p-distance",
    )

    assert report.additive is False
    assert report.tested_quartet_count == 1
    assert report.skipped_quartet_count == 0
    assert report.max_violation == 0.25
    assert [
        (
            row.best_split,
            row.violation_magnitude,
        )
        for row in report.violating_quartets
    ] == [("A,B|C,D", 0.25)]


def test_write_distance_additivity_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = diagnose_imported_distance_matrix_additivity(
        fixture("example_distance_matrix_equal_row_sum_nonultrametric.tsv")
    )

    outputs = write_distance_additivity_artifacts(
        tmp_path / "distance-additivity", report
    )

    table_lines = (
        outputs["four_point_violations"].read_text(encoding="utf-8").splitlines()
    )
    assert table_lines[0] == (
        "quartet\tsplit_ab_cd_sum\tsplit_ac_bd_sum\tsplit_ad_bc_sum\tbest_split\tviolation_magnitude"
    )
    assert table_lines[1] == "A|B|C|D\t2\t8\t14\tA,B|C,D\t6"
    assert outputs["run_json"].exists()
