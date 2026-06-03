from __future__ import annotations

import csv
import json
from pathlib import Path

from bijux_phylogenetics.phylo.dating import (
    solve_dating_calibration_constraints,
    write_dating_calibration_constraint_artifacts,
    write_dating_calibration_constraint_run_json,
    write_dating_calibration_constraint_summary_tsv,
    write_dating_calibration_constraints_tsv,
    write_dating_calibration_issues_tsv,
    write_dating_calibration_node_windows_tsv,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def build_feasible_report():
    return solve_dating_calibration_constraints(
        fixture(
            "trees",
            "penalized_likelihood_cross_validation_substitution_tree_5_taxa.nwk",
        ),
        fixture("metadata", "dating_calibration_constraints_5_taxa.tsv"),
    )


def build_contradictory_report():
    return solve_dating_calibration_constraints(
        fixture(
            "trees",
            "penalized_likelihood_cross_validation_substitution_tree_5_taxa.nwk",
        ),
        fixture(
            "metadata",
            "dating_calibration_constraints_contradictory_5_taxa.tsv",
        ),
    )


def test_write_dating_calibration_constraint_summary_tsv_writes_expected_row(
    tmp_path: Path,
) -> None:
    report = build_feasible_report()
    output_path = tmp_path / "summary.tsv"

    write_dating_calibration_constraint_summary_tsv(output_path, report)

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert rows == [rows[0]]
    row = rows[0]
    assert row["tree_path"] == str(
        fixture(
            "trees",
            "penalized_likelihood_cross_validation_substitution_tree_5_taxa.nwk",
        )
    )
    assert row["calibration_path"] == str(
        fixture("metadata", "dating_calibration_constraints_5_taxa.tsv")
    )
    assert row["valid_calibration_count"] == "4"
    assert row["contradictory_calibration_count"] == "0"
    assert row["feasible"] == "true"


def test_write_dating_calibration_constraints_tsv_marks_contradictory_rows(
    tmp_path: Path,
) -> None:
    report = build_contradictory_report()
    output_path = tmp_path / "constraints.tsv"

    write_dating_calibration_constraints_tsv(output_path, report)

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 3
    rows_by_id = {row["calibration_id"]: row for row in rows}
    assert rows_by_id["cal-root"]["contradictory"] == "true"
    assert rows_by_id["cal-root"]["issue_codes"] == "chronology-conflict"
    assert float(rows_by_id["cal-ab"]["fixed_date"]) == 1991.0


def test_write_dating_calibration_node_windows_tsv_writes_propagated_bounds(
    tmp_path: Path,
) -> None:
    report = build_feasible_report()
    output_path = tmp_path / "node_windows.tsv"

    write_dating_calibration_node_windows_tsv(output_path, report)

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    rows_by_taxa = {row["descendant_taxa"]: row for row in rows}
    assert float(rows_by_taxa["A|B"]["minimum_bound"]) == 1996.0
    assert float(rows_by_taxa["A|B"]["effective_lower_bound"]) == 1996.0
    assert float(rows_by_taxa["D|E"]["maximum_bound"]) == 1997.0
    assert rows_by_taxa["D|E"]["contradictory"] == "false"


def test_write_dating_calibration_issues_tsv_writes_conflict_rows(
    tmp_path: Path,
) -> None:
    report = build_contradictory_report()
    output_path = tmp_path / "issues.tsv"

    write_dating_calibration_issues_tsv(output_path, report)

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 2
    assert [row["code"] for row in rows] == [
        "chronology-conflict",
        "chronology-conflict",
    ]
    assert rows[0]["scope_kind"] == "edge"
    assert rows[0]["related_node_ids"]


def test_write_dating_calibration_constraint_run_json_serializes_report_fields(
    tmp_path: Path,
) -> None:
    report = build_contradictory_report()
    output_path = tmp_path / "run.json"

    write_dating_calibration_constraint_run_json(output_path, report)

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["feasible"] is False
    assert payload["contradictory_calibration_count"] == 3
    assert len(payload["constraint_rows"]) == 3
    assert len(payload["issue_rows"]) == 2


def test_write_dating_calibration_constraint_artifacts_materializes_governed_outputs(
    tmp_path: Path,
) -> None:
    report = build_contradictory_report()

    outputs = write_dating_calibration_constraint_artifacts(tmp_path, report)

    assert sorted(outputs) == [
        "constraints_path",
        "issues_path",
        "node_windows_path",
        "run_json_path",
        "summary_path",
    ]
    assert (
        outputs["summary_path"]
        .read_text(encoding="utf-8")
        .startswith("tree_path\tcalibration_path\ttip_count\t")
    )
    assert (
        outputs["constraints_path"]
        .read_text(encoding="utf-8")
        .startswith("calibration_id\ttarget_kind\ttarget_label\t")
    )
    assert (
        outputs["node_windows_path"]
        .read_text(encoding="utf-8")
        .startswith("node_id\tnode_kind\tnode_label\t")
    )
    assert (
        outputs["issues_path"]
        .read_text(encoding="utf-8")
        .startswith("scope_kind\tscope_id\tcode\t")
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["feasible"] is False
    assert payload["contradictory_node_count"] == 3
