from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from bijux_phylogenetics.phylo.dating import (
    fit_least_squares_dating_from_metadata,
    write_least_squares_branch_residuals_tsv,
    write_least_squares_dating_artifacts,
    write_least_squares_dating_run_json,
    write_least_squares_dating_summary_tsv,
    write_least_squares_node_dates_tsv,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "metadata")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def build_report():
    return fit_least_squares_dating_from_metadata(
        fixture("least_squares_dating_substitution_tree_4_taxa.nwk"),
        fixture("least_squares_dating_tip_dates_4_taxa.tsv"),
    )


def test_write_least_squares_dating_summary_tsv_writes_expected_row(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "summary.tsv"

    write_least_squares_dating_summary_tsv(output_path, report)

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert rows == [rows[0]]
    row = rows[0]
    assert row["tree_path"] == str(
        fixture("least_squares_dating_substitution_tree_4_taxa.nwk")
    )
    assert row["metadata_path"] == str(
        fixture("least_squares_dating_tip_dates_4_taxa.tsv")
    )
    assert row["taxon_column"] == "taxon"
    assert row["date_column"] == "date"
    assert row["tip_count"] == "4"
    assert row["internal_node_count"] == "3"
    assert row["branch_count"] == "6"
    assert row["parameter_count"] == "4"
    assert float(row["minimum_tip_date"]) == 2007.0
    assert float(row["maximum_tip_date"]) == 2009.0
    assert float(row["root_date"]) == pytest.approx(2000.0, abs=1e-6)
    assert float(row["estimated_clock_rate"]) == pytest.approx(0.25, abs=1e-9)
    assert float(row["residual_sum_squares"]) == pytest.approx(0.0, abs=1e-12)
    assert float(row["condition_number"]) > 1.0
    assert row["exact_fit"] == "true"
    assert row["optimizer_name"] == "closed-form-linear-least-squares"
    assert row["converged"] == "true"


def test_write_least_squares_node_dates_tsv_writes_expected_rows(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "node_dates.tsv"

    write_least_squares_node_dates_tsv(output_path, report)

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 7
    rows_by_descendant_taxa = {
        tuple(row["descendant_taxa"].split("|")): row for row in rows
    }
    assert rows_by_descendant_taxa[("A", "B", "C", "D")]["node_kind"] == "root"
    assert float(
        rows_by_descendant_taxa[("A", "B", "C", "D")]["estimated_date"]
    ) == pytest.approx(
        2000.0,
        abs=1e-6,
    )
    assert rows_by_descendant_taxa[("A", "B")]["node_kind"] == "internal"
    assert float(
        rows_by_descendant_taxa[("A", "B")]["estimated_date"]
    ) == pytest.approx(
        2006.0,
        abs=1e-6,
    )
    assert rows_by_descendant_taxa[("A",)]["node_kind"] == "tip"
    assert rows_by_descendant_taxa[("A",)]["node_label"] == "A"
    assert rows_by_descendant_taxa[("A",)]["fixed_tip_date"] == "true"
    assert float(rows_by_descendant_taxa[("D",)]["time_height"]) == pytest.approx(
        0.0,
        abs=1e-12,
    )


def test_write_least_squares_branch_residuals_tsv_writes_expected_rows(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "branch_residuals.tsv"

    write_least_squares_branch_residuals_tsv(output_path, report)

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 6
    rows_by_descendant_taxa = {
        tuple(row["descendant_taxa"].split("|")): row for row in rows
    }
    assert float(
        rows_by_descendant_taxa[("A", "B", "C")]["fitted_time_duration"]
    ) == pytest.approx(4.0, abs=1e-6)
    assert float(
        rows_by_descendant_taxa[("A", "B", "C")]["observed_branch_length"]
    ) == pytest.approx(1.0, abs=1e-12)
    assert float(
        rows_by_descendant_taxa[("A", "B", "C")]["fitted_branch_length"]
    ) == pytest.approx(1.0, abs=1e-8)
    assert float(rows_by_descendant_taxa[("A",)]["residual"]) == pytest.approx(
        0.0,
        abs=1e-8,
    )
    assert rows_by_descendant_taxa[("D",)]["child_name"] == "D"
    assert float(
        rows_by_descendant_taxa[("D",)]["fitted_time_duration"]
    ) == pytest.approx(
        9.0,
        abs=1e-6,
    )


def test_write_least_squares_dating_run_json_serializes_report_fields(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "run.json"

    write_least_squares_dating_run_json(output_path, report)

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["tip_count"] == 4
    assert payload["internal_node_count"] == 3
    assert payload["branch_count"] == 6
    assert payload["parameter_count"] == 4
    assert payload["root_date"] == 1999.9999999931315
    assert payload["estimated_clock_rate"] == 0.2500000000709406
    assert len(payload["node_rows"]) == 7
    assert len(payload["branch_rows"]) == 6


def test_write_least_squares_dating_artifacts_materializes_governed_outputs(
    tmp_path: Path,
) -> None:
    report = build_report()

    outputs = write_least_squares_dating_artifacts(tmp_path, report)

    assert sorted(outputs) == [
        "branch_residuals_path",
        "dated_tree_path",
        "node_dates_path",
        "run_json_path",
        "summary_path",
    ]
    assert (
        outputs["dated_tree_path"].read_text(encoding="utf-8").strip()
        == report.dated_tree_newick
    )
    assert (
        outputs["summary_path"]
        .read_text(encoding="utf-8")
        .startswith("tree_path\tmetadata_path\ttaxon_column\tdate_column\t")
    )
    assert (
        outputs["node_dates_path"]
        .read_text(encoding="utf-8")
        .startswith("node_id\tnode_kind\tnode_label\tdescendant_taxa\t")
    )
    assert (
        outputs["branch_residuals_path"]
        .read_text(encoding="utf-8")
        .startswith("branch_id\tchild_name\tdescendant_taxa\tparent_date\t")
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["estimated_clock_rate"] == 0.2500000000709406
