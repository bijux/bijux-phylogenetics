from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from bijux_phylogenetics.phylo.dating import (
    fit_penalized_likelihood_dating_from_metadata,
    write_penalized_likelihood_branch_rate_tsv,
    write_penalized_likelihood_dating_artifacts,
    write_penalized_likelihood_dating_run_json,
    write_penalized_likelihood_dating_summary_tsv,
    write_penalized_likelihood_node_dates_tsv,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "metadata")
DATE_ABS_TOLERANCE = 1e-5
RATE_ABS_TOLERANCE = 1e-7
SCORE_REL_TOLERANCE = 1e-6
SCORE_ABS_TOLERANCE = 1e-12


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
    return fit_penalized_likelihood_dating_from_metadata(
        fixture("penalized_likelihood_dating_substitution_tree_4_taxa.nwk"),
        fixture("penalized_likelihood_dating_tip_dates_4_taxa.tsv"),
        smoothing_parameter=0.01,
    )


def test_write_penalized_likelihood_dating_summary_tsv_writes_expected_row(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "summary.tsv"

    write_penalized_likelihood_dating_summary_tsv(output_path, report)

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert rows == [rows[0]]
    row = rows[0]
    assert row["tree_path"] == str(
        fixture("penalized_likelihood_dating_substitution_tree_4_taxa.nwk")
    )
    assert row["metadata_path"] == str(
        fixture("penalized_likelihood_dating_tip_dates_4_taxa.tsv")
    )
    assert row["taxon_column"] == "taxon"
    assert row["date_column"] == "date"
    assert row["tip_count"] == "4"
    assert row["internal_node_count"] == "3"
    assert row["branch_count"] == "6"
    assert row["parameter_count"] == "10"
    assert float(row["root_date"]) == pytest.approx(
        1985.738765803845, abs=DATE_ABS_TOLERANCE
    )
    assert float(row["smoothing_parameter"]) == pytest.approx(0.01, abs=1e-12)
    assert float(row["data_score"]) == pytest.approx(
        2.9289583718815336e-06,
        rel=SCORE_REL_TOLERANCE,
        abs=SCORE_ABS_TOLERANCE,
    )
    assert float(row["penalty_score"]) == pytest.approx(
        0.00013482416705344673,
        rel=SCORE_REL_TOLERANCE,
        abs=SCORE_ABS_TOLERANCE,
    )
    assert float(row["total_score"]) == pytest.approx(
        0.00013775312542532825,
        rel=SCORE_REL_TOLERANCE,
        abs=SCORE_ABS_TOLERANCE,
    )
    assert float(row["condition_number"]) > 1.0
    assert row["optimizer_name"] == (
        "bounded-coordinate-search with closed-form penalized log-rate solve"
    )
    assert row["optimization_pass_count"] == "5"
    assert row["function_evaluation_count"] == "771"
    assert row["converged"] == "true"


def test_write_penalized_likelihood_node_dates_tsv_writes_expected_rows(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "node_dates.tsv"

    write_penalized_likelihood_node_dates_tsv(output_path, report)

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
        1985.738765803845,
        abs=DATE_ABS_TOLERANCE,
    )
    assert float(
        rows_by_descendant_taxa[("A", "B", "C")]["estimated_rate"]
    ) == pytest.approx(
        0.08392676134136667,
        abs=RATE_ABS_TOLERANCE,
    )
    assert rows_by_descendant_taxa[("A",)]["node_kind"] == "tip"
    assert rows_by_descendant_taxa[("A",)]["node_label"] == "A"
    assert rows_by_descendant_taxa[("A",)]["fixed_tip_date"] == "true"
    assert float(rows_by_descendant_taxa[("D",)]["time_height"]) == pytest.approx(
        0.0,
        abs=1e-12,
    )


def test_write_penalized_likelihood_branch_rate_tsv_writes_expected_rows(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "branch_rates.tsv"

    write_penalized_likelihood_branch_rate_tsv(output_path, report)

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 6
    rows_by_descendant_taxa = {
        tuple(row["descendant_taxa"].split("|")): row for row in rows
    }
    assert float(
        rows_by_descendant_taxa[("A", "B", "C")]["estimated_branch_rate"]
    ) == pytest.approx(0.08136342580549531, abs=RATE_ABS_TOLERANCE)
    assert float(
        rows_by_descendant_taxa[("A", "B", "C")]["data_score_contribution"]
    ) == pytest.approx(
        3.8334258705446777e-07,
        rel=SCORE_REL_TOLERANCE,
        abs=SCORE_ABS_TOLERANCE,
    )
    assert float(
        rows_by_descendant_taxa[("A",)]["smoothing_penalty_contribution"]
    ) == pytest.approx(
        2.2364498263509753e-06,
        rel=SCORE_REL_TOLERANCE,
        abs=SCORE_ABS_TOLERANCE,
    )
    assert rows_by_descendant_taxa[("D",)]["child_name"] == "D"
    assert float(
        rows_by_descendant_taxa[("D",)]["fitted_time_duration"]
    ) == pytest.approx(23.261234196154992, abs=1e-9)


def test_write_penalized_likelihood_dating_run_json_serializes_report_fields(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "run.json"

    write_penalized_likelihood_dating_run_json(output_path, report)

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["tip_count"] == 4
    assert payload["internal_node_count"] == 3
    assert payload["branch_count"] == 6
    assert payload["parameter_count"] == 10
    assert payload["smoothing_parameter"] == 0.01
    assert payload["data_score"] == pytest.approx(
        2.9289583718815336e-06,
        rel=SCORE_REL_TOLERANCE,
        abs=SCORE_ABS_TOLERANCE,
    )
    assert payload["penalty_score"] == pytest.approx(
        0.00013482416705344673,
        rel=SCORE_REL_TOLERANCE,
        abs=SCORE_ABS_TOLERANCE,
    )
    assert payload["total_score"] == pytest.approx(
        0.00013775312542532825,
        rel=SCORE_REL_TOLERANCE,
        abs=SCORE_ABS_TOLERANCE,
    )
    assert len(payload["node_rows"]) == 7
    assert len(payload["branch_rows"]) == 6


def test_write_penalized_likelihood_dating_artifacts_materializes_governed_outputs(
    tmp_path: Path,
) -> None:
    report = build_report()

    outputs = write_penalized_likelihood_dating_artifacts(tmp_path, report)

    assert sorted(outputs) == [
        "branch_rates_path",
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
        outputs["branch_rates_path"]
        .read_text(encoding="utf-8")
        .startswith("branch_id\tchild_name\tdescendant_taxa\tparent_date\t")
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["smoothing_parameter"] == 0.01
