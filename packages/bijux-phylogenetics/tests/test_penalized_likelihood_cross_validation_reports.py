from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from bijux_phylogenetics.phylo.dating import (
    cross_validate_penalized_likelihood_smoothing_from_metadata,
    write_penalized_likelihood_cross_validation_artifacts,
    write_penalized_likelihood_cross_validation_candidates_tsv,
    write_penalized_likelihood_cross_validation_predictions_tsv,
    write_penalized_likelihood_cross_validation_run_json,
    write_penalized_likelihood_cross_validation_summary_tsv,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def build_report():
    return cross_validate_penalized_likelihood_smoothing_from_metadata(
        fixture(
            "trees",
            "penalized_likelihood_cross_validation_substitution_tree_5_taxa.nwk",
        ),
        fixture(
            "metadata",
            "penalized_likelihood_cross_validation_tip_dates_5_taxa.tsv",
        ),
        fixture(
            "metadata",
            "penalized_likelihood_cross_validation_calibrations_5_taxa.tsv",
        ),
        smoothing_parameters=[0.01, 0.1, 1.0, 10.0, 100.0],
    )


def test_write_penalized_cross_validation_summary_tsv_writes_expected_row(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "summary.tsv"

    write_penalized_likelihood_cross_validation_summary_tsv(output_path, report)

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
    assert row["metadata_path"] == str(
        fixture(
            "metadata",
            "penalized_likelihood_cross_validation_tip_dates_5_taxa.tsv",
        )
    )
    assert row["calibration_path"] == str(
        fixture(
            "metadata",
            "penalized_likelihood_cross_validation_calibrations_5_taxa.tsv",
        )
    )
    assert row["usable_calibration_count"] == "4"
    assert row["candidate_count"] == "5"
    assert float(row["selected_smoothing_parameter"]) == pytest.approx(0.01, abs=1e-12)
    assert float(row["selected_root_mean_squared_error"]) == pytest.approx(
        0.28564007040784534,
        abs=1e-15,
    )


def test_write_penalized_cross_validation_candidates_tsv_writes_expected_rows(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "candidate_scores.tsv"

    write_penalized_likelihood_cross_validation_candidates_tsv(output_path, report)

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 5
    assert rows[0]["selected"] == "true"
    assert float(rows[0]["smoothing_parameter"]) == pytest.approx(0.01, abs=1e-12)
    assert float(rows[-1]["smoothing_parameter"]) == pytest.approx(100.0, abs=1e-12)


def test_write_penalized_cross_validation_predictions_tsv_writes_expected_rows(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "prediction_errors.tsv"

    write_penalized_likelihood_cross_validation_predictions_tsv(output_path, report)

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 20
    held_out_root_rows = [
        row for row in rows if row["held_out_calibration_id"] == "cal-root"
    ]
    assert len(held_out_root_rows) == 5
    assert float(held_out_root_rows[0]["smoothing_parameter"]) == pytest.approx(
        0.01,
        abs=1e-12,
    )
    assert float(held_out_root_rows[0]["absolute_error"]) == pytest.approx(
        0.5633684388644724,
        abs=1e-15,
    )


def test_write_penalized_cross_validation_run_json_serializes_report_fields(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "run.json"

    write_penalized_likelihood_cross_validation_run_json(output_path, report)

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["usable_calibration_count"] == 4
    assert payload["candidate_count"] == 5
    assert payload["selected_smoothing_parameter"] == 0.01
    assert len(payload["candidate_rows"]) == 5
    assert len(payload["prediction_rows"]) == 20
    assert payload["selected_fit"]["smoothing_parameter"] == 0.01


def test_write_penalized_cross_validation_artifacts_materializes_governed_outputs(
    tmp_path: Path,
) -> None:
    report = build_report()

    outputs = write_penalized_likelihood_cross_validation_artifacts(tmp_path, report)

    assert sorted(outputs) == [
        "branch_rates_path",
        "calibrations_path",
        "candidate_scores_path",
        "dated_tree_path",
        "node_dates_path",
        "prediction_errors_path",
        "run_json_path",
        "summary_path",
    ]
    assert (
        outputs["dated_tree_path"].read_text(encoding="utf-8").strip()
        == report.selected_fit.dated_tree_newick
    )
    assert (
        outputs["summary_path"]
        .read_text(encoding="utf-8")
        .startswith("tree_path\tmetadata_path\tcalibration_path\t")
    )
    assert (
        outputs["candidate_scores_path"]
        .read_text(encoding="utf-8")
        .startswith("smoothing_parameter\tfold_count\tmean_absolute_error\t")
    )
    assert (
        outputs["prediction_errors_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "smoothing_parameter\theld_out_calibration_id\theld_out_target_label\t"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["selected_smoothing_parameter"] == 0.01
