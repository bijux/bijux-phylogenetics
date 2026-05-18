from __future__ import annotations

import csv
import json
from pathlib import Path

from bijux_phylogenetics.parity import (
    list_geiger_parity_cases,
    run_geiger_parity_cases,
    write_geiger_parity_observation_table,
    write_geiger_parity_summary_table,
)
from tests.support.fake_geiger_parity import fake_geiger_rscript


def test_list_geiger_parity_cases_returns_governed_registry() -> None:
    cases = list_geiger_parity_cases()

    assert [case.case_id for case in cases] == [
        "fitcontinuous-bm-example-tree",
        "fitcontinuous-bm-brownian-sigma-recovery",
        "fitcontinuous-ou-ou-parameter-recovery",
        "fitcontinuous-eb-early-burst-rate-recovery",
    ]
    assert cases[0].function_name == "geiger::fitContinuous(model='BM')"
    assert cases[2].function_name == "geiger::fitContinuous(model='OU')"
    assert cases[3].function_name == "geiger::fitContinuous(model='EB')"
    assert cases[1].optimizer_settings is not None
    assert cases[2].optimizer_settings["bijux_optimizer_name"] == (
        "governed-two-stage-grid-search"
    )
    assert all(path.is_file() for case in cases for path in case.input_fixtures)


def test_run_geiger_parity_cases_reports_passes_against_fake_runner(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")

    report = run_geiger_parity_cases(
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 4
    assert report.passed_case_count == 4
    assert report.failed_case_count == 0
    assert report.skipped_case_count == 0
    assert report.all_passed is True
    assert len(report.summary_rows) == 3
    observation = next(
        item
        for item in report.observations
        if item.case_id == "fitcontinuous-ou-ou-parameter-recovery"
    )
    assert observation.model_name == "OU"
    assert observation.r_version == "4.4.0"
    assert observation.geiger_version == "2.0.11"
    assert observation.optimizer_settings is not None
    assert observation.bijux_summary["parameter_name"] == "alpha"
    assert observation.reference_rows is not None
    assert any(row["parameter"] == "alpha" for row in observation.reference_rows)


def test_run_geiger_parity_cases_counts_skips_when_geiger_is_unavailable(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        geiger_available=False,
    )

    report = run_geiger_parity_cases(
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )

    assert report.case_count == 4
    assert report.passed_case_count == 0
    assert report.failed_case_count == 0
    assert report.skipped_case_count == 4
    assert report.all_passed is False
    assert all(
        item.mismatch_reason == "geiger_package_unavailable"
        for item in report.observations
    )


def test_run_geiger_parity_cases_persists_failure_artifacts_for_mismatches(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        summary_overrides={
            "fitcontinuous-bm-example-tree": {"root_state": 999.0},
        },
    )
    failure_root = tmp_path / "geiger-parity-failures"

    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-bm-example-tree"],
        rscript_executable=str(rscript),
        failure_root=failure_root,
    )

    assert report.case_count == 1
    assert report.failed_case_count == 1
    observation = report.observations[0]
    assert observation.status == "failed"
    assert observation.mismatch_reason == "summary_field_mismatch:root_state"
    assert observation.reproducible_artifact_root is not None
    artifact_root = observation.reproducible_artifact_root
    assert artifact_root.is_dir()
    assert (artifact_root / "case.json").is_file()
    assert (artifact_root / "reference-summary.json").is_file()
    assert (artifact_root / "bijux-summary.json").is_file()
    assert (
        artifact_root / "mismatch-reason.txt"
    ).read_text(encoding="utf-8") == "summary_field_mismatch:root_state"
    stored_summary = json.loads(
        (artifact_root / "reference-summary.json").read_text(encoding="utf-8")
    )
    assert stored_summary["root_state"] == 999.0


def test_write_geiger_parity_tables_writes_summary_and_observations(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")
    report = run_geiger_parity_cases(
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )
    summary_path = tmp_path / "geiger-parity-summary.tsv"
    observation_path = tmp_path / "geiger-parity-observations.tsv"

    write_geiger_parity_summary_table(summary_path, report)
    write_geiger_parity_observation_table(observation_path, report)

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("function_name\tcase_count")
    assert any("geiger::fitContinuous(model='OU')" in row for row in summary_rows[1:])
    with observation_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 4
    assert rows[0]["model_name"] in {"BM", "OU", "EB"}
    optimizer_settings = json.loads(rows[0]["optimizer_settings"])
    assert "reference_control_policy" in optimizer_settings
