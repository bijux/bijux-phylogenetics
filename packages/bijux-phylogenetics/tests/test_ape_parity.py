from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from bijux_phylogenetics.ape_parity import (
    list_ape_parity_cases,
    run_ape_parity_cases,
    write_ape_parity_observation_table,
    write_ape_parity_summary_table,
)
from tests.support.fake_reference_parity import fake_ape_rscript


def _r_package_available(rscript: str, package_name: str) -> bool:
    repository_root = Path(__file__).resolve().parents[3]
    environment = dict(os.environ)
    r_library = repository_root / "artifacts" / "r-lib"
    if r_library.is_dir():
        environment["R_LIBS_USER"] = str(r_library)
    result = subprocess.run(
        [
            rscript,
            "-e",
            f"cat(requireNamespace('{package_name}', quietly=TRUE), '\\n')",
        ],
        capture_output=True,
        check=False,
        cwd=repository_root,
        env=environment,
        text=True,
    )
    return result.returncode == 0 and result.stdout.strip() == "TRUE"


def test_list_ape_parity_cases_returns_governed_read_tree_registry() -> None:
    cases = list_ape_parity_cases()

    assert [case.case_id for case in cases] == [
        "read-tree-example-rooted",
        "read-tree-example-unrooted",
    ]
    assert {case.function_name for case in cases} == {"ape::read.tree"}
    assert {case.operation for case in cases} == {"read-tree-summary"}


def test_run_ape_parity_cases_passes_against_fake_reference_runner(
    tmp_path: Path,
) -> None:
    rscript = fake_ape_rscript(tmp_path / "fake-ape-rscript")

    report = run_ape_parity_cases(
        rscript_executable=str(rscript),
        failure_root=tmp_path / "ape-parity-failures",
    )

    assert report.all_passed is True
    assert report.case_count == 2
    assert report.passed_case_count == 2
    assert report.failed_case_count == 0
    assert report.skipped_case_count == 0
    assert [row.function_name for row in report.summary_rows] == ["ape::read.tree"]
    assert all(observation.r_version == "4.6.0" for observation in report.observations)
    assert all(observation.ape_version == "5.0.0" for observation in report.observations)
    assert all(observation.reproducible_artifact_root is None for observation in report.observations)


def test_run_ape_parity_cases_records_failure_bundle_for_summary_mismatch(
    tmp_path: Path,
) -> None:
    rscript = fake_ape_rscript(
        tmp_path / "fake-ape-rscript",
        summary_overrides={"rooted": False},
    )

    report = run_ape_parity_cases(
        case_ids=["read-tree-example-rooted"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "ape-parity-failures",
    )

    assert report.all_passed is False
    observation = report.observations[0]
    assert observation.status == "failed"
    assert observation.mismatch_reason == "summary_mismatch"
    assert observation.reproducible_artifact_root is not None
    artifact_root = observation.reproducible_artifact_root
    assert artifact_root.exists()
    comparison_payload = json.loads(
        (artifact_root / "comparison.json").read_text(encoding="utf-8")
    )
    assert comparison_payload["mismatch_reason"] == "summary_mismatch"
    observed_summary = json.loads(
        (artifact_root / "reference-summary.observed.json").read_text(encoding="utf-8")
    )
    assert observed_summary["rooted"] is False
    assert (artifact_root / "bijux-normalized-tree.nwk").exists()


def test_run_ape_parity_cases_marks_missing_rscript_as_skipped(tmp_path: Path) -> None:
    report = run_ape_parity_cases(
        case_ids=["read-tree-example-rooted"],
        rscript_executable=str(tmp_path / "missing-rscript"),
        failure_root=tmp_path / "ape-parity-failures",
    )

    observation = report.observations[0]
    assert observation.status == "skipped"
    assert observation.mismatch_reason == "rscript_unavailable"
    assert observation.reproducible_artifact_root is not None
    assert observation.reproducible_artifact_root.exists()


def test_write_ape_parity_tables_writes_summary_and_observations(tmp_path: Path) -> None:
    rscript = fake_ape_rscript(tmp_path / "fake-ape-rscript")
    report = run_ape_parity_cases(
        rscript_executable=str(rscript),
        failure_root=tmp_path / "ape-parity-failures",
    )
    summary_path = tmp_path / "ape-parity-summary.tsv"
    observation_path = tmp_path / "ape-parity-observations.tsv"

    write_ape_parity_summary_table(summary_path, report)
    write_ape_parity_observation_table(observation_path, report)

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0] == (
        "function_name\tcase_count\tpassed_case_count\tfailed_case_count\tskipped_case_count"
    )
    with observation_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 2
    assert rows[0]["function_name"] == "ape::read.tree"
    assert rows[0]["status"] == "passed"
    assert rows[0]["bijux_version"]


def test_run_ape_parity_cases_records_live_environment_status(tmp_path: Path) -> None:
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is not available")
    if not _r_package_available(rscript, "jsonlite"):
        pytest.skip("jsonlite is required for live ape parity validation")

    report = run_ape_parity_cases(
        case_ids=["read-tree-example-rooted"],
        rscript_executable=rscript,
        failure_root=tmp_path / "ape-parity-failures",
    )

    observation = report.observations[0]
    if _r_package_available(rscript, "ape"):
        assert observation.status == "passed"
        assert observation.ape_version
    else:
        assert observation.status == "skipped"
        assert observation.mismatch_reason == "ape_package_unavailable"
        assert observation.reproducible_artifact_root is not None
