from __future__ import annotations

import json
from pathlib import Path
from typing import cast

import pytest

from bijux_phylogenetics_dev.quality.runtime_budget_governance import (
    build_runtime_budget_report,
    check_runtime_budgets,
    parse_runtime_budget_observations,
)

REPO_ROOT = Path(__file__).resolve().parents[3]


def _write_log(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_parse_runtime_budget_observations_filters_long_call_rows(
    tmp_path: Path,
) -> None:
    log_path = tmp_path / "pytest.log"
    _write_log(
        log_path,
        [
            "============================== slowest durations ===============================",
            "1217.48s call     packages/bijux-phylogenetics/tests/parity/geiger/test_geiger_parity_cli.py::test_parity_cli_runs_live_geiger_harness_and_writes_tables",
            "130.00s setup    packages/bijux-phylogenetics/tests/parity/geiger/test_geiger_parity_cli.py::test_parity_cli_runs_live_geiger_harness_and_writes_tables",
            "99.98s call     packages/bijux-phylogenetics/tests/workflows/pipelines/test_runtime.py::test_cli_evidence_book_rerun_json_output",
        ],
    )

    observations = parse_runtime_budget_observations(log_path)

    assert len(observations) == 1
    assert observations[0].nodeid.endswith(
        "test_parity_cli_runs_live_geiger_harness_and_writes_tables"
    )
    assert observations[0].phase == "call"
    assert observations[0].duration_seconds == 1217.48


def test_build_runtime_budget_report_classifies_governed_long_running_cases(
    tmp_path: Path,
) -> None:
    log_path = tmp_path / "pytest.log"
    _write_log(
        log_path,
        [
            "============================== slowest durations ===============================",
            "1217.48s call     packages/bijux-phylogenetics/tests/parity/geiger/test_geiger_parity_cli.py::test_parity_cli_runs_live_geiger_harness_and_writes_tables",
            "138.42s call     packages/bijux-phylogenetics/tests/workflows/benchmarks/test_native_maximum_likelihood_accuracy_benchmark.py::test_benchmark_native_maximum_likelihood_accuracy_reports_truth_and_wrapper_lanes",
            "107.52s call     packages/bijux-phylogenetics/tests/workflows/pipelines/test_runtime.py::test_cli_evidence_book_build_selected_evidence_json_output",
            "101.00s call     packages/bijux-phylogenetics-dev/tests/test_publish_readiness.py::test_build_publish_readiness_report_exposes_repository_blockers",
        ],
    )

    payload = build_runtime_budget_report(REPO_ROOT, log_path)
    counts_by_classification = cast(
        dict[str, int],
        payload["counts_by_classification"],
    )
    observations = cast(list[dict[str, object]], payload["observations"])

    assert payload["issue_count"] == 0
    assert payload["observation_count"] == 4
    assert payload["classified_observation_count"] == 4
    assert counts_by_classification["split"] == 2
    assert counts_by_classification["optimize"] == 1
    assert counts_by_classification["frozen"] == 1
    observations_by_nodeid = {cast(str, row["nodeid"]): row for row in observations}
    assert (
        observations_by_nodeid[
            "packages/bijux-phylogenetics/tests/parity/geiger/test_geiger_parity_cli.py::test_parity_cli_runs_live_geiger_harness_and_writes_tables"
        ]["classification"]
        == "split"
    )
    assert (
        observations_by_nodeid[
            "packages/bijux-phylogenetics/tests/workflows/benchmarks/test_native_maximum_likelihood_accuracy_benchmark.py::test_benchmark_native_maximum_likelihood_accuracy_reports_truth_and_wrapper_lanes"
        ]["classification"]
        == "optimize"
    )
    assert (
        observations_by_nodeid[
            "packages/bijux-phylogenetics-dev/tests/test_publish_readiness.py::test_build_publish_readiness_report_exposes_repository_blockers"
        ]["classification"]
        == "frozen"
    )
    assert all(row["slow_marked"] is True for row in observations)


def test_build_runtime_budget_report_flags_unclassified_long_running_cases(
    tmp_path: Path,
) -> None:
    log_path = tmp_path / "pytest.log"
    _write_log(
        log_path,
        [
            "============================== slowest durations ===============================",
            "150.00s call     packages/bijux-phylogenetics/tests/test_unknown_surface.py::test_unclassified_case",
        ],
    )

    payload = build_runtime_budget_report(REPO_ROOT, log_path)
    issues = cast(list[dict[str, object]], payload["issues"])

    assert payload["issue_count"] == 1
    issue = issues[0]
    assert issue["code"] == "missing-runtime-budget-classification"
    assert cast(str, issue["nodeid"]).endswith(
        "test_unknown_surface.py::test_unclassified_case"
    )


def test_check_runtime_budgets_writes_json_and_raises_on_drift(
    tmp_path: Path,
) -> None:
    log_path = tmp_path / "pytest.log"
    json_out = tmp_path / "runtime-budget.json"
    _write_log(
        log_path,
        [
            "============================== slowest durations ===============================",
            "150.00s call     packages/bijux-phylogenetics/tests/test_unknown_surface.py::test_unclassified_case",
        ],
    )

    with pytest.raises(SystemExit, match="runtime budget audit failed"):
        check_runtime_budgets(REPO_ROOT, log_path, json_out=json_out)

    payload = json.loads(json_out.read_text(encoding="utf-8"))
    assert payload["issue_count"] == 1
