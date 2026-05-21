from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.validation as validation_api
from bijux_phylogenetics.validation import (
    aggregate_pytest_junit_reports,
    build_release_truth_report,
    parse_pytest_junit_report,
)


def _write_junit_report(
    path: Path,
    *,
    suite_name: str,
    tests: int,
    failures: int,
    skipped: int,
    errors: int = 0,
) -> Path:
    path.write_text(
        (
            '<testsuites name="pytest">'
            f'<testsuite name="{suite_name}" tests="{tests}" failures="{failures}" errors="{errors}" skipped="{skipped}" />'
            "</testsuites>\n"
        ),
        encoding="utf-8",
    )
    return path


def test_parse_pytest_junit_report_counts_pass_fail_skip_and_error(
    tmp_path: Path,
) -> None:
    report_path = _write_junit_report(
        tmp_path / "full.xml",
        suite_name="full-suite",
        tests=12,
        failures=2,
        skipped=3,
        errors=1,
    )

    summary = parse_pytest_junit_report(report_path, suite_kind="total-tests")

    assert summary.suite_kind == "total-tests"
    assert summary.suite_name == "pytest"
    assert summary.total_tests == 12
    assert summary.passed_tests == 6
    assert summary.failed_tests == 2
    assert summary.skipped_tests == 3
    assert summary.error_tests == 1


@pytest.mark.slow
def test_aggregate_and_build_release_truth_report_uses_actual_runtime_surfaces(
    tmp_path: Path,
) -> None:
    total_left = _write_junit_report(
        tmp_path / "unit.xml",
        suite_name="unit",
        tests=10,
        failures=1,
        skipped=2,
    )
    total_right = _write_junit_report(
        tmp_path / "integration.xml",
        suite_name="integration",
        tests=8,
        failures=0,
        skipped=1,
    )
    real_engine = _write_junit_report(
        tmp_path / "real-engine.xml",
        suite_name="real-engine",
        tests=5,
        failures=0,
        skipped=2,
    )

    sessions, aggregate = aggregate_pytest_junit_reports(
        [total_left, total_right],
        suite_kind="total-tests",
    )
    assert len(sessions) == 2
    assert aggregate.total_tests == 18
    assert aggregate.passed_tests == 14
    assert aggregate.failed_tests == 1
    assert aggregate.skipped_tests == 3

    report = build_release_truth_report(
        test_report_paths=[total_left, total_right],
        real_engine_test_report_paths=[real_engine],
        include_extended_parity=False,
        stress_tier="small",
    )

    assert report.total_tests.total_tests == 18
    assert report.total_tests.passed_tests == 14
    assert report.real_engine_tests.total_tests == 5
    assert report.real_engine_tests.passed_tests == 3
    assert any(
        workflow.surface == "fasta-to-tree" for workflow in report.supported_workflows
    )
    assert any(
        workflow.surface == "phylogenetic-logistic"
        for workflow in report.experimental_workflows
    )
    assert any(
        dataset.demo_command == "rabies-cross-host-geography-panel"
        for dataset in report.flagship_datasets
    )
    assert report.reference_parity.case_count > 0
    assert len(report.stress_suite.observations) == 5
    assert report.known_limitations


def test_public_runtime_exports_release_truth_builder() -> None:
    assert validation_api.build_release_truth_report is build_release_truth_report
