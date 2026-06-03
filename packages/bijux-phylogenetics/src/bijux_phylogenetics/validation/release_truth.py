from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from defusedxml import ElementTree

from bijux_phylogenetics.benchmark import (
    LargeDatasetStressSuiteReport,
    benchmark_large_dataset_stress_suite,
)
from bijux_phylogenetics.datasets.catalog import (
    PublicDatasetSurface,
    list_flagship_dataset_surfaces,
)
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    MethodTierAssessment,
    method_tier_warnings,
    release_method_tier_inventory,
)
from bijux_phylogenetics.parity import (
    ReferenceParityReport,
    validate_reference_parity_examples,
)
from bijux_phylogenetics.validation.reference import (
    CoreWorkflowValidationReport,
    LevelOneReleaseGateReport,
    build_core_workflow_validation_report,
    build_level_one_release_gate_report,
)


@dataclass(frozen=True, slots=True)
class PytestSessionSummary:
    """One parsed pytest JUnit session summary."""

    source_path: Path
    suite_kind: str
    suite_name: str
    total_tests: int
    passed_tests: int
    failed_tests: int
    skipped_tests: int
    error_tests: int


@dataclass(frozen=True, slots=True)
class PytestSessionAggregate:
    """Aggregate pytest session counts for one release-report test lane."""

    suite_kind: str
    session_count: int
    source_paths: list[Path]
    total_tests: int
    passed_tests: int
    failed_tests: int
    skipped_tests: int
    error_tests: int


@dataclass(slots=True)
class ReleaseTruthReport:
    """One machine-produced report of what the current release does and does not cover."""

    total_test_sessions: list[PytestSessionSummary]
    total_tests: PytestSessionAggregate
    real_engine_test_sessions: list[PytestSessionSummary]
    real_engine_tests: PytestSessionAggregate
    workflow_validation: CoreWorkflowValidationReport
    release_gate: LevelOneReleaseGateReport
    reference_parity: ReferenceParityReport
    stress_suite: LargeDatasetStressSuiteReport
    method_inventory: list[MethodTierAssessment]
    supported_workflows: list[MethodTierAssessment]
    experimental_workflows: list[MethodTierAssessment]
    advisory_workflows: list[MethodTierAssessment]
    parser_only_workflows: list[MethodTierAssessment]
    flagship_datasets: list[PublicDatasetSurface]
    known_limitations: list[str]


def _iter_junit_testsuites(root: Any) -> list[Any]:
    if root.tag == "testsuite":
        return [root]
    suites = list(root.findall("testsuite"))
    if suites:
        return suites
    nested = list(root.findall(".//testsuite"))
    if nested:
        return nested
    raise ValueError("pytest JUnit report does not contain any testsuite elements")


def parse_pytest_junit_report(path: Path, *, suite_kind: str) -> PytestSessionSummary:
    """Parse one pytest JUnit XML report into stable pass/fail/skip counts."""
    if not path.exists():
        raise FileNotFoundError(path)
    root = ElementTree.fromstring(path.read_text(encoding="utf-8"))
    suites = _iter_junit_testsuites(root)
    total_tests = sum(int(suite.attrib.get("tests", "0")) for suite in suites)
    failed_tests = sum(int(suite.attrib.get("failures", "0")) for suite in suites)
    skipped_tests = sum(int(suite.attrib.get("skipped", "0")) for suite in suites)
    error_tests = sum(int(suite.attrib.get("errors", "0")) for suite in suites)
    passed_tests = max(total_tests - failed_tests - skipped_tests - error_tests, 0)
    suite_name = root.attrib.get("name") or path.stem
    return PytestSessionSummary(
        source_path=path,
        suite_kind=suite_kind,
        suite_name=suite_name,
        total_tests=total_tests,
        passed_tests=passed_tests,
        failed_tests=failed_tests,
        skipped_tests=skipped_tests,
        error_tests=error_tests,
    )


def aggregate_pytest_junit_reports(
    paths: list[Path],
    *,
    suite_kind: str,
) -> tuple[list[PytestSessionSummary], PytestSessionAggregate]:
    """Aggregate one or more pytest JUnit reports for a release test lane."""
    if not paths:
        raise ValueError(f"{suite_kind} requires at least one pytest JUnit report")
    sessions = [
        parse_pytest_junit_report(path, suite_kind=suite_kind) for path in paths
    ]
    return sessions, PytestSessionAggregate(
        suite_kind=suite_kind,
        session_count=len(sessions),
        source_paths=[session.source_path for session in sessions],
        total_tests=sum(session.total_tests for session in sessions),
        passed_tests=sum(session.passed_tests for session in sessions),
        failed_tests=sum(session.failed_tests for session in sessions),
        skipped_tests=sum(session.skipped_tests for session in sessions),
        error_tests=sum(session.error_tests for session in sessions),
    )


def _group_method_inventory(
    inventory: list[MethodTierAssessment],
) -> tuple[
    list[MethodTierAssessment],
    list[MethodTierAssessment],
    list[MethodTierAssessment],
    list[MethodTierAssessment],
]:
    return (
        [item for item in inventory if item.tier == "supported"],
        [item for item in inventory if item.tier == "experimental"],
        [item for item in inventory if item.tier == "advisory"],
        [item for item in inventory if item.tier == "parser-only"],
    )


def _release_limitations(
    release_gate: LevelOneReleaseGateReport,
    parity: ReferenceParityReport,
    stress_suite: LargeDatasetStressSuiteReport,
    method_inventory: list[MethodTierAssessment],
) -> list[str]:
    limitations: list[str] = []
    for item in release_gate.validation.limitations:
        if item not in limitations:
            limitations.append(item)
    for item in release_gate.dataset_blockers:
        if item not in limitations:
            limitations.append(item)
    for item in release_gate.dataset_warnings:
        if item not in limitations:
            limitations.append(item)
    for item in parity.limitations:
        if item not in limitations:
            limitations.append(item)
    for item in stress_suite.limitations:
        if item not in limitations:
            limitations.append(item)
    for assessment in method_inventory:
        for warning in method_tier_warnings(assessment):
            if warning not in limitations:
                limitations.append(warning)
    return limitations


def build_release_truth_report(
    *,
    test_report_paths: list[Path],
    real_engine_test_report_paths: list[Path],
    fixtures_root: Path | None = None,
    include_extended_parity: bool = False,
    stress_tier: str = "small",
) -> ReleaseTruthReport:
    """Build one governed release-truth report from actual test and workflow outputs."""
    total_sessions, total_tests = aggregate_pytest_junit_reports(
        test_report_paths,
        suite_kind="total-tests",
    )
    real_engine_sessions, real_engine_tests = aggregate_pytest_junit_reports(
        real_engine_test_report_paths,
        suite_kind="real-engine-tests",
    )
    workflow_validation = build_core_workflow_validation_report(
        fixtures_root=fixtures_root
    )
    release_gate = build_level_one_release_gate_report(fixtures_root=fixtures_root)
    reference_parity = validate_reference_parity_examples(
        include_extended=include_extended_parity
    )
    stress_suite = benchmark_large_dataset_stress_suite(tier=stress_tier)
    method_inventory = release_method_tier_inventory()
    (
        supported_workflows,
        experimental_workflows,
        advisory_workflows,
        parser_only_workflows,
    ) = _group_method_inventory(method_inventory)
    return ReleaseTruthReport(
        total_test_sessions=total_sessions,
        total_tests=total_tests,
        real_engine_test_sessions=real_engine_sessions,
        real_engine_tests=real_engine_tests,
        workflow_validation=workflow_validation,
        release_gate=release_gate,
        reference_parity=reference_parity,
        stress_suite=stress_suite,
        method_inventory=method_inventory,
        supported_workflows=supported_workflows,
        experimental_workflows=experimental_workflows,
        advisory_workflows=advisory_workflows,
        parser_only_workflows=parser_only_workflows,
        flagship_datasets=list_flagship_dataset_surfaces(),
        known_limitations=_release_limitations(
            release_gate,
            reference_parity,
            stress_suite,
            method_inventory,
        ),
    )
