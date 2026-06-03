from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.validation import build_release_truth_report

from ..artifacts import report_sidecar_path, section, write_machine_manifest
from ..ledger import sha256
from ..models import ReleaseTruthReportBuildResult


def render_release_truth_report(
    *,
    out_path: Path,
    test_report_paths: list[Path],
    real_engine_test_report_paths: list[Path],
    fixtures_root: Path | None = None,
    include_extended_parity: bool = False,
    stress_tier: str = "small",
) -> ReleaseTruthReportBuildResult:
    """Render one machine-produced report of the current release truth surface."""
    release_truth = build_release_truth_report(
        test_report_paths=test_report_paths,
        real_engine_test_report_paths=real_engine_test_report_paths,
        fixtures_root=fixtures_root,
        include_extended_parity=include_extended_parity,
        stress_tier=stress_tier,
    )
    title = "Bijux Release Truth Report"
    reviewer_summary = [
        f"total tests: {release_truth.total_tests.passed_tests} passed, {release_truth.total_tests.failed_tests} failed, {release_truth.total_tests.skipped_tests} skipped",
        f"real-engine tests: {release_truth.real_engine_tests.passed_tests} passed, {release_truth.real_engine_tests.failed_tests} failed, {release_truth.real_engine_tests.skipped_tests} skipped",
        f"supported workflows: {len(release_truth.supported_workflows)}, experimental workflows: {len(release_truth.experimental_workflows)}",
        f"flagship datasets: {len(release_truth.flagship_datasets)}, reference parity cases: {release_truth.reference_parity.case_count}, stress workloads: {len(release_truth.stress_suite.observations)}",
        f"release gate decision: {release_truth.release_gate.gate.decision}",
    ]
    sections = [
        section("reviewer-summary", reviewer_summary),
        section("total-tests", asdict(release_truth.total_tests)),
        section("real-engine-tests", asdict(release_truth.real_engine_tests)),
        section(
            "supported-workflows",
            [asdict(item) for item in release_truth.supported_workflows],
        ),
        section(
            "experimental-workflows",
            [asdict(item) for item in release_truth.experimental_workflows],
        ),
        section(
            "advisory-workflows",
            [asdict(item) for item in release_truth.advisory_workflows],
        ),
        section(
            "parser-only-workflows",
            [asdict(item) for item in release_truth.parser_only_workflows],
        ),
        section(
            "flagship-datasets",
            [asdict(item) for item in release_truth.flagship_datasets],
        ),
        section("workflow-validation", asdict(release_truth.workflow_validation)),
        section("release-gate", asdict(release_truth.release_gate)),
        section("reference-parity", asdict(release_truth.reference_parity)),
        section("stress-suite", asdict(release_truth.stress_suite)),
        section("known-limitations", release_truth.known_limitations),
    ]
    fixture_paths = sorted(
        {
            path
            for suite in release_truth.workflow_validation.suites
            for fixture in suite.fixtures
            for path in fixture.fixture_paths
        }
    )
    input_paths = [*test_report_paths, *real_engine_test_report_paths, *fixture_paths]
    machine_manifest = {
        "report_kind": "release-truth",
        "title": title,
        "input_paths": [str(path) for path in input_paths],
        "input_checksums": {
            str(path): sha256(path) for path in input_paths if path.exists()
        },
        "sections": [name for name, _ in sections],
        "metrics": {
            "total_tests": release_truth.total_tests.total_tests,
            "total_tests_passed": release_truth.total_tests.passed_tests,
            "total_tests_failed": release_truth.total_tests.failed_tests,
            "total_tests_skipped": release_truth.total_tests.skipped_tests,
            "real_engine_tests": release_truth.real_engine_tests.total_tests,
            "real_engine_tests_passed": release_truth.real_engine_tests.passed_tests,
            "real_engine_tests_failed": release_truth.real_engine_tests.failed_tests,
            "real_engine_tests_skipped": release_truth.real_engine_tests.skipped_tests,
            "supported_workflow_count": len(release_truth.supported_workflows),
            "experimental_workflow_count": len(release_truth.experimental_workflows),
            "flagship_dataset_count": len(release_truth.flagship_datasets),
            "reference_parity_case_count": release_truth.reference_parity.case_count,
            "stress_workload_count": len(release_truth.stress_suite.observations),
        },
        "reviewer_summary": reviewer_summary,
        "limitations": release_truth.known_limitations,
    }
    machine_manifest_path = write_machine_manifest(
        report_sidecar_path(out_path), machine_manifest
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return ReleaseTruthReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="release-truth",
        title=title,
        release_truth=release_truth,
        machine_manifest=machine_manifest,
    )
