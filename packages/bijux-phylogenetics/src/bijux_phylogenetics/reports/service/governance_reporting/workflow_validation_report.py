from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.validation import (
    build_core_workflow_validation_report,
)

from ..artifacts import report_sidecar_path, section, write_machine_manifest
from ..ledger import sha256
from ..models import WorkflowValidationReportBuildResult


def render_workflow_validation_report(
    *,
    out_path: Path,
    fixtures_root: Path | None = None,
) -> WorkflowValidationReportBuildResult:
    """Render the Level 1 workflow validation fixture report."""
    validation = build_core_workflow_validation_report(fixtures_root=fixtures_root)
    title = "Bijux Core Workflow Validation Report"
    reviewer_summary = [
        f"fixture checks passed: {validation.passed_fixture_count}/{validation.total_fixture_count}",
        f"validated workflow surfaces: {len(validation.workflows)}",
        f"known failure-gallery cases: {len(validation.failure_gallery)}",
    ]
    sections = [
        section("reviewer-summary", reviewer_summary),
        section(
            "validation-overview",
            {
                "total_fixture_count": validation.total_fixture_count,
                "passed_fixture_count": validation.passed_fixture_count,
                "failed_fixture_count": validation.failed_fixture_count,
            },
        ),
        section("validation-suites", [asdict(suite) for suite in validation.suites]),
        section("workflow-coverage", [asdict(row) for row in validation.workflows]),
        section("failure-gallery", [asdict(row) for row in validation.failure_gallery]),
        section(
            "maturity-classification",
            [asdict(row) for row in validation.maturity_classifications],
        ),
        section("limitations", validation.limitations),
    ]
    fixture_paths = sorted(
        {
            path
            for suite in validation.suites
            for fixture in suite.fixtures
            for path in fixture.fixture_paths
        }
    )
    machine_manifest = {
        "report_kind": "workflow-validation",
        "title": title,
        "input_paths": [str(path) for path in fixture_paths],
        "input_checksums": {
            str(path): sha256(path) for path in fixture_paths if path.exists()
        },
        "sections": [name for name, _ in sections],
        "metrics": {
            "total_fixture_count": validation.total_fixture_count,
            "passed_fixture_count": validation.passed_fixture_count,
            "workflow_count": len(validation.workflows),
        },
        "reviewer_summary": reviewer_summary,
        "limitations": validation.limitations,
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
    return WorkflowValidationReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="workflow-validation",
        title=title,
        validation=validation,
        machine_manifest=machine_manifest,
    )
