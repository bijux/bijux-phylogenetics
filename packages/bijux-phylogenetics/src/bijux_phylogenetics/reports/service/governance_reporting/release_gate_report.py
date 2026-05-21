from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.validation import (
    build_level_one_release_gate_report,
)

from ..artifacts import report_sidecar_path, section, write_machine_manifest
from ..ledger import sha256
from ..models import ReleaseGateReportBuildResult


def render_level_one_release_gate_report(
    *,
    out_path: Path,
    fixtures_root: Path | None = None,
) -> ReleaseGateReportBuildResult:
    """Render the Level 1 release gate for the checked-in workflow fixtures."""
    release_gate = build_level_one_release_gate_report(fixtures_root=fixtures_root)
    title = "Bijux Level 1 Release Gate Report"
    reviewer_summary = [
        f"gate decision: {release_gate.gate.decision}",
        f"dataset readiness: {release_gate.dataset_readiness_decision}",
        f"retained taxa: {len(release_gate.gate.retained_taxa)}, excluded taxa: {len(release_gate.gate.excluded_taxa)}",
    ]
    sections = [
        section("reviewer-summary", reviewer_summary),
        section("gate-decision", asdict(release_gate.gate)),
        section(
            "dataset-readiness",
            {
                "decision": release_gate.dataset_readiness_decision,
                "blockers": release_gate.dataset_blockers,
                "warnings": release_gate.dataset_warnings,
            },
        ),
        section(
            "taxon-loss-traceability",
            {
                "first_loss_stage": release_gate.taxon_first_loss_stage,
                "exclusion_causes": release_gate.exclusion_causes,
            },
        ),
        section("workflow-validation", asdict(release_gate.validation)),
        section("limitations", release_gate.validation.limitations),
    ]
    fixture_paths = sorted(
        {
            path
            for suite in release_gate.validation.suites
            for fixture in suite.fixtures
            for path in fixture.fixture_paths
        }
    )
    machine_manifest = {
        "report_kind": "release-gate",
        "title": title,
        "input_paths": [str(path) for path in fixture_paths],
        "input_checksums": {
            str(path): sha256(path) for path in fixture_paths if path.exists()
        },
        "sections": [name for name, _ in sections],
        "metrics": {
            "retained_taxa": len(release_gate.gate.retained_taxa),
            "excluded_taxa": len(release_gate.gate.excluded_taxa),
            "blocked_analysis_count": len(release_gate.gate.blocked_analyses),
        },
        "reviewer_summary": reviewer_summary,
        "limitations": release_gate.validation.limitations,
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
    return ReleaseGateReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="release-gate",
        title=title,
        release_gate=release_gate,
        machine_manifest=machine_manifest,
    )
