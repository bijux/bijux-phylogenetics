from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.validation import build_shared_discrete_fixture_catalog

from ..artifacts import report_sidecar_path, section, write_machine_manifest
from ..models import SharedDiscreteFixtureCatalogReportBuildResult


def render_shared_discrete_fixture_catalog_report(
    *,
    out_path: Path,
) -> SharedDiscreteFixtureCatalogReportBuildResult:
    """Render the governed shared discrete fixture catalog for live geiger parity."""

    fixture_catalog = build_shared_discrete_fixture_catalog()
    title = "Bijux Shared Discrete Fixture Catalog Report"
    reviewer_summary = [
        f"fixture count: {fixture_catalog.fixture_count}",
        f"live fitdiscrete cases: {fixture_catalog.live_case_count}",
        "covered models: "
        + ", ".join(
            sorted(
                {
                    model_name
                    for fixture in fixture_catalog.fixtures
                    for model_name in fixture.covered_model_names
                }
            )
        ),
    ]
    sections = [
        section("reviewer-summary", reviewer_summary),
        section(
            "fixture-overview",
            {
                "fixture_count": fixture_catalog.fixture_count,
                "live_case_count": fixture_catalog.live_case_count,
            },
        ),
        section("fixtures", [asdict(fixture) for fixture in fixture_catalog.fixtures]),
        section("known-limitations", fixture_catalog.limitations),
    ]
    machine_manifest = {
        "report_kind": "shared-discrete-fixture-catalog",
        "title": title,
        "input_paths": [],
        "input_checksums": {},
        "sections": [name for name, _ in sections],
        "metrics": {
            "goal_id": fixture_catalog.goal_id,
            "fixture_count": fixture_catalog.fixture_count,
            "live_case_count": fixture_catalog.live_case_count,
        },
        "reviewer_summary": reviewer_summary,
        "limitations": fixture_catalog.limitations,
    }
    machine_manifest_path = write_machine_manifest(
        report_sidecar_path(out_path),
        machine_manifest,
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return SharedDiscreteFixtureCatalogReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="shared-discrete-fixture-catalog",
        title=title,
        fixture_catalog=fixture_catalog,
        machine_manifest=machine_manifest,
    )
