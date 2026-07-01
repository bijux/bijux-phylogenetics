from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.validation import build_invariant_mixture_family_coverage

from ..artifacts import report_sidecar_path, section, write_machine_manifest
from ..models import InvariantMixtureFamilyCoverageReportBuildResult


def render_invariant_mixture_family_coverage_report(
    *,
    out_path: Path,
) -> InvariantMixtureFamilyCoverageReportBuildResult:
    """Render the governed invariant-mixture family coverage report."""

    coverage = build_invariant_mixture_family_coverage()
    passed_family_count = sum(1 for family in coverage.families if family.passed)
    title = "Bijux Invariant Mixture Family Coverage Report"
    reviewer_summary = [
        f"family count: {coverage.family_count}",
        f"passed families: {passed_family_count}",
        "covered families: "
        + ", ".join(family.family_name for family in coverage.families),
    ]
    sections = [
        section("reviewer-summary", reviewer_summary),
        section(
            "family-overview",
            {
                "family_count": coverage.family_count,
                "passed": coverage.passed,
                "passed_family_count": passed_family_count,
            },
        ),
        section("families", [asdict(family) for family in coverage.families]),
        section("known-limitations", coverage.limitations),
    ]
    machine_manifest = {
        "report_kind": "invariant-mixture-family-coverage",
        "title": title,
        "input_paths": [],
        "input_checksums": {},
        "sections": [name for name, _ in sections],
        "metrics": {
            "goal_id": coverage.goal_id,
            "family_count": coverage.family_count,
            "passed": coverage.passed,
            "passed_family_count": passed_family_count,
        },
        "reviewer_summary": reviewer_summary,
        "limitations": coverage.limitations,
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
    return InvariantMixtureFamilyCoverageReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="invariant-mixture-family-coverage",
        title=title,
        coverage=coverage,
        machine_manifest=machine_manifest,
    )
