from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.validation import build_gamma_invariant_family_coverage

from ..artifacts import report_sidecar_path, section, write_machine_manifest
from ..models import GammaInvariantMixtureFamilyCoverageReportBuildResult


def render_gamma_invariant_family_coverage_report(
    *,
    out_path: Path,
) -> GammaInvariantMixtureFamilyCoverageReportBuildResult:
    """Render the governed gamma-plus-invariant family coverage report."""

    coverage = build_gamma_invariant_family_coverage()
    passed_family_count = sum(1 for family in coverage.families if family.passed)
    title = "Bijux Gamma Invariant Family Coverage Report"
    reviewer_summary = [
        f"family count: {coverage.family_count}",
        f"model count: {coverage.model_count}",
        f"passed families: {passed_family_count}",
        f"passed models: {coverage.passed_model_count}",
        "covered families: "
        + ", ".join(family.family_name for family in coverage.families),
    ]
    sections = [
        section("reviewer-summary", reviewer_summary),
        section(
            "family-overview",
            {
                "family_count": coverage.family_count,
                "model_count": coverage.model_count,
                "passed": coverage.passed,
                "passed_family_count": passed_family_count,
                "passed_model_count": coverage.passed_model_count,
            },
        ),
        section("families", [asdict(family) for family in coverage.families]),
        section("known-limitations", coverage.limitations),
    ]
    machine_manifest = {
        "report_kind": "gamma-invariant-family-coverage",
        "title": title,
        "input_paths": [],
        "input_checksums": {},
        "sections": [name for name, _ in sections],
        "metrics": {
            "goal_id": coverage.goal_id,
            "family_count": coverage.family_count,
            "model_count": coverage.model_count,
            "passed": coverage.passed,
            "passed_family_count": passed_family_count,
            "passed_model_count": coverage.passed_model_count,
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
    return GammaInvariantMixtureFamilyCoverageReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="gamma-invariant-family-coverage",
        title=title,
        coverage=coverage,
        machine_manifest=machine_manifest,
    )
