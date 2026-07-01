from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.reports.service import (
    render_invariant_mixture_family_coverage_report,
)


def test_render_invariant_mixture_family_coverage_report_writes_family_sections(
    tmp_path: Path,
) -> None:
    output = tmp_path / "invariant-mixture-family-coverage.html"
    result = render_invariant_mixture_family_coverage_report(
        out_path=output,
    )

    text = output.read_text(encoding="utf-8")
    assert result.report_kind == "invariant-mixture-family-coverage"
    assert result.machine_manifest["sections"] == [
        "reviewer-summary",
        "family-overview",
        "families",
        "known-limitations",
    ]
    assert result.machine_manifest["metrics"]["goal_id"] == 29
    assert result.machine_manifest["metrics"]["family_count"] == 5
    assert result.machine_manifest["metrics"]["passed"] is True
    assert result.machine_manifest["metrics"]["passed_family_count"] == 5
    assert result.machine_manifest_path.exists()
    assert "Bijux Invariant Mixture Family Coverage Report" in text
    assert "protein-empirical-wag" in text
