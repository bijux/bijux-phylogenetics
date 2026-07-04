from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.reports.service import (
    render_gamma_invariant_family_coverage_report,
)


def test_render_gamma_invariant_family_coverage_report_writes_family_sections(
    tmp_path: Path,
) -> None:
    output = tmp_path / "gamma-invariant-family-coverage.html"
    result = render_gamma_invariant_family_coverage_report(
        out_path=output,
    )

    text = output.read_text(encoding="utf-8")
    assert result.report_kind == "gamma-invariant-family-coverage"
    assert result.machine_manifest["sections"] == [
        "reviewer-summary",
        "family-overview",
        "families",
        "known-limitations",
    ]
    assert result.machine_manifest["metrics"]["goal_id"] == 30
    assert result.machine_manifest["metrics"]["family_count"] == 7
    assert result.machine_manifest["metrics"]["model_count"] == 17
    assert result.machine_manifest["metrics"]["passed"] is True
    assert result.machine_manifest["metrics"]["passed_family_count"] == 7
    assert result.machine_manifest["metrics"]["passed_model_count"] == 17
    assert result.machine_manifest_path.exists()
    assert "Bijux Gamma Invariant Family Coverage Report" in text
    assert "protein-empirical" in text
    assert "boundary_warnings" in text


def test_cli_report_gamma_invariant_family_coverage_json_output_uses_contract(
    tmp_path: Path,
    capsys,
) -> None:
    output = tmp_path / "gamma-invariant-family-coverage.html"
    exit_code = main(
        [
            "report",
            "gamma-invariant-family-coverage",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    families = payload["data"]["coverage"]["families"]
    nucleotide = next(
        family for family in families if family["family_name"] == "nucleotide"
    )

    assert exit_code == 0
    assert payload["command"] == "report"
    assert payload["outputs"] == [str(output), str(output.with_suffix(".json"))]
    assert payload["data"]["report_kind"] == "gamma-invariant-family-coverage"
    assert payload["metrics"]["goal_id"] == 30
    assert payload["metrics"]["family_count"] == 7
    assert payload["metrics"]["model_count"] == 17
    assert payload["metrics"]["passed"] is True
    assert nucleotide["covered_model_names"] == ["JC69", "K80", "F81", "HKY85", "GTR"]
    assert nucleotide["passed_model_count"] == 5
    assert nucleotide["model_cases"][0]["alpha"] == 0.8
    assert nucleotide["model_cases"][0]["boundary_warnings"] == [
        "invariant proportion hit lower search boundary"
    ]
