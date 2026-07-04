from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.reports.service import (
    render_shared_continuous_fixture_catalog_report,
)


def test_render_shared_continuous_fixture_catalog_report_writes_fixture_sections(
    tmp_path: Path,
) -> None:
    output = tmp_path / "shared-continuous-fixture-catalog.html"
    result = render_shared_continuous_fixture_catalog_report(
        out_path=output,
    )

    text = output.read_text(encoding="utf-8")
    assert result.report_kind == "shared-continuous-fixture-catalog"
    assert result.machine_manifest["sections"] == [
        "reviewer-summary",
        "fixture-overview",
        "fixtures",
        "known-limitations",
    ]
    assert result.machine_manifest["metrics"]["goal_id"] == 252
    assert result.machine_manifest["metrics"]["fixture_count"] == 6
    assert result.machine_manifest["metrics"]["live_case_count"] >= 6
    assert result.machine_manifest_path.exists()
    assert "Bijux Shared Continuous Fixture Catalog Report" in text


def test_cli_report_shared_continuous_fixture_catalog_json_output_uses_catalog_contract(
    tmp_path: Path,
    capsys,
) -> None:
    output = tmp_path / "shared-continuous-fixture-catalog.html"
    exit_code = main(
        [
            "report",
            "shared-continuous-fixture-catalog",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    fixtures = payload["data"]["fixture_catalog"]["fixtures"]
    missing_values = next(
        fixture
        for fixture in fixtures
        if fixture["fixture_id"] == "geiger_continuous_missing_values_twenty_four_taxa"
    )

    assert exit_code == 0
    assert payload["command"] == "report"
    assert payload["outputs"] == [str(output), str(output.with_suffix(".json"))]
    assert payload["data"]["report_kind"] == "shared-continuous-fixture-catalog"
    assert payload["metrics"]["goal_id"] == 252
    assert payload["metrics"]["fixture_count"] == 6
    assert "trend" in missing_values["supported_model_names"]
    assert "trend" not in missing_values["covered_model_names"]
