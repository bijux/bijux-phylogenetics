from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.reports.service import (
    render_shared_discrete_fixture_catalog_report,
)


def test_render_shared_discrete_fixture_catalog_report_writes_fixture_sections(
    tmp_path: Path,
) -> None:
    output = tmp_path / "shared-discrete-fixture-catalog.html"
    result = render_shared_discrete_fixture_catalog_report(
        out_path=output,
    )

    text = output.read_text(encoding="utf-8")
    assert result.report_kind == "shared-discrete-fixture-catalog"
    assert result.machine_manifest["sections"] == [
        "reviewer-summary",
        "fixture-overview",
        "fixtures",
        "known-limitations",
    ]
    assert result.machine_manifest["metrics"]["goal_id"] == 253
    assert result.machine_manifest["metrics"]["fixture_count"] == 17
    assert result.machine_manifest["metrics"]["live_case_count"] >= 17
    assert result.machine_manifest_path.exists()
    assert "Bijux Shared Discrete Fixture Catalog Report" in text


def test_cli_report_shared_discrete_fixture_catalog_json_output_uses_catalog_contract(
    tmp_path: Path,
    capsys,
) -> None:
    output = tmp_path / "shared-discrete-fixture-catalog.html"
    exit_code = main(
        [
            "report",
            "shared-discrete-fixture-catalog",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    fixtures = payload["data"]["fixture_catalog"]["fixtures"]
    sym_fixture = next(
        fixture
        for fixture in fixtures
        if fixture["fixture_id"] == "geiger_discrete_sym_three_state_twenty_four_taxa"
    )

    assert exit_code == 0
    assert payload["command"] == "report"
    assert payload["outputs"] == [str(output), str(output.with_suffix(".json"))]
    assert payload["data"]["report_kind"] == "shared-discrete-fixture-catalog"
    assert payload["metrics"]["goal_id"] == 253
    assert payload["metrics"]["fixture_count"] == 17
    assert sym_fixture["covered_model_names"] == ["SYM"]
    assert sym_fixture["trait_kind"] == "three-state-discrete"
