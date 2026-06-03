from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.reports.service import (
    render_production_scale_readiness_report,
)


@pytest.mark.slow
def test_render_production_scale_readiness_report_writes_scale_sections(
    tmp_path: Path,
) -> None:
    output = tmp_path / "production-scale-readiness.html"
    result = render_production_scale_readiness_report(
        out_path=output,
        replicates=1,
        tree_tip_counts=[8, 16],
        alignment_size_classes=[
            ("sequences-4-sites-16", 4, 16),
            ("sequences-6-sites-24", 6, 24),
        ],
        tree_set_size_classes=[
            ("trees-8-taxa-6", 8, 6),
            ("trees-12-taxa-8", 12, 8),
        ],
        stress_tiers=["small"],
    )

    text = output.read_text(encoding="utf-8")
    assert result.report_kind == "production-scale-readiness"
    assert result.machine_manifest["sections"] == [
        "reviewer-summary",
        "scale-definitions",
        "scale-coverage",
        "production-scale-readiness",
        "known-limitations",
    ]
    assert result.machine_manifest["metrics"]["goal_id"] == 225
    assert result.machine_manifest["metrics"]["small_ready_workflow_count"] >= 1
    assert result.machine_manifest["metrics"]["below_small_workflow_count"] >= 1
    assert result.machine_manifest_path.exists()
    assert "Bijux Production-Scale Readiness Report" in text


@pytest.mark.slow
def test_cli_report_production_scale_readiness_json_output_uses_scale_contract(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    output = tmp_path / "production-scale-readiness.html"
    exit_code = main(
        [
            "report",
            "production-scale-readiness",
            "--replicates",
            "1",
            "--tree-tip-count",
            "8",
            "--tree-tip-count",
            "16",
            "--sequence-count",
            "4",
            "--alignment-length",
            "16",
            "--sequence-count",
            "6",
            "--alignment-length",
            "24",
            "--posterior-tree-count",
            "8",
            "--tree-set-tip-count",
            "6",
            "--posterior-tree-count",
            "12",
            "--tree-set-tip-count",
            "8",
            "--stress-tier",
            "small",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    entries = payload["data"]["production_scale_readiness"]["entries"]
    tree_validation = next(
        entry for entry in entries if entry["workflow"] == "tree-validation"
    )

    assert exit_code == 0
    assert payload["command"] == "report"
    assert payload["outputs"] == [str(output), str(output.with_suffix(".json"))]
    assert payload["data"]["report_kind"] == "production-scale-readiness"
    assert payload["metrics"]["goal_id"] == 225
    assert payload["metrics"]["small_ready_workflow_count"] >= 1
    assert tree_validation["scale_dimensions"] == ["taxa"]
    assert tree_validation["highest_ready_scale"] == "small"
