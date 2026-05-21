from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

from .support.rabies_cross_host_geography_package import (
    build_stub_rabies_cross_host_geography_package,
    refresh_stub_rabies_cross_host_geography_package,
)


def test_cli_report_package_comparison_writes_comparison_bundle(
    tmp_path: Path,
    monkeypatch,
    capsys,
) -> None:
    left = build_stub_rabies_cross_host_geography_package(
        tmp_path / "left", monkeypatch
    )
    right = build_stub_rabies_cross_host_geography_package(
        tmp_path / "right",
        monkeypatch,
    )
    right.workflow_bundle.selected_model = "GTR+F+I"
    right.workflow_bundle.scientific_findings_path.write_text(
        "finding_id\tquestion\tclaim\tevidence\tcaution\tsource_artifact\n"
        "comparative_longitude\tquestion\tchanged claim\tevidence\tcaution\tscientific-findings.tsv\n",
        encoding="utf-8",
    )
    right = refresh_stub_rabies_cross_host_geography_package(right)
    output_root = tmp_path / "comparison"

    exit_code = main(
        [
            "report",
            "package-comparison",
            str(left.package_manifest_path),
            str(right.package_manifest_path),
            "--out-dir",
            str(output_root),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["report_kind"] == "rabies_cross_host_geography_package"
    assert payload["metrics"]["dataset_id"] == "rabies_cross_host_geography_panel"
    assert payload["metrics"]["overall_comparison_status"] == "risk"
    assert payload["metrics"]["changed_artifact_count"] > 0
    assert payload["metrics"]["scientific_finding_difference_count"] == 1
    assert payload["data"]["artifact_table_path"] == str(
        output_root / "publication-package-comparison-artifacts.tsv"
    )
    assert payload["data"]["check_table_path"] == str(
        output_root / "publication-package-comparison-checks.tsv"
    )
    assert payload["data"]["summary_path"] == str(
        output_root / "publication-package-comparison-summary.json"
    )
    assert payload["data"]["report_path"] == str(
        output_root / "publication-package-comparison-report.html"
    )
