from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

from .support.rabies_cross_host_geography_package import (
    build_stub_rabies_cross_host_geography_package,
)


def test_cli_report_package_revalidation_writes_revalidation_bundle(
    tmp_path: Path,
    monkeypatch,
    capsys,
) -> None:
    package_result = build_stub_rabies_cross_host_geography_package(
        tmp_path / "package",
        monkeypatch,
    )
    output_root = tmp_path / "revalidation"

    exit_code = main(
        [
            "report",
            "package-revalidation",
            str(package_result.package_manifest_path),
            "--out-dir",
            str(output_root),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["report_kind"] == "rabies_cross_host_geography_package"
    assert payload["metrics"]["all_original_artifacts_match"] is True
    assert payload["metrics"]["overall_revalidation_status"] == "pass"
    assert payload["metrics"]["unexpected_file_count"] == 0
    assert payload["metrics"]["blocked_check_count"] == 0
    assert payload["data"]["artifact_table_path"] == str(
        output_root / "publication-package-revalidation-artifacts.tsv"
    )
    assert payload["data"]["check_table_path"] == str(
        output_root / "publication-package-revalidation-checks.tsv"
    )
    assert payload["data"]["summary_path"] == str(
        output_root / "publication-package-revalidation-summary.json"
    )
    assert payload["data"]["report_path"] == str(
        output_root / "publication-package-revalidation-report.html"
    )
