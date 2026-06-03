from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.reports.publication.tree import build_tree_report_package


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_cli_report_reviewer_audit_checklist_writes_tsv(tmp_path: Path, capsys) -> None:
    package_dir = tmp_path / "tree-report-package"
    package_result = build_tree_report_package(
        tree_fixture("example_tree_support_left.nwk"),
        out_dir=package_dir,
    )
    output_path = tmp_path / "reviewer-audit-checklist.tsv"

    exit_code = main(
        [
            "report",
            "reviewer-audit-checklist",
            str(package_result.manifest_path),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["report_kind"] == "tree_package"
    assert payload["metrics"]["item_count"] == 5
    assert payload["metrics"]["blocked_item_count"] == 0
    assert output_path.exists()
    lines = output_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "section\tstatus\tsummary\tevidence\tartifact_paths"
    assert any(line.startswith("support_surface\t") for line in lines[1:])
