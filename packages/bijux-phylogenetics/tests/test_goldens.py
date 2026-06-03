from __future__ import annotations

from pathlib import Path
import re

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.reports.service import render_tree_report

from .support.scientific_output_assertions import (
    assert_scientific_outputs_equivalent,
)

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
EXPECTED = Path(__file__).parent / "fixtures" / "expected"


def _normalize_dynamic_output(text: str) -> str:
    return re.sub(r'"version": "[^"]+"', '"version": "<version>"', text)


def test_inspect_cli_matches_golden_output(capsys, monkeypatch) -> None:
    monkeypatch.chdir(PACKAGE_ROOT)
    exit_code = main(["inspect", "tests/fixtures/trees/example_tree.nwk", "--json"])
    captured = capsys.readouterr()
    assert exit_code == 0
    expected = (EXPECTED / "inspect_tree.json").read_text(encoding="utf-8")
    assert _normalize_dynamic_output(captured.out) == _normalize_dynamic_output(
        expected
    )


def test_tree_report_matches_golden_html(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(PACKAGE_ROOT)
    output = tmp_path / "tree-report.html"
    render_tree_report(
        tree_path=Path("tests/fixtures/trees/example_tree.nwk"), out_path=output
    )
    report = assert_scientific_outputs_equivalent(EXPECTED / "tree_report.html", output)
    assert report.compared_file_count == 1
