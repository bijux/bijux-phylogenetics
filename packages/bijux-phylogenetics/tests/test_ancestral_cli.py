from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.cli import main


FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_ancestral_continuous_cli_can_export_table(tmp_path: Path, capsys) -> None:
    table_path = tmp_path / "ancestral.tsv"
    exit_code = main(
        [
            "ancestral",
            "continuous",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--table-out",
            str(table_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["model"] == "brownian"
    assert "estimate\tstandard_error" in table_path.read_text(encoding="utf-8")


def test_ancestral_discrete_cli_reports_sparse_state_warning(capsys) -> None:
    exit_code = main(
        [
            "ancestral",
            "discrete",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_ancestral_sparse.tsv")),
            "--trait",
            "habitat",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["warnings"] == [
        "one or more discrete states are represented by fewer than two taxa and should be interpreted cautiously"
    ]


def test_ancestral_render_cli_writes_svg_with_internal_annotations(tmp_path: Path, capsys) -> None:
    output = tmp_path / "ancestral.svg"
    exit_code = main(
        [
            "ancestral",
            "render",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--kind",
            "continuous",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["rendered_internal_annotation_count"] == 3
    assert 'class="internal-annotation-label"' in output.read_text(encoding="utf-8")


def test_ancestral_report_cli_writes_html_and_svg(tmp_path: Path, capsys) -> None:
    output = tmp_path / "ancestral-report.html"
    exit_code = main(
        [
            "ancestral",
            "report",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--kind",
            "continuous",
            "--compare-model",
            "ou",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["report_kind"] == "ancestral-state"
    assert output.exists()
    assert output.with_suffix(".svg").exists()
