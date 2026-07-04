from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

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


def test_cli_report_ancestral_methods_summary_writes_continuous_metrics(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "ancestral-methods-summary.md"

    exit_code = main(
        [
            "report",
            "ancestral-methods-summary",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--kind",
            "continuous",
            "--model",
            "brownian",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["reconstruction_kind"] == "continuous"
    assert payload["metrics"]["model"] == "brownian"
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    assert output_path.exists()


def test_cli_report_ancestral_methods_summary_writes_discrete_metrics(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "ancestral-discrete-methods-summary.md"

    exit_code = main(
        [
            "report",
            "ancestral-methods-summary",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--kind",
            "discrete",
            "--model",
            "equal-rates",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["reconstruction_kind"] == "discrete"
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    text = output_path.read_text(encoding="utf-8")
    assert "Ancestral Reconstruction Methods Summary" in text
    assert "- discrete model: `equal-rates`" in text
