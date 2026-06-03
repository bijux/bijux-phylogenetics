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


def test_comparative_disparity_cli_writes_summary_and_clade_ledgers(
    tmp_path: Path,
    capsys,
) -> None:
    summary_out = tmp_path / "disparity-summary.tsv"
    clades_out = tmp_path / "disparity-clades.tsv"
    excluded_out = tmp_path / "disparity-excluded.tsv"

    exit_code = main(
        [
            "comparative",
            "disparity",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--traits",
            "response,predictor_one",
            "--summary-out",
            str(summary_out),
            "--clades-out",
            str(clades_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    assert payload["metrics"]["trait_column_count"] == 2
    assert payload["metrics"]["clade_count"] == 3
    assert payload["metrics"]["root_disparity"] == 5.5
    assert payload["data"]["method_formula"]
    assert summary_out.exists()
    assert clades_out.exists()
    assert excluded_out.exists()
