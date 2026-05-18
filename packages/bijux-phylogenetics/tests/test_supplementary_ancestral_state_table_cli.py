from __future__ import annotations

import csv
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


def test_cli_report_supplementary_ancestral_state_table_writes_continuous_metrics(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "supplementary-ancestral-continuous.tsv"

    exit_code = main(
        [
            "report",
            "supplementary-ancestral-state-table",
            "--tree",
            str(fixture("example_tree.nwk")),
            "--traits",
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--reconstruction-kind",
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
    assert payload["metrics"]["row_count"] == 3
    assert payload["metrics"]["reconstruction_kind"] == "continuous"
    assert payload["metrics"]["model"] == "brownian"
    assert payload["metrics"]["analysis_taxon_count"] == 4
    assert payload["metrics"]["excluded_taxon_count"] == 0
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 3
    assert rows[0]["reconstruction_kind"] == "continuous"
    assert rows[0]["estimate_value"] != ""


def test_cli_report_supplementary_ancestral_state_table_writes_discrete_probabilities(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "supplementary-ancestral-discrete.tsv"

    exit_code = main(
        [
            "report",
            "supplementary-ancestral-state-table",
            "--tree",
            str(fixture("example_tree.nwk")),
            "--traits",
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--reconstruction-kind",
            "discrete",
            "--model",
            "equal-rates",
            "--root-prior-mode",
            "equal",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["row_count"] == 3
    assert payload["metrics"]["reconstruction_kind"] == "discrete"
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["analysis_taxon_count"] == 4
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 3
    assert rows[0]["reconstruction_kind"] == "discrete"
    assert rows[0]["most_likely_state"] != ""
    assert rows[0]["state_probabilities"].startswith("{")
