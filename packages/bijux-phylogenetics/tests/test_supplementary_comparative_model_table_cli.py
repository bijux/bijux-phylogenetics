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


def test_cli_report_supplementary_comparative_model_table_writes_metrics(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "supplementary-comparative-model.tsv"

    exit_code = main(
        [
            "report",
            "supplementary-comparative-model-table",
            "--tree",
            str(fixture("example_tree_six_taxa.nwk")),
            "--traits",
            str(fixture("example_traits_comparative_multiple.tsv")),
            "--formula",
            "response_growth ~ predictor_one",
            "--formula",
            "response_growth ~ predictor_two",
            "--formula",
            "response_growth ~ predictor_one + predictor_two",
            "--lambda-value",
            "0.0",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["row_count"] == 7
    assert payload["metrics"]["model_count"] == 3
    assert payload["metrics"]["selected_formula"] == "response_growth ~ predictor_one"
    assert payload["metrics"]["selected_criterion"] == "AICc"
    assert payload["metrics"]["excluded_taxon_count"] == 0
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert any(row["coefficient_name"] == "predictor_one" for row in rows)
