from __future__ import annotations

import csv
import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

RABIES_EXPECTED_BUNDLE = (
    Path(__file__).resolve().parent.parent
    / "src"
    / "bijux_phylogenetics"
    / "resources"
    / "datasets"
    / "pathogens"
    / "rabies_method_sensitivity_panel"
    / "expected"
)


def test_cli_report_supplementary_batch_summary_table_writes_metrics(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "supplementary-batch-summary.tsv"

    exit_code = main(
        [
            "report",
            "supplementary-batch-summary-table",
            "--workflow-bundle-root",
            str(RABIES_EXPECTED_BUNDLE),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["row_count"] == 5
    assert payload["metrics"]["dataset_row_count"] == 1
    assert payload["metrics"]["variant_row_count"] == 4
    assert payload["metrics"]["workflow_status"] == "succeeded"
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 5
    assert rows[0]["row_scope"] == "dataset"
