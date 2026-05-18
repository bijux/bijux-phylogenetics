from __future__ import annotations

import csv
import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def _write_iqtree_model_fixture(report_path: Path, sidecar_path: Path) -> None:
    report_path.write_text(
        "\n".join(
            [
                " No. Model         -LnL         df  AIC          AICc         BIC",
                "  1  JC            23.429       5   56.857       86.857       57.255",
                "  2  K2P+I         21.104       7   56.208       168.208      56.764",
                "Akaike Information Criterion:           K2P+I",
                "Corrected Akaike Information Criterion: JC",
                "Bayesian Information Criterion:         K2P+I",
                "Best-fit model according to BIC: K2P+I",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    sidecar_path.write_text(
        "\n".join(
            [
                "best_model_AIC: K2P+I",
                "best_model_AICc: JC",
                "best_model_BIC: K2P+I",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_cli_report_supplementary_model_selection_table_writes_metrics(
    tmp_path: Path, capsys
) -> None:
    report_path = tmp_path / "example.iqtree"
    sidecar_path = tmp_path / "example.model"
    output_path = tmp_path / "supplementary-model-selection.tsv"
    _write_iqtree_model_fixture(report_path, sidecar_path)

    exit_code = main(
        [
            "report",
            "supplementary-model-selection-table",
            "--iqtree-report",
            str(report_path),
            "--model-sidecar",
            str(sidecar_path),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["row_count"] == 2
    assert payload["metrics"]["candidate_count"] == 2
    assert payload["metrics"]["selected_model"] == "K2P+I"
    assert payload["metrics"]["selected_criterion"] == "BIC"
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[-1]["selected_model"] == "True"
