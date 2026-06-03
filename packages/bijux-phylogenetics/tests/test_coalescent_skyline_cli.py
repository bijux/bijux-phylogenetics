from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def test_cli_simulate_coalescent_writes_skyline_ledgers(
    tmp_path: Path,
    capsys,
) -> None:
    output_path = tmp_path / "simulated-coalescent.trees"
    skyline_path = tmp_path / "coalescent-skyline.tsv"

    exit_code = main(
        [
            "simulate",
            "tree-coalescent",
            "--tree-count",
            "64",
            "--tip-count",
            "5",
            "--population-size",
            "2.5",
            "--waiting-time-tolerance",
            "0.2",
            "--seed",
            "19",
            "--out",
            str(output_path),
            "--skyline-table-out",
            str(skyline_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["skyline_interval_count"] == 4
    assert payload["metrics"]["skyline_high_uncertainty_count"] == 1
    assert payload["outputs"] == [str(output_path), str(skyline_path)]
    assert output_path.read_text(encoding="utf-8").count(";\n") == 64
    assert skyline_path.read_text(encoding="utf-8").startswith(
        "interval\tlineage_count\tduration\teffective_population_size_estimate\t"
    )
