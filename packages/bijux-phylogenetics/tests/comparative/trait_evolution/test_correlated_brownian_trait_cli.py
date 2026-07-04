from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_simulate_correlated_brownian_cli_writes_collection_and_summary(
    tmp_path: Path,
    capsys,
) -> None:
    output_path = tmp_path / "correlated-brownian.tsv"
    summary_path = tmp_path / "correlated-brownian-summary.tsv"

    exit_code = main(
        [
            "simulate",
            "traits-brownian-correlated",
            str(fixture("example_tree.nwk")),
            "--trait",
            "trait_alpha",
            "--trait",
            "trait_beta",
            "--root-state",
            "1.0",
            "--root-state",
            "-0.5",
            "--correlation-row",
            "1.0,0.4",
            "--correlation-row",
            "0.4,1.0",
            "--trait-standard-deviation",
            "2.0",
            "--trait-standard-deviation",
            "1.5",
            "--replicates",
            "8",
            "--seed",
            "7",
            "--out",
            str(output_path),
            "--summary-out",
            str(summary_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["tip_count"] == 4
    assert payload["metrics"]["trait_count"] == 2
    assert payload["metrics"]["replicate_count"] == 8
    assert output_path.read_text(encoding="utf-8").splitlines()[0] == (
        "replicate_index\ttaxon\ttrait\tvalue"
    )
    assert summary_path.read_text(encoding="utf-8").splitlines()[0] == (
        "row_kind\tlabel\tmean_value\tstandard_deviation\tminimum\tmedian\tmaximum\tcovariance\tcorrelation"
    )
