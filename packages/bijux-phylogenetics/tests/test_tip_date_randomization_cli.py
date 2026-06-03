from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "metadata")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_cli_tip_date_randomization_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "tip-date-randomization"

    exit_code = main(
        [
            "diagnose",
            "tip-date-randomization",
            str(fixture("root_to_tip_regression_diagnostic_tree_7_taxa.nwk")),
            "--metadata",
            str(fixture("root_to_tip_regression_dates_7_taxa.tsv")),
            "--permutations",
            "19",
            "--seed",
            "17",
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["tip_count"] == 7
    assert payload["metrics"]["observed_slope"] == pytest.approx(
        2.392857142857143,
        abs=1e-12,
    )
    assert payload["metrics"]["observed_r_squared"] == pytest.approx(
        0.6390945330296127,
        abs=1e-12,
    )
    assert payload["metrics"]["permutations"] == 19
    assert payload["metrics"]["seed"] == 17
    assert payload["metrics"]["p_value"] == pytest.approx(0.05, abs=1e-12)
    assert payload["metrics"]["null_distribution_mean"] == pytest.approx(
        0.15226741397913907,
        abs=1e-15,
    )
    assert payload["data"]["observed_regression"]["outliers"][0]["tip"] == "G"
    assert len(payload["data"]["permutation_rows"]) == 19
    assert sorted(Path(path).name for path in payload["outputs"]) == [
        "observed_outliers.tsv",
        "observed_residuals.tsv",
        "permutations.tsv",
        "run.json",
        "summary.tsv",
    ]
    assert (out_dir / "summary.tsv").is_file()
    assert (out_dir / "permutations.tsv").is_file()
    assert (out_dir / "observed_residuals.tsv").is_file()
    assert (out_dir / "observed_outliers.tsv").is_file()
    assert (out_dir / "run.json").is_file()
