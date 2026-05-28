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


def test_cli_root_to_tip_regression_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "root-to-tip-regression"

    exit_code = main(
        [
            "diagnose",
            "root-to-tip-regression",
            str(fixture("root_to_tip_regression_diagnostic_tree_7_taxa.nwk")),
            "--metadata",
            str(fixture("root_to_tip_regression_dates_7_taxa.tsv")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["tip_count"] == 7
    assert payload["metrics"]["slope"] == pytest.approx(2.392857142857143, abs=1e-12)
    assert payload["metrics"]["intercept"] == pytest.approx(
        -1.821428571428572,
        abs=1e-12,
    )
    assert payload["metrics"]["r_squared"] == pytest.approx(
        0.6390945330296127,
        abs=1e-12,
    )
    assert payload["metrics"]["outlier_count"] == 1
    assert payload["data"]["outliers"][0]["tip"] == "G"
    assert sorted(Path(path).name for path in payload["outputs"]) == [
        "outliers.tsv",
        "residuals.tsv",
        "run.json",
        "summary.tsv",
    ]
    assert (out_dir / "summary.tsv").is_file()
    assert (out_dir / "residuals.tsv").is_file()
    assert (out_dir / "outliers.tsv").is_file()
    assert (out_dir / "run.json").is_file()
