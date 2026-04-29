from __future__ import annotations

import json
import math
from pathlib import Path

from bijux_phylogenetics.cli import main


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


def test_comparative_readiness_cli_reports_analysis_taxa(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "readiness",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_validate.tsv")),
            "--trait",
            "height_cm",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["analysis_taxa"] == 4
    assert payload["data"]["ready"] is True


def test_comparative_signal_cli_reports_lambda_and_p_value(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "signal",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--permutations",
            "19",
            "--seed",
            "7",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert 0.0 <= payload["metrics"]["pagels_lambda"] <= 1.0
    assert 0.0 < payload["metrics"]["signal_p_value"] <= 1.0


def test_comparative_pgls_cli_fits_multiple_predictors(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "pgls",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--response",
            "response",
            "--predictors",
            "predictor_one",
            "predictor_two",
            "--lambda-value",
            "0.0",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    coefficients = {
        coefficient["name"]: coefficient["estimate"]
        for coefficient in payload["data"]["model"]["coefficients"]
    }
    assert exit_code == 0
    assert math.isclose(coefficients["intercept"], 1.0)
    assert math.isclose(coefficients["predictor_one"], 0.5)
    assert math.isclose(coefficients["predictor_two"], 1.0)


def test_comparative_pgls_cli_reports_categorical_predictor_error(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "pgls",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--response",
            "response",
            "--predictors",
            "predictor_one",
            "habitat",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["errors"][0]["code"] == "comparative_method_error"
