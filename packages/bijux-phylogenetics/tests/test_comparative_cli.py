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


def test_comparative_pgls_cli_encodes_categorical_predictor(capsys) -> None:
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
    assert exit_code == 0
    assert payload["data"]["inputs"]["predictors"][1]["reference_level"] == "forest"
    coefficients = {
        coefficient["name"]: coefficient["estimate"]
        for coefficient in payload["data"]["model"]["coefficients"]
    }
    assert math.isclose(coefficients["habitat[tundra]"], -2.0)


def test_comparative_pgls_cli_accepts_formula_interactions(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "pgls",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_comparative_interaction.tsv")),
            "--formula",
            "response ~ predictor_one * habitat",
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
    assert payload["data"]["inputs"]["formula"]["interaction_terms"] == ["predictor_one:habitat"]
    assert math.isclose(coefficients["predictor_one:habitat[tundra]"], 0.5)


def test_comparative_multiple_testing_cli_reports_adjusted_counts(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "multiple-testing",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_comparative_multiple.tsv")),
            "--responses",
            "response_growth",
            "response_range",
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
    assert exit_code == 0
    assert payload["metrics"]["response_count"] == 2
    assert payload["metrics"]["test_count"] == 4


def test_comparative_report_cli_reports_audit_and_limitations(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "report",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--response",
            "response",
            "--predictors",
            "predictor_one",
            "--lambda-value",
            "0.0",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["audit_row_count"] == 3
    assert payload["metrics"]["limitation_count"] >= 2


def test_comparative_influence_cli_reports_taxa(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "influence",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--response",
            "response",
            "--predictors",
            "predictor_one",
            "--lambda-value",
            "0.0",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["taxon_count"] == 4


def test_comparative_compare_trees_cli_reports_deltas(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "compare-trees",
            str(fixture("example_tree.nwk")),
            str(fixture("example_tree_topology_diff.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--response",
            "response",
            "--predictors",
            "predictor_one",
            "--lambda-value",
            "0.0",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["coefficient_delta_count"] >= 1


def test_comparative_compare_pruning_cli_reports_dropped_taxa(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "compare-pruning",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_comparative_interaction.tsv")),
            "--formula",
            "response ~ predictor_one + habitat",
            "--drop-taxa",
            "F",
            "--lambda-value",
            "0.0",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["dropped_taxa"] == 1
