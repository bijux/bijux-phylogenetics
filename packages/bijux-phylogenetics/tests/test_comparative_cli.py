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
        coefficient["name"]: coefficient
        for coefficient in payload["data"]["model"]["coefficients"]
    }
    assert exit_code == 0
    assert payload["metrics"]["coefficient_count"] == 3
    assert payload["metrics"]["confidence_interval_count"] == 3
    assert payload["metrics"]["residual_degrees_of_freedom"] == 1
    assert payload["metrics"]["coefficient_inference_distribution"] == "student-t"
    assert math.isclose(coefficients["intercept"]["estimate"], 1.0)
    assert math.isclose(coefficients["predictor_one"]["estimate"], 0.5)
    assert math.isclose(coefficients["predictor_two"]["estimate"], 1.0)
    assert coefficients["predictor_one"]["inference_distribution"] == "student-t"
    assert coefficients["predictor_one"]["degrees_of_freedom"] == 1
    assert payload["metrics"]["lambda_estimation_mode"] == "fixed"
    assert payload["metrics"]["lambda_profile_point_count"] == 1
    assert payload["metrics"]["lambda_lower_95_confidence_interval"] is None
    assert payload["metrics"]["lambda_upper_95_confidence_interval"] is None
    assert (
        coefficients["predictor_one"]["lower_95_confidence_interval"]
        < coefficients["predictor_one"]["estimate"]
        < coefficients["predictor_one"]["upper_95_confidence_interval"]
    )
    assert "test_statistic" in coefficients["predictor_one"]


def test_comparative_pgls_cli_reports_estimated_lambda_profile(capsys) -> None:
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
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["lambda_estimation_mode"] == "estimated"
    assert payload["metrics"]["lambda_profile_point_count"] == 101
    assert payload["metrics"]["lambda_lower_95_confidence_interval"] == 0.0
    assert payload["metrics"]["lambda_upper_95_confidence_interval"] == 1.0
    assert payload["data"]["model"]["lambda_fit"]["mode"] == "estimated"
    assert len(payload["data"]["model"]["lambda_fit"]["profile_rows"]) == 101


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
    assert payload["data"]["inputs"]["formula_audit"]["parameter_count"] == 3
    assert payload["metrics"]["categorical_contrast_predictor_count"] == 1
    assert payload["metrics"]["categorical_contrast_row_count"] == 2
    coefficients = {
        coefficient["name"]: coefficient["estimate"]
        for coefficient in payload["data"]["model"]["coefficients"]
    }
    assert math.isclose(coefficients["habitat[tundra]"], -2.0)
    contrast_rows = payload["data"]["categorical_contrasts"]["rows"]
    assert contrast_rows[0]["is_reference_level"] is True
    assert contrast_rows[0]["level"] == "forest"
    assert contrast_rows[1]["coefficient_name"] == "habitat[tundra]"


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
    assert payload["metrics"]["interaction_term_count"] == 1
    assert payload["metrics"]["interaction_coefficient_row_count"] == 1
    assert payload["data"]["inputs"]["formula"]["interaction_terms"] == [
        "predictor_one:habitat"
    ]
    assert payload["data"]["inputs"]["formula_audit"]["interaction_terms"][0][
        "encoded_columns"
    ] == ["predictor_one:habitat[tundra]"]
    assert math.isclose(coefficients["predictor_one:habitat[tundra]"], 0.5)
    interaction_row = payload["data"]["interaction_coefficients"]["rows"][0]
    assert interaction_row["interaction_kind"] == "continuous-by-categorical"
    assert interaction_row["component_levels"] == [None, "tundra"]


def test_comparative_pgls_cli_reports_transformed_terms(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "pgls",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--formula",
            "response ~ log(predictor_one) + habitat",
            "--lambda-value",
            "0.0",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["transformed_term_count"] == 1
    assert payload["data"]["inputs"]["formula_audit"]["transformed_terms"] == [
        "log(predictor_one)"
    ]


def test_comparative_pgls_cli_writes_interceptless_model_matrix(
    tmp_path: Path, capsys
) -> None:
    matrix_out = tmp_path / "comparative-model-matrix.tsv"
    exit_code = main(
        [
            "comparative",
            "pgls",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--formula",
            "response ~ 0 + habitat",
            "--model-matrix-out",
            str(matrix_out),
            "--lambda-value",
            "0.0",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["intercept_included"] is False
    assert payload["metrics"]["model_matrix_row_count"] == 4
    assert payload["metrics"]["model_matrix_column_count"] == 2
    assert payload["data"]["inputs"]["formula"]["include_intercept"] is False
    assert payload["data"]["inputs"]["model_matrix"]["encoded_columns"] == [
        "habitat[forest]",
        "habitat[tundra]",
    ]
    assert matrix_out.exists()
    written_rows = matrix_out.read_text(encoding="utf-8").splitlines()
    assert written_rows[0] == "taxon\tresponse_value\thabitat[forest]\thabitat[tundra]"


def test_comparative_pgls_cli_writes_categorical_contrast_table(
    tmp_path: Path, capsys
) -> None:
    contrasts_out = tmp_path / "categorical-contrasts.tsv"
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
            "--categorical-contrasts-out",
            str(contrasts_out),
            "--lambda-value",
            "0.0",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["categorical_contrast_predictor_count"] == 1
    assert payload["metrics"]["categorical_contrast_row_count"] == 2
    assert contrasts_out.exists()
    written_rows = contrasts_out.read_text(encoding="utf-8").splitlines()
    assert written_rows[0].startswith("predictor\tsource_column\tencoding_scheme")
    assert any("habitat\thabitat\treference-level\tforest\tforest\ttrue" in row for row in written_rows[1:])


def test_comparative_pgls_cli_writes_interaction_coefficient_table(
    tmp_path: Path, capsys
) -> None:
    interaction_out = tmp_path / "interaction-coefficients.tsv"
    exit_code = main(
        [
            "comparative",
            "pgls",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_comparative_interaction.tsv")),
            "--formula",
            "response ~ predictor_one * habitat",
            "--interaction-coefficients-out",
            str(interaction_out),
            "--lambda-value",
            "0.0",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["interaction_term_count"] == 1
    assert payload["metrics"]["interaction_coefficient_row_count"] == 1
    assert interaction_out.exists()
    written_rows = interaction_out.read_text(encoding="utf-8").splitlines()
    assert written_rows[0].startswith(
        "interaction_term\tinteraction_kind\tcoefficient_name"
    )
    assert any(
        "predictor_one:habitat\tcontinuous-by-categorical\tpredictor_one:habitat[tundra]"
        in row
        for row in written_rows[1:]
    )


def test_comparative_pgls_cli_writes_lambda_profile_table(
    tmp_path: Path, capsys
) -> None:
    profile_out = tmp_path / "lambda-profile.tsv"
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
            "--lambda-profile-out",
            str(profile_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["lambda_estimation_mode"] == "estimated"
    assert payload["metrics"]["lambda_profile_point_count"] == 101
    assert profile_out.exists()
    written_rows = profile_out.read_text(encoding="utf-8").splitlines()
    assert written_rows[0].startswith("mode\tlambda_value\tlog_likelihood")
    assert len(written_rows) == 102


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
    assert payload["metrics"]["family_size"] == 4
    assert (
        payload["metrics"]["raw_significant_count"]
        >= payload["metrics"]["significant_count"]
    )


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
    assert "formula_audit" in payload["data"]["snapshot"]["pgls_inputs"]
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
    assert payload["metrics"]["top_taxa"] == 3


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
    assert "conclusion_changed" in payload["metrics"]


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
    assert "conclusion_changed" in payload["metrics"]
