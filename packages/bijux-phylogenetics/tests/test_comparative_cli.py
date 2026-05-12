from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

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
    assert 0.0 <= payload["metrics"]["lambda_likelihood_ratio_p_value"] <= 1.0
    assert payload["metrics"]["permutation_row_count"] == 19
    assert len(payload["data"]["signal_test"]["permutation_rows"]) == 19


def test_comparative_signal_cli_writes_summary_and_permutation_ledgers(
    tmp_path: Path, capsys
) -> None:
    summary_out = tmp_path / "phylogenetic-signal-summary.tsv"
    permutations_out = tmp_path / "phylogenetic-signal-permutations.tsv"
    exit_code = main(
        [
            "comparative",
            "signal",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--permutations",
            "7",
            "--seed",
            "5",
            "--summary-out",
            str(summary_out),
            "--permutations-out",
            str(permutations_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["permutation_row_count"] == 7
    assert summary_out.exists()
    assert permutations_out.exists()
    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    permutation_rows = permutations_out.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("trait\ttaxon_count\tblombergs_k")
    assert permutation_rows[0].startswith(
        "trait\tobserved_k\testimated_lambda\tpermutations"
    )
    assert len(summary_rows) == 2
    assert len(permutation_rows) == 8


def test_comparative_contrasts_cli_reports_regression_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "contrasts",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--predictor-trait",
            "predictor_one",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["taxon_count"] == 4
    assert payload["metrics"]["contrast_count"] == 3
    assert payload["metrics"]["regression_row_count"] == 3
    assert math.isclose(
        payload["metrics"]["regression_slope"],
        0.9576271186440678,
        abs_tol=1e-12,
    )
    assert 0.0 < payload["metrics"]["regression_p_value"] < 0.5
    assert payload["data"]["contrast_report"]["trait"] == "response"
    assert payload["data"]["regression"]["predictor_trait"] == "predictor_one"
    assert len(payload["data"]["regression"]["rows"]) == 3


def test_comparative_contrasts_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    contrasts_out = tmp_path / "independent-contrasts.tsv"
    regression_out = tmp_path / "independent-contrast-regression.tsv"
    exit_code = main(
        [
            "comparative",
            "contrasts",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--predictor-trait",
            "predictor_one",
            "--contrasts-out",
            str(contrasts_out),
            "--regression-out",
            str(regression_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["contrast_count"] == 3
    assert payload["metrics"]["regression_row_count"] == 3
    assert contrasts_out.exists()
    assert regression_out.exists()
    contrast_rows = contrasts_out.read_text(encoding="utf-8").splitlines()
    regression_rows = regression_out.read_text(encoding="utf-8").splitlines()
    assert contrast_rows[0].startswith("trait\tnode\tleft_taxa\tright_taxa")
    assert regression_rows[0].startswith(
        "response_trait\tpredictor_trait\tnode\tpredictor_contrast"
    )
    assert len(contrast_rows) == 4
    assert len(regression_rows) == 4


def test_comparative_contrasts_cli_requires_predictor_for_regression_output() -> None:
    with pytest.raises(SystemExit) as excinfo:
        main(
            [
                "comparative",
                "contrasts",
                str(fixture("example_tree.nwk")),
                str(fixture("example_traits_comparative.tsv")),
                "--trait",
                "response",
                "--regression-out",
                "artifacts/independent-contrast-regression.tsv",
            ]
        )
    assert excinfo.value.code == 2


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


def test_comparative_brownian_pgls_cli_reports_covariance_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "brownian-pgls",
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
    coefficients = {
        coefficient["name"]: coefficient["estimate"]
        for coefficient in payload["data"]["model"]["coefficients"]
    }
    assert exit_code == 0
    assert payload["metrics"]["covariance_model"] == "brownian-shared-path"
    assert payload["metrics"]["lambda_value"] == 1.0
    assert payload["metrics"]["tree_is_ultrametric"] is True
    assert payload["metrics"]["covariance_row_count"] == 16
    assert payload["metrics"]["positive_definite_before_stabilization"] is True
    assert math.isclose(coefficients["intercept"], 0.305084745762712, abs_tol=1e-6)
    assert math.isclose(coefficients["predictor_one"], 0.957627118644068, abs_tol=1e-6)


def test_comparative_brownian_pgls_cli_writes_covariance_table(
    tmp_path: Path, capsys
) -> None:
    covariance_out = tmp_path / "brownian-covariance.tsv"
    exit_code = main(
        [
            "comparative",
            "brownian-pgls",
            str(fixture("example_tree_internal_long_branch.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--response",
            "response",
            "--predictors",
            "predictor_one",
            "--covariance-out",
            str(covariance_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["tree_is_ultrametric"] is False
    assert payload["metrics"]["covariance_row_count"] == 16
    assert covariance_out.exists()
    written_rows = covariance_out.read_text(encoding="utf-8").splitlines()
    assert written_rows[0].startswith("left_taxon\tright_taxon\tis_diagonal")
    assert len(written_rows) == 17


def test_comparative_ou_pgls_cli_reports_alpha_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "ou-pgls",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--response",
            "response",
            "--predictors",
            "predictor_one",
            "--alpha",
            "1.0",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    coefficients = {
        coefficient["name"]: coefficient["estimate"]
        for coefficient in payload["data"]["model"]["coefficients"]
    }
    assert exit_code == 0
    assert payload["metrics"]["covariance_model"] == "ou-stationary-root"
    assert payload["metrics"]["alpha"] == 1.0
    assert payload["metrics"]["alpha_estimation_mode"] == "fixed"
    assert payload["metrics"]["alpha_profile_point_count"] == 1
    assert payload["metrics"]["covariance_row_count"] == 16
    assert payload["metrics"]["aic"] > 0.0
    assert math.isclose(coefficients["intercept"], 0.43120431282304317, abs_tol=1e-6)
    assert math.isclose(
        coefficients["predictor_one"], 0.9087628025668772, abs_tol=1e-6
    )


def test_comparative_ou_pgls_cli_writes_covariance_and_alpha_profile_tables(
    tmp_path: Path, capsys
) -> None:
    covariance_out = tmp_path / "ou-covariance.tsv"
    profile_out = tmp_path / "ou-alpha-profile.tsv"
    exit_code = main(
        [
            "comparative",
            "ou-pgls",
            str(fixture("example_tree_internal_long_branch.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--response",
            "response",
            "--predictors",
            "predictor_one",
            "--alpha",
            "estimate",
            "--covariance-out",
            str(covariance_out),
            "--alpha-profile-out",
            str(profile_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["alpha_estimation_mode"] == "estimated"
    assert payload["metrics"]["alpha_profile_point_count"] == 8
    assert covariance_out.exists()
    assert profile_out.exists()
    covariance_rows = covariance_out.read_text(encoding="utf-8").splitlines()
    profile_rows = profile_out.read_text(encoding="utf-8").splitlines()
    assert covariance_rows[0].startswith("left_taxon\tright_taxon\tis_diagonal")
    assert profile_rows[0].startswith("alpha_estimation_mode\talpha\tlog_likelihood")
    assert len(covariance_rows) == 17
    assert len(profile_rows) == 9


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


def test_comparative_logistic_cli_reports_binary_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "logistic",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_phylogenetic_logistic.tsv")),
            "--response",
            "presence",
            "--predictors",
            "body_size",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    coefficients = {
        coefficient["name"]: coefficient["estimate"]
        for coefficient in payload["data"]["coefficients"]
    }
    assert exit_code == 0
    assert payload["metrics"]["taxon_count"] == 6
    assert payload["metrics"]["success_count"] == 3
    assert payload["metrics"]["failure_count"] == 3
    assert payload["metrics"]["coefficient_count"] == 2
    assert payload["metrics"]["fitted_row_count"] == 6
    assert payload["metrics"]["lambda_value"] == 1.0
    assert payload["metrics"]["approximation_method"] == (
        "phylogenetic-working-correlation-gee"
    )
    assert payload["metrics"]["converged"] is True
    assert payload["metrics"]["warning_count"] == 0
    assert payload["metrics"]["coefficient_inference_distribution"] == "wald-normal"
    assert coefficients["body_size"] > 0.0


def test_comparative_logistic_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    coefficients_out = tmp_path / "phylogenetic-logistic-coefficients.tsv"
    fitted_out = tmp_path / "phylogenetic-logistic-fitted.tsv"
    excluded_out = tmp_path / "phylogenetic-logistic-excluded.tsv"
    exit_code = main(
        [
            "comparative",
            "logistic",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_phylogenetic_logistic_separated.tsv")),
            "--formula",
            "presence ~ habitat",
            "--coefficients-out",
            str(coefficients_out),
            "--fitted-out",
            str(fitted_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["separation_detected"] is True
    assert payload["metrics"]["warning_count"] >= 1
    assert coefficients_out.exists()
    assert fitted_out.exists()
    assert excluded_out.exists()
    coefficient_rows = coefficients_out.read_text(encoding="utf-8").splitlines()
    fitted_rows = fitted_out.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    assert coefficient_rows[0].startswith("response\tterm\testimate")
    assert fitted_rows[0].startswith(
        "taxon\tobserved_response\tfitted_probability\tlinear_predictor"
    )
    assert excluded_rows == ["taxon\treason\tdetails"]


def test_comparative_model_selection_cli_reports_ranked_formula_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "model-selection",
            str(fixture("example_tree_eight_taxa.nwk")),
            str(fixture("example_traits_phylogenetic_logistic_model_selection.tsv")),
            "--formula",
            "presence ~ body_size",
            "--formula",
            "presence ~ habitat",
            "--formula",
            "presence ~ body_size + habitat",
            "--lambda-value",
            "1.0",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    ranked_formulas = [
        row["formula"]
        for row in sorted(payload["data"]["rows"], key=lambda row: row["rank"])
    ]
    assert exit_code == 0
    assert payload["metrics"]["model_family"] == "logistic"
    assert payload["metrics"]["model_count"] == 3
    assert payload["metrics"]["analysis_taxon_count"] == 8
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["pairwise_comparison_count"] == 3
    assert payload["metrics"]["best_formula"] == "presence ~ body_size"
    assert ranked_formulas == [
        "presence ~ body_size",
        "presence ~ body_size + habitat",
        "presence ~ habitat",
    ]


def test_comparative_model_selection_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    ranking_out = tmp_path / "comparative-model-ranking.tsv"
    pairwise_out = tmp_path / "comparative-model-pairwise.tsv"
    excluded_out = tmp_path / "comparative-model-excluded.tsv"
    exit_code = main(
        [
            "comparative",
            "model-selection",
            str(fixture("example_tree_eight_taxa.nwk")),
            str(fixture("example_traits_phylogenetic_logistic_model_selection.tsv")),
            "--formula",
            "presence ~ body_size",
            "--formula",
            "presence ~ habitat",
            "--formula",
            "presence ~ body_size + habitat",
            "--lambda-value",
            "1.0",
            "--ranking-out",
            str(ranking_out),
            "--pairwise-out",
            str(pairwise_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["selected_criterion"] == "AICc"
    assert ranking_out.exists()
    assert pairwise_out.exists()
    assert excluded_out.exists()
    ranking_rows = ranking_out.read_text(encoding="utf-8").splitlines()
    pairwise_rows = pairwise_out.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    assert ranking_rows[0].startswith("formula\tmodel_family\tparameter_count")
    assert pairwise_rows[0].startswith("left_formula\tright_formula\tcomparison_kind")
    assert excluded_rows == ["taxon\treason\tmissing_columns"]
    assert len(ranking_rows) == 4
    assert len(pairwise_rows) == 4


def test_comparative_clade_residuals_cli_reports_heavy_clade_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "clade-residuals",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_comparative_multiple.tsv")),
            "--formula",
            "response_growth ~ predictor_two",
            "--lambda-value",
            "0.0",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["model_family"] == "pgls"
    assert payload["metrics"]["taxon_count"] == 6
    assert payload["metrics"]["clade_count"] == 4
    assert payload["metrics"]["residual_heavy_clade_count"] == 1
    assert payload["metrics"]["top_influential_clade"] == "E|F"
    assert payload["metrics"]["standardized_residual_method"] == (
        "leveraged-gls-residual"
    )


def test_comparative_clade_residuals_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    taxa_out = tmp_path / "comparative-residual-taxa.tsv"
    clades_out = tmp_path / "comparative-residual-clades.tsv"
    exit_code = main(
        [
            "comparative",
            "clade-residuals",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_comparative_multiple.tsv")),
            "--formula",
            "response_growth ~ predictor_two",
            "--lambda-value",
            "0.0",
            "--taxa-out",
            str(taxa_out),
            "--clades-out",
            str(clades_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["top_influential_clade"] == "E|F"
    assert taxa_out.exists()
    assert clades_out.exists()
    taxa_rows = taxa_out.read_text(encoding="utf-8").splitlines()
    clade_rows = clades_out.read_text(encoding="utf-8").splitlines()
    assert taxa_rows[0].startswith("taxon\tobserved_value\tfitted_value\tresidual")
    assert clade_rows[0].startswith("clade_id\tnode_label\ttaxon_count\ttaxa")
    assert len(taxa_rows) == 7
    assert len(clade_rows) == 5


def test_comparative_clade_stability_cli_reports_influential_clades(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "clade-stability",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_comparative_multiple.tsv")),
            "--formula",
            "response_growth ~ predictor_two",
            "--lambda-value",
            "0.0",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["model_family"] == "pgls"
    assert payload["metrics"]["baseline_taxon_count"] == 6
    assert payload["metrics"]["baseline_term_count"] == 2
    assert payload["metrics"]["candidate_clade_count"] == 4
    assert payload["metrics"]["blocked_clade_count"] == 1
    assert payload["metrics"]["coefficient_change_row_count"] == 6
    assert payload["metrics"]["top_influential_clade"] == "A|B"
    assert payload["metrics"]["minimum_major_clade_size"] == 2


def test_comparative_clade_stability_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    clades_out = tmp_path / "comparative-clade-stability.tsv"
    terms_out = tmp_path / "comparative-clade-coefficients.tsv"
    exit_code = main(
        [
            "comparative",
            "clade-stability",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_comparative_multiple.tsv")),
            "--formula",
            "response_growth ~ predictor_two",
            "--lambda-value",
            "0.0",
            "--clades-out",
            str(clades_out),
            "--terms-out",
            str(terms_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["top_influential_clade"] == "A|B"
    assert clades_out.exists()
    assert terms_out.exists()
    clade_rows = clades_out.read_text(encoding="utf-8").splitlines()
    term_rows = terms_out.read_text(encoding="utf-8").splitlines()
    assert clade_rows[0].startswith("clade_id\tnode_label\tdropped_taxon_count")
    assert term_rows[0].startswith("clade_id\tnode_label\tterm\tbaseline_estimate")
    assert len(clade_rows) == 5
    assert len(term_rows) == 7


def test_comparative_posterior_pgls_cli_reports_stability_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "posterior-pgls",
            str(fixture("example_posterior_tree_set_six_taxa.nwk")),
            str(fixture("example_traits_comparative_multiple.tsv")),
            "--formula",
            "response_growth ~ predictor_two",
            "--lambda-value",
            "estimate",
            "--significance-threshold",
            "0.1",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["total_tree_count"] == 5
    assert payload["metrics"]["burnin_tree_count"] == 0
    assert payload["metrics"]["kept_tree_count"] == 5
    assert payload["metrics"]["analysis_taxon_count"] == 6
    assert payload["metrics"]["rooted_topology_count"] == 5
    assert payload["metrics"]["unrooted_topology_count"] == 4
    assert payload["metrics"]["tree_fit_row_count"] == 5
    assert payload["metrics"]["coefficient_row_count"] == 10
    assert payload["metrics"]["coefficient_summary_count"] == 2
    assert payload["metrics"]["stable_supported_term_count"] == 1
    assert payload["metrics"]["direction_conflict_term_count"] == 0
    assert payload["metrics"]["lambda_mode"] == "estimate"


def test_comparative_posterior_pgls_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    trees_out = tmp_path / "posterior-tree-pgls-trees.tsv"
    coefficients_out = tmp_path / "posterior-tree-pgls-coefficients.tsv"
    summary_out = tmp_path / "posterior-tree-pgls-summary.tsv"
    exit_code = main(
        [
            "comparative",
            "posterior-pgls",
            str(fixture("example_posterior_tree_set_six_taxa.nwk")),
            str(fixture("example_traits_comparative_multiple.tsv")),
            "--formula",
            "response_growth ~ predictor_two",
            "--lambda-value",
            "estimate",
            "--significance-threshold",
            "0.1",
            "--trees-out",
            str(trees_out),
            "--coefficients-out",
            str(coefficients_out),
            "--summary-out",
            str(summary_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["coefficient_summary_count"] == 2
    assert trees_out.exists()
    assert coefficients_out.exists()
    assert summary_out.exists()
    tree_rows = trees_out.read_text(encoding="utf-8").splitlines()
    coefficient_rows = coefficients_out.read_text(encoding="utf-8").splitlines()
    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    assert tree_rows[0].startswith("source_tree_index\tpost_burnin_index\trooted_topology_id")
    assert coefficient_rows[0].startswith("source_tree_index\tpost_burnin_index\trooted_topology_id\tterm")
    assert summary_rows[0].startswith("term\ttree_fit_count\tpositive_tree_count")
    assert len(tree_rows) == 6
    assert len(coefficient_rows) == 11
    assert len(summary_rows) == 3


def test_comparative_multivariate_cli_reports_shared_taxa_and_associations(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "multivariate",
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
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["response_count"] == 2
    assert payload["metrics"]["predictor_count"] == 2
    assert payload["metrics"]["analysis_taxa"] == 6
    assert payload["metrics"]["excluded_taxa"] == 0
    assert payload["metrics"]["residual_covariance_row_count"] == 4
    assert payload["metrics"]["residual_association_count"] == 1
    assert len(payload["data"]["response_models"]) == 2


def test_comparative_multivariate_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    covariance_out = tmp_path / "multivariate-residual-covariance.tsv"
    associations_out = tmp_path / "multivariate-residual-associations.tsv"
    excluded_out = tmp_path / "multivariate-excluded-taxa.tsv"
    exit_code = main(
        [
            "comparative",
            "multivariate",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_comparative_multivariate_missing.tsv")),
            "--responses",
            "response_growth",
            "response_range",
            "--predictors",
            "predictor_one",
            "predictor_two",
            "--lambda-value",
            "0.0",
            "--covariance-out",
            str(covariance_out),
            "--associations-out",
            str(associations_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["analysis_taxa"] == 5
    assert payload["metrics"]["excluded_taxa"] == 1
    assert covariance_out.exists()
    assert associations_out.exists()
    assert excluded_out.exists()
    covariance_rows = covariance_out.read_text(encoding="utf-8").splitlines()
    association_rows = associations_out.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    assert covariance_rows[0].startswith(
        "left_response\tright_response\tpair_count\tis_diagonal"
    )
    assert association_rows[0].startswith(
        "left_response\tright_response\tpair_count\tcovariance\tcorrelation"
    )
    assert excluded_rows[0] == "taxon\treason\tmissing_columns"


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
