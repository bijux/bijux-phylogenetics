from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main

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
    assert payload["metrics"]["tree_is_ultrametric"] is True
    assert payload["metrics"]["ultrametric_policy"] == (
        "accept-rooted-trees-and-report-ultrametricity"
    )
    assert payload["metrics"]["missing_value_policy"] == (
        "prune-overlapping-missing-values"
    )
    assert payload["metrics"]["pruned_missing_value_taxon_count"] == 0
    assert payload["metrics"]["signal_seed"] == 7
    assert 0.0 <= payload["metrics"]["lambda_likelihood_ratio_p_value"] <= 1.0
    assert payload["metrics"]["permutation_row_count"] == 19
    assert len(payload["data"]["signal_test"]["permutation_rows"]) == 19


def test_comparative_signal_cli_reports_pruned_missing_value_policy(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "signal",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_brownian_missing.tsv")),
            "--trait",
            "response_growth",
            "--permutations",
            "11",
            "--seed",
            "17",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["pruned_missing_value_taxon_count"] == 1
    assert payload["metrics"]["missing_value_policy"] == (
        "prune-overlapping-missing-values"
    )
    assert payload["data"]["input_audit"]["pruned_missing_value_taxa"] == ["B"]


def test_comparative_signal_cli_reports_non_ultrametric_acceptance(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "signal",
            str(fixture("example_tree_internal_long_branch.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--permutations",
            "11",
            "--seed",
            "19",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["tree_is_ultrametric"] is False
    assert payload["metrics"]["ultrametric_policy"] == (
        "accept-rooted-trees-and-report-ultrametricity"
    )


def test_comparative_signal_cli_reports_lambda_optimizer_diagnostics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "signal",
            str(fixture("example_tree_phytools_non_ultrametric_twenty_four_taxa.nwk")),
            str(
                fixture(
                    "example_traits_phytools_signal_non_ultrametric_twenty_four_taxa.tsv"
                )
            ),
            "--trait",
            "signal_strong",
            "--permutations",
            "11",
            "--seed",
            "19",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["lambda_optimizer_name"] == "two-stage-grid-search"
    assert payload["metrics"]["lambda_optimizer_function_evaluation_count"] > 0
    assert (
        payload["metrics"]["lambda_optimizer_function_evaluation_count"]
        == payload["data"]["pagels_lambda"]["optimizer_diagnostics"][
            "function_evaluation_count"
        ]
    )
    assert payload["metrics"]["lambda_optimizer_hit_upper_boundary"] is True
    assert payload["metrics"]["lambda_likelihood_ratio_statistic"] > 0.0
    assert (
        payload["metrics"]["lambda_log_likelihood"]
        >= payload["data"]["pagels_lambda"]["null_log_likelihood"]
    )
    assert (
        payload["metrics"]["signal_null_k_minimum"]
        <= payload["metrics"]["signal_null_k_mean"]
        <= payload["metrics"]["signal_null_k_maximum"]
    )


def test_comparative_signal_cli_rejects_constant_trait_values(
    tmp_path: Path, capsys
) -> None:
    traits_path = tmp_path / "constant-traits.tsv"
    traits_path.write_text(
        "taxon\tresponse\nA\t2.0\nB\t2.0\nC\t2.0\nD\t2.0\n",
        encoding="utf-8",
    )
    exit_code = main(
        [
            "comparative",
            "signal",
            str(fixture("example_tree.nwk")),
            str(traits_path),
            "--trait",
            "response",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"] == [
        {
            "code": "comparative_method_error",
            "message": (
                "phylogenetic signal requires at least two distinct numeric trait values after pruning"
            ),
        }
    ]


def test_comparative_brownian_pgls_cli_reports_negative_branch_length_failure(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "brownian-pgls",
            str(fixture("example_tree_negative_length.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--response",
            "response",
            "--predictors",
            "predictor_one",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 2
    assert payload["errors"][0]["code"] == "comparative_method_error"
    assert (
        payload["errors"][0]["details"]["failure_reason"]
        == "brownian_covariance_negative_branch_lengths"
    )
    assert payload["errors"][0]["details"]["evidence"]["minimum_branch_length"] < 0.0


def test_comparative_ou_pgls_cli_reports_negative_branch_length_failure(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "ou-pgls",
            str(fixture("example_tree_negative_length.nwk")),
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
    assert exit_code == 2
    assert payload["errors"][0]["code"] == "comparative_method_error"
    assert (
        payload["errors"][0]["details"]["failure_reason"]
        == "ou_covariance_negative_branch_lengths"
    )
    assert payload["errors"][0]["details"]["evidence"]["minimum_branch_length"] < 0.0


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
    assert (
        payload["metrics"]["signal_null_k_minimum"]
        <= payload["metrics"]["signal_null_k_mean"]
        <= payload["metrics"]["signal_null_k_maximum"]
    )
    assert summary_out.exists()
    assert permutations_out.exists()
    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    permutation_rows = permutations_out.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("trait\ttaxon_count\tblombergs_k")
    assert "signal_null_k_mean" in summary_rows[0]
    assert permutation_rows[0].startswith(
        "trait\tobserved_k\testimated_lambda\tpermutations"
    )
    assert len(summary_rows) == 2
    assert len(permutation_rows) == 8


def test_comparative_discrete_mk_cli_reports_er_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "discrete-mk",
            str(fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk")),
            str(fixture("example_traits_phytools_signal_twenty_four_taxa.tsv")),
            "--trait",
            "binary_state",
            "--taxon-column",
            "taxon",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["taxon_count"] == 24
    assert payload["metrics"]["observed_state_count"] == 2
    assert payload["metrics"]["parameter_count"] == 1
    assert payload["metrics"]["aic"] < payload["metrics"]["aicc"]
    assert payload["metrics"]["optimizer_name"] == "golden-section-search"
    assert payload["metrics"]["optimizer_converged"] is True
    assert payload["metrics"]["transition_rate_count"] == 2
    assert payload["metrics"]["baseline_model"] is None


@pytest.mark.slow
def test_comparative_discrete_mk_cli_reports_symmetric_baseline_metrics(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "discrete-mk",
            str(fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk")),
            str(fixture("example_traits_phytools_signal_twenty_four_taxa.tsv")),
            "--trait",
            "region_state",
            "--taxon-column",
            "taxon",
            "--model",
            "symmetric",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "symmetric"
    assert payload["metrics"]["parameter_count"] == 3
    assert payload["metrics"]["transition_rate_count"] == 6
    assert payload["metrics"]["optimizer_converged"] is True
    assert payload["metrics"]["overparameterized"] is False
    assert payload["metrics"]["baseline_model"] == "equal-rates"
    assert payload["metrics"]["preferred_model_by_aic"] == "equal-rates"
    assert payload["metrics"]["delta_aic"] > 0.0


@pytest.mark.slow
def test_comparative_discrete_mk_cli_reports_lambda_transform_metrics(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "discrete-mk",
            str(fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk")),
            str(
                fixture(
                    "example_traits_geiger_discrete_model_panel_twenty_four_taxa.tsv"
                )
            ),
            "--trait",
            "er_binary_transform_weak_signal",
            "--taxon-column",
            "taxon",
            "--transform",
            "lambda",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["transform"] == "lambda"
    assert payload["metrics"]["parameter_count"] == 2
    assert payload["metrics"]["transform_parameter_name"] == "lambda"
    assert payload["metrics"]["transform_parameter_value"] == 0.0
    assert payload["metrics"]["transform_function_evaluation_count"] == 26
    assert payload["metrics"]["transform_warning_count"] >= 1


@pytest.mark.slow
def test_comparative_discrete_mk_cli_reports_kappa_transform_metrics(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "discrete-mk",
            str(fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk")),
            str(
                fixture(
                    "example_traits_geiger_discrete_model_panel_twenty_four_taxa.tsv"
                )
            ),
            "--trait",
            "er_binary_truth",
            "--taxon-column",
            "taxon",
            "--transform",
            "kappa",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["transform"] == "kappa"
    assert payload["metrics"]["parameter_count"] == 2
    assert payload["metrics"]["transform_parameter_name"] == "kappa"
    assert math.isclose(
        payload["metrics"]["transform_parameter_value"],
        0.9011763252454394,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert payload["metrics"]["transform_function_evaluation_count"] == 26
    assert payload["metrics"]["transform_warning_count"] >= 1


def test_comparative_discrete_mk_cli_reports_delta_transform_metrics(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "discrete-mk",
            str(fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk")),
            str(
                fixture(
                    "example_traits_geiger_discrete_model_panel_twenty_four_taxa.tsv"
                )
            ),
            "--trait",
            "er_binary_delta_time_sensitive",
            "--taxon-column",
            "taxon",
            "--transform",
            "delta",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["transform"] == "delta"
    assert payload["metrics"]["parameter_count"] == 2
    assert payload["metrics"]["transform_parameter_name"] == "delta"
    assert 0.0 < payload["metrics"]["transform_parameter_value"] < 3.0
    assert payload["metrics"]["transform_function_evaluation_count"] == 26
    assert payload["metrics"]["transform_warning_count"] >= 0


def test_comparative_discrete_mk_cli_reports_early_burst_transform_metrics(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "discrete-mk",
            str(fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk")),
            str(
                fixture(
                    "example_traits_geiger_discrete_model_panel_twenty_four_taxa.tsv"
                )
            ),
            "--trait",
            "er_binary_delta_time_sensitive",
            "--taxon-column",
            "taxon",
            "--transform",
            "early-burst",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["transform"] == "early-burst"
    assert payload["metrics"]["parameter_count"] == 2
    assert payload["metrics"]["transform_parameter_name"] == "a"
    assert payload["metrics"]["transform_parameter_value"] < 0.0
    assert payload["metrics"]["transform_function_evaluation_count"] == 26
    assert payload["metrics"]["transform_warning_count"] >= 1


def test_comparative_discrete_mk_cli_rejects_meristic_parity_claim(
    capsys,
) -> None:
    with pytest.raises(SystemExit) as excinfo:
        main(
            [
                "comparative",
                "discrete-mk",
                str(fixture("example_tree.nwk")),
                str(fixture("example_traits_geography.tsv")),
                "--trait",
                "region",
                "--taxon-column",
                "taxon",
                "--model",
                "meristic",
                "--json",
            ]
        )

    assert excinfo.value.code == 2
    captured = capsys.readouterr()
    assert "explicitly excluded this round" in captured.err
    assert "ordered-state Mk support is not claimed as meristic parity" in (
        captured.err
    )


@pytest.mark.slow
def test_comparative_discrete_mk_cli_reports_ard_boundary_metrics(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "discrete-mk",
            str(fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk")),
            str(fixture("example_traits_phytools_discrete_ard_twenty_four_taxa.tsv")),
            "--trait",
            "region_state",
            "--taxon-column",
            "taxon",
            "--model",
            "all-rates-different",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "all-rates-different"
    assert payload["metrics"]["parameter_count"] == 12
    assert payload["metrics"]["transition_rate_count"] == 12
    assert payload["metrics"]["optimizer_converged"] is False
    assert payload["metrics"]["optimizer_hit_lower_parameter_bound"] is True
    assert payload["metrics"]["optimizer_hit_upper_parameter_bound"] is False
    assert payload["metrics"]["overparameterized"] is False
    assert payload["metrics"]["baseline_model"] == "equal-rates"
    assert payload["metrics"]["preferred_model_by_aic"] == "equal-rates"
    assert payload["metrics"]["delta_aic"] > 0.0


def test_comparative_discrete_mk_cli_writes_summary_and_rate_ledgers(
    tmp_path: Path, capsys
) -> None:
    summary_out = tmp_path / "discrete-mk-summary.tsv"
    rates_out = tmp_path / "discrete-mk-rates.tsv"

    exit_code = main(
        [
            "comparative",
            "discrete-mk",
            str(fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk")),
            str(
                fixture("example_traits_phytools_discrete_missing_twenty_four_taxa.tsv")
            ),
            "--trait",
            "region_state",
            "--taxon-column",
            "taxon",
            "--summary-out",
            str(summary_out),
            "--rates-out",
            str(rates_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["pruned_missing_value_taxon_count"] == 1
    assert summary_out.exists()
    assert rates_out.exists()
    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    rate_rows = rates_out.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("trait\ttaxon_column\tmodel")
    assert "aicc" in summary_rows[0]
    assert rate_rows[0].startswith(
        "source_state\ttarget_state\ttransition_allowed\tstep_distance\trate"
    )
    assert len(summary_rows) == 2
    assert len(rate_rows) == 7


def test_comparative_correlated_traits_cli_reports_coupling_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "correlated-traits",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--left-trait",
            "response",
            "--right-trait",
            "predictor_one",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["analysis_kind"] == "continuous-brownian-contrasts"
    assert payload["metrics"]["tree_taxon_count"] == 4
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["observation_row_count"] == 3
    assert payload["metrics"]["comparison_row_count"] == 2
    assert payload["metrics"]["association_measure_name"] == "evolutionary_correlation"
    assert math.isclose(
        payload["metrics"]["association_measure_value"], 0.8871275993361114
    )
    assert payload["metrics"]["better_model"] == "correlated"


def test_comparative_correlated_traits_cli_rejects_meristic_parity_claim(
    capsys,
) -> None:
    try:
        main(
            [
                "comparative",
                "correlated-traits",
                str(fixture("example_tree_eight_taxa.nwk")),
                str(fixture("example_traits_correlated_binary_missing.tsv")),
                "--left-trait",
                "trait_left",
                "--right-trait",
                "trait_right",
                "--analysis-kind",
                "binary",
                "--binary-model",
                "meristic",
                "--json",
            ]
        )
    except SystemExit as error:
        assert error.code == 2
    captured = capsys.readouterr()
    assert "explicitly excluded this round" in captured.err
    assert "integer-state meristic contract" in captured.err


def test_comparative_correlated_traits_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    summary_out = tmp_path / "correlated-traits-summary.tsv"
    comparison_out = tmp_path / "correlated-traits-comparison.tsv"
    observations_out = tmp_path / "correlated-traits-observations.tsv"
    excluded_out = tmp_path / "correlated-traits-excluded.tsv"
    exit_code = main(
        [
            "comparative",
            "correlated-traits",
            str(fixture("example_tree_eight_taxa.nwk")),
            str(fixture("example_traits_correlated_binary_missing.tsv")),
            "--left-trait",
            "trait_left",
            "--right-trait",
            "trait_right",
            "--summary-out",
            str(summary_out),
            "--comparison-out",
            str(comparison_out),
            "--observations-out",
            str(observations_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["analysis_kind"] == "binary-joint-state"
    assert payload["metrics"]["analyzed_taxon_count"] == 5
    assert payload["metrics"]["excluded_taxon_count"] == 4
    assert payload["metrics"]["joint_state_count"] == 3
    assert summary_out.exists()
    assert comparison_out.exists()
    assert observations_out.exists()
    assert excluded_out.exists()
    assert (
        summary_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("analysis_kind\tleft_trait\tright_trait")
    )
    assert (
        comparison_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("model_kind\tmodel_description\tparameter_count")
    )
    assert (
        observations_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("row_kind\tlabel\ttaxon\tleft_taxa")
    )
    assert excluded_out.read_text(encoding="utf-8").splitlines() == [
        "taxon\treason\tmissing_traits",
        "B\tmissing_trait_value\ttrait_right",
        "G\tmissing_from_trait_table\ttrait_left,trait_right",
        "H\tmissing_from_trait_table\ttrait_left,trait_right",
        "I\tmissing_from_tree\t",
    ]


def test_comparative_brownian_cli_reports_model_fit_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "brownian",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["tree_taxon_count"] == 4
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert math.isclose(payload["metrics"]["sigma_squared"], 4.774305191647407)
    assert math.isclose(payload["metrics"]["log_likelihood"], -6.0415206788726)
    assert "aic" in payload["metrics"]
    assert "aicc" in payload["metrics"]


def test_comparative_brownian_cli_writes_summary_and_exclusion_ledgers(
    tmp_path: Path, capsys
) -> None:
    summary_out = tmp_path / "brownian-trait-summary.tsv"
    excluded_out = tmp_path / "brownian-trait-excluded.tsv"
    exit_code = main(
        [
            "comparative",
            "brownian",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_brownian_missing.tsv")),
            "--trait",
            "response_growth",
            "--summary-out",
            str(summary_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    assert payload["metrics"]["excluded_taxon_count"] == 3
    assert summary_out.exists()
    assert excluded_out.exists()
    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("trait\ttaxon_column\ttree_taxon_count")
    assert excluded_rows == [
        "taxon\treason",
        "B\tmissing_trait_value",
        "C\tnon_numeric_trait_value",
        "G\tabsent_from_tree",
    ]


def test_comparative_brownian_regimes_cli_reports_model_fit_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "brownian-regimes",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            str(fixture("example_branch_regimes.tsv")),
            "--trait",
            "response",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["tree_taxon_count"] == 4
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["regime_count"] == 2
    assert payload["metrics"]["profile_row_count"] >= 162
    assert payload["metrics"]["better_model"] in {"brownian", "brownian-regimes"}
    assert 0.0 <= payload["metrics"]["likelihood_ratio_p_value"] <= 1.0


def test_comparative_brownian_regimes_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    summary_out = tmp_path / "brownian-regimes-summary.tsv"
    rates_out = tmp_path / "brownian-regimes-rates.tsv"
    comparison_out = tmp_path / "brownian-regimes-comparison.tsv"
    profile_out = tmp_path / "brownian-regimes-profile.tsv"
    branch_out = tmp_path / "brownian-regimes-branches.tsv"
    excluded_out = tmp_path / "brownian-regimes-excluded.tsv"
    exit_code = main(
        [
            "comparative",
            "brownian-regimes",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_brownian_missing.tsv")),
            str(fixture("example_branch_regimes_six_taxa.tsv")),
            "--trait",
            "response_growth",
            "--summary-out",
            str(summary_out),
            "--rates-out",
            str(rates_out),
            "--comparison-out",
            str(comparison_out),
            "--profile-out",
            str(profile_out),
            "--branches-out",
            str(branch_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    assert payload["metrics"]["excluded_taxon_count"] == 3
    assert summary_out.exists()
    assert rates_out.exists()
    assert comparison_out.exists()
    assert profile_out.exists()
    assert branch_out.exists()
    assert excluded_out.exists()
    assert (
        summary_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("trait\ttaxon_column\tbranch_id_column\tregime_column")
    )
    assert (
        rates_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("regime\tbranch_count\tcontributing_branch_count")
    )
    comparison_rows = comparison_out.read_text(encoding="utf-8").splitlines()
    assert comparison_rows[0].startswith("row_kind\tmodel\tcomparison_id")
    assert any(
        row.startswith("likelihood_ratio_test\t\tbrownian-vs-brownian-regimes\t")
        for row in comparison_rows[1:]
    )
    assert (
        profile_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("regime\tsigma_squared\tlog_likelihood")
    )
    assert (
        branch_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("branch_id\tregime\tbranch_length")
    )
    assert excluded_out.read_text(encoding="utf-8").splitlines() == [
        "taxon\treason",
        "B\tmissing_trait_value",
        "C\tnon_numeric_trait_value",
        "G\tabsent_from_tree",
    ]


def test_comparative_regime_map_cli_reconstructs_and_renders_regimes(
    tmp_path: Path, capsys
) -> None:
    summary_out = tmp_path / "trait-regime-summary.tsv"
    branch_out = tmp_path / "trait-regime-branches.tsv"
    node_out = tmp_path / "trait-regime-nodes.tsv"
    excluded_out = tmp_path / "trait-regime-excluded.tsv"
    svg_out = tmp_path / "trait-regime.svg"
    exit_code = main(
        [
            "comparative",
            "regime-map",
            str(fixture("example_tree.nwk")),
            "--table",
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--summary-out",
            str(summary_out),
            "--branches-out",
            str(branch_out),
            "--nodes-out",
            str(node_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--svg-out",
            str(svg_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["source_kind"] == "tip-state-reconstruction"
    assert payload["metrics"]["tree_taxon_count"] == 4
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["regime_count"] == 3
    assert payload["metrics"]["branch_count"] == 6
    assert payload["metrics"]["node_count"] == 7
    assert payload["metrics"]["ambiguous_branch_count"] == 1
    assert payload["metrics"]["rendered_internal_annotation_count"] >= 2
    assert payload["metrics"]["rendered_categorical_trait_count"] == 4
    assert (
        summary_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("source_kind\ttrait\ttaxon_column\treconstruction_model")
    )
    assert (
        branch_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("branch_id\tchild_node_name\tis_tip_branch")
    )
    assert (
        node_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("node_id\tnode_name\tis_tip")
    )
    assert excluded_out.read_text(encoding="utf-8").splitlines() == ["taxon\treason"]
    assert svg_out.exists()


def test_comparative_regime_map_cli_normalizes_user_map(tmp_path: Path, capsys) -> None:
    branch_out = tmp_path / "trait-regime-branches.tsv"
    svg_out = tmp_path / "trait-regime.svg"
    exit_code = main(
        [
            "comparative",
            "regime-map",
            str(fixture("example_tree.nwk")),
            "--regime-map",
            str(fixture("example_branch_regimes.tsv")),
            "--branches-out",
            str(branch_out),
            "--svg-out",
            str(svg_out),
            "--layout",
            "circular",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["source_kind"] == "user-provided-map"
    assert payload["metrics"]["regime_count"] == 2
    assert payload["metrics"]["branch_count"] == 6
    assert payload["metrics"]["node_count"] == 0
    assert payload["metrics"]["ambiguous_branch_count"] == 0
    assert payload["metrics"]["rendered_internal_annotation_count"] >= 2
    assert payload["metrics"]["rendered_categorical_trait_count"] == 4
    assert (
        branch_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("branch_id\tchild_node_name\tis_tip_branch")
    )
    assert svg_out.exists()


def test_comparative_ou_cli_reports_model_fit_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "ou",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["tree_taxon_count"] == 4
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert math.isclose(payload["metrics"]["alpha"], 33.333333)
    assert math.isclose(payload["metrics"]["theta"], 2.7503175510416416)
    assert math.isclose(payload["metrics"]["log_likelihood"], -5.260966005271784)
    assert payload["metrics"]["sigma_squared"] > 0.0
    assert "aic" in payload["metrics"]
    assert "aicc" in payload["metrics"]


def test_comparative_ou_cli_writes_summary_and_exclusion_ledgers(
    tmp_path: Path, capsys
) -> None:
    summary_out = tmp_path / "ou-trait-summary.tsv"
    excluded_out = tmp_path / "ou-trait-excluded.tsv"
    exit_code = main(
        [
            "comparative",
            "ou",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_continuous_evolution_missing.tsv")),
            "--trait",
            "response_growth",
            "--summary-out",
            str(summary_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    assert payload["metrics"]["excluded_taxon_count"] == 3
    assert summary_out.exists()
    assert excluded_out.exists()
    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("trait\ttaxon_column\ttree_taxon_count")
    assert excluded_rows == [
        "taxon\treason",
        "B\tmissing_trait_value",
        "C\tnon_numeric_trait_value",
        "G\tabsent_from_tree",
    ]


def test_comparative_early_burst_cli_reports_model_fit_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "early-burst",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["tree_taxon_count"] == 4
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["rate_change"] == 0.0
    assert payload["metrics"]["better_model"] == "brownian"
    assert payload["metrics"]["identifiability_warning_count"] == 3
    assert payload["metrics"]["profile_row_count"] == 161
    assert "aic" in payload["metrics"]
    assert "aicc" in payload["metrics"]


def test_comparative_early_burst_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    summary_out = tmp_path / "early-burst-summary.tsv"
    excluded_out = tmp_path / "early-burst-excluded.tsv"
    comparison_out = tmp_path / "early-burst-comparison.tsv"
    profile_out = tmp_path / "early-burst-profile.tsv"
    exit_code = main(
        [
            "comparative",
            "early-burst",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_continuous_evolution_missing.tsv")),
            "--trait",
            "response_growth",
            "--summary-out",
            str(summary_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--comparison-out",
            str(comparison_out),
            "--profile-out",
            str(profile_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    assert payload["metrics"]["excluded_taxon_count"] == 3
    assert summary_out.exists()
    assert excluded_out.exists()
    assert comparison_out.exists()
    assert profile_out.exists()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    comparison_rows = comparison_out.read_text(encoding="utf-8").splitlines()
    profile_rows = profile_out.read_text(encoding="utf-8").splitlines()
    assert excluded_rows == [
        "taxon\treason",
        "B\tmissing_trait_value",
        "C\tnon_numeric_trait_value",
        "G\tabsent_from_tree",
    ]
    assert comparison_rows[0].startswith("row_kind\tmodel\tcomparison_id")
    assert profile_rows[0].startswith("trait\trate_change\tlog_likelihood")


def test_comparative_rate_through_time_cli_reports_interval_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "rate-through-time",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--interval-count",
            "4",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["tree_taxon_count"] == 4
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["interval_count"] == 4
    assert payload["metrics"]["nonempty_interval_count"] >= 2
    assert payload["metrics"]["tree_depth"] > 0.0
    assert payload["metrics"]["trend_direction"] in {
        "slowdown",
        "acceleration",
        "stable",
        "insufficient_data",
    }


def test_comparative_rate_through_time_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    summary_out = tmp_path / "trait-rate-through-time-summary.tsv"
    intervals_out = tmp_path / "trait-rate-through-time-intervals.tsv"
    excluded_out = tmp_path / "trait-rate-through-time-excluded.tsv"
    exit_code = main(
        [
            "comparative",
            "rate-through-time",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_continuous_evolution_missing.tsv")),
            "--trait",
            "response_growth",
            "--interval-count",
            "4",
            "--summary-out",
            str(summary_out),
            "--intervals-out",
            str(intervals_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    assert payload["metrics"]["excluded_taxon_count"] == 3
    assert summary_out.exists()
    assert intervals_out.exists()
    assert excluded_out.exists()
    assert (
        summary_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("trait\ttaxon_column\ttree_taxon_count")
    )
    assert (
        intervals_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("interval_index\tstart_depth\tend_depth")
    )
    assert excluded_out.read_text(encoding="utf-8").splitlines() == [
        "taxon\treason",
        "B\tmissing_trait_value",
        "C\tnon_numeric_trait_value",
        "G\tabsent_from_tree",
    ]


def test_comparative_clade_traits_cli_reports_summary_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "clade-traits",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_clade_summary.tsv")),
            "--trait",
            "body_mass",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["tree_taxon_count"] == 6
    assert payload["metrics"]["analyzed_taxon_count"] == 6
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["trait_kind"] == "continuous"
    assert payload["metrics"]["clade_count"] == 4
    assert payload["metrics"]["exceptional_clade_count"] == 1
    assert payload["metrics"]["top_exceptional_clade"] == "E|F"
    assert payload["metrics"]["top_exceptionality_score"] > 0.0


def test_comparative_clade_traits_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    summary_out = tmp_path / "clade-traits-summary.tsv"
    clades_out = tmp_path / "clade-traits.tsv"
    excluded_out = tmp_path / "clade-traits-excluded.tsv"
    exit_code = main(
        [
            "comparative",
            "clade-traits",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_continuous_evolution_missing.tsv")),
            "--trait",
            "response_growth",
            "--trait-kind",
            "continuous",
            "--summary-out",
            str(summary_out),
            "--clades-out",
            str(clades_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    assert payload["metrics"]["excluded_taxon_count"] == 3
    assert payload["metrics"]["trait_kind"] == "continuous"
    assert summary_out.exists()
    assert clades_out.exists()
    assert excluded_out.exists()
    assert (
        summary_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("trait\ttaxon_column\ttrait_kind")
    )
    assert (
        clades_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("clade_id\tnode_label\ttrait_kind")
    )
    assert excluded_out.read_text(encoding="utf-8").splitlines() == [
        "taxon\treason",
        "B\tmissing_trait_value",
        "C\tnon_numeric_trait_value",
        "G\tabsent_from_tree",
    ]


def test_comparative_phylogenetic_residuals_cli_reports_outlier_metrics(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "phylogenetic-residuals",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_phylogenetic_residuals.tsv")),
            "--response",
            "brain_mass",
            "--predictor",
            "body_mass",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["tree_taxon_count"] == 6
    assert payload["metrics"]["analyzed_taxon_count"] == 6
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["method"] == "lambda"
    assert payload["metrics"]["outlier_count"] >= 1
    assert payload["metrics"]["top_outlier_taxon"] == "F"


def test_comparative_phylogenetic_residuals_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    summary_out = tmp_path / "phylogenetic-residual-summary.tsv"
    residuals_out = tmp_path / "phylogenetic-residuals.tsv"
    coefficients_out = tmp_path / "phylogenetic-residual-coefficients.tsv"
    excluded_out = tmp_path / "phylogenetic-residual-excluded.tsv"
    exit_code = main(
        [
            "comparative",
            "phylogenetic-residuals",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_phylogenetic_residuals_missing.tsv")),
            "--response",
            "brain_mass",
            "--predictor",
            "body_mass",
            "--method",
            "brownian",
            "--summary-out",
            str(summary_out),
            "--residuals-out",
            str(residuals_out),
            "--coefficients-out",
            str(coefficients_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["analyzed_taxon_count"] == 5
    assert payload["metrics"]["excluded_taxon_count"] == 2
    assert summary_out.exists()
    assert residuals_out.exists()
    assert coefficients_out.exists()
    assert excluded_out.exists()
    assert (
        summary_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("response\tpredictor\tmethod")
    )
    assert (
        residuals_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("taxon\tinput_order\ttree_tip_label\tobserved_value")
    )
    assert (
        coefficients_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("name\testimate\tstandard_error\tp_value")
    )
    assert excluded_out.read_text(encoding="utf-8").splitlines() == [
        "taxon\treason\tdetails",
        "E\tmissing_value\ttaxon is missing required value(s): brain_mass",
        "G\tabsent_from_tree\ttaxon is present in the trait table but absent from the tree",
    ]


def test_comparative_phylogenetic_anova_cli_reports_group_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "phylogenetic-anova",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_phylogenetic_anova.tsv")),
            "--response",
            "trait_value",
            "--group",
            "habitat",
            "--simulations",
            "16",
            "--seed",
            "7",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["tree_taxon_count"] == 6
    assert payload["metrics"]["analyzed_taxon_count"] == 6
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["group_count"] == 2
    assert payload["metrics"]["simulation_count"] == 16
    assert 0.0 <= payload["metrics"]["p_value"] <= 1.0
    assert payload["metrics"]["low_sample_group_count"] == 1


def test_comparative_phylogenetic_anova_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    summary_out = tmp_path / "phylogenetic-anova-summary.tsv"
    groups_out = tmp_path / "phylogenetic-anova-groups.tsv"
    pairwise_out = tmp_path / "phylogenetic-anova-pairwise.tsv"
    simulations_out = tmp_path / "phylogenetic-anova-simulations.tsv"
    excluded_out = tmp_path / "phylogenetic-anova-excluded.tsv"
    exit_code = main(
        [
            "comparative",
            "phylogenetic-anova",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_phylogenetic_anova_missing.tsv")),
            "--response",
            "trait_value",
            "--group",
            "habitat",
            "--simulations",
            "16",
            "--seed",
            "7",
            "--summary-out",
            str(summary_out),
            "--groups-out",
            str(groups_out),
            "--pairwise-out",
            str(pairwise_out),
            "--simulations-out",
            str(simulations_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["analyzed_taxon_count"] == 5
    assert payload["metrics"]["excluded_taxon_count"] == 2
    assert summary_out.exists()
    assert groups_out.exists()
    assert pairwise_out.exists()
    assert simulations_out.exists()
    assert excluded_out.exists()
    assert (
        summary_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("response\tgroup\ttaxon_column")
    )
    assert (
        groups_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("group\ttaxon_count\ttaxa\tmean")
    )
    assert (
        pairwise_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("left_group\tright_group\tleft_taxon_count")
    )
    assert (
        simulations_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("simulation_index\tf_statistic\tat_or_above_observed")
    )
    assert excluded_out.read_text(encoding="utf-8").splitlines() == [
        "taxon\treason\tdetails",
        "F\tmissing_value\ttaxon is missing required value(s): habitat",
        "G\tabsent_from_tree\ttaxon is present in the trait table but absent from the tree",
    ]


def test_comparative_trait_outliers_cli_reports_top_outlier_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "trait-outliers",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_trait_outliers.tsv")),
            "--trait",
            "body_mass",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["tree_taxon_count"] == 6
    assert payload["metrics"]["analyzed_taxon_count"] == 6
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["selected_model"] in {"brownian", "ou"}
    assert payload["metrics"]["outlier_count"] == 1
    assert payload["metrics"]["top_outlier_taxon"] == "F"
    assert payload["metrics"]["top_abs_standardized_residual"] >= 2.0


def test_comparative_trait_outliers_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    summary_out = tmp_path / "trait-outlier-summary.tsv"
    outliers_out = tmp_path / "trait-outliers.tsv"
    excluded_out = tmp_path / "trait-outlier-excluded.tsv"
    exit_code = main(
        [
            "comparative",
            "trait-outliers",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_continuous_evolution_missing.tsv")),
            "--trait",
            "response_growth",
            "--summary-out",
            str(summary_out),
            "--outliers-out",
            str(outliers_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["analyzed_taxon_count"] == 4
    assert payload["metrics"]["excluded_taxon_count"] == 3
    assert summary_out.exists()
    assert outliers_out.exists()
    assert excluded_out.exists()
    assert (
        summary_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("trait\ttaxon_column\ttree_taxon_count")
    )
    assert (
        outliers_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("taxon\tobserved_value\tconditional_expected_value\tresidual")
    )
    assert excluded_out.read_text(encoding="utf-8").splitlines() == [
        "taxon\treason",
        "B\tmissing_trait_value",
        "C\tnon_numeric_trait_value",
        "G\tabsent_from_tree",
    ]


def test_comparative_trait_imputation_cli_reports_holdout_metrics(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "trait-imputation",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_trait_imputation.tsv")),
            "--trait",
            "body_mass",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["tree_taxon_count"] == 6
    assert payload["metrics"]["observed_taxon_count"] == 5
    assert payload["metrics"]["imputed_taxon_count"] == 1
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["holdout_validation_status"] == "performed"
    assert payload["metrics"]["holdout_count"] == 5
    assert payload["metrics"]["holdout_mean_absolute_error"] is not None
    assert 0.0 <= payload["metrics"]["holdout_interval_coverage"] <= 1.0


def test_comparative_trait_imputation_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    summary_out = tmp_path / "trait-imputation-summary.tsv"
    imputations_out = tmp_path / "trait-imputations.tsv"
    holdout_out = tmp_path / "trait-imputation-holdout.tsv"
    excluded_out = tmp_path / "trait-imputation-excluded.tsv"
    exit_code = main(
        [
            "comparative",
            "trait-imputation",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_continuous_evolution_missing.tsv")),
            "--trait",
            "response_growth",
            "--summary-out",
            str(summary_out),
            "--imputations-out",
            str(imputations_out),
            "--holdout-out",
            str(holdout_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["observed_taxon_count"] == 4
    assert payload["metrics"]["imputed_taxon_count"] == 1
    assert payload["metrics"]["excluded_taxon_count"] == 2
    assert summary_out.exists()
    assert imputations_out.exists()
    assert holdout_out.exists()
    assert excluded_out.exists()
    assert (
        summary_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("trait\ttaxon_column\tmodel")
    )
    assert (
        imputations_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("taxon\tmissing_reason\tobserved_support_taxon_count")
    )
    assert (
        holdout_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("taxon\tobserved_value\tpredicted_value\tresidual")
    )
    assert excluded_out.read_text(encoding="utf-8").splitlines() == [
        "taxon\treason",
        "C\tnon_numeric_trait_value",
        "G\tabsent_from_tree",
    ]


def test_comparative_covariance_audit_cli_reports_pgls_profile_metrics(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "covariance-audit",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--analysis",
            "pgls",
            "--response",
            "response",
            "--predictors",
            "predictor_one",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["analysis"] == "pgls"
    assert payload["metrics"]["covariance_model"] == "pagel-lambda"
    assert payload["metrics"]["matrix_dimension"] == 4
    assert payload["metrics"]["matrix_rank"] == 4
    assert payload["metrics"]["fit_strategy"] == "exact"
    assert payload["metrics"]["singular"] is False
    assert payload["metrics"]["near_singular"] is False
    assert payload["metrics"]["matched_taxon_count"] == 4
    assert payload["metrics"]["missing_from_traits_count"] == 0
    assert payload["metrics"]["extra_trait_taxon_count"] == 0
    assert payload["metrics"]["analysis_taxon_count"] == 4
    assert payload["metrics"]["zero_length_branch_count"] == 0
    assert payload["metrics"]["negative_branch_length_count"] == 0
    assert payload["metrics"]["candidate_row_count"] == 101
    assert payload["metrics"]["blocker_count"] == 0


def test_comparative_covariance_audit_cli_reports_duplicate_trait_blockers(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "covariance-audit",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative_duplicate.tsv")),
            "--analysis",
            "pgls",
            "--response",
            "response",
            "--predictors",
            "predictor_one",
            "--lambda-value",
            "1.0",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["fit_strategy"] == "failure"
    assert payload["metrics"]["duplicate_trait_taxon_count"] == 1
    assert payload["metrics"]["candidate_row_count"] == 0
    assert payload["metrics"]["blocker_count"] == 1
    assert payload["data"]["duplicate_trait_taxa"] == ["A"]
    assert payload["data"]["blockers"] == ["trait table contains duplicate taxon keys"]


def test_comparative_covariance_audit_cli_reports_taxon_overlap_metrics(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "covariance-audit",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits.tsv")),
            "--analysis",
            "brownian-trait",
            "--trait",
            "value",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["fit_strategy"] == "exact"
    assert payload["metrics"]["matched_taxon_count"] == 3
    assert payload["metrics"]["missing_from_traits_count"] == 1
    assert payload["metrics"]["extra_trait_taxon_count"] == 1
    assert payload["metrics"]["analysis_taxon_count"] == 3
    assert payload["data"]["missing_from_traits"] == ["D"]
    assert payload["data"]["extra_trait_taxa"] == ["E"]


def test_comparative_covariance_audit_cli_reports_zero_length_regularization(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "covariance-audit",
            str(fixture("example_tree_zero_lengths.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--analysis",
            "pgls",
            "--response",
            "response",
            "--predictors",
            "predictor_one",
            "--lambda-value",
            "1.0",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["fit_strategy"] == "regularization"
    assert payload["metrics"]["zero_length_branch_count"] == 3
    assert payload["metrics"]["negative_branch_length_count"] == 0
    assert payload["metrics"]["matrix_rank"] == 3
    assert payload["metrics"]["singular"] is True
    assert payload["metrics"]["near_singular"] is True
    assert "tree contains zero-length branches" in payload["warnings"]


def test_comparative_covariance_audit_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    summary_out = tmp_path / "covariance-audit-summary.tsv"
    candidates_out = tmp_path / "covariance-audit-candidates.tsv"
    excluded_out = tmp_path / "covariance-audit-excluded.tsv"
    exit_code = main(
        [
            "comparative",
            "covariance-audit",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative_missing_predictor.tsv")),
            "--analysis",
            "pgls",
            "--response",
            "response",
            "--predictors",
            "predictor_one",
            "--lambda-value",
            "1.0",
            "--summary-out",
            str(summary_out),
            "--candidates-out",
            str(candidates_out),
            "--excluded-taxa-out",
            str(excluded_out),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["candidate_row_count"] == 1
    assert payload["metrics"]["blocker_count"] == 0
    assert summary_out.exists()
    assert candidates_out.exists()
    assert excluded_out.exists()
    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    candidate_rows = candidates_out.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("analysis\tcovariance_model\tanalysis_label")
    assert candidate_rows[0].startswith("candidate_label\tparameter_name")
    assert excluded_rows[0] == "taxon\treason\tdetails"
    assert len(summary_rows) == 2
    assert len(candidate_rows) == 2
    assert len(excluded_rows) == 2


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
    assert contrast_rows[0].startswith("trait\tnode_id\tnode\tleft_taxa\tright_taxa")
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
    assert payload["metrics"]["aic"] == payload["data"]["model"]["aic"]
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
    assert payload["metrics"]["aic"] == payload["data"]["model"]["aic"]
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
    assert any(
        "habitat\thabitat\treference-level\tforest\tforest\ttrue" in row
        for row in written_rows[1:]
    )


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
    assert math.isclose(coefficients["predictor_one"], 0.9087628025668772, abs_tol=1e-6)


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
    assert payload["metrics"]["warning_count"] == 2
    assert payload["metrics"]["method_tier"] == "experimental"
    assert payload["metrics"]["method_approximation"] == (
        "phylogenetic-working-correlation-gee"
    )
    assert payload["metrics"]["method_excluded_reference_surfaces"] == [
        "ape::compar.gee"
    ]
    assert payload["warnings"][0].startswith("experimental method tier:")
    assert any("ape::compar.gee parity" in warning for warning in payload["warnings"])
    assert payload["metrics"]["coefficient_inference_distribution"] == "wald-normal"
    assert coefficients["body_size"] > 0.0


def test_comparative_logistic_cli_writes_review_ledgers(tmp_path: Path, capsys) -> None:
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
    assert payload["metrics"]["method_tier"] == "experimental"
    assert any(
        warning.startswith("experimental method tier:")
        for warning in payload["warnings"]
    )
    assert any("ape::compar.gee parity" in warning for warning in payload["warnings"])
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
    assert tree_rows[0].startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id"
    )
    assert coefficient_rows[0].startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id\tterm"
    )
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
    assert payload["metrics"]["residual_covariance_response_count"] == 2
    assert payload["metrics"]["residual_covariance_matrix_rank"] == 1
    assert math.isinf(payload["metrics"]["residual_covariance_condition_number"])
    assert payload["metrics"]["residual_covariance_singular"] is True
    assert payload["metrics"]["residual_covariance_near_singular"] is True
    assert payload["metrics"]["residual_covariance_row_count"] == 4
    assert payload["metrics"]["residual_correlation_row_count"] == 4
    assert payload["metrics"]["residual_association_count"] == 1
    assert payload["metrics"]["response_model_count"] == 2
    assert payload["metrics"]["coefficient_row_count"] == 6
    assert len(payload["data"]["response_models"]) == 2
    assert payload["data"]["missing_value_policy"] == (
        "shared_complete_case_across_responses_and_predictor_terms"
    )
    assert payload["data"]["numerical_tolerance"] == 1e-12


def test_comparative_multivariate_cli_reports_full_rank_covariance_metrics(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "multivariate",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_comparative_multivariate_full_rank.tsv")),
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
    assert payload["metrics"]["residual_covariance_response_count"] == 2
    assert payload["metrics"]["residual_covariance_matrix_rank"] == 2
    assert math.isfinite(payload["metrics"]["residual_covariance_condition_number"])
    assert payload["metrics"]["residual_covariance_singular"] is False
    assert payload["metrics"]["residual_covariance_near_singular"] is False


def test_comparative_multivariate_cli_reports_heterogeneous_lambda_warning(
    capsys,
) -> None:
    exit_code = main(
        [
            "comparative",
            "multivariate",
            str(fixture("example_tree_six_taxa.nwk")),
            str(
                fixture(
                    "example_traits_comparative_multivariate_heterogeneous_lambda.tsv"
                )
            ),
            "--responses",
            "response_growth",
            "response_range",
            "--predictors",
            "predictor_one",
            "predictor_two",
            "--lambda-value",
            "estimate",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["warning_count"] == 1
    assert payload["data"]["response_model_rows"][0]["lambda_value"] == 1.0
    assert payload["data"]["response_model_rows"][1]["lambda_value"] == 0.0
    assert payload["data"]["warnings"] == [
        "response models resolved materially different Pagel lambda values (0 to 1), so shared residual covariance and correlation compare residuals fit under different phylogenetic error assumptions"
    ]


def test_comparative_multivariate_cli_writes_review_ledgers(
    tmp_path: Path, capsys
) -> None:
    models_out = tmp_path / "multivariate-response-models.tsv"
    coefficients_out = tmp_path / "multivariate-response-coefficients.tsv"
    covariance_out = tmp_path / "multivariate-residual-covariance.tsv"
    correlation_out = tmp_path / "multivariate-residual-correlation.tsv"
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
            "--response-models-out",
            str(models_out),
            "--coefficients-out",
            str(coefficients_out),
            "--covariance-out",
            str(covariance_out),
            "--correlation-out",
            str(correlation_out),
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
    assert payload["metrics"]["response_model_count"] == 2
    assert payload["metrics"]["coefficient_row_count"] == 6
    assert payload["metrics"]["residual_correlation_row_count"] == 4
    assert models_out.exists()
    assert coefficients_out.exists()
    assert covariance_out.exists()
    assert correlation_out.exists()
    assert associations_out.exists()
    assert excluded_out.exists()
    model_rows = models_out.read_text(encoding="utf-8").splitlines()
    coefficient_rows = coefficients_out.read_text(encoding="utf-8").splitlines()
    covariance_rows = covariance_out.read_text(encoding="utf-8").splitlines()
    correlation_rows = correlation_out.read_text(encoding="utf-8").splitlines()
    association_rows = associations_out.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    assert model_rows[0].startswith("response\tformula\tpredictor_term_count")
    assert coefficient_rows[0].startswith("response\tformula\tterm\testimate")
    assert covariance_rows[0].startswith(
        "left_response\tright_response\tpair_count\tis_diagonal"
    )
    assert correlation_rows[0].startswith(
        "left_response\tright_response\tpair_count\tis_diagonal\tcorrelation"
    )
    assert association_rows[0].startswith(
        "left_response\tright_response\tpair_count\tcovariance\tcorrelation"
    )
    assert excluded_rows[0] == (
        "taxon\treason\tmissing_columns\tblocking_responses\tdetails"
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


def test_comparative_report_cli_can_export_full_review_package(
    tmp_path: Path, capsys
) -> None:
    out_dir = tmp_path / "comparative-report-package"
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
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["taxon_count"] == 4
    assert payload["metrics"]["coefficient_count"] >= 2
    assert payload["metrics"]["package_output_count"] == 12
    assert payload["data"]["output_dir"] == str(out_dir)
    assert (out_dir / "comparative-report.html").exists()
    assert (out_dir / "comparative-methods-summary.md").exists()
    assert (out_dir / "reviewer-audit-checklist.tsv").exists()
    assert (out_dir / "comparative-summary.tsv").exists()
    assert (out_dir / "coefficient-table.tsv").exists()
    assert (out_dir / "residual-summary.tsv").exists()
    assert (out_dir / "signal-summary.tsv").exists()
    assert (out_dir / "model-comparison.tsv").exists()
    assert (out_dir / "interpretation-table.tsv").exists()
    assert (out_dir / "audit-table.tsv").exists()
    assert (out_dir / "contrast-table.tsv").exists()
    assert (out_dir / "comparative-report.manifest.json").exists()
    summary_rows = (
        (out_dir / "comparative-summary.tsv").read_text(encoding="utf-8").splitlines()
    )
    assert summary_rows[0].startswith("response\tformula\tpredictor_count")


def test_comparative_report_cli_can_export_methods_summary_markdown(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "comparative-methods-summary.md"
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
            "--methods-summary-out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["package_output_count"] == 1
    assert output_path.exists()
    text = output_path.read_text(encoding="utf-8")
    assert "Comparative Analysis Methods Summary" in text
    assert "- predictor terms: `predictor_one`" in text


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
