from __future__ import annotations

import math
from pathlib import Path
from types import SimpleNamespace

import pytest

from bijux_phylogenetics.ancestral import (
    reconstruct_continuous_evolutionary_mode_states,
)
from bijux_phylogenetics.comparative import (
    ContinuousModeSearchControls,
    compare_continuous_evolutionary_modes,
    compare_fitcontinuous_model_ranking,
    fit_continuous_evolutionary_mode,
    rescale_tree_early_burst,
    rescale_tree_ornstein_uhlenbeck,
    rescale_tree_pagel_delta,
    rescale_tree_pagel_kappa,
    rescale_tree_pagel_lambda,
    rescale_tree_white_noise,
)
from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    log_determinant,
    quadratic_form,
)
from bijux_phylogenetics.comparative.common import (
    build_brownian_covariance_matrix,
    load_comparative_dataset,
    stable_covariance,
)
from bijux_phylogenetics.comparative.continuous.model_fitting import (
    ComparativeResidualSummary,
)
import bijux_phylogenetics.comparative.evolutionary_modes as evolutionary_modes_module
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_geiger_continuous_fixture,
)
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures"
EXAMPLE_TREE = FIXTURE_ROOT / "trees" / "example_tree.nwk"
EXAMPLE_TREE_NEGATIVE = FIXTURE_ROOT / "trees" / "example_tree_negative_length.nwk"
EXAMPLE_TRAITS = FIXTURE_ROOT / "metadata" / "example_traits_comparative.tsv"


def _analytical_gaussian_intercept_fit(
    trait_values: list[float],
    covariance: list[list[float]],
) -> tuple[float, float, float]:
    inverse_covariance = invert_matrix(covariance)
    ones = [1.0] * len(trait_values)
    denominator = quadratic_form(ones, inverse_covariance)
    theta = (
        sum(
            ones[row_index]
            * sum(
                inverse_covariance[row_index][column_index] * trait_values[column_index]
                for column_index in range(len(trait_values))
            )
            for row_index in range(len(trait_values))
        )
        / denominator
    )
    residuals = [value - theta for value in trait_values]
    sigma_squared = quadratic_form(residuals, inverse_covariance) / len(trait_values)
    log_likelihood = -0.5 * (
        len(trait_values) * math.log(2.0 * math.pi * sigma_squared)
        + log_determinant(covariance)
        + len(trait_values)
    )
    return theta, sigma_squared, log_likelihood


def test_rescale_tree_ornstein_uhlenbeck_reports_deterministic_branch_lengths() -> None:
    report = rescale_tree_ornstein_uhlenbeck(EXAMPLE_TREE, alpha=1.0)

    assert report.mode == "ornstein-uhlenbeck"
    assert report.parameter_name == "alpha"
    assert report.tip_count == 4
    assert math.isclose(report.transformed_total_branch_length, 0.706662964349163)
    assert report.branch_rows[0].node == "A"
    assert math.isclose(report.branch_rows[0].parent_depth, 0.2)
    assert math.isclose(report.branch_rows[0].child_depth, 0.3)
    assert math.isclose(
        report.branch_rows[0].transformed_branch_length, 0.090634623461009
    )


def test_rescale_tree_early_burst_zero_matches_original_tree_length() -> None:
    baseline = rescale_tree_early_burst(EXAMPLE_TREE, rate_change=0.0)

    assert baseline.mode == "early-burst"
    assert math.isclose(
        baseline.original_total_branch_length,
        baseline.transformed_total_branch_length,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_rescale_tree_early_burst_supports_negative_rate_change() -> None:
    report = rescale_tree_early_burst(EXAMPLE_TREE, rate_change=-2.0)

    assert report.mode == "early-burst"
    assert report.parameter_name == "rate_change"
    assert math.isclose(report.parameter_value, -2.0)
    assert report.transformed_total_branch_length > 0.0


def test_rescale_tree_pagel_lambda_reports_deterministic_branch_lengths() -> None:
    report = rescale_tree_pagel_lambda(EXAMPLE_TREE, lambda_value=0.5)

    assert report.mode == "pagel-lambda"
    assert report.parameter_name == "lambda"
    assert report.tip_count == 4
    assert math.isclose(report.parameter_value, 0.5)
    assert math.isclose(report.transformed_total_branch_length, 1.05)
    assert report.branch_rows[0].node == "A"
    assert math.isclose(report.branch_rows[0].original_branch_length, 0.1)
    assert math.isclose(report.branch_rows[0].transformed_branch_length, 0.2)


def test_rescale_tree_pagel_kappa_reports_deterministic_branch_lengths() -> None:
    report = rescale_tree_pagel_kappa(EXAMPLE_TREE, kappa=0.5)

    assert report.mode == "pagel-kappa"
    assert report.parameter_name == "kappa"
    assert report.tip_count == 4
    assert math.isclose(report.parameter_value, 0.5)
    assert math.isclose(report.transformed_total_branch_length, 2.290324084550388)
    assert report.branch_rows[0].node == "A"
    assert math.isclose(report.branch_rows[0].original_branch_length, 0.1)
    assert math.isclose(
        report.branch_rows[0].transformed_branch_length, 0.316227766016838
    )


def test_rescale_tree_pagel_delta_reports_deterministic_branch_lengths() -> None:
    report = rescale_tree_pagel_delta(EXAMPLE_TREE, delta=0.5)

    assert report.mode == "pagel-delta"
    assert report.parameter_name == "delta"
    assert report.tip_count == 4
    assert math.isclose(report.parameter_value, 0.5)
    assert math.isclose(report.transformed_total_branch_length, 0.781845944964795)
    assert report.branch_rows[0].node == "A"
    assert math.isclose(report.branch_rows[0].original_branch_length, 0.1)
    assert math.isclose(
        report.branch_rows[0].transformed_branch_length, 0.055051025721682
    )


def test_rescale_tree_white_noise_reports_deterministic_branch_lengths() -> None:
    report = rescale_tree_white_noise(EXAMPLE_TREE, sigsq=2.5)

    assert report.mode == "white-noise"
    assert report.parameter_name == "sigsq"
    assert report.tip_count == 4
    assert math.isclose(report.parameter_value, 2.5)
    assert math.isclose(report.transformed_total_branch_length, 10.0)
    assert report.branch_rows[0].node == "A"
    assert math.isclose(report.branch_rows[0].original_branch_length, 0.1)
    assert math.isclose(report.branch_rows[0].transformed_branch_length, 2.5)


def test_rescale_tree_ornstein_uhlenbeck_rejects_negative_alpha() -> None:
    with pytest.raises(ComparativeMethodError, match="OU alpha must be non-negative"):
        rescale_tree_ornstein_uhlenbeck(EXAMPLE_TREE, alpha=-0.5)


def test_rescale_tree_pagel_lambda_rejects_negative_branch_lengths() -> None:
    with pytest.raises(
        ComparativeMethodError,
        match="Pagel lambda cannot transform negative branch lengths",
    ):
        rescale_tree_pagel_lambda(EXAMPLE_TREE_NEGATIVE, lambda_value=0.5)


def test_rescale_tree_pagel_kappa_rejects_negative_branch_lengths() -> None:
    with pytest.raises(
        ComparativeMethodError,
        match="Pagel kappa cannot transform negative branch lengths",
    ):
        rescale_tree_pagel_kappa(EXAMPLE_TREE_NEGATIVE, kappa=0.5)


def test_rescale_tree_pagel_delta_rejects_negative_branch_lengths() -> None:
    with pytest.raises(
        ComparativeMethodError,
        match="Pagel delta cannot transform negative branch lengths",
    ):
        rescale_tree_pagel_delta(EXAMPLE_TREE_NEGATIVE, delta=0.5)


def test_rescale_tree_early_burst_rejects_negative_branch_lengths() -> None:
    with pytest.raises(
        ComparativeMethodError,
        match="Early-burst rescaling cannot transform negative branch lengths",
    ):
        rescale_tree_early_burst(EXAMPLE_TREE_NEGATIVE, rate_change=-2.0)


def test_rescale_tree_white_noise_rejects_negative_branch_lengths() -> None:
    with pytest.raises(
        ComparativeMethodError,
        match="White-noise rescaling cannot transform negative branch lengths",
    ):
        rescale_tree_white_noise(EXAMPLE_TREE_NEGATIVE, sigsq=1.0)


def test_rescale_tree_white_noise_rejects_negative_sigsq() -> None:
    with pytest.raises(
        ComparativeMethodError,
        match="White-noise sigsq must be non-negative",
    ):
        rescale_tree_white_noise(EXAMPLE_TREE, sigsq=-1.0)


def test_fit_continuous_evolutionary_mode_supports_pagel_delta_strong_signal() -> None:
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_brownian_signal_twenty_four_taxa"
    )

    fit = fit_continuous_evolutionary_mode(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        mode="pagel-delta",
        taxon_column=fixture.taxon_column,
    )

    assert fit.mode == "pagel-delta"
    assert fit.parameter_name == "delta"
    assert fit.parameter_value is not None
    assert 1.3 < fit.parameter_value < 1.7
    assert fit.optimizer_diagnostics is not None
    assert fit.optimizer_diagnostics.hit_lower_boundary is False
    assert fit.optimizer_diagnostics.hit_upper_boundary is False
    assert fit.aicc >= fit.aic
    assert "depth transformation" in fit.assumptions[0]
    assert [warning.kind for warning in fit.identifiability_warnings] == [
        "flat_likelihood"
    ]


def test_fit_continuous_evolutionary_mode_supports_pagel_delta_weak_signal() -> None:
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_white_noise_twenty_four_taxa"
    )

    fit = fit_continuous_evolutionary_mode(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        mode="pagel-delta",
        taxon_column=fixture.taxon_column,
    )

    assert fit.mode == "pagel-delta"
    assert fit.parameter_name == "delta"
    assert fit.parameter_value is not None
    assert math.isclose(fit.parameter_value, 3.0, abs_tol=1e-12)
    assert fit.optimizer_diagnostics is not None
    assert fit.optimizer_diagnostics.hit_lower_boundary is False
    assert fit.optimizer_diagnostics.hit_upper_boundary is True
    assert [warning.kind for warning in fit.identifiability_warnings] == [
        "boundary_delta",
        "flat_likelihood",
        "late_change_limit",
    ]
    assert fit.boundary_assessment is not None
    assert fit.boundary_assessment.affected_parameter == "delta"
    assert fit.boundary_assessment.hit_upper_boundary is True
    assert fit.boundary_assessment.near_upper_boundary is True
    assert fit.boundary_assessment.flat_likelihood_near_boundary is True
    assert fit.boundary_assessment.stable_conclusion_supported is False


def test_fit_continuous_evolutionary_mode_supports_white_noise_baseline() -> None:
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_white_noise_twenty_four_taxa"
    )

    fit = fit_continuous_evolutionary_mode(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        mode="white-noise",
        taxon_column=fixture.taxon_column,
    )

    assert fit.mode == "white-noise"
    assert fit.parameter_name is None
    assert fit.parameter_value is None
    assert fit.optimizer_diagnostics is None
    assert fit.aicc >= fit.aic
    assert fit.identifiability_warnings[0].kind == "no_phylogenetic_correlation"
    assert "identity covariance" in fit.assumptions[0]


def test_fit_continuous_evolutionary_mode_white_noise_fits_high_signal_worse_than_brownian() -> (
    None
):
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_brownian_signal_twenty_four_taxa"
    )

    brownian = fit_continuous_evolutionary_mode(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        mode="brownian",
        taxon_column=fixture.taxon_column,
    )
    white = fit_continuous_evolutionary_mode(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        mode="white-noise",
        taxon_column=fixture.taxon_column,
    )

    assert white.log_likelihood < brownian.log_likelihood
    assert white.aic > brownian.aic


def test_fit_continuous_evolutionary_mode_brownian_log_likelihood_matches_analytical_known_answer() -> (
    None
):
    dataset = load_comparative_dataset(
        EXAMPLE_TREE,
        EXAMPLE_TRAITS,
        trait="response",
        taxon_column="taxon",
        minimum_taxa=3,
        require_rooted=True,
        require_binary=False,
    )
    covariance = stable_covariance(
        build_brownian_covariance_matrix(dataset.tree, dataset.taxa)
    )
    theta, sigma_squared, log_likelihood = _analytical_gaussian_intercept_fit(
        dataset.trait_values,
        covariance,
    )

    fit = fit_continuous_evolutionary_mode(
        EXAMPLE_TREE,
        EXAMPLE_TRAITS,
        trait="response",
        mode="brownian",
        taxon_column="taxon",
    )

    assert fit.likelihood_constant_policy == (
        "full-gaussian-loglikelihood-includes-normalizing-constant"
    )
    assert fit.likelihood_comparison_policy == (
        "raw-loglikelihood-and-derived-aic-are-directly-comparable-when-the-shared-gaussian-constant-policy-matches"
    )
    assert math.isclose(fit.root_state, theta, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(fit.rate, sigma_squared, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(
        fit.log_likelihood,
        log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_fit_continuous_evolutionary_mode_white_noise_log_likelihood_matches_analytical_known_answer() -> (
    None
):
    dataset = load_comparative_dataset(
        EXAMPLE_TREE,
        EXAMPLE_TRAITS,
        trait="response",
        taxon_column="taxon",
        minimum_taxa=3,
        require_rooted=True,
        require_binary=False,
    )
    identity_covariance = stable_covariance(
        [
            [
                1.0 if row_index == column_index else 0.0
                for column_index in range(len(dataset.taxa))
            ]
            for row_index in range(len(dataset.taxa))
        ]
    )
    theta, sigma_squared, log_likelihood = _analytical_gaussian_intercept_fit(
        dataset.trait_values,
        identity_covariance,
    )

    fit = fit_continuous_evolutionary_mode(
        EXAMPLE_TREE,
        EXAMPLE_TRAITS,
        trait="response",
        mode="white-noise",
        taxon_column="taxon",
    )

    assert fit.likelihood_constant_policy == (
        "full-gaussian-loglikelihood-includes-normalizing-constant"
    )
    assert fit.likelihood_comparison_policy == (
        "raw-loglikelihood-and-derived-aic-are-directly-comparable-when-the-shared-gaussian-constant-policy-matches"
    )
    assert math.isclose(fit.root_state, theta, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(fit.rate, sigma_squared, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(
        fit.log_likelihood,
        log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_fit_continuous_evolutionary_mode_white_noise_handles_missing_values() -> None:
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_missing_values_twenty_four_taxa"
    )

    fit = fit_continuous_evolutionary_mode(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        mode="white-noise",
        taxon_column=fixture.taxon_column,
    )

    assert fit.mode == "white-noise"
    assert fit.taxon_count == 22
    assert fit.optimizer_diagnostics is None
    assert fit.identifiability_warnings[0].kind == "no_phylogenetic_correlation"


def test_fit_continuous_evolutionary_mode_records_bounded_search_controls() -> None:
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_brownian_signal_twenty_four_taxa"
    )

    fit = fit_continuous_evolutionary_mode(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        mode="pagel-lambda",
        taxon_column=fixture.taxon_column,
        search_controls=ContinuousModeSearchControls(
            coarse_grid_point_count=41,
            fine_grid_point_count=31,
            initial_parameter_value=0.35,
        ),
        lambda_bounds=(0.2, 0.6),
    )

    assert fit.parameter_value is not None
    assert math.isclose(fit.parameter_value, 0.6, abs_tol=1e-12)
    assert fit.optimizer_diagnostics is not None
    assert fit.optimizer_diagnostics.parameter_search_strategy == (
        "bounded-two-stage-grid-search"
    )
    assert fit.optimizer_diagnostics.starting_parameter_policy == (
        "user-provided-first-evaluation"
    )
    assert math.isclose(
        fit.optimizer_diagnostics.starting_parameter_value,
        0.35,
        abs_tol=1e-12,
    )
    assert fit.optimizer_diagnostics.coarse_grid_point_count == 41
    assert fit.optimizer_diagnostics.fine_grid_point_count == 31
    assert fit.optimizer_diagnostics.function_evaluation_count >= 71
    assert fit.optimizer_diagnostics.hit_upper_boundary is True


def test_fit_continuous_evolutionary_mode_rejects_out_of_bounds_initial_parameter_value() -> (
    None
):
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_brownian_signal_twenty_four_taxa"
    )

    with pytest.raises(
        ComparativeMethodError,
        match="initial_parameter_value must fall within the declared bounded search interval",
    ):
        fit_continuous_evolutionary_mode(
            fixture.tree_path,
            fixture.traits_path,
            trait=fixture.trait_name,
            mode="pagel-lambda",
            taxon_column=fixture.taxon_column,
            search_controls=ContinuousModeSearchControls(
                initial_parameter_value=0.8,
            ),
            lambda_bounds=(0.2, 0.6),
        )


def test_fit_continuous_evolutionary_mode_rejects_invalid_grid_control_counts() -> None:
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_brownian_signal_twenty_four_taxa"
    )

    with pytest.raises(
        ComparativeMethodError,
        match="coarse_grid_point_count must be at least 2",
    ):
        fit_continuous_evolutionary_mode(
            fixture.tree_path,
            fixture.traits_path,
            trait=fixture.trait_name,
            mode="pagel-kappa",
            taxon_column=fixture.taxon_column,
            search_controls=ContinuousModeSearchControls(coarse_grid_point_count=1),
        )


def test_fit_continuous_evolutionary_mode_rejects_search_controls_for_white_noise() -> (
    None
):
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_white_noise_twenty_four_taxa"
    )

    with pytest.raises(
        ComparativeMethodError,
        match="white-noise mode does not expose bounded parameter-search controls",
    ):
        fit_continuous_evolutionary_mode(
            fixture.tree_path,
            fixture.traits_path,
            trait=fixture.trait_name,
            mode="white-noise",
            taxon_column=fixture.taxon_column,
            search_controls=ContinuousModeSearchControls(),
        )


def test_fit_continuous_evolutionary_mode_explicitly_excludes_standard_error_review() -> (
    None
):
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_standard_error_review_twenty_four_taxa"
    )

    with pytest.raises(
        ComparativeMethodError,
        match="standard-error parity is explicitly excluded in this round",
    ):
        fit_continuous_evolutionary_mode(
            fixture.tree_path,
            fixture.traits_path,
            trait=fixture.trait_name,
            mode="ornstein-uhlenbeck",
            taxon_column=fixture.taxon_column,
            standard_error_trait=fixture.standard_error_trait_name,
        )


def test_compare_continuous_evolutionary_modes_explicitly_excludes_standard_error_review() -> (
    None
):
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_standard_error_review_twenty_four_taxa"
    )

    with pytest.raises(
        ComparativeMethodError,
        match="standard-error parity is explicitly excluded in this round",
    ):
        compare_continuous_evolutionary_modes(
            fixture.tree_path,
            fixture.traits_path,
            trait=fixture.trait_name,
            taxon_column=fixture.taxon_column,
            standard_error_trait=fixture.standard_error_trait_name,
        )


def test_compare_fitcontinuous_model_ranking_blocks_mixed_likelihood_constant_policies(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_load_dataset(*args, **kwargs):
        return SimpleNamespace(
            tree_path=EXAMPLE_TREE,
            traits_path=EXAMPLE_TRAITS,
            trait="response",
            taxon_column="taxon",
            taxa=["A", "B", "C", "D"],
            readiness=SimpleNamespace(tree_taxa=4),
        )

    def fake_fit(dataset, *, mode, **kwargs):
        policy = (
            "full-gaussian-loglikelihood-includes-normalizing-constant"
            if mode == "brownian"
            else "shifted-gaussian-loglikelihood-without-normalizing-constant"
        )
        return evolutionary_modes_module.ContinuousEvolutionaryModeFitReport(
            tree_path=EXAMPLE_TREE,
            traits_path=EXAMPLE_TRAITS,
            taxon_column="taxon",
            trait="response",
            taxon_count=4,
            taxa=["A", "B", "C", "D"],
            mode=mode,
            parameter_name=None,
            parameter_value=None,
            root_state=0.0,
            rate=1.0,
            log_likelihood=-4.0 if mode == "brownian" else -3.0,
            aic=12.0 if mode == "brownian" else 10.0,
            aicc=24.0 if mode == "brownian" else 22.0,
            likelihood_constant_policy=policy,
            likelihood_comparison_policy="synthetic-test-policy",
            fitted_values=[0.0, 0.0, 0.0, 0.0],
            residuals=[0.0, 0.0, 0.0, 0.0],
            transformed_tree_newick="(A:1,B:1,C:1,D:1);",
            confidence_intervals=[],
            residual_diagnostics=ComparativeResidualSummary(
                residual_mean=0.0,
                residual_variance=0.0,
                residual_skewness=0.0,
                max_abs_standardized_residual=0.0,
                phylogenetic_residual_lambda=0.0,
                outlier_taxa=[],
                warnings=[],
            ),
            optimizer_diagnostics=None,
            optimizer_profile_rows=None,
            identifiability_warnings=[],
            assumptions=[],
        )

    monkeypatch.setattr(
        evolutionary_modes_module,
        "load_comparative_dataset",
        fake_load_dataset,
    )
    monkeypatch.setattr(
        evolutionary_modes_module,
        "_fit_evolutionary_mode_from_dataset",
        fake_fit,
    )

    with pytest.raises(
        ComparativeMethodError,
        match="mixed likelihood constant policies prevent ranking incompatible continuous-mode models",
    ):
        compare_fitcontinuous_model_ranking(
            EXAMPLE_TREE,
            EXAMPLE_TRAITS,
            trait="response",
            taxon_column="taxon",
            modes=("brownian", "white-noise"),
        )


def test_compare_fitcontinuous_model_ranking_records_likelihood_policy() -> None:
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_brownian_signal_twenty_four_taxa"
    )

    report = compare_fitcontinuous_model_ranking(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        taxon_column=fixture.taxon_column,
        modes=("brownian", "white-noise"),
    )

    assert report.likelihood_constant_policy == (
        "full-gaussian-loglikelihood-includes-normalizing-constant"
    )
    assert report.likelihood_comparison_policy == (
        "relative-aic-and-aicc-ranking-is-permitted-only-when-all-candidate-modes-share-one-gaussian-likelihood-constant-policy"
    )
    assert report.noncomparable_likelihood_models == []


def test_compare_fitcontinuous_model_ranking_withholds_stable_conclusion_for_boundary_dominated_selection() -> (
    None
):
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_white_noise_twenty_four_taxa"
    )

    report = compare_fitcontinuous_model_ranking(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        taxon_column=fixture.taxon_column,
        modes=("pagel-lambda",),
    )

    assert report.better_model == "pagel-lambda"
    assert report.selected_model_boundary_assessment is not None
    assert report.selected_model_boundary_assessment.affected_parameter == "lambda"
    assert report.selected_model_boundary_assessment.hit_lower_boundary is True
    assert (
        report.selected_model_boundary_assessment.boundary_dominates_interpretation
        is True
    )
    assert report.stable_conclusion_supported is False
    assert any(
        "stable conclusion support is withheld" in warning
        for warning in report.warnings
    )


@pytest.mark.parametrize("mode", ["trend", "mean_trend", "rate_trend"])
def test_fit_continuous_evolutionary_mode_explicitly_excludes_geiger_trend_aliases(
    mode: str,
) -> None:
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_trend_proxy_twenty_four_taxa"
    )

    with pytest.raises(
        ComparativeMethodError,
        match="trend-mode parity is explicitly excluded in this round",
    ):
        fit_continuous_evolutionary_mode(
            fixture.tree_path,
            fixture.traits_path,
            trait=fixture.trait_name,
            mode=mode,
            taxon_column=fixture.taxon_column,
        )


def test_fit_continuous_evolutionary_mode_supports_early_burst() -> None:
    fit = fit_continuous_evolutionary_mode(
        EXAMPLE_TREE,
        EXAMPLE_TRAITS,
        trait="response",
        mode="early-burst",
    )

    assert fit.mode == "early-burst"
    assert fit.parameter_name == "rate_change"
    assert fit.parameter_value is not None
    assert fit.transformed_tree_newick.endswith(";")
    assert fit.log_likelihood < 0.0
    assert fit.aic > 0.0
    assert fit.optimizer_diagnostics is not None
    assert fit.optimizer_diagnostics.optimizer_name == "governed-two-stage-grid-search"
    assert [warning.kind for warning in fit.identifiability_warnings] == [
        "boundary_rate_change",
        "flat_likelihood_profile",
        "brownian_like_rate_change",
    ]
    assert fit.boundary_assessment is not None
    assert fit.boundary_assessment.affected_parameter == "rate_change"
    assert fit.boundary_assessment.hit_lower_boundary is True
    assert fit.boundary_assessment.near_lower_boundary is True
    assert fit.boundary_assessment.flat_likelihood_near_boundary is True
    assert fit.boundary_assessment.stable_conclusion_supported is False


def test_fit_continuous_evolutionary_mode_supports_pagel_kappa_strong_signal() -> None:
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_brownian_signal_twenty_four_taxa"
    )

    fit = fit_continuous_evolutionary_mode(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        mode="pagel-kappa",
        taxon_column=fixture.taxon_column,
    )

    assert fit.mode == "pagel-kappa"
    assert fit.parameter_name == "kappa"
    assert fit.parameter_value is not None
    assert 1.0 < fit.parameter_value < 1.5
    assert fit.optimizer_diagnostics is not None
    assert fit.optimizer_diagnostics.hit_lower_boundary is False
    assert fit.optimizer_diagnostics.hit_upper_boundary is False
    assert fit.aicc >= fit.aic
    assert "branch-length power transformation" in fit.assumptions[0]
    assert [warning.kind for warning in fit.identifiability_warnings] == [
        "flat_likelihood"
    ]


def test_fit_continuous_evolutionary_mode_supports_pagel_kappa_weak_signal() -> None:
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_white_noise_twenty_four_taxa"
    )

    fit = fit_continuous_evolutionary_mode(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        mode="pagel-kappa",
        taxon_column=fixture.taxon_column,
    )

    assert fit.mode == "pagel-kappa"
    assert fit.parameter_name == "kappa"
    assert fit.parameter_value is not None
    assert math.isclose(fit.parameter_value, 0.0, abs_tol=1e-12)
    assert fit.optimizer_diagnostics is not None
    assert fit.optimizer_diagnostics.hit_lower_boundary is True
    assert fit.optimizer_diagnostics.hit_upper_boundary is False
    assert [warning.kind for warning in fit.identifiability_warnings] == [
        "boundary_kappa",
        "flat_likelihood",
        "punctuational_limit",
    ]
    assert fit.boundary_assessment is not None
    assert fit.boundary_assessment.affected_parameter == "kappa"
    assert fit.boundary_assessment.hit_lower_boundary is True
    assert fit.boundary_assessment.near_lower_boundary is True
    assert fit.boundary_assessment.flat_likelihood_near_boundary is True
    assert fit.boundary_assessment.stable_conclusion_supported is False


def test_fit_continuous_evolutionary_mode_recovers_shared_early_burst_fixture() -> None:
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_early_burst_known_truth_twenty_four_taxa"
    )

    fit = fit_continuous_evolutionary_mode(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        mode="early-burst",
        taxon_column=fixture.taxon_column,
        early_burst_bounds=(0.0, 10.0),
    )

    assert fit.parameter_name == "rate_change"
    assert fit.parameter_value is not None
    assert 4.0 <= fit.parameter_value <= 5.0
    assert fit.optimizer_diagnostics is not None
    assert fit.optimizer_diagnostics.hit_lower_boundary is False
    assert fit.optimizer_diagnostics.hit_upper_boundary is False
    assert [warning.kind for warning in fit.identifiability_warnings] == [
        "flat_likelihood_profile"
    ]
    assert fit.boundary_assessment is not None
    assert fit.boundary_assessment.affected_parameter == "rate_change"
    assert fit.boundary_assessment.boundary_dominates_interpretation is False
    assert fit.boundary_assessment.stable_conclusion_supported is True


def test_fit_continuous_evolutionary_mode_reports_ou_identifiability_and_bounds() -> (
    None
):
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_nonultrametric_control_twenty_four_taxa"
    )

    fit = fit_continuous_evolutionary_mode(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        mode="ornstein-uhlenbeck",
        taxon_column=fixture.taxon_column,
    )

    assert fit.mode == "ornstein-uhlenbeck"
    assert fit.parameter_name == "alpha"
    assert fit.parameter_value is not None
    assert fit.optimizer_diagnostics is not None
    assert fit.optimizer_diagnostics.hit_lower_boundary is True
    assert fit.optimizer_diagnostics.hit_upper_boundary is False
    assert [warning.kind for warning in fit.identifiability_warnings] == [
        "boundary_alpha",
        "flat_likelihood",
        "weak_pull_to_optimum",
    ]
    assert fit.boundary_assessment is not None
    assert fit.boundary_assessment.affected_parameter == "alpha"
    assert fit.boundary_assessment.hit_lower_boundary is True
    assert fit.boundary_assessment.near_lower_boundary is True
    assert fit.boundary_assessment.flat_likelihood_near_boundary is True
    assert fit.boundary_assessment.stable_conclusion_supported is False


def test_fit_continuous_evolutionary_mode_reports_ou_aicc_and_missing_value_context() -> (
    None
):
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_missing_values_twenty_four_taxa"
    )

    fit = fit_continuous_evolutionary_mode(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        mode="ornstein-uhlenbeck",
        taxon_column=fixture.taxon_column,
    )

    assert fit.aicc >= fit.aic
    assert fit.optimizer_diagnostics is not None
    assert fit.optimizer_diagnostics.lower_bound > 0.0
    assert [warning.kind for warning in fit.identifiability_warnings] == [
        "flat_likelihood",
        "weak_pull_to_optimum",
    ]


def test_fit_continuous_evolutionary_mode_supports_pagel_lambda_strong_signal() -> None:
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_brownian_signal_twenty_four_taxa"
    )

    fit = fit_continuous_evolutionary_mode(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        mode="pagel-lambda",
        taxon_column=fixture.taxon_column,
    )

    assert fit.mode == "pagel-lambda"
    assert fit.parameter_name == "lambda"
    assert fit.parameter_value is not None
    assert math.isclose(fit.parameter_value, 1.0, abs_tol=1e-12)
    assert fit.optimizer_diagnostics is not None
    assert fit.optimizer_diagnostics.hit_lower_boundary is False
    assert fit.optimizer_diagnostics.hit_upper_boundary is True
    assert fit.aicc >= fit.aic
    assert "phytools::phylosig" in fit.assumptions[1]
    assert [warning.kind for warning in fit.identifiability_warnings] == [
        "boundary_lambda",
        "flat_likelihood",
        "brownian_limit",
    ]
    assert fit.boundary_assessment is not None
    assert fit.boundary_assessment.affected_parameter == "lambda"
    assert fit.boundary_assessment.hit_upper_boundary is True
    assert fit.boundary_assessment.near_upper_boundary is True
    assert fit.boundary_assessment.boundary_dominates_interpretation is True


def test_fit_continuous_evolutionary_mode_supports_pagel_lambda_weak_signal() -> None:
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_white_noise_twenty_four_taxa"
    )

    fit = fit_continuous_evolutionary_mode(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        mode="pagel-lambda",
        taxon_column=fixture.taxon_column,
    )

    assert fit.mode == "pagel-lambda"
    assert fit.parameter_name == "lambda"
    assert fit.parameter_value is not None
    assert math.isclose(fit.parameter_value, 0.0, abs_tol=1e-12)
    assert fit.optimizer_diagnostics is not None
    assert fit.optimizer_diagnostics.hit_lower_boundary is True
    assert fit.optimizer_diagnostics.hit_upper_boundary is False
    assert [warning.kind for warning in fit.identifiability_warnings] == [
        "boundary_lambda",
        "flat_likelihood",
        "weak_phylogenetic_signal",
    ]
    assert fit.boundary_assessment is not None
    assert fit.boundary_assessment.affected_parameter == "lambda"
    assert fit.boundary_assessment.hit_lower_boundary is True
    assert fit.boundary_assessment.near_lower_boundary is True
    assert fit.boundary_assessment.flat_likelihood_near_boundary is True
    assert fit.boundary_assessment.stable_conclusion_supported is False


def test_fit_continuous_evolutionary_mode_reports_pagel_lambda_missing_value_context() -> (
    None
):
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_missing_values_twenty_four_taxa"
    )

    fit = fit_continuous_evolutionary_mode(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        mode="pagel-lambda",
        taxon_column=fixture.taxon_column,
    )

    assert fit.parameter_name == "lambda"
    assert fit.parameter_value is not None
    assert math.isclose(fit.parameter_value, 1.0, abs_tol=1e-12)
    assert fit.aicc >= fit.aic
    assert fit.optimizer_diagnostics is not None
    assert fit.optimizer_diagnostics.hit_upper_boundary is True
    assert [warning.kind for warning in fit.identifiability_warnings] == [
        "boundary_lambda",
        "flat_likelihood",
        "brownian_limit",
    ]


def test_compare_continuous_evolutionary_modes_reports_likelihood_ratios() -> None:
    report = compare_continuous_evolutionary_modes(
        EXAMPLE_TREE,
        EXAMPLE_TRAITS,
        trait="response",
    )

    assert report.better_model == "brownian"
    assert [row.model for row in report.rows] == [
        "brownian",
        "early-burst",
        "ornstein-uhlenbeck",
    ]
    assert sum(1 for row in report.rows if row.selected) == 1
    assert [row.comparison_id for row in report.likelihood_ratio_tests] == [
        "brownian-vs-ornstein-uhlenbeck",
        "brownian-vs-early-burst",
        "ornstein-uhlenbeck-vs-early-burst",
    ]
    assert all(row.degrees_of_freedom == 1 for row in report.likelihood_ratio_tests)
    assert all(row.statistic >= 0.0 for row in report.likelihood_ratio_tests)


def test_compare_fitcontinuous_model_ranking_reports_ranked_all_model_surface() -> None:
    report = compare_fitcontinuous_model_ranking(
        EXAMPLE_TREE,
        EXAMPLE_TRAITS,
        trait="response",
    )

    assert report.better_model == "white-noise"
    assert [row.model for row in report.rows] == [
        "white-noise",
        "brownian",
        "early-burst",
        "ornstein-uhlenbeck",
        "pagel-delta",
        "pagel-kappa",
        "pagel-lambda",
    ]
    assert report.rows[0].rank == 1
    assert report.rows[0].delta_aicc == 0.0
    assert report.rows[1].rank == 2
    assert report.rows[2].comparable is False
    assert (
        report.rows[2].comparability_note
        == "sample size is too small to compute finite AICc for this parameter count"
    )
    assert any(
        "likelihood-ratio tests remain reported only for the nested brownian/ornstein-uhlenbeck/early-burst subset"
        in warning
        for warning in report.warnings
    )


def test_compare_fitcontinuous_model_ranking_records_model_confidence_surface() -> None:
    report = compare_fitcontinuous_model_ranking(
        EXAMPLE_TREE,
        EXAMPLE_TRAITS,
        trait="response",
    )

    comparable_rows = [row for row in report.rows if row.comparable]
    noncomparable_rows = [row for row in report.rows if not row.comparable]

    assert report.model_confidence_weight_basis == "AICc"
    assert report.model_confidence_delta_threshold == 2.0
    assert comparable_rows
    assert math.isclose(
        sum(row.akaike_weight or 0.0 for row in comparable_rows),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert all(row.akaike_weight is not None for row in comparable_rows)
    assert all(row.within_delta_aic_threshold is not None for row in comparable_rows)
    assert all(row.within_delta_aicc_threshold is not None for row in comparable_rows)
    assert all(row.akaike_weight is None for row in noncomparable_rows)
    assert all(row.within_delta_aic_threshold is None for row in noncomparable_rows)
    assert all(row.within_delta_aicc_threshold is None for row in noncomparable_rows)
    assert report.selected_model_akaike_weight == report.rows[0].akaike_weight
    assert report.models_within_delta_aicc_threshold == [
        row.model for row in comparable_rows if row.within_delta_aicc_threshold
    ]
    assert "model confidence" in report.uncertainty_language


def test_reconstruct_continuous_evolutionary_mode_states_supports_early_burst() -> None:
    report = reconstruct_continuous_evolutionary_mode_states(
        EXAMPLE_TREE,
        EXAMPLE_TRAITS,
        trait="response",
        mode="early-burst",
        rate_change=0.5,
    )

    assert report.mode == "early-burst"
    assert report.parameter_name == "rate_change"
    assert math.isclose(report.parameter_value or 0.0, 0.5)
    assert report.transformed_tree_newick.endswith(";")
    assert report.reconstruction.analysis_tree_newick == report.transformed_tree_newick
    assert report.reconstruction.model == "brownian"
    assert len(report.reconstruction.estimates) >= 4
