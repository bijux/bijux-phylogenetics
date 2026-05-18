from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.ancestral import (
    reconstruct_continuous_evolutionary_mode_states,
)
from bijux_phylogenetics.comparative import (
    compare_continuous_evolutionary_modes,
    fit_continuous_evolutionary_mode,
    rescale_tree_early_burst,
    rescale_tree_ornstein_uhlenbeck,
    rescale_tree_pagel_delta,
    rescale_tree_pagel_kappa,
)
from bijux_phylogenetics.fixtures import get_shared_geiger_continuous_fixture
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures"
EXAMPLE_TREE = FIXTURE_ROOT / "trees" / "example_tree.nwk"
EXAMPLE_TREE_NEGATIVE = FIXTURE_ROOT / "trees" / "example_tree_negative_length.nwk"
EXAMPLE_TRAITS = FIXTURE_ROOT / "metadata" / "example_traits_comparative.tsv"


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


def test_rescale_tree_ornstein_uhlenbeck_rejects_negative_alpha() -> None:
    with pytest.raises(ComparativeMethodError, match="OU alpha must be non-negative"):
        rescale_tree_ornstein_uhlenbeck(EXAMPLE_TREE, alpha=-0.5)


def test_rescale_tree_pagel_kappa_rejects_negative_branch_lengths() -> None:
    with pytest.raises(
        ComparativeMethodError,
        match="Pagel kappa cannot transform negative branch lengths",
    ):
        rescale_tree_pagel_kappa(EXAMPLE_TREE_NEGATIVE, kappa=0.5)


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


def test_fit_continuous_evolutionary_mode_white_noise_fits_high_signal_worse_than_brownian() -> None:
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


def test_fit_continuous_evolutionary_mode_supports_pagel_lambda_strong_signal() -> (
    None
):
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
        "ornstein-uhlenbeck",
        "early-burst",
    ]
    assert sum(1 for row in report.rows if row.selected) == 1
    assert [row.comparison_id for row in report.likelihood_ratio_tests] == [
        "brownian-vs-ornstein-uhlenbeck",
        "brownian-vs-early-burst",
        "ornstein-uhlenbeck-vs-early-burst",
    ]
    assert all(row.degrees_of_freedom == 1 for row in report.likelihood_ratio_tests)
    assert all(row.statistic >= 0.0 for row in report.likelihood_ratio_tests)


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
