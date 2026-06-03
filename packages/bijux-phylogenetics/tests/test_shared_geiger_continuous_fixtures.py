from __future__ import annotations

import csv
from statistics import median

import pytest

from bijux_phylogenetics.comparative.common import summarize_numeric_trait_readiness
from bijux_phylogenetics.comparative.evolutionary_modes import (
    compare_fitcontinuous_model_ranking,
)
from bijux_phylogenetics.comparative.signal import (
    compute_blombergs_k,
    estimate_pagels_lambda,
)
import bijux_phylogenetics.datasets.shared_fixtures as fixtures_api
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_geiger_continuous_fixture,
    list_shared_geiger_continuous_fixtures,
)
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    assess_tree_ultrametricity,
)
from bijux_phylogenetics.runtime.errors import ComparativeMethodError


def _read_numeric_trait_values(fixture_id: str) -> list[float]:
    fixture = get_shared_geiger_continuous_fixture(fixture_id)
    with fixture.traits_path.open(encoding="utf-8", newline="") as handle:
        return [
            float(row[fixture.trait_name])
            for row in csv.DictReader(handle, delimiter="\t")
            if row[fixture.trait_name]
        ]


def test_shared_geiger_continuous_fixture_catalog_covers_required_goal_cases() -> None:
    fixtures = list_shared_geiger_continuous_fixtures()
    feature_tags = {tag for fixture in fixtures for tag in fixture.feature_tags}
    supported_models = {
        model for fixture in fixtures for model in fixture.supported_model_names
    }

    assert {
        "twenty-plus-taxa-ultrametric-tree",
        "one-hundred-plus-taxa-ultrametric-tree",
        "non-ultrametric-control",
        "brownian-known-truth",
        "ou-known-truth",
        "early-burst-known-truth",
        "low-signal",
        "high-signal",
        "missing-trait-values",
        "constant-trait-negative-case",
        "extreme-outlier",
        "trait-standard-error",
        "known-truth-simulation",
    } <= feature_tags
    assert {"BM", "OU", "EB", "lambda", "kappa", "delta", "white", "trend"} <= (
        supported_models
    )


def test_shared_geiger_continuous_fixture_lookup_resolves_linked_shared_surfaces() -> (
    None
):
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_ou_known_truth_twenty_four_taxa"
    )

    assert fixture.tree_fixture.tip_count == 24
    assert fixture.trait_table_fixture.row_count == 24
    assert fixture.trait_name == "ou_truth"
    assert fixture.supported_model_names == ("OU",)
    assert fixture.simulation_metadata == {
        "model": "ornstein-uhlenbeck",
        "root_state": 0.0,
        "sigma": 0.2,
        "alpha": 0.5,
        "theta": 1.0,
        "seed": 101,
    }


def test_shared_geiger_continuous_fixture_catalog_includes_large_and_nonultrametric_cases() -> (
    None
):
    medium_fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_brownian_signal_twenty_four_taxa"
    )
    large_fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_brownian_signal_one_hundred_twenty_eight_taxa"
    )
    nonultrametric_fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_nonultrametric_control_twenty_four_taxa"
    )

    medium_ultrametric = assess_tree_ultrametricity(medium_fixture.tree_path)
    large_ultrametric = assess_tree_ultrametricity(large_fixture.tree_path)
    nonultrametric = assess_tree_ultrametricity(nonultrametric_fixture.tree_path)

    assert medium_fixture.tree_fixture.tip_count >= 20
    assert medium_ultrametric.ultrametric is True
    assert large_fixture.tree_fixture.tip_count >= 100
    assert large_ultrametric.ultrametric is True
    assert nonultrametric_fixture.tree_fixture.tip_count >= 20
    assert nonultrametric.ultrametric is False


@pytest.mark.slow
def test_shared_geiger_continuous_fixture_catalog_supports_ou_and_early_burst_truth_cases() -> (
    None
):
    brownian_fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_brownian_signal_twenty_four_taxa"
    )
    ou_fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_ou_known_truth_twenty_four_taxa"
    )
    early_burst_fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_early_burst_known_truth_twenty_four_taxa"
    )
    white_fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_white_noise_twenty_four_taxa"
    )

    brownian_comparison = compare_fitcontinuous_model_ranking(
        brownian_fixture.tree_path,
        brownian_fixture.traits_path,
        trait=brownian_fixture.trait_name,
        taxon_column=brownian_fixture.taxon_column,
    )
    ou_comparison = compare_fitcontinuous_model_ranking(
        ou_fixture.tree_path,
        ou_fixture.traits_path,
        trait=ou_fixture.trait_name,
        taxon_column=ou_fixture.taxon_column,
    )
    early_burst_comparison = compare_fitcontinuous_model_ranking(
        early_burst_fixture.tree_path,
        early_burst_fixture.traits_path,
        trait=early_burst_fixture.trait_name,
        taxon_column=early_burst_fixture.taxon_column,
    )
    white_comparison = compare_fitcontinuous_model_ranking(
        white_fixture.tree_path,
        white_fixture.traits_path,
        trait=white_fixture.trait_name,
        taxon_column=white_fixture.taxon_column,
    )

    assert brownian_comparison.better_model == "brownian"
    assert ou_comparison.better_model == "ornstein-uhlenbeck"
    assert early_burst_comparison.better_model == "early-burst"
    assert white_comparison.better_model == "white-noise"
    assert len(brownian_comparison.rows) == 7
    assert brownian_comparison.rows[0].rank == 1
    assert brownian_comparison.rows[0].delta_aicc == 0.0


@pytest.mark.slow
def test_shared_geiger_continuous_fixture_catalog_supports_signal_strength_boundaries() -> (
    None
):
    strong_small = get_shared_geiger_continuous_fixture(
        "geiger_continuous_brownian_signal_twenty_four_taxa"
    )
    weak_small = get_shared_geiger_continuous_fixture(
        "geiger_continuous_white_noise_twenty_four_taxa"
    )
    strong_large = get_shared_geiger_continuous_fixture(
        "geiger_continuous_brownian_signal_one_hundred_twenty_eight_taxa"
    )
    weak_large = get_shared_geiger_continuous_fixture(
        "geiger_continuous_white_noise_one_hundred_twenty_eight_taxa"
    )

    strong_small_k = compute_blombergs_k(
        strong_small.tree_path,
        strong_small.traits_path,
        trait=strong_small.trait_name,
    )
    weak_small_k = compute_blombergs_k(
        weak_small.tree_path,
        weak_small.traits_path,
        trait=weak_small.trait_name,
    )
    strong_large_lambda = estimate_pagels_lambda(
        strong_large.tree_path,
        strong_large.traits_path,
        trait=strong_large.trait_name,
    )
    weak_large_lambda = estimate_pagels_lambda(
        weak_large.tree_path,
        weak_large.traits_path,
        trait=weak_large.trait_name,
    )

    assert strong_small_k.k > weak_small_k.k
    assert strong_large_lambda.lambda_value > weak_large_lambda.lambda_value


def test_shared_geiger_continuous_fixture_catalog_handles_missing_constant_outlier_and_trend_cases() -> (
    None
):
    missing_fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_missing_values_twenty_four_taxa"
    )
    constant_fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_constant_negative_twenty_four_taxa"
    )

    missing_readiness = summarize_numeric_trait_readiness(
        missing_fixture.tree_path,
        missing_fixture.traits_path,
        trait=missing_fixture.trait_name,
        taxon_column=missing_fixture.taxon_column,
    )
    outlier_values = _read_numeric_trait_values(
        "geiger_continuous_outlier_signal_twenty_four_taxa"
    )
    trend_fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_trend_proxy_twenty_four_taxa"
    )
    trend_values = _read_numeric_trait_values(trend_fixture.fixture_id)

    assert len(missing_readiness.analysis_taxa) == 22
    assert missing_readiness.pruned_missing_value_taxa == ["Phy10"]
    assert missing_readiness.pruned_non_numeric_taxa == ["Phy14"]
    with pytest.raises(ComparativeMethodError):
        compute_blombergs_k(
            constant_fixture.tree_path,
            constant_fixture.traits_path,
            trait=constant_fixture.trait_name,
        )
    assert trend_fixture.geiger_reference_expectation == (
        "fitcontinuous-trend-explicitly-excluded-this-round"
    )
    assert "explicitly excludes geiger trend parity" in trend_fixture.notes
    assert max(outlier_values) - median(outlier_values) > 4.0
    assert trend_values[-1] > trend_values[0]


def test_shared_geiger_continuous_fixture_catalog_tracks_standard_error_review_surface() -> (
    None
):
    fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_standard_error_review_twenty_four_taxa"
    )

    assert fixture.trait_name == "ou_truth"
    assert fixture.standard_error_trait_name == "ou_truth_standard_error"
    assert fixture.geiger_reference_expectation == (
        "fitcontinuous-standard-error-explicitly-excluded-this-round"
    )
    assert "trait-standard-error" in fixture.feature_tags
    assert "explicitly excludes measurement-error variance handling" in fixture.notes


def test_public_runtime_exports_include_shared_geiger_continuous_fixture_surface() -> (
    None
):
    assert (
        fixtures_api.get_shared_geiger_continuous_fixture
        is get_shared_geiger_continuous_fixture
    )
    assert (
        fixtures_api.list_shared_geiger_continuous_fixtures
        is list_shared_geiger_continuous_fixtures
    )
