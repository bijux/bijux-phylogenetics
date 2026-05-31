from __future__ import annotations

import csv

import pytest

from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.comparative.common import summarize_numeric_trait_readiness
from bijux_phylogenetics.comparative.discrete_mk import fit_discrete_mk_model
from bijux_phylogenetics.comparative.pgls import (
    build_pgls_model_matrix,
    run_pgls,
)
from bijux_phylogenetics.comparative.regression import (
    summarize_phylogenetic_anova,
    summarize_phylogenetic_residuals,
)
from bijux_phylogenetics.comparative.signal import (
    compute_blombergs_k,
    compute_phylogenetic_signal_test,
    estimate_pagels_lambda,
)
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_phytools_comparative_fixture,
    list_shared_phytools_comparative_fixtures,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError


def _read_trait_values(fixture_id: str, trait_name: str) -> list[str]:
    fixture = get_shared_phytools_comparative_fixture(fixture_id)
    with fixture.traits_path.open(encoding="utf-8", newline="") as handle:
        return [row[trait_name] for row in csv.DictReader(handle, delimiter="\t")]


def test_shared_phytools_comparative_fixture_catalog_covers_required_goal_cases() -> (
    None
):
    fixtures = list_shared_phytools_comparative_fixtures()
    feature_tags = {tag for fixture in fixtures for tag in fixture.feature_tags}

    assert {
        "twenty-plus-taxa-ultrametric-tree",
        "one-hundred-plus-taxa-ultrametric-tree",
        "non-ultrametric-tree",
        "continuous-trait",
        "strong-phylogenetic-signal",
        "weak-phylogenetic-signal",
        "binary-discrete-trait",
        "multistate-discrete-trait",
        "missing-trait-values",
        "constant-trait-negative-case",
        "tree-trait-taxon-mismatch",
        "branch-length-edge-case",
        "known-truth-simulation",
    } <= feature_tags


def test_shared_phytools_comparative_fixture_lookup_resolves_linked_shared_surfaces() -> (
    None
):
    fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_strong_signal_twenty_four_taxa"
    )

    assert fixture.tree_fixture.tip_count == 24
    assert fixture.trait_table_fixture.row_count == 24
    assert fixture.trait_name == "signal_strong"
    assert fixture.trait_kind == "continuous"
    assert fixture.simulation_metadata == {
        "model": "brownian-motion",
        "root_state": 0.5,
        "sigma": 0.8,
        "seed": 32024,
    }


def test_shared_phytools_comparative_fixture_catalog_includes_large_ultrametric_cases() -> (
    None
):
    small_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_strong_signal_twenty_four_taxa"
    )
    large_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_strong_signal_one_hundred_twenty_eight_taxa"
    )

    assert small_fixture.tree_fixture.tip_count >= 20
    assert "ultrametric" in small_fixture.tree_fixture.feature_tags
    assert large_fixture.tree_fixture.tip_count >= 100
    assert "ultrametric" in large_fixture.tree_fixture.feature_tags


@pytest.mark.slow
def test_shared_phytools_comparative_fixture_catalog_strong_and_weak_signal_cases_behave_as_labeled() -> (
    None
):
    strong_small = get_shared_phytools_comparative_fixture(
        "phytools_continuous_strong_signal_twenty_four_taxa"
    )
    weak_small = get_shared_phytools_comparative_fixture(
        "phytools_continuous_weak_signal_twenty_four_taxa"
    )
    strong_large = get_shared_phytools_comparative_fixture(
        "phytools_continuous_strong_signal_one_hundred_twenty_eight_taxa"
    )
    weak_large = get_shared_phytools_comparative_fixture(
        "phytools_continuous_weak_signal_one_hundred_twenty_eight_taxa"
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
    strong_small_signal_test = compute_phylogenetic_signal_test(
        strong_small.tree_path,
        strong_small.traits_path,
        trait=strong_small.trait_name,
        permutations=199,
        seed=17,
    )
    weak_small_signal_test = compute_phylogenetic_signal_test(
        weak_small.tree_path,
        weak_small.traits_path,
        trait=weak_small.trait_name,
        permutations=199,
        seed=17,
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

    assert strong_small_k.k > 0.5
    assert weak_small_k.k < 0.1
    assert strong_small_signal_test.p_value < weak_small_signal_test.p_value
    assert strong_large_lambda.lambda_value >= 0.95
    assert weak_large_lambda.lambda_value <= 0.1


def test_shared_phytools_comparative_fixture_catalog_handles_missing_constant_and_mismatch_cases() -> (
    None
):
    missing_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_missing_values_twenty_four_taxa"
    )
    mismatch_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_mismatch_twenty_four_taxa"
    )
    constant_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_constant_negative_twenty_four_taxa"
    )

    missing_readiness = summarize_numeric_trait_readiness(
        missing_fixture.tree_path,
        missing_fixture.traits_path,
        trait=missing_fixture.trait_name,
    )
    missing_signal = compute_blombergs_k(
        missing_fixture.tree_path,
        missing_fixture.traits_path,
        trait=missing_fixture.trait_name,
    )
    mismatch_readiness = summarize_numeric_trait_readiness(
        mismatch_fixture.tree_path,
        mismatch_fixture.traits_path,
        trait=mismatch_fixture.trait_name,
    )

    assert missing_signal.taxon_count == 22
    assert missing_signal.input_audit.pruned_missing_value_taxa == ["Phy10"]
    assert missing_readiness.pruned_non_numeric_taxa == ["Phy14"]
    assert mismatch_readiness.missing_from_traits == ["Phy9"]
    assert mismatch_readiness.extra_trait_taxa == ["PhyExtra"]
    with pytest.raises(ComparativeMethodError):
        compute_blombergs_k(
            constant_fixture.tree_path,
            constant_fixture.traits_path,
            trait=constant_fixture.trait_name,
        )


def test_shared_phytools_comparative_fixture_catalog_supports_fast_anc_signal_cases() -> (
    None
):
    strong_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_strong_signal_twenty_four_taxa"
    )
    weak_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_weak_signal_twenty_four_taxa"
    )
    missing_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_missing_values_twenty_four_taxa"
    )

    strong_report = reconstruct_continuous_ancestral_states(
        strong_fixture.tree_path,
        strong_fixture.traits_path,
        trait=strong_fixture.trait_name,
        taxon_column=strong_fixture.taxon_column,
        model="brownian",
        estimator="fast-anc",
    )
    weak_report = reconstruct_continuous_ancestral_states(
        weak_fixture.tree_path,
        weak_fixture.traits_path,
        trait=weak_fixture.trait_name,
        taxon_column=weak_fixture.taxon_column,
        model="brownian",
        estimator="fast-anc",
    )
    missing_report = reconstruct_continuous_ancestral_states(
        missing_fixture.tree_path,
        missing_fixture.traits_path,
        trait=missing_fixture.trait_name,
        taxon_column=missing_fixture.taxon_column,
        model="brownian",
        estimator="fast-anc",
    )

    assert strong_report.estimator == "fast-anc"
    assert weak_report.estimator == "fast-anc"
    assert len([row for row in strong_report.estimates if not row.is_tip]) == 23
    assert len([row for row in weak_report.estimates if not row.is_tip]) == 23
    assert missing_report.taxon_count == 22
    assert sorted(
        missing_report.dropped_missing_taxa + missing_report.dropped_non_numeric_taxa
    ) == ["Phy10", "Phy14"]


def test_shared_phytools_comparative_fixture_catalog_supports_anc_ml_signal_cases() -> (
    None
):
    strong_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_strong_signal_twenty_four_taxa"
    )
    nonultrametric_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_strong_signal_non_ultrametric_twenty_four_taxa"
    )
    missing_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_missing_values_twenty_four_taxa"
    )

    strong_report = reconstruct_continuous_ancestral_states(
        strong_fixture.tree_path,
        strong_fixture.traits_path,
        trait=strong_fixture.trait_name,
        taxon_column=strong_fixture.taxon_column,
        model="brownian",
        estimator="anc-ml",
    )
    nonultrametric_report = reconstruct_continuous_ancestral_states(
        nonultrametric_fixture.tree_path,
        nonultrametric_fixture.traits_path,
        trait=nonultrametric_fixture.trait_name,
        taxon_column=nonultrametric_fixture.taxon_column,
        model="brownian",
        estimator="anc-ml",
    )
    missing_report = reconstruct_continuous_ancestral_states(
        missing_fixture.tree_path,
        missing_fixture.traits_path,
        trait=missing_fixture.trait_name,
        taxon_column=missing_fixture.taxon_column,
        model="brownian",
        estimator="anc-ml",
    )

    assert strong_report.estimator == "anc-ml"
    assert nonultrametric_report.estimator == "anc-ml"
    assert strong_report.optimizer_diagnostics is not None
    assert strong_report.optimizer_diagnostics.converged is True
    assert len([row for row in strong_report.estimates if not row.is_tip]) == 23
    assert len([row for row in nonultrametric_report.estimates if not row.is_tip]) == 23
    assert missing_report.taxon_count == 22
    assert sorted(
        missing_report.dropped_missing_taxa + missing_report.dropped_non_numeric_taxa
    ) == ["Phy10", "Phy14"]


def test_shared_phytools_comparative_fixture_catalog_supports_fitmk_er_cases() -> None:
    binary_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_binary_twenty_four_taxa"
    )
    multistate_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_multistate_twenty_four_taxa"
    )
    binary_missing_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_binary_missing_twenty_four_taxa"
    )
    multistate_missing_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_multistate_missing_twenty_four_taxa"
    )

    binary_report = fit_discrete_mk_model(
        binary_fixture.tree_path,
        binary_fixture.traits_path,
        trait=binary_fixture.trait_name,
        taxon_column=binary_fixture.taxon_column,
        model="equal-rates",
    )
    multistate_report = fit_discrete_mk_model(
        multistate_fixture.tree_path,
        multistate_fixture.traits_path,
        trait=multistate_fixture.trait_name,
        taxon_column=multistate_fixture.taxon_column,
        model="equal-rates",
    )
    binary_missing_report = fit_discrete_mk_model(
        binary_missing_fixture.tree_path,
        binary_missing_fixture.traits_path,
        trait=binary_missing_fixture.trait_name,
        taxon_column=binary_missing_fixture.taxon_column,
        model="equal-rates",
    )
    multistate_missing_report = fit_discrete_mk_model(
        multistate_missing_fixture.tree_path,
        multistate_missing_fixture.traits_path,
        trait=multistate_missing_fixture.trait_name,
        taxon_column=multistate_missing_fixture.taxon_column,
        model="equal-rates",
    )

    assert binary_report.parameter_count == 1
    assert multistate_report.parameter_count == 1
    assert binary_missing_report.input_audit.pruned_missing_value_taxa == ["Phy10"]
    assert multistate_missing_report.input_audit.pruned_missing_value_taxa == ["Phy14"]
    assert len(binary_report.transition_rate_rows) == 2
    assert len(multistate_report.transition_rate_rows) == 6


@pytest.mark.slow
def test_shared_phytools_comparative_fixture_catalog_supports_fitmk_sym_cases() -> None:
    multistate_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_multistate_twenty_four_taxa"
    )
    multistate_missing_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_multistate_missing_twenty_four_taxa"
    )

    multistate_report = fit_discrete_mk_model(
        multistate_fixture.tree_path,
        multistate_fixture.traits_path,
        trait=multistate_fixture.trait_name,
        taxon_column=multistate_fixture.taxon_column,
        model="symmetric",
    )
    multistate_missing_report = fit_discrete_mk_model(
        multistate_missing_fixture.tree_path,
        multistate_missing_fixture.traits_path,
        trait=multistate_missing_fixture.trait_name,
        taxon_column=multistate_missing_fixture.taxon_column,
        model="symmetric",
    )

    assert "symmetric-discrete-model" in multistate_fixture.feature_tags
    assert "symmetric-discrete-model" in multistate_missing_fixture.feature_tags
    assert multistate_report.model == "symmetric"
    assert multistate_report.parameter_count == 3
    assert multistate_report.baseline_comparison is not None
    assert multistate_missing_report.input_audit.pruned_missing_value_taxa == ["Phy14"]


@pytest.mark.slow
def test_shared_phytools_comparative_fixture_catalog_supports_fitmk_ard_cases() -> None:
    binary_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_ard_binary_twenty_four_taxa"
    )
    multistate_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_ard_multistate_twenty_four_taxa"
    )
    binary_missing_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_ard_binary_missing_twenty_four_taxa"
    )
    multistate_missing_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_ard_multistate_missing_twenty_four_taxa"
    )

    binary_report = fit_discrete_mk_model(
        binary_fixture.tree_path,
        binary_fixture.traits_path,
        trait=binary_fixture.trait_name,
        taxon_column=binary_fixture.taxon_column,
        model="all-rates-different",
    )
    multistate_report = fit_discrete_mk_model(
        multistate_fixture.tree_path,
        multistate_fixture.traits_path,
        trait=multistate_fixture.trait_name,
        taxon_column=multistate_fixture.taxon_column,
        model="all-rates-different",
    )
    binary_missing_report = fit_discrete_mk_model(
        binary_missing_fixture.tree_path,
        binary_missing_fixture.traits_path,
        trait=binary_missing_fixture.trait_name,
        taxon_column=binary_missing_fixture.taxon_column,
        model="all-rates-different",
    )
    multistate_missing_report = fit_discrete_mk_model(
        multistate_missing_fixture.tree_path,
        multistate_missing_fixture.traits_path,
        trait=multistate_missing_fixture.trait_name,
        taxon_column=multistate_missing_fixture.taxon_column,
        model="all-rates-different",
    )

    binary_rate_lookup = {
        (row.source_state, row.target_state): row.rate
        for row in binary_report.transition_rate_rows
        if row.transition_allowed
    }

    assert "all-rates-different-model" in binary_fixture.feature_tags
    assert "all-rates-different-model" in multistate_fixture.feature_tags
    assert "all-rates-different-model" in binary_missing_fixture.feature_tags
    assert "all-rates-different-model" in multistate_missing_fixture.feature_tags
    assert binary_report.model == "all-rates-different"
    assert binary_report.parameter_count == 2
    assert multistate_report.parameter_count == 12
    assert binary_rate_lookup[("0", "1")] > 0.0
    assert binary_rate_lookup[("1", "0")] > 0.0
    assert multistate_report.optimizer_diagnostics.hit_lower_parameter_bound is True
    assert binary_missing_report.input_audit.pruned_missing_value_taxa == ["Phy10"]
    assert multistate_missing_report.input_audit.pruned_missing_value_taxa == ["Phy14"]


def test_shared_phytools_comparative_fixture_catalog_covers_discrete_nonultrametric_and_branch_edge_surfaces() -> (
    None
):
    binary_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_binary_twenty_four_taxa"
    )
    multistate_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_multistate_twenty_four_taxa"
    )
    nonultrametric_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_strong_signal_non_ultrametric_twenty_four_taxa"
    )
    branch_edge_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_strong_signal_branch_edge_twenty_four_taxa"
    )

    binary_values = set(
        _read_trait_values(binary_fixture.fixture_id, binary_fixture.trait_name)
    )
    multistate_values = set(
        _read_trait_values(multistate_fixture.fixture_id, multistate_fixture.trait_name)
    )
    nonultrametric_signal = estimate_pagels_lambda(
        nonultrametric_fixture.tree_path,
        nonultrametric_fixture.traits_path,
        trait=nonultrametric_fixture.trait_name,
    )
    branch_edge_tree = load_tree(branch_edge_fixture.tree_path)
    branch_lengths = [
        node.branch_length
        for node in branch_edge_tree.iter_nodes()
        if node.branch_length is not None
    ]

    assert binary_values == {"0", "1"}
    assert len(multistate_values) >= 3
    assert nonultrametric_signal.input_audit.tree_is_ultrametric is False
    assert branch_edge_tree.tip_count == 24
    assert min(branch_lengths) == 0.0
    assert max(branch_lengths) > 2.0


def test_shared_phytools_comparative_fixture_catalog_supports_brownian_pgls_cases() -> (
    None
):
    continuous_fixture = get_shared_phytools_comparative_fixture(
        "phytools_pgls_brownian_continuous_four_taxa"
    )
    categorical_fixture = get_shared_phytools_comparative_fixture(
        "phytools_pgls_brownian_categorical_eight_taxa"
    )
    interaction_fixture = get_shared_phytools_comparative_fixture(
        "phytools_pgls_brownian_interaction_eight_taxa"
    )

    continuous_report = run_pgls(
        continuous_fixture.tree_path,
        continuous_fixture.traits_path,
        response=continuous_fixture.trait_name,
        predictors=["predictor_one"],
        taxon_column=continuous_fixture.taxon_column,
        lambda_value=1.0,
    )
    categorical_report = run_pgls(
        categorical_fixture.tree_path,
        categorical_fixture.traits_path,
        formula="response ~ habitat + diet",
        taxon_column=categorical_fixture.taxon_column,
        lambda_value=1.0,
    )
    interaction_report = run_pgls(
        interaction_fixture.tree_path,
        interaction_fixture.traits_path,
        formula="response ~ habitat * diet",
        taxon_column=interaction_fixture.taxon_column,
        lambda_value=1.0,
    )
    categorical_model_matrix = build_pgls_model_matrix(
        categorical_fixture.tree_path,
        categorical_fixture.traits_path,
        formula="response ~ habitat + diet",
        taxon_column=categorical_fixture.taxon_column,
    )
    interaction_model_matrix = build_pgls_model_matrix(
        interaction_fixture.tree_path,
        interaction_fixture.traits_path,
        formula="response ~ habitat * diet",
        taxon_column=interaction_fixture.taxon_column,
    )

    assert continuous_report.lambda_fit.mode == "fixed"
    assert continuous_report.lambda_value == 1.0
    assert continuous_report.encoded_columns == ["intercept", "predictor_one"]
    assert continuous_report.coefficients[1].name == "predictor_one"
    assert "continuous-predictor-case" in continuous_fixture.feature_tags
    assert categorical_report.lambda_fit.mode == "fixed"
    assert categorical_report.lambda_value == 1.0
    assert categorical_report.encoded_columns == [
        "intercept",
        "habitat[tundra]",
        "diet[herbivore]",
    ]
    assert "categorical-predictor-case" in categorical_fixture.feature_tags
    assert categorical_model_matrix.encoded_columns == [
        "intercept",
        "habitat[tundra]",
        "diet[herbivore]",
    ]
    assert interaction_report.lambda_fit.mode == "fixed"
    assert interaction_report.lambda_value == 1.0
    assert interaction_report.encoded_columns == [
        "intercept",
        "habitat[tundra]",
        "diet[herbivore]",
        "habitat[tundra]:diet[herbivore]",
    ]
    assert "interaction-term-case" in interaction_fixture.feature_tags
    assert interaction_model_matrix.encoded_columns == [
        "intercept",
        "habitat[tundra]",
        "diet[herbivore]",
        "habitat[tundra]:diet[herbivore]",
    ]


def test_shared_phytools_comparative_fixture_catalog_supports_phylogenetic_residual_cases() -> (
    None
):
    brownian_fixture = get_shared_phytools_comparative_fixture(
        "phytools_phyl_resid_bm_allometry_six_taxa"
    )
    lambda_fixture = get_shared_phytools_comparative_fixture(
        "phytools_phyl_resid_lambda_allometry_six_taxa"
    )
    missing_fixture = get_shared_phytools_comparative_fixture(
        "phytools_phyl_resid_lambda_missing_six_taxa"
    )

    brownian_report = summarize_phylogenetic_residuals(
        brownian_fixture.tree_path,
        brownian_fixture.traits_path,
        response=brownian_fixture.trait_name,
        predictor="body_mass",
        taxon_column=brownian_fixture.taxon_column,
        method="brownian",
    )
    lambda_report = summarize_phylogenetic_residuals(
        lambda_fixture.tree_path,
        lambda_fixture.traits_path,
        response=lambda_fixture.trait_name,
        predictor="body_mass",
        taxon_column=lambda_fixture.taxon_column,
        method="lambda",
    )
    missing_report = summarize_phylogenetic_residuals(
        missing_fixture.tree_path,
        missing_fixture.traits_path,
        response=missing_fixture.trait_name,
        predictor="body_mass",
        taxon_column=missing_fixture.taxon_column,
        method="lambda",
    )

    assert brownian_report.lambda_estimation_mode == "fixed"
    assert brownian_report.lambda_value == 1.0
    assert lambda_report.lambda_estimation_mode == "estimated"
    assert lambda_report.analyzed_taxon_count == 6
    assert "F" in lambda_report.outlier_taxa
    assert {row.taxon: row.reason for row in missing_report.excluded_taxa} == {
        "E": "missing_value",
        "G": "absent_from_tree",
    }


def test_shared_phytools_comparative_fixture_catalog_supports_phylogenetic_anova_cases() -> (
    None
):
    clean_fixture = get_shared_phytools_comparative_fixture(
        "phytools_phyl_anova_group_effect_six_taxa"
    )
    missing_fixture = get_shared_phytools_comparative_fixture(
        "phytools_phyl_anova_group_effect_missing_six_taxa"
    )

    clean_report = summarize_phylogenetic_anova(
        clean_fixture.tree_path,
        clean_fixture.traits_path,
        response=clean_fixture.trait_name,
        group="habitat",
        taxon_column=clean_fixture.taxon_column,
        simulations=16,
        seed=9,
    )
    missing_report = summarize_phylogenetic_anova(
        missing_fixture.tree_path,
        missing_fixture.traits_path,
        response=missing_fixture.trait_name,
        group="habitat",
        taxon_column=missing_fixture.taxon_column,
        simulations=16,
        seed=9,
    )

    assert clean_report.group_count == 2
    assert {row.group: row.taxon_count for row in clean_report.group_rows} == {
        "desert": 2,
        "forest": 4,
    }
    assert clean_report.pairwise_rows[0].left_group == "desert"
    assert clean_report.pairwise_rows[0].right_group == "forest"
    assert {row.taxon: row.reason for row in missing_report.excluded_taxa} == {
        "F": "missing_value",
        "G": "absent_from_tree",
    }
