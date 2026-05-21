from __future__ import annotations

import csv

import pytest

from bijux_phylogenetics.ancestral.common import load_discrete_dataset
from bijux_phylogenetics.comparative.discrete_mk import fit_discrete_mk_model
import bijux_phylogenetics.datasets.shared_fixtures as fixtures_api
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_geiger_discrete_fixture,
    list_shared_geiger_discrete_fixtures,
)
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError


def _read_taxa(fixture_id: str) -> list[str]:
    fixture = get_shared_geiger_discrete_fixture(fixture_id)
    with fixture.traits_path.open(encoding="utf-8", newline="") as handle:
        return [
            row[fixture.taxon_column] for row in csv.DictReader(handle, delimiter="\t")
        ]


def test_shared_geiger_discrete_fixture_catalog_covers_required_goal_cases() -> None:
    fixtures = list_shared_geiger_discrete_fixtures()
    feature_tags = {tag for fixture in fixtures for tag in fixture.feature_tags}
    supported_models = {
        model for fixture in fixtures for model in fixture.supported_model_names
    }

    assert {
        "binary-discrete-trait",
        "three-state-discrete-trait",
        "four-plus-state-discrete-trait",
        "equal-rates-known-truth",
        "symmetric-known-truth",
        "all-rates-different-known-truth",
        "weak-phylogenetic-signal-review",
        "rare-state-case",
        "missing-states",
        "constant-state-negative-case",
        "overparameterized-multistate-case",
        "tree-trait-taxon-mismatch",
        "known-transition-matrix-metadata",
    } <= feature_tags
    assert {"ER", "SYM", "ARD"} <= supported_models


def test_shared_geiger_discrete_fixture_lookup_resolves_linked_shared_surfaces() -> (
    None
):
    fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_sym_three_state_twenty_four_taxa"
    )

    assert fixture.tree_fixture.tip_count == 24
    assert fixture.trait_table_fixture.row_count == 24
    assert fixture.trait_kind == "three-state-discrete"
    assert fixture.supported_model_names == ("SYM",)
    assert fixture.transition_matrix_metadata == {
        "states": ["north", "south", "central"],
        "matrix_shape": "3x3",
        "parameterization": "symmetric",
        "pairwise_rate": 0.45,
    }


@pytest.mark.slow
def test_shared_geiger_discrete_fixture_catalog_supports_er_sym_and_ard_fits() -> None:
    er_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_er_binary_twenty_four_taxa"
    )
    sym_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_sym_three_state_twenty_four_taxa"
    )
    ard_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_ard_four_state_twenty_four_taxa"
    )

    er_report = fit_discrete_mk_model(
        er_fixture.tree_path,
        er_fixture.traits_path,
        trait=er_fixture.trait_name,
        taxon_column=er_fixture.taxon_column,
        model="equal-rates",
    )
    sym_report = fit_discrete_mk_model(
        sym_fixture.tree_path,
        sym_fixture.traits_path,
        trait=sym_fixture.trait_name,
        taxon_column=sym_fixture.taxon_column,
        model="symmetric",
    )
    ard_report = fit_discrete_mk_model(
        ard_fixture.tree_path,
        ard_fixture.traits_path,
        trait=ard_fixture.trait_name,
        taxon_column=ard_fixture.taxon_column,
        model="all-rates-different",
    )

    assert er_report.parameter_count == 1
    assert er_report.input_audit.observed_states == ["0", "1"]
    assert sym_report.parameter_count == 3
    assert sym_report.input_audit.observed_states == ["central", "north", "south"]
    assert ard_report.parameter_count == 12
    assert ard_report.input_audit.observed_states == ["east", "north", "south", "west"]
    assert ard_report.baseline_comparison is not None


@pytest.mark.slow
def test_shared_geiger_discrete_fixture_catalog_supports_lambda_transform_review() -> (
    None
):
    strong_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_er_binary_twenty_four_taxa"
    )
    weak_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_transform_weak_signal_twenty_four_taxa"
    )
    missing_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_lambda_missing_binary_twenty_four_taxa"
    )

    strong_report = fit_discrete_mk_model(
        strong_fixture.tree_path,
        strong_fixture.traits_path,
        trait=strong_fixture.trait_name,
        taxon_column=strong_fixture.taxon_column,
        model="equal-rates",
        transform="lambda",
    )
    weak_report = fit_discrete_mk_model(
        weak_fixture.tree_path,
        weak_fixture.traits_path,
        trait=weak_fixture.trait_name,
        taxon_column=weak_fixture.taxon_column,
        model="equal-rates",
        transform="lambda",
    )
    missing_report = fit_discrete_mk_model(
        missing_fixture.tree_path,
        missing_fixture.traits_path,
        trait=missing_fixture.trait_name,
        taxon_column=missing_fixture.taxon_column,
        model="equal-rates",
        transform="lambda",
    )

    assert strong_report.transform_fit is not None
    assert strong_report.transform_fit.parameter_name == "lambda"
    assert strong_report.parameter_count == 2
    assert weak_report.transform_fit is not None
    assert weak_report.transform_fit.parameter_value <= 1e-6
    assert any(
        warning.kind == "weak_phylogenetic_signal"
        for warning in weak_report.transform_fit.warnings
    )
    assert missing_report.transform_fit is not None
    assert missing_report.input_audit.pruned_missing_value_taxa == ["Phy10"]
    assert missing_report.parameter_count == 2


@pytest.mark.slow
def test_shared_geiger_discrete_fixture_catalog_supports_kappa_transform_review() -> (
    None
):
    strong_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_kappa_branch_sensitive_twenty_four_taxa"
    )
    weak_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_kappa_weak_signal_twenty_four_taxa"
    )
    missing_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_kappa_missing_three_state_twenty_four_taxa"
    )

    strong_report = fit_discrete_mk_model(
        strong_fixture.tree_path,
        strong_fixture.traits_path,
        trait=strong_fixture.trait_name,
        taxon_column=strong_fixture.taxon_column,
        model="equal-rates",
        transform="kappa",
    )
    weak_report = fit_discrete_mk_model(
        weak_fixture.tree_path,
        weak_fixture.traits_path,
        trait=weak_fixture.trait_name,
        taxon_column=weak_fixture.taxon_column,
        model="equal-rates",
        transform="kappa",
    )
    missing_report = fit_discrete_mk_model(
        missing_fixture.tree_path,
        missing_fixture.traits_path,
        trait=missing_fixture.trait_name,
        taxon_column=missing_fixture.taxon_column,
        model="symmetric",
        transform="kappa",
    )

    assert strong_report.transform_fit is not None
    assert strong_report.transform_fit.parameter_name == "kappa"
    assert strong_report.parameter_count == 2
    assert weak_report.transform_fit is not None
    assert weak_report.transform_fit.parameter_value <= 1e-6
    assert any(
        warning.kind == "branch_length_flattening_limit"
        for warning in weak_report.transform_fit.warnings
    )
    assert missing_report.transform_fit is not None
    assert missing_report.input_audit.pruned_missing_value_taxa == ["Phy14"]
    assert missing_report.parameter_count == 4


@pytest.mark.slow
def test_shared_geiger_discrete_fixture_catalog_supports_delta_transform_review() -> (
    None
):
    boundary_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_delta_late_change_binary_twenty_four_taxa"
    )
    time_sensitive_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_delta_time_sensitive_twenty_four_taxa"
    )
    missing_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_delta_missing_binary_twenty_four_taxa"
    )

    boundary_report = fit_discrete_mk_model(
        boundary_fixture.tree_path,
        boundary_fixture.traits_path,
        trait=boundary_fixture.trait_name,
        taxon_column=boundary_fixture.taxon_column,
        model="equal-rates",
        transform="delta",
    )
    time_sensitive_report = fit_discrete_mk_model(
        time_sensitive_fixture.tree_path,
        time_sensitive_fixture.traits_path,
        trait=time_sensitive_fixture.trait_name,
        taxon_column=time_sensitive_fixture.taxon_column,
        model="equal-rates",
        transform="delta",
    )
    missing_report = fit_discrete_mk_model(
        missing_fixture.tree_path,
        missing_fixture.traits_path,
        trait=missing_fixture.trait_name,
        taxon_column=missing_fixture.taxon_column,
        model="equal-rates",
        transform="delta",
    )

    assert boundary_report.transform_fit is not None
    assert boundary_report.transform_fit.parameter_name == "delta"
    assert boundary_report.parameter_count == 2
    assert boundary_report.transform_baseline_comparison is not None
    assert time_sensitive_report.transform_fit is not None
    assert time_sensitive_report.transform_fit.parameter_value > 0.0
    assert time_sensitive_report.transform_fit.parameter_value < 3.0
    assert missing_report.transform_fit is not None
    assert missing_report.input_audit.pruned_missing_value_taxa == ["Phy10"]
    assert missing_report.parameter_count == 2


@pytest.mark.slow
def test_shared_geiger_discrete_fixture_catalog_supports_early_burst_review() -> None:
    early_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_early_burst_early_change_twenty_four_taxa"
    )
    weak_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_early_burst_weak_signal_twenty_four_taxa"
    )
    late_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_early_burst_late_change_twenty_four_taxa"
    )
    missing_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_early_burst_missing_binary_twenty_four_taxa"
    )

    early_report = fit_discrete_mk_model(
        early_fixture.tree_path,
        early_fixture.traits_path,
        trait=early_fixture.trait_name,
        taxon_column=early_fixture.taxon_column,
        model="equal-rates",
        transform="early-burst",
    )
    weak_report = fit_discrete_mk_model(
        weak_fixture.tree_path,
        weak_fixture.traits_path,
        trait=weak_fixture.trait_name,
        taxon_column=weak_fixture.taxon_column,
        model="equal-rates",
        transform="early-burst",
    )
    late_report = fit_discrete_mk_model(
        late_fixture.tree_path,
        late_fixture.traits_path,
        trait=late_fixture.trait_name,
        taxon_column=late_fixture.taxon_column,
        model="equal-rates",
        transform="early-burst",
    )
    missing_report = fit_discrete_mk_model(
        missing_fixture.tree_path,
        missing_fixture.traits_path,
        trait=missing_fixture.trait_name,
        taxon_column=missing_fixture.taxon_column,
        model="equal-rates",
        transform="early-burst",
    )

    assert early_report.transform_fit is not None
    assert early_report.transform_fit.parameter_name == "a"
    assert early_report.transform_fit.parameter_value > 0.0
    assert early_report.parameter_count == 2
    assert weak_report.transform_fit is not None
    assert weak_report.transform_fit.parameter_value == 0.0
    assert any(
        warning.kind == "brownian_like_rate_change"
        for warning in weak_report.transform_fit.warnings
    )
    assert late_report.transform_fit is not None
    assert late_report.transform_fit.parameter_value < 0.0
    assert late_report.parameter_count == 2
    assert missing_report.transform_fit is not None
    assert missing_report.input_audit.pruned_missing_value_taxa == ["Phy10"]
    assert missing_report.parameter_count == 2


@pytest.mark.slow
def test_shared_geiger_discrete_fixture_catalog_handles_missing_sparse_and_mismatch_cases() -> (
    None
):
    missing_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_missing_three_state_twenty_four_taxa"
    )
    sparse_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_sparse_six_state_twenty_four_taxa"
    )
    mismatch_fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_mismatch_four_state_twenty_four_taxa"
    )

    missing_report = fit_discrete_mk_model(
        missing_fixture.tree_path,
        missing_fixture.traits_path,
        trait=missing_fixture.trait_name,
        taxon_column=missing_fixture.taxon_column,
        model="symmetric",
    )
    sparse_report = fit_discrete_mk_model(
        sparse_fixture.tree_path,
        sparse_fixture.traits_path,
        trait=sparse_fixture.trait_name,
        taxon_column=sparse_fixture.taxon_column,
        model="all-rates-different",
    )
    mismatch_dataset = load_discrete_dataset(
        mismatch_fixture.tree_path,
        mismatch_fixture.traits_path,
        trait=mismatch_fixture.trait_name,
        taxon_column=mismatch_fixture.taxon_column,
    )

    assert missing_report.input_audit.pruned_missing_value_taxa == ["Phy14"]
    assert sparse_report.parameter_count == 30
    assert sparse_report.overparameterized is True
    assert sparse_report.input_audit.sparse_states == ["f"]
    assert any(
        "overparameterized" in warning for warning in sparse_report.input_audit.warnings
    )
    assert len(mismatch_dataset.taxa) == 23
    assert mismatch_dataset.dropped_missing_taxa == []
    assert mismatch_dataset.alignment_report.dropped_tree_taxa == ["Phy9"]
    assert mismatch_dataset.alignment_report.dropped_trait_taxa == ["PhyExtra"]
    assert mismatch_dataset.alignment_report.dropped_missing_value_taxa == []
    assert "PhyExtra" in _read_taxa(
        "geiger_discrete_mismatch_four_state_twenty_four_taxa"
    )
    assert "Phy9" not in _read_taxa(
        "geiger_discrete_mismatch_four_state_twenty_four_taxa"
    )


def test_shared_geiger_discrete_fixture_catalog_blocks_constant_state_inputs() -> None:
    fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_constant_negative_twenty_four_taxa"
    )

    with pytest.raises(AncestralReconstructionError):
        fit_discrete_mk_model(
            fixture.tree_path,
            fixture.traits_path,
            trait=fixture.trait_name,
            taxon_column=fixture.taxon_column,
            model="equal-rates",
        )


def test_public_runtime_exports_include_shared_geiger_discrete_fixture_surface() -> (
    None
):
    assert (
        fixtures_api.get_shared_geiger_discrete_fixture
        is get_shared_geiger_discrete_fixture
    )
    assert (
        fixtures_api.list_shared_geiger_discrete_fixtures
        is list_shared_geiger_discrete_fixtures
    )
