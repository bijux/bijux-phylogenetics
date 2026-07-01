from __future__ import annotations

import bijux_phylogenetics.comparative.discrete_mk.transforms as discrete_transform_module
import bijux_phylogenetics.comparative.evolutionary_modes.search as continuous_search_module
from bijux_phylogenetics.comparative import fit_continuous_evolutionary_mode
from bijux_phylogenetics.comparative.discrete_mk import fit_discrete_mk_model
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_geiger_continuous_fixture,
    get_shared_geiger_discrete_fixture,
)


def test_pagel_lambda_routes_through_shared_bounded_scalar_search(monkeypatch) -> None:
    observed_bounds: list[tuple[float, float]] = []
    original_search = continuous_search_module.run_bounded_maximization

    def observe_search(**kwargs):
        observed_bounds.append((kwargs["lower_bound"], kwargs["upper_bound"]))
        return original_search(**kwargs)

    monkeypatch.setattr(
        continuous_search_module,
        "run_bounded_maximization",
        observe_search,
    )

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

    assert fit.optimizer_diagnostics is not None
    assert observed_bounds == [(0.0, 1.0)]


def test_discrete_kappa_transform_routes_through_shared_bounded_scalar_search(
    monkeypatch,
) -> None:
    observed_bounds: list[tuple[float, float]] = []
    original_search = discrete_transform_module.run_bounded_maximization

    def observe_search(**kwargs):
        observed_bounds.append((kwargs["lower_bound"], kwargs["upper_bound"]))
        return original_search(**kwargs)

    monkeypatch.setattr(
        discrete_transform_module,
        "run_bounded_maximization",
        observe_search,
    )

    fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_kappa_branch_sensitive_twenty_four_taxa"
    )
    report = fit_discrete_mk_model(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        taxon_column=fixture.taxon_column,
        model="equal-rates",
        transform="kappa",
    )

    assert report.transform_fit is not None
    assert observed_bounds == [(0.0, 1.0)]


def test_discrete_delta_transform_routes_through_shared_bounded_scalar_search(
    monkeypatch,
) -> None:
    observed_bounds: list[tuple[float, float]] = []
    original_search = discrete_transform_module.run_bounded_maximization

    def observe_search(**kwargs):
        observed_bounds.append((kwargs["lower_bound"], kwargs["upper_bound"]))
        return original_search(**kwargs)

    monkeypatch.setattr(
        discrete_transform_module,
        "run_bounded_maximization",
        observe_search,
    )

    fixture = get_shared_geiger_discrete_fixture(
        "geiger_discrete_delta_late_change_binary_twenty_four_taxa"
    )
    report = fit_discrete_mk_model(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        taxon_column=fixture.taxon_column,
        model="equal-rates",
        transform="delta",
    )

    assert report.transform_fit is not None
    assert observed_bounds == [(0.006737946999085467, 3.0)]
