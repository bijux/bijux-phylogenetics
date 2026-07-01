from __future__ import annotations

from pathlib import Path

import bijux_phylogenetics.comparative.continuous.model_fitting as continuous_model_fitting_module
import bijux_phylogenetics.comparative.discrete_mk.transforms as discrete_transform_module
import bijux_phylogenetics.comparative.evolutionary_modes.search as continuous_search_module
import bijux_phylogenetics.comparative.pgls.fitting as pgls_fitting_module
import bijux_phylogenetics.comparative.pgls.ou_covariance as pgls_ou_covariance_module
import bijux_phylogenetics.comparative.signal.core as signal_module
from bijux_phylogenetics.comparative import fit_continuous_evolutionary_mode
from bijux_phylogenetics.comparative.continuous import fit_ornstein_uhlenbeck_model
from bijux_phylogenetics.comparative.discrete_mk import fit_discrete_mk_model
from bijux_phylogenetics.comparative.pgls import run_pgls
from bijux_phylogenetics.comparative.pgls.ou_covariance import (
    summarize_ou_covariance_pgls,
)
from bijux_phylogenetics.comparative.signal import estimate_pagels_lambda
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_geiger_continuous_fixture,
    get_shared_geiger_discrete_fixture,
)

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


def test_pagel_kappa_routes_through_shared_bounded_scalar_search(monkeypatch) -> None:
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
        mode="pagel-kappa",
        taxon_column=fixture.taxon_column,
    )

    assert fit.optimizer_diagnostics is not None
    assert observed_bounds == [(0.0, 3.0)]


def test_pagel_delta_routes_through_shared_bounded_scalar_search(monkeypatch) -> None:
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
        mode="pagel-delta",
        taxon_column=fixture.taxon_column,
    )

    assert fit.optimizer_diagnostics is not None
    assert observed_bounds == [(0.0, 3.0)]


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


def test_discrete_lambda_transform_routes_through_shared_bounded_scalar_search(
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
        "geiger_discrete_lambda_missing_binary_twenty_four_taxa"
    )
    report = fit_discrete_mk_model(
        fixture.tree_path,
        fixture.traits_path,
        trait=fixture.trait_name,
        taxon_column=fixture.taxon_column,
        model="equal-rates",
        transform="lambda",
    )

    assert report.transform_fit is not None
    assert observed_bounds == [(0.0, 1.0)]


def test_signal_pagel_lambda_routes_through_shared_bounded_scalar_fit(
    monkeypatch,
) -> None:
    observed_calls: list[tuple[tuple[float, float], str, str]] = []
    original_fit = signal_module.fit_bounded_scalar_parameter

    def observe_fit(**kwargs):
        observed_calls.append(
            (
                (kwargs["lower_bound"], kwargs["upper_bound"]),
                kwargs["parameter_name"],
                kwargs["owner_name"],
            )
        )
        return original_fit(**kwargs)

    monkeypatch.setattr(signal_module, "fit_bounded_scalar_parameter", observe_fit)

    report = estimate_pagels_lambda(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )

    assert report.optimizer_diagnostics is not None
    assert observed_calls == [
        ((0.0, 1.0), "lambda", "phylogenetic signal Pagel lambda fit")
    ]


def test_pgls_pagel_lambda_routes_through_shared_bounded_scalar_fit(
    monkeypatch,
) -> None:
    observed_calls: list[tuple[tuple[float, float], str, str]] = []
    original_fit = pgls_fitting_module.fit_bounded_scalar_parameter

    def observe_fit(**kwargs):
        observed_calls.append(
            (
                (kwargs["lower_bound"], kwargs["upper_bound"]),
                kwargs["parameter_name"],
                kwargs["owner_name"],
            )
        )
        return original_fit(**kwargs)

    monkeypatch.setattr(pgls_fitting_module, "fit_bounded_scalar_parameter", observe_fit)

    report = run_pgls(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        lambda_value="estimate",
    )

    assert report.lambda_fit.optimizer_diagnostics is not None
    assert observed_calls == [((0.0, 1.0), "lambda", "PGLS lambda fit")]


def test_ou_covariance_alpha_routes_through_shared_bounded_scalar_fit(
    monkeypatch,
) -> None:
    observed_calls: list[tuple[tuple[float, float], str, str]] = []
    original_fit = pgls_ou_covariance_module.fit_bounded_scalar_parameter

    def observe_fit(**kwargs):
        observed_calls.append(
            (
                (kwargs["lower_bound"], kwargs["upper_bound"]),
                kwargs["parameter_name"],
                kwargs["owner_name"],
            )
        )
        return original_fit(**kwargs)

    monkeypatch.setattr(
        pgls_ou_covariance_module,
        "fit_bounded_scalar_parameter",
        observe_fit,
    )

    report = summarize_ou_covariance_pgls(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        alpha="estimate",
    )

    assert report.alpha_optimizer_diagnostics is not None
    assert len(observed_calls) == 1
    bounds, parameter_name, owner_name = observed_calls[0]
    assert bounds[0] > 0.0
    assert bounds[1] > bounds[0]
    assert parameter_name == "alpha"
    assert owner_name == "OU covariance PGLS alpha fit"


def test_ou_trait_alpha_routes_through_shared_bounded_scalar_fit(monkeypatch) -> None:
    observed_calls: list[tuple[tuple[float, float], str, str]] = []
    original_fit = continuous_model_fitting_module.fit_bounded_scalar_parameter

    def observe_fit(**kwargs):
        observed_calls.append(
            (
                (kwargs["lower_bound"], kwargs["upper_bound"]),
                kwargs["parameter_name"],
                kwargs["owner_name"],
            )
        )
        return original_fit(**kwargs)

    monkeypatch.setattr(
        continuous_model_fitting_module,
        "fit_bounded_scalar_parameter",
        observe_fit,
    )

    report = fit_ornstein_uhlenbeck_model(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )

    assert report.optimizer_diagnostics is not None
    first_bounds, first_parameter_name, first_owner_name = observed_calls[0]
    assert first_bounds[0] > 0.0
    assert first_bounds[1] > first_bounds[0]
    assert first_parameter_name == "alpha"
    assert first_owner_name == "OU trait model alpha fit"
    assert observed_calls[1] == (
        (0.0, 1.0),
        "lambda",
        "continuous residual lambda fit",
    )
