from __future__ import annotations

import math

import bijux_phylogenetics.comparative.search.bounded as comparative_bounded_module
import bijux_phylogenetics.comparative.search.bounded_scalar_fit as comparative_fit_module
import bijux_phylogenetics.phylo.likelihood.parameter_search as phylo_parameter_search_module
from bijux_phylogenetics.comparative.search import BoundedSearchControls


def test_phylo_bounded_likelihood_search_delegates_to_shared_golden_section_owner(
    monkeypatch,
) -> None:
    observed_owner_names: list[str] = []
    observed_parameter_names: list[str] = []
    original_search = (
        phylo_parameter_search_module.run_bounded_parameter_golden_section_maximization
    )

    def observe_search(**kwargs):
        observed_owner_names.append(kwargs["owner_name"])
        observed_parameter_names.append(kwargs["parameter_name"])
        return original_search(**kwargs)

    monkeypatch.setattr(
        phylo_parameter_search_module,
        "run_bounded_parameter_golden_section_maximization",
        observe_search,
    )

    result = phylo_parameter_search_module.run_bounded_likelihood_search(
        lower_bound=0.0,
        upper_bound=4.0,
        evaluate=lambda parameter_value: (
            parameter_value,
            -((parameter_value - 2.0) ** 2),
        ),
        parameter_name="branch length",
        owner_name="JC69 branch-length fit",
    )

    assert observed_owner_names == ["JC69 branch-length fit"]
    assert observed_parameter_names == ["branch length"]
    assert math.isclose(result.parameter_value, 2.0, abs_tol=1e-6)


def test_phylo_coordinate_search_delegates_to_shared_coordinate_owner(
    monkeypatch,
) -> None:
    observed_owner_names: list[str] = []
    original_search = phylo_parameter_search_module.run_bounded_coordinate_optimization

    def observe_search(**kwargs):
        observed_owner_names.append(kwargs["owner_name"])
        return original_search(**kwargs)

    monkeypatch.setattr(
        phylo_parameter_search_module,
        "run_bounded_coordinate_optimization",
        observe_search,
    )

    result = phylo_parameter_search_module.run_bounded_coordinate_likelihood_search(
        initial_values={"alpha": 0.8, "kappa": 1.0},
        bounds_by_name={"alpha": (0.1, 4.0), "kappa": (0.1, 5.0)},
        evaluate=lambda values: (
            dict(values),
            -((values["alpha"] - 1.5) ** 2) - ((values["kappa"] - 2.0) ** 2),
        ),
        owner_name="joint substitution fit",
    )

    assert observed_owner_names == ["joint substitution fit"]
    assert math.isclose(result.parameter_values["alpha"], 1.5, abs_tol=1e-6)
    assert math.isclose(result.parameter_values["kappa"], 2.0, abs_tol=1e-6)


def test_comparative_bounded_search_delegates_to_shared_parameter_owner(
    monkeypatch,
) -> None:
    observed_owner_names: list[str] = []
    observed_parameter_names: list[str] = []
    original_search = comparative_bounded_module.run_bounded_parameter_maximization

    def observe_search(**kwargs):
        observed_owner_names.append(kwargs["owner_name"])
        observed_parameter_names.append(kwargs["parameter_name"])
        return original_search(**kwargs)

    monkeypatch.setattr(
        comparative_bounded_module,
        "run_bounded_parameter_maximization",
        observe_search,
    )

    result = comparative_bounded_module.run_bounded_maximization(
        lower_bound=0.0,
        upper_bound=4.0,
        controls=BoundedSearchControls(
            coarse_grid_point_count=9,
            fine_grid_point_count=11,
            initial_parameter_value=3.0,
            refinement_start_count=2,
        ),
        evaluate=lambda parameter_value: (
            parameter_value,
            -((parameter_value - 1.5) ** 2),
        ),
        parameter_name="pagel lambda",
        owner_name="comparative lambda search",
    )

    assert observed_owner_names == ["comparative lambda search"]
    assert observed_parameter_names == ["pagel lambda"]
    assert math.isclose(result.parameter_value, 1.5, abs_tol=1e-6)


def test_comparative_profiled_fit_delegates_to_shared_profile_owner(
    monkeypatch,
) -> None:
    observed_owner_names: list[str] = []
    observed_parameter_names: list[str] = []
    original_fit = comparative_fit_module.fit_profiled_bounded_parameter

    def observe_fit(**kwargs):
        observed_owner_names.append(kwargs["owner_name"])
        observed_parameter_names.append(kwargs["parameter_name"])
        return original_fit(**kwargs)

    monkeypatch.setattr(
        comparative_fit_module,
        "fit_profiled_bounded_parameter",
        observe_fit,
    )

    result = comparative_fit_module.fit_bounded_scalar_parameter(
        lower_bound=0.0,
        upper_bound=4.0,
        controls=BoundedSearchControls(
            coarse_grid_point_count=9,
            fine_grid_point_count=11,
            initial_parameter_value=3.0,
            refinement_start_count=2,
        ),
        evaluate=lambda parameter_value: (
            parameter_value,
            -((parameter_value - 1.5) ** 2),
        ),
        profile_parameter_values=[0.0, 1.0, 1.5, 2.0, 4.0],
        confidence_interval_drop=0.5,
        parameter_name="OU alpha",
        owner_name="comparative OU alpha fit",
    )

    assert observed_owner_names == ["comparative OU alpha fit"]
    assert observed_parameter_names == ["OU alpha"]
    assert math.isclose(result.parameter_value, 1.5, abs_tol=1e-6)
    assert result.lower_95_confidence_interval == 1.0
    assert result.upper_95_confidence_interval == 2.0
