from __future__ import annotations

import math

import pytest

from bijux_phylogenetics.core.bounded_parameter_optimization import (
    BoundedParameterSearchControls,
    fit_profiled_bounded_parameter,
    run_bounded_coordinate_optimization,
    run_bounded_parameter_golden_section_maximization,
    run_bounded_parameter_maximization,
)


def test_run_bounded_parameter_maximization_reports_shared_optimizer_audit_surface() -> (
    None
):
    result = run_bounded_parameter_maximization(
        lower_bound=0.0,
        upper_bound=4.0,
        controls=BoundedParameterSearchControls(
            coarse_grid_point_count=9,
            fine_grid_point_count=11,
            initial_parameter_value=3.0,
            refinement_start_count=2,
        ),
        evaluate=lambda parameter_value: (
            parameter_value,
            -((parameter_value - 1.5) ** 2),
        ),
        parameter_name="lambda",
        owner_name="shared comparative search",
    )

    assert math.isclose(result.parameter_value, 1.5, abs_tol=1e-6)
    assert result.payload == result.parameter_value
    assert result.diagnostics.lower_bound == 0.0
    assert result.diagnostics.upper_bound == 4.0
    assert (
        result.diagnostics.starting_parameter_policy == "user-provided-first-evaluation"
    )
    assert result.diagnostics.starting_parameter_value == 3.0
    assert result.converged is True
    assert result.diagnostics.failure_reason is None
    assert result.profile_rows


def test_run_bounded_parameter_golden_section_maximization_reports_failure_reason() -> (
    None
):
    result = run_bounded_parameter_golden_section_maximization(
        lower_bound=0.0,
        upper_bound=4.0,
        evaluate=lambda parameter_value: (
            parameter_value,
            -((parameter_value - 2.0) ** 2),
        ),
        parameter_name="branch length",
        owner_name="shared phylo search",
        tolerance=1e-30,
        max_iterations=1,
    )

    assert result.converged is False
    assert result.diagnostics.failure_reason == "max-iterations-exhausted"


def test_run_bounded_coordinate_optimization_reports_all_parameter_diagnostics() -> (
    None
):
    result = run_bounded_coordinate_optimization(
        initial_values={"kappa": 1.0, "alpha": 0.8},
        bounds_by_name={"kappa": (0.1, 5.0), "alpha": (0.1, 4.0)},
        evaluate=lambda values: (
            dict(values),
            -((values["kappa"] - 2.0) ** 2) - ((values["alpha"] - 1.5) ** 2),
        ),
        owner_name="joint substitution search",
    )

    assert math.isclose(result.parameter_values["kappa"], 2.0, abs_tol=1e-6)
    assert math.isclose(result.parameter_values["alpha"], 1.5, abs_tol=1e-6)
    assert result.payload == result.parameter_values
    assert result.function_evaluation_count >= 3
    assert result.optimization_pass_count >= 1
    assert result.converged is True
    assert set(result.parameter_diagnostics) == {"kappa", "alpha"}


def test_run_bounded_coordinate_optimization_rejects_mismatched_parameter_names() -> (
    None
):
    with pytest.raises(
        ValueError,
        match="identical parameter names in starts and bounds",
    ):
        run_bounded_coordinate_optimization(
            initial_values={"kappa": 1.0},
            bounds_by_name={"alpha": (0.1, 4.0)},
            evaluate=lambda values: (dict(values), 0.0),
        )


def test_fit_profiled_bounded_parameter_reports_supported_interval() -> None:
    fit = fit_profiled_bounded_parameter(
        lower_bound=0.0,
        upper_bound=4.0,
        controls=BoundedParameterSearchControls(
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
        parameter_name="lambda",
        owner_name="shared profile fit",
    )

    assert math.isclose(fit.parameter_value, 1.5, abs_tol=1e-6)
    assert fit.payload == fit.parameter_value
    assert fit.optimizer_diagnostics.starting_parameter_policy == (
        "user-provided-first-evaluation"
    )
    assert fit.lower_confidence_interval == 1.0
    assert fit.upper_confidence_interval == 2.0
    assert [row.parameter_value for row in fit.profile_rows] == [
        0.0,
        1.0,
        1.5,
        2.0,
        4.0,
    ]
    assert fit.profile_rows[2].within_confidence_interval is True
    assert fit.profile_rows[0].within_confidence_interval is False


def test_run_bounded_parameter_maximization_preserves_exact_interval_endpoints() -> (
    None
):
    visited: list[float] = []

    result = run_bounded_parameter_maximization(
        lower_bound=0.016667,
        upper_bound=3.333333,
        controls=BoundedParameterSearchControls(
            coarse_grid_point_count=81,
            fine_grid_point_count=11,
            initial_parameter_value=0.333333,
            refinement_start_count=3,
        ),
        evaluate=lambda parameter_value: (
            visited.append(parameter_value) or parameter_value,
            -((parameter_value - 1.0) ** 2),
        ),
        parameter_name="alpha",
        owner_name="OU trait model alpha fit",
    )

    assert result.converged is True
    assert min(visited) >= 0.016667
    assert max(visited) <= 3.333333
