from __future__ import annotations

import math

from bijux_phylogenetics.comparative.search import (
    BoundedSearchControls,
    fit_bounded_scalar_parameter,
)


def test_comparative_bounded_scalar_fit_reports_shared_optimizer_diagnostics() -> None:
    fit = fit_bounded_scalar_parameter(
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
    )

    assert math.isclose(fit.parameter_value, 1.5, abs_tol=1e-6)
    assert fit.payload == fit.parameter_value
    assert fit.optimizer_diagnostics.lower_bound == 0.0
    assert fit.optimizer_diagnostics.upper_bound == 4.0
    assert (
        fit.optimizer_diagnostics.starting_parameter_policy
        == "user-provided-first-evaluation"
    )
    assert fit.optimizer_diagnostics.starting_parameter_value == 3.0
    assert fit.optimizer_diagnostics.converged is True
    assert fit.optimizer_diagnostics.failure_reason is None
    assert fit.lower_95_confidence_interval == 1.0
    assert fit.upper_95_confidence_interval == 2.0
    assert [row.parameter_value for row in fit.profile_rows] == [
        0.0,
        1.0,
        1.5,
        2.0,
        4.0,
    ]
    assert fit.profile_rows[2].within_95_confidence_interval is True
    assert fit.profile_rows[0].within_95_confidence_interval is False
