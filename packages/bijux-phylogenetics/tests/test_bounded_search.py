from __future__ import annotations

import math

from bijux_phylogenetics.comparative.bounded_search import (
    BoundedSearchControls,
    run_bounded_maximization,
)


def test_run_bounded_maximization_records_user_start_and_multi_start_count() -> None:
    result = run_bounded_maximization(
        lower_bound=0.0,
        upper_bound=4.0,
        controls=BoundedSearchControls(
            coarse_grid_point_count=9,
            fine_grid_point_count=11,
            initial_parameter_value=3.0,
            refinement_start_count=2,
        ),
        evaluate=lambda parameter: (parameter, -((parameter - 1.75) ** 2)),
    )

    assert math.isclose(result.parameter_value, 1.75, rel_tol=0.0, abs_tol=5e-3)
    assert result.payload == result.parameter_value
    assert result.diagnostics.starting_parameter_policy == "user-provided-first-evaluation"
    assert math.isclose(
        result.diagnostics.starting_parameter_value,
        3.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert result.diagnostics.refinement_start_count == 2
    assert result.diagnostics.function_evaluation_count >= 9


def test_run_bounded_maximization_can_escape_weaker_coarse_candidate_with_multi_start() -> (
    None
):
    def objective(parameter: float) -> tuple[float, float]:
        left_peak = 9.95 - ((parameter - 0.95) ** 2)
        right_peak = 10.0 - (20.0 * ((parameter - 3.05) ** 2))
        return parameter, max(left_peak, right_peak)

    result = run_bounded_maximization(
        lower_bound=0.0,
        upper_bound=4.0,
        controls=BoundedSearchControls(
            coarse_grid_point_count=5,
            fine_grid_point_count=17,
            refinement_start_count=2,
        ),
        evaluate=objective,
    )

    assert math.isclose(result.parameter_value, 3.05, rel_tol=0.0, abs_tol=1e-2)
    assert result.diagnostics.refinement_start_count == 2
    assert not result.diagnostics.hit_lower_boundary
    assert not result.diagnostics.hit_upper_boundary
