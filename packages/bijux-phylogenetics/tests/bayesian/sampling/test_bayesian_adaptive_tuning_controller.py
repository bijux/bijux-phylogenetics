from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.adaptive_tuning import (
    AdaptiveTuningController,
    build_adaptive_tuning_controller,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_build_adaptive_tuning_controller_records_validated_policy() -> None:
    controller = build_adaptive_tuning_controller(
        proposal_name="normal-random-walk",
        scale_parameter_name="step-size",
        initial_scale=1.5,
        target_acceptance_rate=0.3,
        burnin_iteration_count=12,
        adaptation_window_size=4,
        decrease_factor=0.6,
        increase_factor=1.5,
        minimum_scale=0.1,
        maximum_scale=10.0,
    )

    assert isinstance(controller, AdaptiveTuningController)
    assert controller.proposal_name == "normal-random-walk"
    assert controller.scale_parameter_name == "step-size"
    assert controller.initial_scale == 1.5
    assert controller.target_acceptance_rate == 0.3
    assert controller.burnin_iteration_count == 12
    assert controller.adaptation_window_size == 4
    assert controller.decrease_factor == 0.6
    assert controller.increase_factor == 1.5
    assert controller.minimum_scale == 0.1
    assert controller.maximum_scale == 10.0


def test_build_adaptive_tuning_controller_accepts_zero_burnin_for_pre_frozen_policy() -> (
    None
):
    controller = build_adaptive_tuning_controller(
        proposal_name="clock-rate-scaling",
        scale_parameter_name="log-scale-standard-deviation",
        initial_scale=0.25,
        target_acceptance_rate=0.44,
        burnin_iteration_count=0,
        adaptation_window_size=5,
    )

    assert controller.burnin_iteration_count == 0


@pytest.mark.parametrize(
    ("field_name", "kwargs", "match"),
    [
        (
            "target_acceptance_rate",
            {"target_acceptance_rate": 1.0},
            "lie strictly between 0 and 1",
        ),
        (
            "decrease_factor",
            {"decrease_factor": 1.0},
            "lie strictly between 0.0 and 1.0",
        ),
        (
            "increase_factor",
            {"increase_factor": 1.0},
            "be greater than 1.0",
        ),
        (
            "minimum_scale",
            {"minimum_scale": 5.0, "maximum_scale": 4.0},
            "minimum_scale' to be less than or equal to 'maximum_scale",
        ),
        (
            "initial_scale",
            {"initial_scale": 12.0, "maximum_scale": 10.0},
            "initial_scale' to lie within the configured scale bounds",
        ),
        (
            "burnin_iteration_count",
            {"burnin_iteration_count": -1},
            "to be nonnegative",
        ),
    ],
)
def test_build_adaptive_tuning_controller_rejects_invalid_policy_values(
    field_name: str,
    kwargs: dict[str, float | int],
    match: str,
) -> None:
    base_kwargs: dict[str, object] = {
        "proposal_name": "normal-random-walk",
        "scale_parameter_name": "step-size",
        "initial_scale": 1.5,
        "target_acceptance_rate": 0.3,
        "burnin_iteration_count": 12,
        "adaptation_window_size": 4,
        "decrease_factor": 0.6,
        "increase_factor": 1.5,
        "minimum_scale": 0.1,
        "maximum_scale": 10.0,
    }
    base_kwargs.update(kwargs)

    with pytest.raises(PhylogeneticsError, match=match):
        build_adaptive_tuning_controller(**base_kwargs)
