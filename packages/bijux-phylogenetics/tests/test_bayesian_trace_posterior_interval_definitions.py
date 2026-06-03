from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.trace_posterior_intervals import (
    HighestPosteriorDensityInterval,
    TracePosteriorIntervalRow,
    compute_equal_tail_interval,
    compute_highest_posterior_density_interval,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_compute_highest_posterior_density_interval_differs_from_equal_tail_on_skewed_fixture() -> (
    None
):
    values = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 8.0, 9.0, 10.0, 11.0]

    hpd_interval = compute_highest_posterior_density_interval(
        values,
        mass_fraction=0.6,
    )
    equal_tail_interval = compute_equal_tail_interval(
        values,
        mass_fraction=0.6,
    )

    assert isinstance(hpd_interval, HighestPosteriorDensityInterval)
    assert isinstance(
        TracePosteriorIntervalRow(
            parameter_name="x",
            sample_count=len(values),
            mass_fraction=0.6,
            hpd_lower_bound=hpd_interval.lower_bound,
            hpd_upper_bound=hpd_interval.upper_bound,
            equal_tail_lower_bound=equal_tail_interval[0],
            equal_tail_upper_bound=equal_tail_interval[1],
        ),
        TracePosteriorIntervalRow,
    )
    assert (hpd_interval.lower_bound, hpd_interval.upper_bound) == (0.0, 2.0)
    assert equal_tail_interval == pytest.approx((0.0, 9.2), rel=0, abs=1e-12)
    assert (hpd_interval.lower_bound, hpd_interval.upper_bound) != equal_tail_interval


@pytest.mark.parametrize(
    ("values", "mass_fraction", "expected_code"),
    [
        ([], 0.95, "trace_posterior_interval_series_empty"),
        ([1.0, 2.0], 0.0, "trace_posterior_interval_mass_fraction_out_of_range"),
        ([1.0, 2.0], 1.1, "trace_posterior_interval_mass_fraction_out_of_range"),
    ],
)
def test_compute_highest_posterior_density_interval_rejects_invalid_inputs(
    values: list[float],
    mass_fraction: float,
    expected_code: str,
) -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        compute_highest_posterior_density_interval(
            values,
            mass_fraction=mass_fraction,
        )

    assert error_info.value.code == expected_code
