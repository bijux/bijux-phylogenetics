from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.trace_autocorrelation import (
    TraceAutocorrelationLagRow,
    compute_trace_autocorrelation,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_compute_trace_autocorrelation_distinguishes_independent_and_autocorrelated_traces() -> (
    None
):
    independent_lag_one = compute_trace_autocorrelation(
        [1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0],
        lag=1,
    )
    autocorrelated_lag_one = compute_trace_autocorrelation(
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
        lag=1,
    )

    assert isinstance(
        TraceAutocorrelationLagRow(
            parameter_name="x",
            lag=1,
            autocorrelation=autocorrelated_lag_one,
        ),
        TraceAutocorrelationLagRow,
    )
    assert independent_lag_one == 0.0
    assert autocorrelated_lag_one == 0.625


@pytest.mark.parametrize(
    ("values", "lag", "expected_code"),
    [
        ([1.0], 1, "trace_autocorrelation_series_too_short"),
        ([1.0, 2.0, 3.0], 0, "trace_autocorrelation_integer_not_positive"),
        ([1.0, 2.0, 3.0], 3, "trace_autocorrelation_lag_out_of_range"),
    ],
)
def test_compute_trace_autocorrelation_rejects_invalid_inputs(
    values: list[float],
    lag: int,
    expected_code: str,
) -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        compute_trace_autocorrelation(values, lag=lag)

    assert error_info.value.code == expected_code
