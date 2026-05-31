from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.trace_effective_sample_size import (
    TraceEffectiveSampleSizeRow,
    compute_trace_effective_sample_size,
    compute_trace_integrated_autocorrelation_time,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_compute_trace_effective_sample_size_uses_autocorrelation_time_logic() -> None:
    independent_series = [1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0]
    autocorrelated_series = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

    independent_tau, independent_last_positive_lag = (
        compute_trace_integrated_autocorrelation_time(independent_series)
    )
    autocorrelated_tau, autocorrelated_last_positive_lag = (
        compute_trace_integrated_autocorrelation_time(autocorrelated_series)
    )
    independent_ess = compute_trace_effective_sample_size(independent_series)
    autocorrelated_ess = compute_trace_effective_sample_size(autocorrelated_series)

    assert isinstance(
        TraceEffectiveSampleSizeRow(
            parameter_name="x",
            sample_count=8,
            integrated_autocorrelation_time=autocorrelated_tau,
            effective_sample_size=autocorrelated_ess,
            last_positive_lag=autocorrelated_last_positive_lag,
        ),
        TraceEffectiveSampleSizeRow,
    )
    assert independent_tau == 1.0
    assert independent_last_positive_lag is None
    assert independent_ess == 8.0
    assert autocorrelated_tau == pytest.approx(2.797619047619047, rel=0, abs=1e-12)
    assert autocorrelated_last_positive_lag == 2
    assert autocorrelated_ess == pytest.approx(2.859574468085106, rel=0, abs=1e-12)
    assert autocorrelated_ess < independent_ess


@pytest.mark.parametrize(
    ("values", "maximum_lag", "expected_code"),
    [
        ([], None, "trace_effective_sample_size_series_empty"),
        ([1.0], 1, "trace_effective_sample_size_maximum_lag_singleton_invalid"),
        ([1.0, 2.0, 3.0], 3, "trace_effective_sample_size_maximum_lag_out_of_range"),
    ],
)
def test_compute_trace_effective_sample_size_rejects_invalid_inputs(
    values: list[float],
    maximum_lag: int | None,
    expected_code: str,
) -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        compute_trace_effective_sample_size(values, maximum_lag=maximum_lag)

    assert error_info.value.code == expected_code
