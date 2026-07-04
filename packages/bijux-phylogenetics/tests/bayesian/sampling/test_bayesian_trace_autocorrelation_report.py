from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.trace_autocorrelation import (
    MetropolisHastingsTraceAutocorrelationReport,
    summarize_metropolis_hastings_trace_autocorrelation,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from .test_bayesian_burnin_policy_application import _build_chain_report


def test_summarize_metropolis_hastings_trace_autocorrelation_reports_scalar_parameter_lags() -> (
    None
):
    independent_report = summarize_metropolis_hastings_trace_autocorrelation(
        chain_report=_build_chain_report([1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0]),
        maximum_lag=2,
    )
    autocorrelated_report = summarize_metropolis_hastings_trace_autocorrelation(
        chain_report=_build_chain_report([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]),
        maximum_lag=2,
    )

    assert isinstance(independent_report, MetropolisHastingsTraceAutocorrelationReport)
    independent_x = independent_report.parameter_reports[0]
    autocorrelated_x = autocorrelated_report.parameter_reports[0]

    assert independent_x.parameter_name == "x"
    assert [row.lag for row in independent_x.lag_rows] == [1, 2]
    assert independent_x.lag_rows[0].autocorrelation == 0.0
    assert autocorrelated_x.lag_rows[0].autocorrelation == 0.625
    assert (
        autocorrelated_x.lag_rows[0].autocorrelation
        > independent_x.lag_rows[0].autocorrelation
    )


def test_summarize_metropolis_hastings_trace_autocorrelation_rejects_out_of_range_maximum_lag() -> (
    None
):
    with pytest.raises(PhylogeneticsError) as error_info:
        summarize_metropolis_hastings_trace_autocorrelation(
            chain_report=_build_chain_report([1.0, 2.0, 3.0, 4.0]),
            maximum_lag=4,
        )

    assert error_info.value.code == "trace_autocorrelation_maximum_lag_out_of_range"
