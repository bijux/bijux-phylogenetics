from __future__ import annotations

from bijux_phylogenetics.bayesian.trace_effective_sample_size import (
    MetropolisHastingsTraceEffectiveSampleSizeReport,
    summarize_metropolis_hastings_trace_effective_sample_size,
)

from .test_bayesian_burnin_policy_application import _build_chain_report


def test_summarize_metropolis_hastings_trace_effective_sample_size_reports_per_parameter_ess() -> (
    None
):
    independent_report = summarize_metropolis_hastings_trace_effective_sample_size(
        chain_report=_build_chain_report([1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0]),
    )
    autocorrelated_report = summarize_metropolis_hastings_trace_effective_sample_size(
        chain_report=_build_chain_report([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]),
    )

    assert isinstance(
        independent_report, MetropolisHastingsTraceEffectiveSampleSizeReport
    )
    independent_x = independent_report.parameter_rows[0]
    autocorrelated_x = autocorrelated_report.parameter_rows[0]

    assert independent_x.parameter_name == "x"
    assert independent_x.sample_count == 8
    assert independent_x.effective_sample_size == 8.0
    assert independent_x.integrated_autocorrelation_time == 1.0
    assert autocorrelated_x.last_positive_lag == 2
    assert autocorrelated_x.effective_sample_size < independent_x.effective_sample_size
