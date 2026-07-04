from __future__ import annotations

from bijux_phylogenetics.bayesian.trace_posterior_intervals import (
    MetropolisHastingsTracePosteriorIntervalReport,
    summarize_metropolis_hastings_trace_posterior_intervals,
)

from .test_bayesian_burnin_policy_application import _build_chain_report


def test_summarize_metropolis_hastings_trace_posterior_intervals_reports_per_parameter_hpd() -> (
    None
):
    report = summarize_metropolis_hastings_trace_posterior_intervals(
        chain_report=_build_chain_report(
            [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 8.0, 9.0, 10.0, 11.0]
        ),
        mass_fraction=0.6,
    )

    assert isinstance(report, MetropolisHastingsTracePosteriorIntervalReport)
    row = report.parameter_rows[0]
    assert row.parameter_name == "x"
    assert row.sample_count == 10
    assert row.mass_fraction == 0.6
    assert (row.hpd_lower_bound, row.hpd_upper_bound) == (0.0, 2.0)
    assert (row.equal_tail_lower_bound, row.equal_tail_upper_bound) == (0.0, 9.2)
