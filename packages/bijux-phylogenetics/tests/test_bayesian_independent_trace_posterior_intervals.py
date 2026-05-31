from __future__ import annotations

from bijux_phylogenetics.bayesian.independent_chains import (
    IndependentMetropolisHastingsDiagnosticsReport,
    IndependentMetropolisHastingsRunReport,
    build_independent_metropolis_hastings_chain_report,
)
from bijux_phylogenetics.bayesian.trace_posterior_intervals import (
    IndependentMetropolisHastingsTracePosteriorIntervalReport,
    summarize_independent_metropolis_hastings_trace_posterior_intervals,
)

from .test_bayesian_burnin_policy_application import _build_chain_report


def test_summarize_independent_metropolis_hastings_trace_posterior_intervals_preserves_chain_names() -> (
    None
):
    run_report = IndependentMetropolisHastingsRunReport(
        chain_reports=[
            build_independent_metropolis_hastings_chain_report(
                chain_name="left",
                chain_report=_build_chain_report(
                    [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 8.0, 9.0, 10.0, 11.0]
                ),
            ),
            build_independent_metropolis_hastings_chain_report(
                chain_name="right",
                chain_report=_build_chain_report(
                    [0.0, 0.0, 0.0, 1.0, 2.0, 8.0, 9.0, 10.0, 11.0, 12.0]
                ),
            ),
        ],
        diagnostics=IndependentMetropolisHastingsDiagnosticsReport(
            chain_count=2,
            mean_acceptance_rate=1.0,
            minimum_acceptance_rate=1.0,
            maximum_acceptance_rate=1.0,
            comparison_rows=[],
            warnings=[],
        ),
    )

    report = summarize_independent_metropolis_hastings_trace_posterior_intervals(
        run_report=run_report,
        mass_fraction=0.6,
    )

    assert isinstance(report, IndependentMetropolisHastingsTracePosteriorIntervalReport)
    assert [chain_report.chain_name for chain_report in report.chain_reports] == [
        "left",
        "right",
    ]
    assert (
        report.chain_reports[0]
        .posterior_interval_report.parameter_rows[0]
        .hpd_upper_bound
        == 2.0
    )
