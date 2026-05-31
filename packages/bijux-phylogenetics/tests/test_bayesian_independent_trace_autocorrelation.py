from __future__ import annotations

from bijux_phylogenetics.bayesian.independent_chains import (
    IndependentMetropolisHastingsDiagnosticsReport,
    IndependentMetropolisHastingsRunReport,
    build_independent_metropolis_hastings_chain_report,
)
from bijux_phylogenetics.bayesian.trace_autocorrelation import (
    IndependentMetropolisHastingsTraceAutocorrelationReport,
    summarize_independent_metropolis_hastings_trace_autocorrelation,
)

from .test_bayesian_burnin_policy_application import _build_chain_report


def test_summarize_independent_metropolis_hastings_trace_autocorrelation_preserves_chain_names() -> (
    None
):
    run_report = IndependentMetropolisHastingsRunReport(
        chain_reports=[
            build_independent_metropolis_hastings_chain_report(
                chain_name="independent",
                chain_report=_build_chain_report(
                    [1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0]
                ),
            ),
            build_independent_metropolis_hastings_chain_report(
                chain_name="autocorrelated",
                chain_report=_build_chain_report(
                    [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
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

    report = summarize_independent_metropolis_hastings_trace_autocorrelation(
        run_report=run_report,
        maximum_lag=1,
    )

    assert isinstance(report, IndependentMetropolisHastingsTraceAutocorrelationReport)
    assert [chain_report.chain_name for chain_report in report.chain_reports] == [
        "independent",
        "autocorrelated",
    ]
    independent_lag_one = (
        report.chain_reports[0]
        .autocorrelation_report.parameter_reports[0]
        .lag_rows[0]
        .autocorrelation
    )
    autocorrelated_lag_one = (
        report.chain_reports[1]
        .autocorrelation_report.parameter_reports[0]
        .lag_rows[0]
        .autocorrelation
    )
    assert independent_lag_one == 0.0
    assert autocorrelated_lag_one == 0.625
