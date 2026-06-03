from __future__ import annotations

from bijux_phylogenetics.bayesian.burnin_policies import (
    IndependentMetropolisHastingsBurninReport,
    apply_independent_metropolis_hastings_burnin_policy,
    build_metropolis_hastings_burnin_policy,
)
from bijux_phylogenetics.bayesian.independent_chains import (
    IndependentMetropolisHastingsDiagnosticsReport,
    IndependentMetropolisHastingsRunReport,
    build_independent_metropolis_hastings_chain_report,
)

from .test_bayesian_burnin_policy_application import _build_chain_report


def test_apply_independent_metropolis_hastings_burnin_policy_keeps_named_chain_reports() -> (
    None
):
    run_report = IndependentMetropolisHastingsRunReport(
        chain_reports=[
            build_independent_metropolis_hastings_chain_report(
                chain_name="alpha",
                chain_report=_build_chain_report([8.0, 6.0, 4.0, 2.0, 1.0, 1.0]),
            ),
            build_independent_metropolis_hastings_chain_report(
                chain_name="beta",
                chain_report=_build_chain_report([7.0, 5.0, 3.0, 2.0, 2.0, 2.0]),
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

    burnin_report = apply_independent_metropolis_hastings_burnin_policy(
        run_report=run_report,
        policy=build_metropolis_hastings_burnin_policy(
            policy_name="fixed-count",
            discarded_sample_count=2,
        ),
    )

    assert isinstance(burnin_report, IndependentMetropolisHastingsBurninReport)
    assert [row.chain_name for row in burnin_report.chain_reports] == ["alpha", "beta"]
    assert [
        row.burnin_report.retained_sample_count for row in burnin_report.chain_reports
    ] == [
        4,
        4,
    ]
    assert burnin_report.minimum_retained_sample_count == 4
    assert burnin_report.maximum_retained_sample_count == 4
