from __future__ import annotations

from bijux_phylogenetics.bayesian.burnin_policies import (
    apply_metropolis_hastings_burnin_policy,
    build_metropolis_hastings_burnin_policy,
    diagnose_metropolis_hastings_burnin,
)

from .test_bayesian_burnin_policy_application import _build_chain_report


def test_diagnostic_burnin_policy_selects_first_stable_tail() -> None:
    chain_report = _build_chain_report([8.0, 6.0, 4.0, 2.0, 1.0, 1.0, 1.0, 1.0])

    diagnostic_report = diagnose_metropolis_hastings_burnin(
        chain_report=chain_report,
        mean_shift_threshold=1.25,
        minimum_retained_sample_count=4,
    )

    assert diagnostic_report.selected_discarded_sample_count == 2
    assert diagnostic_report.stabilized_tail_found is True
    assert [row.acceptable for row in diagnostic_report.candidate_rows] == [
        False,
        False,
        True,
        True,
        True,
    ]


def test_apply_diagnostic_burnin_policy_records_diagnostic_ledger() -> None:
    chain_report = _build_chain_report([8.0, 6.0, 4.0, 2.0, 1.0, 1.0, 1.0, 1.0])

    burnin_report = apply_metropolis_hastings_burnin_policy(
        chain_report=chain_report,
        policy=build_metropolis_hastings_burnin_policy(
            policy_name="diagnostic-suggested",
            mean_shift_threshold=1.25,
            minimum_retained_sample_count=4,
        ),
    )

    assert burnin_report.discarded_sample_count == 2
    assert burnin_report.retained_sample_count == 6
    assert burnin_report.diagnostic_report is not None
    assert (
        burnin_report.diagnostic_report.selected_discarded_sample_count
        == burnin_report.discarded_sample_count
    )
