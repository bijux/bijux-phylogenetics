from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.burnin_policies import (
    METROPOLIS_HASTINGS_BURNIN_POLICY_NAMES,
    MetropolisHastingsBurninPolicy,
    build_metropolis_hastings_burnin_policy,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_build_metropolis_hastings_burnin_policy_supports_all_goal_policy_names() -> (
    None
):
    assert METROPOLIS_HASTINGS_BURNIN_POLICY_NAMES == (
        "none",
        "fixed-count",
        "fixed-fraction",
        "diagnostic-suggested",
    )

    none_policy = build_metropolis_hastings_burnin_policy(policy_name="none")
    fixed_count_policy = build_metropolis_hastings_burnin_policy(
        policy_name="fixed-count",
        discarded_sample_count=3,
    )
    fixed_fraction_policy = build_metropolis_hastings_burnin_policy(
        policy_name="fixed-fraction",
        discarded_fraction=0.25,
    )
    diagnostic_policy = build_metropolis_hastings_burnin_policy(
        policy_name="diagnostic-suggested",
    )

    assert isinstance(none_policy, MetropolisHastingsBurninPolicy)
    assert fixed_count_policy.discarded_sample_count == 3
    assert fixed_fraction_policy.discarded_fraction == 0.25
    assert diagnostic_policy.mean_shift_threshold == 0.5
    assert diagnostic_policy.minimum_retained_sample_count == 4


@pytest.mark.parametrize(
    ("kwargs", "expected_code"),
    [
        (
            {"policy_name": "fixed-count"},
            "metropolis_hastings_burnin_policy_discarded_sample_count_missing",
        ),
        (
            {"policy_name": "fixed-fraction"},
            "metropolis_hastings_burnin_policy_discarded_fraction_missing",
        ),
        (
            {"policy_name": "none", "discarded_sample_count": 1},
            "metropolis_hastings_burnin_policy_discarded_sample_count_unexpected",
        ),
        (
            {"policy_name": "diagnostic-suggested", "discarded_fraction": 0.1},
            "metropolis_hastings_burnin_policy_discarded_fraction_unexpected",
        ),
        (
            {"policy_name": "unsupported-policy"},
            "metropolis_hastings_burnin_policy_name_invalid",
        ),
    ],
)
def test_build_metropolis_hastings_burnin_policy_rejects_invalid_policy_payloads(
    kwargs: dict[str, object],
    expected_code: str,
) -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        build_metropolis_hastings_burnin_policy(**kwargs)

    assert error_info.value.code == expected_code
