from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.discrete_trait_mk import (
    DISCRETE_TRAIT_MK_MODELS,
    DISCRETE_TRAIT_MK_ROOT_PRIOR_MODES,
    build_discrete_trait_mk_model_definition,
    build_discrete_trait_mk_proposal_schedule,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    build_exponential_discrete_trait_rate_prior,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


@pytest.mark.parametrize(
    ("requested_model_name", "expected_model_name"),
    [
        ("ER", "equal-rates"),
        ("SYM", "symmetric"),
        ("ARD", "all-rates-different"),
    ],
)
def test_discrete_trait_mk_model_definition_resolves_supported_aliases(
    requested_model_name: str,
    expected_model_name: str,
) -> None:
    model_definition = build_discrete_trait_mk_model_definition(
        transition_model_name=requested_model_name,
        rate_prior=build_exponential_discrete_trait_rate_prior(rate=0.75),
    )

    assert model_definition.transition_model_name == expected_model_name
    assert model_definition.root_prior_mode == "equal"
    assert model_definition.initial_rate == 1.0


def test_discrete_trait_mk_model_definition_support_constants_match_contract() -> None:
    assert DISCRETE_TRAIT_MK_MODELS == (
        "equal-rates",
        "symmetric",
        "all-rates-different",
    )
    assert DISCRETE_TRAIT_MK_ROOT_PRIOR_MODES == ("equal", "empirical", "fixed")


def test_discrete_trait_mk_model_definition_requires_fixed_root_state_only_in_fixed_mode() -> (
    None
):
    with pytest.raises(
        PhylogeneticsError,
        match="accepts 'fixed_root_state' only when root_prior_mode is 'fixed'",
    ):
        build_discrete_trait_mk_model_definition(
            transition_model_name="ER",
            rate_prior=build_exponential_discrete_trait_rate_prior(rate=0.75),
            root_prior_mode="equal",
            fixed_root_state="north",
        )

    with pytest.raises(
        PhylogeneticsError,
        match="requires 'fixed_root_state' when root_prior_mode is 'fixed'",
    ):
        build_discrete_trait_mk_model_definition(
            transition_model_name="ER",
            rate_prior=build_exponential_discrete_trait_rate_prior(rate=0.75),
            root_prior_mode="fixed",
        )


def test_discrete_trait_mk_proposal_schedule_records_model_and_rate_kernel() -> None:
    model_definition = build_discrete_trait_mk_model_definition(
        transition_model_name="SYM",
        rate_prior=build_exponential_discrete_trait_rate_prior(rate=1.25),
        root_prior_mode="empirical",
        initial_rate=0.5,
    )

    proposal_schedule = build_discrete_trait_mk_proposal_schedule(
        model_definition=model_definition,
        rate_log_scale_standard_deviation=0.35,
    )

    assert proposal_schedule.transition_model_name == "symmetric"
    assert proposal_schedule.rate_log_scale_standard_deviation == 0.35
