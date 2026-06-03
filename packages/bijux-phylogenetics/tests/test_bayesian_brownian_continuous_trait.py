from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.brownian_continuous_trait import (
    build_brownian_continuous_trait_model_definition,
    build_brownian_continuous_trait_proposal_schedule,
    run_brownian_continuous_trait_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.continuous_trait_location_priors import (
    build_fixed_continuous_trait_location_prior,
    build_normal_continuous_trait_location_prior,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    build_exponential_continuous_trait_scalar_prior,
    build_fixed_continuous_trait_scalar_prior,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_brownian_continuous_trait_model_requires_sampled_root_state_prior() -> None:
    with pytest.raises(PhylogeneticsError, match="non-fixed root-state prior"):
        build_brownian_continuous_trait_model_definition(
            root_state_prior=build_fixed_continuous_trait_location_prior(
                fixed_value=0.0
            ),
            sigma_squared_prior=build_exponential_continuous_trait_scalar_prior(
                rate=1.0
            ),
        )


def test_brownian_continuous_trait_model_requires_sampled_sigma_squared_prior() -> None:
    with pytest.raises(PhylogeneticsError, match="non-fixed sigma-squared prior"):
        build_brownian_continuous_trait_model_definition(
            root_state_prior=build_normal_continuous_trait_location_prior(
                mean=0.0,
                standard_deviation=2.0,
            ),
            sigma_squared_prior=build_fixed_continuous_trait_scalar_prior(
                fixed_value=0.5
            ),
        )


def test_brownian_continuous_trait_runner_emits_parameter_summaries() -> None:
    report = run_brownian_continuous_trait_metropolis_hastings(
        tree=_load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk"),
        tip_values={"A": 0.2, "B": 0.5, "C": 1.0, "D": 1.4},
        model_definition=build_brownian_continuous_trait_model_definition(
            root_state_prior=build_normal_continuous_trait_location_prior(
                mean=0.0,
                standard_deviation=2.0,
            ),
            sigma_squared_prior=build_exponential_continuous_trait_scalar_prior(
                rate=1.0
            ),
        ),
        proposal_schedule=build_brownian_continuous_trait_proposal_schedule(
            model_definition=build_brownian_continuous_trait_model_definition(
                root_state_prior=build_normal_continuous_trait_location_prior(
                    mean=0.0,
                    standard_deviation=2.0,
                ),
                sigma_squared_prior=build_exponential_continuous_trait_scalar_prior(
                    rate=1.0
                ),
            ),
            root_state_move_weight=1.0,
            root_state_standard_deviation=0.25,
            sigma_squared_move_weight=1.0,
            sigma_squared_log_scale_standard_deviation=0.25,
        ),
        iteration_count=20,
        sample_every=1,
        seed=5,
    )

    assert report.taxa == ["A", "B", "C", "D"]
    assert len(report.posterior_rows) == 21
    assert all(
        row.topology_id == report.posterior_rows[0].topology_id
        for row in report.posterior_rows
    )
    assert all(math.isfinite(row.log_likelihood) for row in report.posterior_rows)
    assert all(math.isfinite(row.posterior_log_score) for row in report.posterior_rows)
    assert all(row.sigma_squared > 0.0 for row in report.posterior_rows)
    assert {
        step_row.proposal_changed_fields for step_row in report.chain_report.step_rows
    } == {
        ("scalar_parameters.root-state",),
        ("scalar_parameters.sigma-squared",),
    }

    summary_by_parameter = {
        summary.parameter_name: summary for summary in report.parameter_summaries
    }
    assert set(summary_by_parameter) == {"root-state", "sigma-squared"}
    assert summary_by_parameter["root-state"].sample_count == len(report.posterior_rows)
    assert summary_by_parameter["sigma-squared"].sample_count == len(
        report.posterior_rows
    )
    assert math.isfinite(summary_by_parameter["root-state"].posterior_mean)
    assert math.isfinite(summary_by_parameter["sigma-squared"].posterior_mean)
    assert (
        summary_by_parameter["root-state"].hpd_95_lower
        <= summary_by_parameter["root-state"].hpd_95_upper
    )
    assert (
        summary_by_parameter["sigma-squared"].hpd_95_lower
        <= summary_by_parameter["sigma-squared"].hpd_95_upper
    )


def test_brownian_continuous_trait_runner_requires_exact_tip_taxa() -> None:
    model_definition = build_brownian_continuous_trait_model_definition(
        root_state_prior=build_normal_continuous_trait_location_prior(
            mean=0.0,
            standard_deviation=2.0,
        ),
        sigma_squared_prior=build_exponential_continuous_trait_scalar_prior(rate=1.0),
    )
    proposal_schedule = build_brownian_continuous_trait_proposal_schedule(
        model_definition=model_definition,
        root_state_move_weight=1.0,
        root_state_standard_deviation=0.25,
        sigma_squared_move_weight=1.0,
        sigma_squared_log_scale_standard_deviation=0.25,
    )

    with pytest.raises(PhylogeneticsError, match="match the tree tip set exactly"):
        run_brownian_continuous_trait_metropolis_hastings(
            tree=_load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk"),
            tip_values={"A": 0.2, "B": 0.5, "C": 1.0, "X": 1.4},
            model_definition=model_definition,
            proposal_schedule=proposal_schedule,
            iteration_count=4,
            sample_every=1,
            seed=0,
        )


def test_brownian_continuous_trait_runner_rejects_constant_tip_values() -> None:
    model_definition = build_brownian_continuous_trait_model_definition(
        root_state_prior=build_normal_continuous_trait_location_prior(
            mean=0.0,
            standard_deviation=2.0,
        ),
        sigma_squared_prior=build_exponential_continuous_trait_scalar_prior(rate=1.0),
    )
    proposal_schedule = build_brownian_continuous_trait_proposal_schedule(
        model_definition=model_definition,
        root_state_move_weight=1.0,
        root_state_standard_deviation=0.25,
        sigma_squared_move_weight=1.0,
        sigma_squared_log_scale_standard_deviation=0.25,
    )

    with pytest.raises(PhylogeneticsError, match="requires non-constant tip values"):
        run_brownian_continuous_trait_metropolis_hastings(
            tree=_load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk"),
            tip_values={"A": 0.5, "B": 0.5, "C": 0.5, "D": 0.5},
            model_definition=model_definition,
            proposal_schedule=proposal_schedule,
            iteration_count=4,
            sample_every=1,
            seed=0,
        )


def _load_rooted_tree_fixture(name: str):
    tree = load_tree(fixture("trees", name))
    tree.rooted = True
    return tree
