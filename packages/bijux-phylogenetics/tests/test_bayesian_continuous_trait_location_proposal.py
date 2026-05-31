from __future__ import annotations

import math
from pathlib import Path
from random import Random

from bijux_phylogenetics.bayesian.continuous_trait_location_priors import (
    build_normal_continuous_trait_location_prior,
    evaluate_continuous_trait_location_log_prior,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_continuous_trait_location_move,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.clades import rooted_topology_fingerprint

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_continuous_trait_location_proposal_keeps_finite_support() -> None:
    current_state = _build_scored_root_state()
    rng = Random(17)

    for _ in range(100):
        proposal = propose_continuous_trait_location_move(
            current_state,
            rng,
            standard_deviation=0.4,
        )
        assert proposal.is_valid is True
        assert proposal.invalid_reason is None
        assert proposal.changed_fields == ("scalar_parameters.root-state",)
        assert proposal.proposed_model_parameters is not None
        assert math.isfinite(
            proposal.proposed_model_parameters.scalar_parameters["root-state"]
        )


def test_continuous_trait_location_proposal_changes_real_valued_posterior_surface() -> (
    None
):
    current_state = _build_scored_root_state()

    proposal = propose_continuous_trait_location_move(
        current_state,
        Random(2),
        standard_deviation=0.4,
    )

    assert proposal.is_valid is True
    assert proposal.invalid_reason is None
    assert proposal.changed_fields == ("scalar_parameters.root-state",)
    proposed_tree = proposal.proposed_tree
    proposed_model_parameters = proposal.proposed_model_parameters
    assert proposed_tree is not None
    assert proposed_model_parameters is not None
    assert rooted_topology_fingerprint(proposed_tree) == rooted_topology_fingerprint(
        current_state.tree.to_tree()
    )
    assert proposal.log_hastings_ratio == 0.0
    proposed_state = score_bayesian_phylogenetic_state(
        tree=proposed_tree,
        model_parameters=proposed_model_parameters,
        update_prior_components=_root_state_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )
    assert not math.isclose(
        proposed_state.posterior_log_score,
        current_state.posterior_log_score,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_sampler_uses_continuous_trait_location_proposal() -> None:
    initial_state = _build_scored_root_state()

    report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=lambda current_state, rng: propose_continuous_trait_location_move(
            current_state,
            rng,
            standard_deviation=0.4,
        ),
        update_prior_components=_root_state_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        iteration_count=6,
        sample_every=1,
        seed=0,
    )

    assert all(step_row.proposal_valid is True for step_row in report.step_rows)
    assert all(
        step_row.proposal_changed_fields == ("scalar_parameters.root-state",)
        for step_row in report.step_rows
    )
    assert any(
        not math.isclose(
            state.model_parameters.scalar_parameters["root-state"],
            initial_state.model_parameters.scalar_parameters["root-state"],
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for state in report.sampled_states
    )


def test_continuous_trait_location_proposal_requires_scalar_parameter() -> None:
    current_state = score_bayesian_phylogenetic_state(
        tree=_load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk"),
        model_parameters=build_bayesian_model_parameter_state(),
        update_prior_components=_flat_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )

    proposal = propose_continuous_trait_location_move(
        current_state,
        Random(0),
        standard_deviation=0.5,
    )

    assert proposal.is_valid is False
    assert proposal.invalid_reason == (
        "continuous-trait location proposal requires one 'root-state' scalar parameter"
    )
    assert proposal.proposed_tree is None
    assert proposal.proposed_model_parameters is None


def _build_scored_root_state() -> BayesianPhylogeneticState:
    return score_bayesian_phylogenetic_state(
        tree=_load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk"),
        model_parameters=build_bayesian_model_parameter_state(
            scalar_parameters={"root-state": 0.5}
        ),
        update_prior_components=_root_state_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )


def _root_state_prior_components(
    state: BayesianPhylogeneticState,
) -> list[BayesianPriorComponentState]:
    log_prior = evaluate_continuous_trait_location_log_prior(
        parameter_value=state.model_parameters.scalar_parameters["root-state"],
        prior_model=build_normal_continuous_trait_location_prior(
            mean=0.0,
            standard_deviation=1.0,
        ),
        parameter_name="root-state",
    )
    return [
        BayesianPriorComponentState(
            component_name="continuous-trait:root-state",
            family="normal",
            log_prior=log_prior,
        )
    ]


def _flat_prior_components(
    state: BayesianPhylogeneticState,
) -> list[BayesianPriorComponentState]:
    return [
        BayesianPriorComponentState(
            component_name="flat",
            family=None,
            log_prior=0.0,
        )
    ]


def _zero_log_likelihood(_state: BayesianPhylogeneticState) -> float:
    return 0.0


def _load_rooted_tree_fixture(name: str):
    tree = load_tree(fixture("trees", name))
    tree.rooted = True
    return tree
