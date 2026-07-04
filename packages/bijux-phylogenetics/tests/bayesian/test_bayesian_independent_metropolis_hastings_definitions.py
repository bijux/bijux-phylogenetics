from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.independent_chains import (
    build_independent_metropolis_hastings_chain_definition,
    build_independent_metropolis_hastings_chain_report,
    build_independent_metropolis_hastings_diagnostics_report,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    build_metropolis_hastings_proposal,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_independent_chain_definition_preserves_chain_name_and_seed() -> None:
    definition = build_independent_metropolis_hastings_chain_definition(
        chain_name=" chain_alpha ",
        seed=13,
    )

    assert definition.chain_name == "chain_alpha"
    assert definition.seed == 13


def test_independent_chain_diagnostics_reject_duplicated_traces() -> None:
    duplicated_chain_reports = [
        build_independent_metropolis_hastings_chain_report(
            chain_name="chain_alpha",
            chain_report=_run_deterministic_chain(seed=3),
        ),
        build_independent_metropolis_hastings_chain_report(
            chain_name="chain_beta",
            chain_report=_run_deterministic_chain(seed=17),
        ),
    ]

    with pytest.raises(
        PhylogeneticsError,
        match="duplicated chain traces",
    ):
        build_independent_metropolis_hastings_diagnostics_report(
            chain_reports=duplicated_chain_reports
        )


def _run_deterministic_chain(seed: int):
    return run_metropolis_hastings_sampler(
        initial_state=_build_scored_normal_target_state(0.0),
        propose_state=_propose_deterministic_zero_shift,
        update_prior_components=_normal_target_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        iteration_count=4,
        sample_every=2,
        seed=seed,
    )


def _propose_deterministic_zero_shift(
    current_state: BayesianPhylogeneticState,
    _rng,
):
    return build_metropolis_hastings_proposal(
        changed_fields=("scalar_parameters.x",),
        log_forward_density=0.0,
        log_reverse_density=0.0,
        is_valid=True,
        proposed_tree=current_state.tree.to_tree(),
        proposed_model_parameters=current_state.model_parameters,
    )


def _build_scored_normal_target_state(x_value: float) -> BayesianPhylogeneticState:
    return score_bayesian_phylogenetic_state(
        tree=PhyloTree(
            TreeNode(
                children=[
                    TreeNode(name="A", branch_length=0.1),
                    TreeNode(name="B", branch_length=0.2),
                ]
            ),
            rooted=True,
        ),
        model_parameters=build_bayesian_model_parameter_state(
            scalar_parameters={"x": x_value}
        ),
        update_prior_components=_normal_target_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )


def _normal_target_prior_components(
    state: BayesianPhylogeneticState,
) -> list[BayesianPriorComponentState]:
    x_value = state.model_parameters.scalar_parameters["x"]
    return [
        BayesianPriorComponentState(
            component_name="normal-target",
            family="gaussian",
            log_prior=-0.5 * (x_value**2),
            parameter_values={"mean": 0.0, "variance": 1.0},
        )
    ]


def _zero_log_likelihood(_state: BayesianPhylogeneticState) -> float:
    return 0.0
