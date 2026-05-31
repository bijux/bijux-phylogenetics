from __future__ import annotations

from bijux_phylogenetics.bayesian.burnin_policies import (
    MetropolisHastingsBurninReport,
    apply_metropolis_hastings_burnin_policy,
    build_metropolis_hastings_burnin_policy,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsRunReport,
    MetropolisHastingsStepRow,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


def test_apply_metropolis_hastings_burnin_policy_changes_retained_sample_count() -> (
    None
):
    chain_report = _build_chain_report([8.0, 6.0, 4.0, 2.0, 1.0, 1.0, 1.0, 1.0])

    none_report = apply_metropolis_hastings_burnin_policy(
        chain_report=chain_report,
        policy=build_metropolis_hastings_burnin_policy(policy_name="none"),
    )
    fixed_count_report = apply_metropolis_hastings_burnin_policy(
        chain_report=chain_report,
        policy=build_metropolis_hastings_burnin_policy(
            policy_name="fixed-count",
            discarded_sample_count=1,
        ),
    )
    fixed_fraction_report = apply_metropolis_hastings_burnin_policy(
        chain_report=chain_report,
        policy=build_metropolis_hastings_burnin_policy(
            policy_name="fixed-fraction",
            discarded_fraction=0.5,
        ),
    )

    assert isinstance(none_report, MetropolisHastingsBurninReport)
    assert none_report.retained_sample_count == 8
    assert fixed_count_report.retained_sample_count == 7
    assert fixed_fraction_report.retained_sample_count == 4
    assert [row.iteration_index for row in fixed_count_report.discarded_rows] == [0]
    assert [row.iteration_index for row in fixed_fraction_report.retained_rows] == [
        4,
        5,
        6,
        7,
    ]


def _build_chain_report(sampled_x_values: list[float]) -> MetropolisHastingsRunReport:
    sampled_states = [_build_state(x_value) for x_value in sampled_x_values]
    iteration_count = len(sampled_states) - 1
    step_rows = [
        MetropolisHastingsStepRow(
            iteration_index=iteration_index,
            proposal_changed_fields=("scalar_parameters.x",),
            proposal_valid=True,
            proposal_invalid_reason=None,
            log_forward_density=0.0,
            log_reverse_density=0.0,
            accepted=True,
            log_hastings_ratio=0.0,
            current_posterior_log_score=sampled_states[
                iteration_index - 1
            ].posterior_log_score,
            proposed_posterior_log_score=sampled_states[
                iteration_index
            ].posterior_log_score,
            log_acceptance_ratio=0.0,
            recorded_posterior_log_score=sampled_states[
                iteration_index
            ].posterior_log_score,
        )
        for iteration_index in range(1, len(sampled_states))
    ]
    return MetropolisHastingsRunReport(
        iteration_count=iteration_count,
        sample_every=1,
        seed=11,
        accepted_count=iteration_count,
        rejected_count=0,
        acceptance_rate=1.0,
        initial_state=sampled_states[0],
        final_state=sampled_states[-1],
        sampled_states=sampled_states,
        step_rows=step_rows,
    )


def _build_state(x_value: float) -> BayesianPhylogeneticState:
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
        update_prior_components=_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )


def _prior_components(
    state: BayesianPhylogeneticState,
) -> list[BayesianPriorComponentState]:
    x_value = state.model_parameters.scalar_parameters["x"]
    return [
        BayesianPriorComponentState(
            component_name="quadratic-x",
            family="gaussian",
            log_prior=-(x_value**2),
            parameter_values={"mean": 0.0, "variance": 0.5},
        )
    ]


def _zero_log_likelihood(_state: BayesianPhylogeneticState) -> float:
    return 0.0
