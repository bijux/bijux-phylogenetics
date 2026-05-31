from __future__ import annotations

import math
from pathlib import Path
from random import Random

from bijux_phylogenetics.bayesian.clock_models import (
    build_relaxed_lognormal_clock_model,
    build_strict_clock_rate_model,
    evaluate_relaxed_lognormal_clock_tree_log_prior,
    evaluate_strict_clock_tree_log_prior,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_clock_rate_move,
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
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_clock_rate_proposal_keeps_rates_on_positive_support() -> None:
    current_state = _build_scored_strict_clock_rate_state()
    rng = Random(17)

    for _ in range(100):
        proposal = propose_clock_rate_move(
            current_state,
            rng,
            log_scale_standard_deviation=0.6,
        )
        assert proposal.is_valid is True
        assert proposal.invalid_reason is None
        assert proposal.changed_fields == ("scalar_parameters.clock-rate",)
        assert proposal.proposed_model_parameters is not None
        assert proposal.proposed_model_parameters.scalar_parameters["clock-rate"] > 0.0


def test_clock_rate_proposal_changes_strict_clock_posterior() -> None:
    current_state = _build_scored_strict_clock_rate_state()

    proposal = propose_clock_rate_move(
        current_state,
        Random(2),
        log_scale_standard_deviation=0.6,
    )

    assert proposal.is_valid is True
    assert proposal.invalid_reason is None
    assert proposal.changed_fields == ("scalar_parameters.clock-rate",)
    proposed_tree = proposal.proposed_tree
    proposed_model_parameters = proposal.proposed_model_parameters
    assert proposed_tree is not None
    assert proposed_model_parameters is not None
    assert rooted_topology_fingerprint(proposed_tree) == rooted_topology_fingerprint(
        current_state.tree.to_tree()
    )
    current_clock_rate = current_state.model_parameters.scalar_parameters["clock-rate"]
    proposed_clock_rate = proposed_model_parameters.scalar_parameters["clock-rate"]
    assert proposed_clock_rate > 0.0
    assert math.isclose(
        proposal.log_hastings_ratio,
        math.log(proposed_clock_rate / current_clock_rate),
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    proposed_state = score_bayesian_phylogenetic_state(
        tree=proposed_tree,
        model_parameters=proposed_model_parameters,
        update_prior_components=_strict_clock_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )

    assert not math.isclose(
        proposed_state.posterior_log_score,
        current_state.posterior_log_score,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_clock_rate_proposal_supports_relaxed_clock_rate_parameters() -> None:
    current_state = _build_scored_relaxed_clock_rate_state()

    proposal = propose_clock_rate_move(
        current_state,
        Random(11),
        log_scale_standard_deviation=0.4,
        parameter_name="mean-clock-rate",
    )

    assert proposal.is_valid is True
    assert proposal.invalid_reason is None
    assert proposal.changed_fields == ("scalar_parameters.mean-clock-rate",)
    proposed_model_parameters = proposal.proposed_model_parameters
    assert proposed_model_parameters is not None

    proposed_state = score_bayesian_phylogenetic_state(
        tree=proposal.proposed_tree,
        model_parameters=proposed_model_parameters,
        update_prior_components=_relaxed_clock_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )

    assert not math.isclose(
        proposed_state.posterior_log_score,
        current_state.posterior_log_score,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_sampler_uses_clock_rate_proposal_on_real_dated_tree_surface() -> None:
    initial_state = _build_scored_strict_clock_rate_state()

    report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=lambda current_state, rng: propose_clock_rate_move(
            current_state,
            rng,
            log_scale_standard_deviation=0.6,
        ),
        update_prior_components=_strict_clock_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        iteration_count=6,
        sample_every=1,
        seed=0,
    )

    assert all(step_row.proposal_valid is True for step_row in report.step_rows)
    assert all(
        step_row.proposal_changed_fields == ("scalar_parameters.clock-rate",)
        for step_row in report.step_rows
    )
    assert any(
        step_row.log_acceptance_ratio is not None
        and not math.isclose(
            step_row.log_acceptance_ratio,
            step_row.log_hastings_ratio,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for step_row in report.step_rows
    )
    assert report.final_state.model_parameters.scalar_parameters["clock-rate"] > 0.0


def test_clock_rate_proposal_requires_scalar_parameter() -> None:
    current_state = score_bayesian_phylogenetic_state(
        tree=_load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk"),
        model_parameters=build_bayesian_model_parameter_state(),
        update_prior_components=_flat_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )

    proposal = propose_clock_rate_move(
        current_state,
        Random(0),
        log_scale_standard_deviation=0.5,
    )

    assert proposal.is_valid is False
    assert proposal.invalid_reason == (
        "clock-rate proposal requires one 'clock-rate' scalar parameter"
    )
    assert proposal.proposed_tree is None
    assert proposal.proposed_model_parameters is None


def _build_scored_strict_clock_rate_state() -> BayesianPhylogeneticState:
    return score_bayesian_phylogenetic_state(
        tree=_scale_tree_branch_lengths(
            _load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk"),
            clock_rate=0.5,
        ),
        model_parameters=build_bayesian_model_parameter_state(
            scalar_parameters={"clock-rate": 0.5}
        ),
        update_prior_components=_strict_clock_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )


def _build_scored_relaxed_clock_rate_state() -> BayesianPhylogeneticState:
    return score_bayesian_phylogenetic_state(
        tree=_load_rooted_tree_fixture(
            "relaxed_rate_summary_substitution_tree_4_taxa.nwk"
        ),
        model_parameters=build_bayesian_model_parameter_state(
            scalar_parameters={"mean-clock-rate": 0.2}
        ),
        update_prior_components=_relaxed_clock_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )


def _strict_clock_prior_components(
    state: BayesianPhylogeneticState,
) -> list[BayesianPriorComponentState]:
    report = evaluate_strict_clock_tree_log_prior(
        state.tree.to_tree(),
        _load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk"),
        build_strict_clock_rate_model(
            global_clock_rate=state.model_parameters.scalar_parameters["clock-rate"]
        ),
    )
    return [
        BayesianPriorComponentState(
            component_name="strict-clock-rate",
            family=report.family,
            log_prior=report.total_log_prior,
            parameter_values={"clock-rate": report.global_clock_rate},
        )
    ]


def _relaxed_clock_prior_components(
    state: BayesianPhylogeneticState,
) -> list[BayesianPriorComponentState]:
    report = evaluate_relaxed_lognormal_clock_tree_log_prior(
        state.tree.to_tree(),
        _load_rooted_tree_fixture("relaxed_rate_summary_dated_tree_4_taxa.nwk"),
        build_relaxed_lognormal_clock_model(
            rate_policy="independent",
            mean_clock_rate=state.model_parameters.scalar_parameters["mean-clock-rate"],
            log_standard_deviation=0.5,
        ),
    )
    return [
        BayesianPriorComponentState(
            component_name="relaxed-lognormal-clock-rate",
            family=report.family,
            log_prior=report.total_log_prior,
            parameter_values={"mean-clock-rate": report.mean_clock_rate},
        )
    ]


def _flat_prior_components(
    _state: BayesianPhylogeneticState,
) -> list[BayesianPriorComponentState]:
    return [
        BayesianPriorComponentState(
            component_name="flat-prior",
            family="constant",
            log_prior=0.0,
        )
    ]


def _zero_log_likelihood(_state: BayesianPhylogeneticState) -> float:
    return 0.0


def _load_rooted_tree_fixture(name: str) -> PhyloTree:
    tree = load_tree(fixture("trees", name))
    tree.rooted = True
    return tree


def _scale_tree_branch_lengths(tree: PhyloTree, *, clock_rate: float) -> PhyloTree:
    scaled_tree = tree.copy()
    for _parent, child in scaled_tree.iter_edges():
        child.branch_length = float(child.branch_length or 0.0) * clock_rate
    return scaled_tree
