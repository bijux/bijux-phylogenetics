from __future__ import annotations

from functools import cache
import math
from pathlib import Path
from random import Random

import numpy
import pytest

from bijux_phylogenetics.ancestral.common import load_discrete_dataset
from bijux_phylogenetics.bayesian.discrete_trait_rate_parameters import (
    parameterize_discrete_trait_rate_rows,
    resolve_discrete_trait_rate_rows,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_discrete_trait_rate_move,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
)
from bijux_phylogenetics.comparative.discrete_mk import fit_discrete_mk_model
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_phytools_comparative_fixture,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood.pruning import (
    build_transition_matrix_evaluator,
    postorder_conditional_likelihoods,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


@pytest.mark.parametrize(
    ("surface_key", "expected_parameter_names"),
    [
        ("equal-rates", {"shared-rate"}),
        (
            "symmetric",
            {"north<->south", "north<->west", "south<->west"},
        ),
        ("all-rates-different", {"0->1", "1->0"}),
    ],
)
@pytest.mark.slow
def test_discrete_trait_rate_proposal_changes_real_mk_likelihood(
    surface_key: str,
    expected_parameter_names: set[str],
) -> None:
    surface = _load_surface(surface_key)
    current_state = _build_scored_discrete_trait_rate_state(surface)

    proposal = propose_discrete_trait_rate_move(
        current_state,
        Random(2),
        log_scale_standard_deviation=0.6,
    )

    assert proposal.is_valid is True
    assert proposal.invalid_reason is None
    assert proposal.proposed_tree is not None
    proposed_model_parameters = proposal.proposed_model_parameters
    assert proposed_model_parameters is not None
    assert (
        set(proposed_model_parameters.vector_parameters["discrete-trait-rates"])
        == expected_parameter_names
    )
    assert proposal.changed_fields[0].startswith(
        "vector_parameters.discrete-trait-rates."
    )

    proposed_state = score_bayesian_phylogenetic_state(
        tree=proposal.proposed_tree,
        model_parameters=proposed_model_parameters,
        update_prior_components=_zero_prior_components,
        update_log_likelihood=lambda state: _discrete_mk_log_likelihood(state, surface),
    )

    assert not math.isclose(
        proposed_state.log_likelihood,
        current_state.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_discrete_trait_rate_proposal_keeps_transition_rates_positive() -> None:
    current_state = _build_scored_discrete_trait_rate_state(_load_ard_surface())
    rng = Random(17)

    for _ in range(100):
        proposal = propose_discrete_trait_rate_move(
            current_state,
            rng,
            log_scale_standard_deviation=0.6,
        )
        assert proposal.is_valid is True
        assert proposal.invalid_reason is None
        assert proposal.proposed_model_parameters is not None
        assert all(
            value > 0.0
            for value in proposal.proposed_model_parameters.vector_parameters[
                "discrete-trait-rates"
            ].values()
        )


def test_sampler_uses_discrete_trait_rate_proposal_on_mk_surface() -> None:
    surface = _load_symmetric_surface()
    initial_state = _build_scored_discrete_trait_rate_state(surface)

    report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=lambda current_state, rng: propose_discrete_trait_rate_move(
            current_state,
            rng,
            log_scale_standard_deviation=0.6,
        ),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=lambda state: _discrete_mk_log_likelihood(state, surface),
        iteration_count=6,
        sample_every=1,
        seed=0,
    )

    assert all(step_row.proposal_valid is True for step_row in report.step_rows)
    assert all(
        step_row.proposal_changed_fields[0].startswith(
            "vector_parameters.discrete-trait-rates."
        )
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


def test_discrete_trait_rate_proposal_requires_rate_vector() -> None:
    current_state = score_bayesian_phylogenetic_state(
        tree=load_tree(
            fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk")
        ),
        model_parameters=build_bayesian_model_parameter_state(),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )

    proposal = propose_discrete_trait_rate_move(
        current_state,
        Random(0),
        log_scale_standard_deviation=0.5,
    )

    assert proposal.is_valid is False
    assert proposal.invalid_reason == (
        "discrete-trait rate proposal requires one 'discrete-trait-rates' vector parameter"
    )
    assert proposal.proposed_tree is None
    assert proposal.proposed_model_parameters is None


class _DiscreteMkSurface:
    def __init__(
        self,
        *,
        fit_report,
        tree_path: Path,
        traits_path: Path,
        trait: str,
        taxon_column: str,
    ) -> None:
        self.fit_report = fit_report
        self.tree_path = tree_path
        self.traits_path = traits_path
        self.trait = trait
        self.taxon_column = taxon_column
        self.dataset = load_discrete_dataset(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
        )
        self.parameterization = parameterize_discrete_trait_rate_rows(
            model=fit_report.model,
            transition_rate_rows=fit_report.transition_rate_rows,
        )


def _build_scored_discrete_trait_rate_state(
    surface: _DiscreteMkSurface,
) -> BayesianPhylogeneticState:
    return score_bayesian_phylogenetic_state(
        tree=load_tree(surface.tree_path),
        model_parameters=build_bayesian_model_parameter_state(
            vector_parameters={
                "discrete-trait-rates": surface.parameterization.parameter_values
            }
        ),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=lambda state: _discrete_mk_log_likelihood(state, surface),
    )


def _discrete_mk_log_likelihood(
    state: BayesianPhylogeneticState,
    surface: _DiscreteMkSurface,
) -> float:
    transition_rate_rows = resolve_discrete_trait_rate_rows(
        model=surface.fit_report.model,
        transition_rate_rows=surface.fit_report.transition_rate_rows,
        parameter_values=state.model_parameters.vector_parameters[
            "discrete-trait-rates"
        ],
    )
    rate_matrix = _build_rate_matrix_from_transition_rows(
        state_order=surface.fit_report.state_order,
        transition_rate_rows=transition_rate_rows,
    )
    tree = load_tree(surface.tree_path)
    state_index = {
        state: index for index, state in enumerate(surface.fit_report.state_order)
    }
    transition_evaluator = build_transition_matrix_evaluator(rate_matrix)
    pruning_pass = postorder_conditional_likelihoods(
        tree,
        state_count=len(surface.fit_report.state_order),
        leaf_likelihood=lambda node: _leaf_likelihood_vector(
            states_by_taxon=surface.dataset.states_by_taxon,
            state_index=state_index,
            state_count=len(surface.fit_report.state_order),
            node_name=node.name,
        ),
        transition_matrix_for_child=lambda child: (
            transition_evaluator.transition_probability_matrix(
                _branch_length_for_mk_likelihood(child.branch_length)
            )
        ),
    )
    root_partial = pruning_pass.conditional_for_node(tree.root)
    observed_root_total = float(root_partial.sum())
    if observed_root_total <= 0.0:
        return float("-inf")
    root_scale = float((root_partial @ root_partial) / observed_root_total)
    if root_scale <= 0.0:
        return float("-inf")
    subtree_log_scale = pruning_pass.subtree_log_scaling_for_node(tree.root)
    if math.isinf(subtree_log_scale) and subtree_log_scale < 0.0:
        return float("-inf")
    return subtree_log_scale + math.log(root_scale)


def _build_rate_matrix_from_transition_rows(
    *,
    state_order: list[str],
    transition_rate_rows,
):
    rate_matrix = numpy.zeros((len(state_order), len(state_order)), dtype=float)
    state_index = {state: index for index, state in enumerate(state_order)}
    for row in transition_rate_rows:
        if row.transition_allowed:
            rate_matrix[
                state_index[row.source_state],
                state_index[row.target_state],
            ] = row.rate
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(rate_matrix[row_index, :].sum())
    return rate_matrix


def _zero_prior_components(
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


@cache
def _load_surface(model: str) -> _DiscreteMkSurface:
    if model == "equal-rates":
        return _load_equal_rates_surface()
    if model == "symmetric":
        return _load_symmetric_surface()
    if model == "all-rates-different":
        return _load_ard_surface()
    raise ValueError(model)


@cache
def _load_equal_rates_surface() -> _DiscreteMkSurface:
    tree_path = fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk")
    traits_path = fixture("example_traits_phytools_signal_twenty_four_taxa.tsv")
    return _DiscreteMkSurface(
        fit_report=fit_discrete_mk_model(
            tree_path,
            traits_path,
            trait="binary_state",
            taxon_column="taxon",
            model="equal-rates",
        ),
        tree_path=tree_path,
        traits_path=traits_path,
        trait="binary_state",
        taxon_column="taxon",
    )


@cache
def _load_symmetric_surface() -> _DiscreteMkSurface:
    tree_path = fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk")
    traits_path = fixture("example_traits_phytools_signal_twenty_four_taxa.tsv")
    return _DiscreteMkSurface(
        fit_report=fit_discrete_mk_model(
            tree_path,
            traits_path,
            trait="region_state",
            taxon_column="taxon",
            model="symmetric",
        ),
        tree_path=tree_path,
        traits_path=traits_path,
        trait="region_state",
        taxon_column="taxon",
    )


@cache
def _load_ard_surface() -> _DiscreteMkSurface:
    fixture_entry = get_shared_phytools_comparative_fixture(
        "phytools_discrete_ard_binary_twenty_four_taxa"
    )
    return _DiscreteMkSurface(
        fit_report=fit_discrete_mk_model(
            fixture_entry.tree_path,
            fixture_entry.traits_path,
            trait=fixture_entry.trait_name,
            taxon_column=fixture_entry.taxon_column,
            model="all-rates-different",
        ),
        tree_path=fixture_entry.tree_path,
        traits_path=fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
    )


def _branch_length_for_mk_likelihood(branch_length: float | None) -> float:
    if branch_length is None:
        return 1.0
    return max(float(branch_length), 0.0)


def _leaf_likelihood_vector(
    *,
    states_by_taxon: dict[str, str],
    state_index: dict[str, int],
    state_count: int,
    node_name: str | None,
) -> numpy.ndarray:
    if node_name is None:
        raise ValueError("discrete Mk likelihood test trees require named leaves")
    likelihood = numpy.zeros(state_count, dtype=float)
    likelihood[state_index[states_by_taxon[node_name]]] = 1.0
    return likelihood
