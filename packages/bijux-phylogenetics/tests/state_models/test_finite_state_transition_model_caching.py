from __future__ import annotations

from pathlib import Path

import numpy

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    FiniteStateLikelihoodSurface,
    clear_shared_finite_state_transition_cache,
    compile_finite_state_likelihood_context,
    evaluate_finite_state_tree_log_likelihood_from_patterns,
    shared_finite_state_transition_cache_entry_count,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    AlignmentSitePattern,
    CompressedAlignmentSitePatterns,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    return FIXTURES / "trees" / name


def test_compiled_context_reuses_shared_transition_cache_for_identical_model_state() -> (
    None
):
    clear_shared_finite_state_transition_cache()
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    evaluation_count = 0
    rate_matrix = numpy.array([[-2.0, 2.0], [1.0, -1.0]], dtype=float)

    def transition(branch_length: float) -> numpy.ndarray:
        nonlocal evaluation_count
        evaluation_count += 1
        return numpy.array(
            [
                [1.0 - branch_length, branch_length],
                [branch_length / 2.0, 1.0 - (branch_length / 2.0)],
            ],
            dtype=float,
        )

    surface = FiniteStateLikelihoodSurface(
        model_name="binary synthetic",
        state_labels=("0", "1"),
        root_prior=numpy.array([0.5, 0.5], dtype=float),
        root_prior_source="provided",
        observation_policy="reject",
        transition_matrix_for_branch_length=transition,
        leaf_likelihood_vector_for_state=_binary_leaf_vector,
        transition_rate_matrix=rate_matrix,
        parameter_values={"exchangeability": 1.0},
    )

    first_context = compile_finite_state_likelihood_context(tree, surface)
    second_context = compile_finite_state_likelihood_context(tree, surface)

    assert evaluation_count == 1
    assert shared_finite_state_transition_cache_entry_count() == 1
    assert first_context.transition_model_state == second_context.transition_model_state
    assert (
        first_context.transition_cache_key_by_node_id
        == second_context.transition_cache_key_by_node_id
    )


def test_compiled_context_separates_shared_transition_cache_when_parameters_change() -> (
    None
):
    clear_shared_finite_state_transition_cache()
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    evaluation_count = 0
    rate_matrix = numpy.array([[-2.0, 2.0], [1.0, -1.0]], dtype=float)

    def transition(branch_length: float) -> numpy.ndarray:
        nonlocal evaluation_count
        evaluation_count += 1
        return numpy.array(
            [
                [1.0 - branch_length, branch_length],
                [branch_length / 2.0, 1.0 - (branch_length / 2.0)],
            ],
            dtype=float,
        )

    first_context = compile_finite_state_likelihood_context(
        tree,
        FiniteStateLikelihoodSurface(
            model_name="binary synthetic",
            state_labels=("0", "1"),
            root_prior=numpy.array([0.5, 0.5], dtype=float),
            root_prior_source="provided",
            observation_policy="reject",
            transition_matrix_for_branch_length=transition,
            leaf_likelihood_vector_for_state=_binary_leaf_vector,
            transition_rate_matrix=rate_matrix,
            parameter_values={"kappa": 2.0},
        ),
    )
    second_context = compile_finite_state_likelihood_context(
        tree,
        FiniteStateLikelihoodSurface(
            model_name="binary synthetic",
            state_labels=("0", "1"),
            root_prior=numpy.array([0.5, 0.5], dtype=float),
            root_prior_source="provided",
            observation_policy="reject",
            transition_matrix_for_branch_length=transition,
            leaf_likelihood_vector_for_state=_binary_leaf_vector,
            transition_rate_matrix=rate_matrix,
            parameter_values={"kappa": 3.0},
        ),
    )

    assert evaluation_count == 2
    assert shared_finite_state_transition_cache_entry_count() == 2
    assert first_context.transition_model_state != second_context.transition_model_state


def test_compiled_context_separates_shared_transition_cache_when_rate_matrix_or_root_policy_changes() -> (
    None
):
    clear_shared_finite_state_transition_cache()
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    evaluation_count = 0
    first_rate_matrix = numpy.array([[-2.0, 2.0], [1.0, -1.0]], dtype=float)
    second_rate_matrix = numpy.array([[-1.0, 1.0], [0.25, -0.25]], dtype=float)

    def transition(branch_length: float) -> numpy.ndarray:
        nonlocal evaluation_count
        evaluation_count += 1
        return numpy.array(
            [
                [1.0 - branch_length, branch_length],
                [branch_length / 2.0, 1.0 - (branch_length / 2.0)],
            ],
            dtype=float,
        )

    first_context = compile_finite_state_likelihood_context(
        tree,
        FiniteStateLikelihoodSurface(
            model_name="binary synthetic",
            state_labels=("0", "1"),
            root_prior=numpy.array([0.5, 0.5], dtype=float),
            root_prior_source="provided",
            observation_policy="reject",
            transition_matrix_for_branch_length=transition,
            leaf_likelihood_vector_for_state=_binary_leaf_vector,
            transition_rate_matrix=first_rate_matrix,
        ),
    )
    second_context = compile_finite_state_likelihood_context(
        tree,
        FiniteStateLikelihoodSurface(
            model_name="binary synthetic",
            state_labels=("0", "1"),
            root_prior=numpy.array([0.7, 0.3], dtype=float),
            root_prior_source="stationary",
            observation_policy="reject",
            transition_matrix_for_branch_length=transition,
            leaf_likelihood_vector_for_state=_binary_leaf_vector,
            transition_rate_matrix=second_rate_matrix,
        ),
    )

    assert evaluation_count == 2
    assert shared_finite_state_transition_cache_entry_count() == 2
    assert first_context.transition_model_state != second_context.transition_model_state


def test_shared_transition_cache_matches_uncached_context_likelihood() -> None:
    clear_shared_finite_state_transition_cache()
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    rate_matrix = numpy.array([[-2.0, 2.0], [1.0, -1.0]], dtype=float)

    def transition(branch_length: float) -> numpy.ndarray:
        return numpy.array(
            [
                [1.0 - branch_length, branch_length],
                [branch_length / 2.0, 1.0 - (branch_length / 2.0)],
            ],
            dtype=float,
        )

    cached_context = compile_finite_state_likelihood_context(
        tree,
        FiniteStateLikelihoodSurface(
            model_name="binary synthetic",
            state_labels=("0", "1"),
            root_prior=numpy.array([0.5, 0.5], dtype=float),
            root_prior_source="provided",
            observation_policy="reject",
            transition_matrix_for_branch_length=transition,
            leaf_likelihood_vector_for_state=_binary_leaf_vector,
            transition_rate_matrix=rate_matrix,
            parameter_values={"exchangeability": 1.0},
        ),
    )
    uncached_context = compile_finite_state_likelihood_context(
        tree,
        FiniteStateLikelihoodSurface(
            model_name="binary synthetic",
            state_labels=("0", "1"),
            root_prior=numpy.array([0.5, 0.5], dtype=float),
            root_prior_source="provided",
            observation_policy="reject",
            transition_matrix_for_branch_length=transition,
            leaf_likelihood_vector_for_state=_binary_leaf_vector,
            parameter_values={"exchangeability": 1.0},
        ),
    )
    compressed_patterns = CompressedAlignmentSitePatterns(
        source_path=None,
        taxon_order=["A", "B"],
        alignment_length=2,
        pattern_count=2,
        patterns=[
            AlignmentSitePattern(
                pattern_id="pattern-1",
                states=("0", "1"),
                weight=1,
                site_positions=[1],
            ),
            AlignmentSitePattern(
                pattern_id="pattern-2",
                states=("1", "0"),
                weight=1,
                site_positions=[2],
            ),
        ],
    )

    cached_log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        cached_context,
        compressed_patterns,
    )
    uncached_log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        uncached_context,
        compressed_patterns,
    )

    assert numpy.isclose(
        cached_log_likelihood,
        uncached_log_likelihood,
        rtol=0.0,
        atol=1e-12,
    )


def _binary_leaf_vector(_taxon_name: str, observed_state: str) -> numpy.ndarray:
    if observed_state == "0":
        return numpy.array([1.0, 0.0], dtype=float)
    if observed_state == "1":
        return numpy.array([0.0, 1.0], dtype=float)
    raise ValueError(f"unexpected synthetic state '{observed_state}'")
