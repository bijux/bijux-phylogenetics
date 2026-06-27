from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    FiniteStateLikelihoodSurface,
    clear_shared_finite_state_transition_cache,
    compile_finite_state_likelihood_context,
    evaluate_finite_state_tree_log_likelihood_from_patterns,
    shared_finite_state_transition_cache_entry_count,
    update_compiled_finite_state_context_branch_length,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    AlignmentSitePattern,
    CompressedAlignmentSitePatterns,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    return FIXTURES / "trees" / name


def test_compiled_branch_update_reuses_shared_transition_cache_entries() -> None:
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
    )
    first_context = compile_finite_state_likelihood_context(tree, surface)
    second_context = compile_finite_state_likelihood_context(tree, surface)
    compressed_patterns = _compressed_patterns()

    first_log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        first_context,
        compressed_patterns,
    )
    changed_leaf_node_id = next(
        node.node_id for node in tree.iter_leaves() if node.name == "A"
    )
    update_compiled_finite_state_context_branch_length(
        first_context,
        child_node_id=changed_leaf_node_id,
        branch_length=0.2,
    )
    updated_log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        first_context,
        compressed_patterns,
    )
    rebuilt_tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    rebuilt_tree.root.children[0].branch_length = 0.2
    rebuilt_tree.refresh()
    rebuilt_context = compile_finite_state_likelihood_context(rebuilt_tree, surface)
    rebuilt_log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        rebuilt_context,
        compressed_patterns,
    )

    assert second_context.transition_model_state == first_context.transition_model_state
    assert evaluation_count == 2
    assert shared_finite_state_transition_cache_entry_count() == 2
    assert not math.isclose(
        first_log_likelihood,
        updated_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        updated_log_likelihood,
        rebuilt_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_compiled_branch_update_noop_preserves_cached_site_evaluation() -> None:
    clear_shared_finite_state_transition_cache()
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    evaluation_count = 0

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

    context = compile_finite_state_likelihood_context(
        tree,
        FiniteStateLikelihoodSurface(
            model_name="binary synthetic",
            state_labels=("0", "1"),
            root_prior=numpy.array([0.5, 0.5], dtype=float),
            root_prior_source="provided",
            observation_policy="reject",
            transition_matrix_for_branch_length=transition,
            leaf_likelihood_vector_for_state=_binary_leaf_vector,
        ),
    )
    compressed_patterns = _compressed_patterns()
    evaluate_finite_state_tree_log_likelihood_from_patterns(context, compressed_patterns)
    cached_site_evaluation = next(iter(context.site_evaluation_by_key.values()))
    changed_leaf_node_id = next(
        node.node_id for node in tree.iter_leaves() if node.name == "A"
    )

    update_compiled_finite_state_context_branch_length(
        context,
        child_node_id=changed_leaf_node_id,
        branch_length=1.0,
    )

    assert next(iter(context.site_evaluation_by_key.values())) is cached_site_evaluation
    assert evaluation_count == 2


def _compressed_patterns() -> CompressedAlignmentSitePatterns:
    return CompressedAlignmentSitePatterns(
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


def _binary_leaf_vector(_taxon_name: str, observed_state: str) -> numpy.ndarray:
    if observed_state == "0":
        return numpy.array([1.0, 0.0], dtype=float)
    return numpy.array([0.0, 1.0], dtype=float)
