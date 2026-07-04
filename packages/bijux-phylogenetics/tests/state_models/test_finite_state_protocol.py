from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    CompressedAlignmentSitePatterns,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    FiniteStateLikelihoodSurface,
    compile_finite_state_likelihood_context,
    evaluate_finite_state_site_log_likelihood,
    evaluate_finite_state_site_log_likelihood_rows_from_patterns,
    evaluate_finite_state_tree_log_likelihood_from_patterns,
)
from bijux_phylogenetics.phylo.likelihood.patterns import AlignmentSitePattern

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    return FIXTURES / "trees" / name


def test_compiled_finite_state_context_reuses_cached_branch_matrices() -> None:
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

    surface = FiniteStateLikelihoodSurface(
        model_name="binary synthetic",
        state_labels=("0", "1"),
        root_prior=numpy.array([0.5, 0.5], dtype=float),
        root_prior_source="provided",
        observation_policy="reject",
        transition_matrix_for_branch_length=transition,
        leaf_likelihood_vector_for_state=_binary_leaf_vector,
    )
    context = compile_finite_state_likelihood_context(tree, surface)

    assert evaluation_count == 2
    assert len(context.transition_matrix_by_node_id) == 2
    assert all(
        node_id.startswith("root:clade:A|B/")
        for node_id in context.transition_matrix_by_node_id
    )


def test_finite_state_surface_matches_analytical_binary_two_tip_likelihood() -> None:
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    surface = FiniteStateLikelihoodSurface(
        model_name="binary analytical",
        state_labels=("0", "1"),
        root_prior=numpy.array([0.5, 0.5], dtype=float),
        root_prior_source="provided",
        observation_policy="reject",
        transition_matrix_for_branch_length=lambda _branch_length: numpy.array(
            [[0.9, 0.1], [0.2, 0.8]],
            dtype=float,
        ),
        leaf_likelihood_vector_for_state=_binary_leaf_vector,
    )
    context = compile_finite_state_likelihood_context(tree, surface)
    compressed_patterns = CompressedAlignmentSitePatterns(
        source_path=None,
        taxon_order=["A", "B"],
        alignment_length=1,
        pattern_count=1,
        patterns=[
            AlignmentSitePattern(
                pattern_id="pattern-1",
                states=("0", "1"),
                weight=1,
                site_positions=[1],
            )
        ],
    )

    log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )
    site_log_likelihood = evaluate_finite_state_site_log_likelihood(
        context,
        ("0", "1"),
        taxon_order=["A", "B"],
    )

    assert math.isclose(log_likelihood, math.log(0.125), abs_tol=1e-12)
    assert math.isclose(site_log_likelihood, math.log(0.125), abs_tol=1e-12)


def test_finite_state_site_rows_reconstruct_pattern_total() -> None:
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    surface = FiniteStateLikelihoodSurface(
        model_name="binary analytical",
        state_labels=("0", "1"),
        root_prior=numpy.array([0.5, 0.5], dtype=float),
        root_prior_source="provided",
        observation_policy="reject",
        transition_matrix_for_branch_length=lambda _branch_length: numpy.array(
            [[0.9, 0.1], [0.2, 0.8]],
            dtype=float,
        ),
        leaf_likelihood_vector_for_state=_binary_leaf_vector,
    )
    context = compile_finite_state_likelihood_context(tree, surface)
    compressed_patterns = CompressedAlignmentSitePatterns(
        source_path=None,
        taxon_order=["A", "B"],
        alignment_length=2,
        pattern_count=1,
        patterns=[
            AlignmentSitePattern(
                pattern_id="pattern-1",
                states=("0", "1"),
                weight=2,
                site_positions=[1, 2],
            )
        ],
    )

    rows, log_likelihood = evaluate_finite_state_site_log_likelihood_rows_from_patterns(
        context,
        compressed_patterns,
    )

    assert len(rows) == 2
    assert rows[0].site_position == 1
    assert rows[1].site_position == 2
    assert rows[0].pattern_id == rows[1].pattern_id == "pattern-1"
    assert math.isclose(log_likelihood, 2.0 * math.log(0.125), abs_tol=1e-12)
    assert math.isclose(rows[0].log_likelihood, math.log(0.125), abs_tol=1e-12)
    assert math.isclose(rows[1].log_likelihood, math.log(0.125), abs_tol=1e-12)


def _binary_leaf_vector(taxon_name: str, observed_state: str) -> numpy.ndarray:
    if observed_state == "0":
        return numpy.array([1.0, 0.0], dtype=float)
    if observed_state == "1":
        return numpy.array([0.0, 1.0], dtype=float)
    raise ValueError(f"unexpected synthetic state '{observed_state}' at {taxon_name}")
