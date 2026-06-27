from __future__ import annotations

import math

import numpy

from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    FiniteStateLikelihoodSurface,
    clear_shared_finite_state_transition_cache,
    compile_finite_state_likelihood_context,
    evaluate_finite_state_tree_log_likelihood_from_patterns,
    update_compiled_finite_state_context_topology,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    AlignmentSitePattern,
    CompressedAlignmentSitePatterns,
)
from bijux_phylogenetics.phylo.topology.rooted_nni import (
    apply_rooted_nni_move,
    resolve_rooted_nni_move_candidate,
)
from bijux_phylogenetics.phylo.topology.rooted_spr import (
    apply_rooted_spr_move,
    resolve_rooted_spr_move_candidate,
)
from bijux_phylogenetics.phylo.topology.rooted_tbr import (
    apply_rooted_tbr_move,
    resolve_rooted_tbr_move_candidate,
)


def test_compiled_topology_refresh_matches_nni_recompilation() -> None:
    original_tree = loads_newick("(((A:0.2,C:0.3):0.4,B:0.5):0.6,D:0.7);")
    moved_tree = apply_rooted_nni_move(
        original_tree,
        resolve_rooted_nni_move_candidate(original_tree, 1)[0],
    )
    _assert_topology_refresh_matches_recompilation(
        original_tree=original_tree,
        moved_tree=moved_tree,
        reusable_clade_id="A|C",
        expect_transition_reuse=True,
    )


def test_compiled_topology_refresh_matches_spr_recompilation() -> None:
    original_tree = loads_newick("(((A:0.2,B:0.3):0.4,C:0.5):0.6,D:0.7);")
    moved_tree = apply_rooted_spr_move(
        original_tree,
        resolve_rooted_spr_move_candidate(original_tree, 1)[0],
    )
    _assert_topology_refresh_matches_recompilation(
        original_tree=original_tree,
        moved_tree=moved_tree,
        reusable_clade_id="A",
        expect_transition_reuse=True,
    )


def test_compiled_topology_refresh_matches_tbr_recompilation() -> None:
    original_tree = loads_newick("((((A:0.2,B:0.3):0.4,C:0.5):0.6,D:0.7):0.8,E:0.9);")
    moved_tree = apply_rooted_tbr_move(
        original_tree,
        resolve_rooted_tbr_move_candidate(original_tree, 1)[0],
    )
    _assert_topology_refresh_matches_recompilation(
        original_tree=original_tree,
        moved_tree=moved_tree,
        reusable_clade_id="E",
        expect_transition_reuse=False,
    )


def _assert_topology_refresh_matches_recompilation(
    *,
    original_tree,
    moved_tree,
    reusable_clade_id: str,
    expect_transition_reuse: bool,
) -> None:
    clear_shared_finite_state_transition_cache()
    surface = _binary_surface()
    compressed_patterns = _compressed_patterns(original_tree.tip_names)
    cached_context = compile_finite_state_likelihood_context(original_tree, surface)
    rebuilt_context = compile_finite_state_likelihood_context(moved_tree, surface)

    original_total = evaluate_finite_state_tree_log_likelihood_from_patterns(
        cached_context,
        compressed_patterns,
    )
    cached_site_evaluation = next(iter(cached_context.site_evaluation_by_key.values()))
    reusable_conditional_before = _conditional_for_clade(
        cached_context,
        cached_site_evaluation.pruning_pass.conditional_by_node_id,
        reusable_clade_id,
    )
    reusable_transition_before = _transition_for_clade(cached_context, reusable_clade_id)

    update_compiled_finite_state_context_topology(
        cached_context,
        tree=moved_tree,
    )
    refreshed_total = evaluate_finite_state_tree_log_likelihood_from_patterns(
        cached_context,
        compressed_patterns,
    )
    rebuilt_total = evaluate_finite_state_tree_log_likelihood_from_patterns(
        rebuilt_context,
        compressed_patterns,
    )

    assert not math.isclose(
        original_total,
        refreshed_total,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        refreshed_total,
        rebuilt_total,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert (
        next(iter(cached_context.site_evaluation_by_key.values()))
        is cached_site_evaluation
    )
    assert (
        _conditional_for_clade(
            cached_context,
            cached_site_evaluation.pruning_pass.conditional_by_node_id,
            reusable_clade_id,
        )
        is reusable_conditional_before
    )
    refreshed_transition = _transition_for_clade(cached_context, reusable_clade_id)
    if expect_transition_reuse:
        assert refreshed_transition is reusable_transition_before
    else:
        assert numpy.allclose(
            refreshed_transition,
            reusable_transition_before,
            rtol=0.0,
            atol=1e-12,
        )


def _conditional_for_clade(
    context,
    conditional_by_node_id: dict[str, numpy.ndarray],
    clade_id: str,
):
    node_id = next(
        node.node_id
        for node in context.tree.iter_nodes(order="preorder")
        if node is not context.tree.root and "|".join(node.descendant_taxa) == clade_id
    )
    return conditional_by_node_id[node_id]


def _transition_for_clade(context, clade_id: str) -> numpy.ndarray:
    node_id = next(
        node.node_id
        for node in context.tree.iter_nodes(order="preorder")
        if node is not context.tree.root and "|".join(node.descendant_taxa) == clade_id
    )
    return context.transition_matrix_by_node_id[node_id]


def _compressed_patterns(taxa: list[str]) -> CompressedAlignmentSitePatterns:
    pattern_states = {
        4: [
            ("0", "1", "0", "1"),
            ("1", "1", "0", "0"),
            ("0", "0", "1", "1"),
        ],
        5: [
            ("0", "1", "0", "1", "0"),
            ("1", "1", "0", "0", "1"),
            ("0", "0", "1", "1", "0"),
        ],
    }
    return CompressedAlignmentSitePatterns(
        source_path=None,
        taxon_order=taxa,
        alignment_length=3,
        pattern_count=3,
        patterns=[
            AlignmentSitePattern(
                pattern_id="pattern-1",
                states=pattern_states[len(taxa)][0],
                weight=1,
                site_positions=[1],
            ),
            AlignmentSitePattern(
                pattern_id="pattern-2",
                states=pattern_states[len(taxa)][1],
                weight=1,
                site_positions=[2],
            ),
            AlignmentSitePattern(
                pattern_id="pattern-3",
                states=pattern_states[len(taxa)][2],
                weight=1,
                site_positions=[3],
            ),
        ],
    )


def _binary_surface() -> FiniteStateLikelihoodSurface:
    return FiniteStateLikelihoodSurface(
        model_name="binary synthetic",
        state_labels=("0", "1"),
        root_prior=numpy.array([0.5, 0.5], dtype=float),
        root_prior_source="provided",
        observation_policy="reject",
        transition_matrix_for_branch_length=lambda branch_length: numpy.array(
            [
                [1.0 - min(branch_length, 0.45), min(branch_length, 0.45)],
                [
                    min(branch_length, 0.45) / 2.0,
                    1.0 - (min(branch_length, 0.45) / 2.0),
                ],
            ],
            dtype=float,
        ),
        leaf_likelihood_vector_for_state=_binary_leaf_vector,
    )


def _binary_leaf_vector(_taxon_name: str, observed_state: str) -> numpy.ndarray:
    if observed_state == "0":
        return numpy.array([1.0, 0.0], dtype=float)
    return numpy.array([0.0, 1.0], dtype=float)
