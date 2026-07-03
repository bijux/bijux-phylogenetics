from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    FiniteStateLikelihoodSurface,
    compile_finite_state_likelihood_analysis,
    evaluate_compiled_finite_state_likelihood_analysis,
    update_compiled_finite_state_analysis_branch_length,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


def test_compiled_analysis_reuses_cached_site_state_between_evaluations() -> None:
    analysis, evaluation_counter = _compiled_four_tip_analysis()

    first_log_likelihood = evaluate_compiled_finite_state_likelihood_analysis(analysis)
    second_log_likelihood = evaluate_compiled_finite_state_likelihood_analysis(analysis)

    assert math.isclose(
        first_log_likelihood,
        second_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert len(analysis.context.site_evaluation_by_key) == 4
    assert evaluation_counter["transition"] == 6


def test_compiled_branch_update_refreshes_only_affected_cached_path() -> None:
    analysis, _evaluation_counter = _compiled_four_tip_analysis()
    evaluate_compiled_finite_state_likelihood_analysis(analysis)

    cached_site_evaluation = next(
        iter(analysis.context.site_evaluation_by_key.values())
    )
    root = analysis.context.tree.root
    left_internal = root.children[0]
    right_internal = root.children[1]
    changed_leaf = left_internal.children[0]

    previous_right_conditional = (
        cached_site_evaluation.pruning_pass.conditional_by_node_id[
            right_internal.node_id
        ]
    )
    previous_left_conditional = (
        cached_site_evaluation.pruning_pass.conditional_by_node_id[
            left_internal.node_id
        ]
    )
    previous_root_conditional = (
        cached_site_evaluation.pruning_pass.conditional_by_node_id[root.node_id]
    )

    update_compiled_finite_state_analysis_branch_length(
        analysis,
        child_node_id=changed_leaf.node_id,
        branch_length=0.37,
    )

    assert (
        cached_site_evaluation.pruning_pass.conditional_by_node_id[
            right_internal.node_id
        ]
        is previous_right_conditional
    )
    assert (
        cached_site_evaluation.pruning_pass.conditional_by_node_id[
            left_internal.node_id
        ]
        is not previous_left_conditional
    )
    assert (
        cached_site_evaluation.pruning_pass.conditional_by_node_id[root.node_id]
        is not previous_root_conditional
    )


def test_compiled_branch_update_matches_full_recompilation() -> None:
    analysis, evaluation_counter = _compiled_four_tip_analysis()
    original_total = evaluate_compiled_finite_state_likelihood_analysis(analysis)

    update_compiled_finite_state_analysis_branch_length(
        analysis,
        child_node_id=analysis.context.tree.root.children[0].children[0].node_id,
        branch_length=0.37,
    )
    refreshed_total = evaluate_compiled_finite_state_likelihood_analysis(analysis)

    rebuilt_tree = _four_tip_tree()
    rebuilt_tree.root.children[0].children[0].branch_length = 0.37
    rebuilt_tree.refresh()
    rebuilt_analysis, rebuilt_counter = _compiled_four_tip_analysis(tree=rebuilt_tree)
    rebuilt_total = evaluate_compiled_finite_state_likelihood_analysis(rebuilt_analysis)

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
    assert evaluation_counter["transition"] == 7
    assert rebuilt_counter["transition"] == 6


def _compiled_four_tip_analysis(
    *,
    tree: PhyloTree | None = None,
) -> tuple[object, dict[str, int]]:
    resolved_tree = _four_tip_tree() if tree is None else tree
    records = [
        AlignmentRecord(identifier="A", sequence="0011"),
        AlignmentRecord(identifier="B", sequence="0111"),
        AlignmentRecord(identifier="C", sequence="1001"),
        AlignmentRecord(identifier="D", sequence="1010"),
    ]
    evaluation_counter = {"transition": 0}

    def transition(branch_length: float) -> numpy.ndarray:
        evaluation_counter["transition"] += 1
        stay = 1.0 - min(branch_length, 0.45)
        switch = 1.0 - stay
        return numpy.array(
            [[stay, switch], [switch / 2.0, 1.0 - (switch / 2.0)]],
            dtype=float,
        )

    analysis = compile_finite_state_likelihood_analysis(
        resolved_tree,
        FiniteStateLikelihoodSurface(
            model_name="binary synthetic",
            state_labels=("0", "1"),
            root_prior=numpy.array([0.6, 0.4], dtype=float),
            root_prior_source="provided",
            observation_policy="reject",
            transition_matrix_for_branch_length=transition,
            leaf_likelihood_vector_for_state=_binary_leaf_vector,
        ),
        records=records,
        pattern_compiler=_pattern_compiler,
        missing_data_policy="reject",
    )
    return analysis, evaluation_counter


def _four_tip_tree() -> PhyloTree:
    return PhyloTree(
        TreeNode(
            children=[
                TreeNode(
                    children=[
                        TreeNode(name="A", branch_length=0.11),
                        TreeNode(name="B", branch_length=0.23),
                    ],
                    branch_length=0.41,
                ),
                TreeNode(
                    children=[
                        TreeNode(name="C", branch_length=0.17),
                        TreeNode(name="D", branch_length=0.29),
                    ],
                    branch_length=0.53,
                ),
            ]
        )
    )


def _pattern_compiler(
    records: list[AlignmentRecord],
    source_path: Path | None,
) -> CompressedAlignmentSitePatterns:
    return compress_alignment_site_patterns_from_records(
        records,
        source_path=source_path,
    )


def _binary_leaf_vector(_taxon_name: str, observed_state: str) -> numpy.ndarray:
    if observed_state == "0":
        return numpy.array([1.0, 0.0], dtype=float)
    return numpy.array([0.0, 1.0], dtype=float)
