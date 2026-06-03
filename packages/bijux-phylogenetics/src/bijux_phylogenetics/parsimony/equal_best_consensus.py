from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError
from bijux_phylogenetics.trees.tree_sets.consensus import (
    _build_consensus_tree_with_threshold_from_trees,
)

from .exact_resampling import (
    DEFAULT_MAX_EXACT_TAXA,
    resolve_resampled_parsimony_context,
    resolve_resampled_parsimony_matrix,
    select_equal_best_trees,
)
from .models import (
    FitchCharacterMatrix,
    ParsimonyCharacterWeights,
    ParsimonyConsensusSummary,
    ParsimonyEqualBestConsensusReport,
    ParsimonyEqualBestTree,
    SankoffCostMatrix,
)

DEFAULT_MAX_RETAINED_EQUAL_BEST_TREES = 128


def summarize_equal_best_parsimony_trees(
    matrix: FitchCharacterMatrix | Path,
    *,
    method: str,
    taxon_column: str | None = None,
    state_order: list[str] | None = None,
    cost_matrix: SankoffCostMatrix | Path | None = None,
    allow_asymmetric_costs: bool = False,
    character_weights: (
        ParsimonyCharacterWeights | Mapping[str, float] | Path | None
    ) = None,
    max_exact_taxa: int = DEFAULT_MAX_EXACT_TAXA,
    max_retained_equal_best_trees: int = DEFAULT_MAX_RETAINED_EQUAL_BEST_TREES,
) -> ParsimonyEqualBestConsensusReport:
    """Enumerate the exact equal-best parsimony tree set and summarize its consensus."""
    if max_retained_equal_best_trees <= 0:
        raise ParsimonyAnalysisError(
            "equal-best parsimony consensus requires a retained-tree cap of at least one",
            code="parsimony_equal_best_consensus_invalid_tree_cap",
            details={
                "max_retained_equal_best_trees": max_retained_equal_best_trees,
            },
        )
    resolved_matrix = resolve_resampled_parsimony_matrix(
        matrix,
        taxon_column=taxon_column,
    )
    (
        resolved_matrix,
        resolved_cost_matrix,
        resolved_weights,
        candidate_trees,
    ) = resolve_resampled_parsimony_context(
        resolved_matrix,
        method=method,
        cost_matrix=cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
        character_weights=character_weights,
        max_exact_taxa=max_exact_taxa,
        workflow_name="equal-best parsimony consensus",
        unsupported_method_code="parsimony_equal_best_consensus_method_unsupported",
        invalid_taxon_limit_code="parsimony_equal_best_consensus_invalid_exact_taxon_limit",
        taxon_limit_exceeded_code="parsimony_equal_best_consensus_taxon_limit_exceeded",
        cost_matrix_required_code="parsimony_equal_best_consensus_cost_matrix_required",
    )
    equal_best_trees, best_score = select_equal_best_trees(
        candidate_trees,
        resolved_matrix,
        method=method,
        state_order=state_order,
        cost_matrix=resolved_cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
        character_weights=resolved_weights,
    )
    retained_equal_best_trees = equal_best_trees[:max_retained_equal_best_trees]
    retained_all_equal_best_trees = len(retained_equal_best_trees) == len(
        equal_best_trees
    )
    strict_consensus = None
    majority_consensus = None
    if retained_all_equal_best_trees:
        strict_tree, strict_included_clade_count = (
            _build_consensus_tree_with_threshold_from_trees(
                list(retained_equal_best_trees),
                threshold=1.0,
            )
        )
        majority_tree, majority_included_clade_count = (
            _build_consensus_tree_with_threshold_from_trees(
                list(retained_equal_best_trees),
                threshold=0.5,
            )
        )
        strict_consensus = ParsimonyConsensusSummary(
            consensus_method="strict",
            consensus_threshold=1.0,
            tree_count=len(retained_equal_best_trees),
            included_clade_count=strict_included_clade_count,
            consensus_newick=dumps_newick(strict_tree),
        )
        majority_consensus = ParsimonyConsensusSummary(
            consensus_method="majority-rule",
            consensus_threshold=0.5,
            tree_count=len(retained_equal_best_trees),
            included_clade_count=majority_included_clade_count,
            consensus_newick=dumps_newick(majority_tree),
        )
    return ParsimonyEqualBestConsensusReport(
        algorithm="parsimony-equal-best-consensus",
        method=method,
        matrix_path=resolved_matrix.matrix_path,
        cost_matrix_path=None
        if resolved_cost_matrix is None
        else resolved_cost_matrix.matrix_path,
        weights_path=resolved_weights.weights_path,
        taxon_column=resolved_matrix.taxon_column,
        taxon_count=resolved_matrix.taxon_count,
        character_count=resolved_matrix.character_count,
        candidate_tree_count=len(candidate_trees),
        max_exact_taxa=max_exact_taxa,
        max_retained_equal_best_trees=max_retained_equal_best_trees,
        best_score=best_score,
        equal_best_tree_count=len(equal_best_trees),
        retained_equal_best_tree_count=len(retained_equal_best_trees),
        retained_all_equal_best_trees=retained_all_equal_best_trees,
        strict_consensus=strict_consensus,
        majority_consensus=majority_consensus,
        equal_best_tree_rows=[
            ParsimonyEqualBestTree(
                tree_index=tree_index,
                total_score=best_score,
                tree_newick=dumps_newick(tree),
            )
            for tree_index, tree in enumerate(retained_equal_best_trees, start=1)
        ],
    )
