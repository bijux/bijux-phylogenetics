from __future__ import annotations

from collections.abc import Mapping
from functools import cache
import itertools
import math

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_clade_id,
    informative_rooted_clade_nodes,
    split_sort_key,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode, descendant_taxa
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .cost_matrix import load_sankoff_cost_matrix
from .matrix import load_parsimony_character_matrix
from .models import (
    FitchCharacterMatrix,
    ParsimonyCharacterWeights,
    SankoffCostMatrix,
)
from .tree_length import tree_length
from .weights import resolve_parsimony_character_weights

SUPPORTED_RESAMPLED_PARSIMONY_METHODS = frozenset(
    {
        "fitch",
        "wagner",
        "sankoff",
        "dollo",
        "camin-sokal",
        "acctran",
        "deltran",
    }
)
DEFAULT_MAX_EXACT_TAXA = 6


def resolve_resampled_parsimony_context(
    matrix: FitchCharacterMatrix,
    *,
    method: str,
    cost_matrix: SankoffCostMatrix | str | None,
    allow_asymmetric_costs: bool,
    character_weights: ParsimonyCharacterWeights | Mapping[str, float] | None,
    max_exact_taxa: int,
    workflow_name: str,
    unsupported_method_code: str,
    invalid_taxon_limit_code: str,
    taxon_limit_exceeded_code: str,
    cost_matrix_required_code: str,
) -> tuple[
    FitchCharacterMatrix,
    SankoffCostMatrix | None,
    ParsimonyCharacterWeights,
    tuple[PhyloTree, ...],
]:
    if method not in SUPPORTED_RESAMPLED_PARSIMONY_METHODS:
        raise ParsimonyAnalysisError(
            f"{workflow_name} requires one supported parsimony method",
            code=unsupported_method_code,
            details={
                "method": method,
                "supported_methods": sorted(SUPPORTED_RESAMPLED_PARSIMONY_METHODS),
            },
        )
    if max_exact_taxa < 4:
        raise ParsimonyAnalysisError(
            f"{workflow_name} exact search requires an exact-search taxon limit of at least four",
            code=invalid_taxon_limit_code,
            details={"max_exact_taxa": max_exact_taxa},
        )
    resolved_cost_matrix = _resolve_cost_matrix(
        method=method,
        cost_matrix=cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
        workflow_name=workflow_name,
        cost_matrix_required_code=cost_matrix_required_code,
    )
    resolved_weights = resolve_parsimony_character_weights(
        matrix.character_ids,
        character_weights,
    )
    taxa = sorted(matrix.states_by_taxon)
    if len(taxa) > max_exact_taxa:
        raise ParsimonyAnalysisError(
            f"{workflow_name} exact search is currently limited to small taxon sets",
            code=taxon_limit_exceeded_code,
            details={
                "taxon_count": len(taxa),
                "max_exact_taxa": max_exact_taxa,
            },
        )
    candidate_trees = enumerate_rooted_binary_trees(tuple(taxa))
    return matrix, resolved_cost_matrix, resolved_weights, candidate_trees


def resolve_resampled_parsimony_matrix(
    matrix: FitchCharacterMatrix | str | None,
    *,
    taxon_column: str | None,
) -> FitchCharacterMatrix:
    if isinstance(matrix, FitchCharacterMatrix):
        return matrix
    if matrix is None:
        raise ValueError("matrix path or matrix object is required")
    return load_parsimony_character_matrix(matrix, taxon_column=taxon_column)


def _resolve_cost_matrix(
    *,
    method: str,
    cost_matrix: SankoffCostMatrix | str | None,
    allow_asymmetric_costs: bool,
    workflow_name: str,
    cost_matrix_required_code: str,
) -> SankoffCostMatrix | None:
    if method != "sankoff":
        return None
    if cost_matrix is None:
        raise ParsimonyAnalysisError(
            f"{workflow_name} requires a Sankoff cost matrix when method='sankoff'",
            code=cost_matrix_required_code,
            details={"method": method},
        )
    if isinstance(cost_matrix, SankoffCostMatrix):
        return cost_matrix
    return load_sankoff_cost_matrix(
        cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
    )


@cache
def enumerate_rooted_binary_trees(taxa: tuple[str, ...]) -> tuple[PhyloTree, ...]:
    if len(taxa) < 2:
        raise ValueError("rooted binary tree enumeration requires at least two taxa")
    return tuple(
        PhyloTree(root=root, source_format="newick", rooted=True)
        for root in _enumerate_rooted_binary_subtrees(taxa)
    )


@cache
def _enumerate_rooted_binary_subtrees(taxa: tuple[str, ...]) -> tuple[TreeNode, ...]:
    if len(taxa) == 1:
        return (TreeNode(name=taxa[0]),)
    anchor = taxa[0]
    trailing_taxa = taxa[1:]
    roots_by_newick: dict[str, TreeNode] = {}
    for left_size in range(len(trailing_taxa)):
        for left_subset in itertools.combinations(trailing_taxa, left_size):
            left_taxa = tuple(sorted((anchor, *left_subset)))
            right_taxa = tuple(taxon for taxon in taxa if taxon not in set(left_taxa))
            if not right_taxa:
                continue
            for left_root in _enumerate_rooted_binary_subtrees(left_taxa):
                for right_root in _enumerate_rooted_binary_subtrees(right_taxa):
                    children = _order_child_subtrees(
                        left_root.copy(), right_root.copy()
                    )
                    root = TreeNode(children=children)
                    key = dumps_newick(
                        PhyloTree(root=root.copy(), source_format="newick", rooted=True)
                    )
                    roots_by_newick.setdefault(key, root)
    return tuple(roots_by_newick[key] for key in sorted(roots_by_newick))


def _order_child_subtrees(left: TreeNode, right: TreeNode) -> list[TreeNode]:
    left_key = tuple(descendant_taxa(left))
    right_key = tuple(descendant_taxa(right))
    if right_key < left_key:
        return [right, left]
    return [left, right]


def select_best_tree(
    candidate_trees: tuple[PhyloTree, ...],
    matrix: FitchCharacterMatrix,
    *,
    method: str,
    state_order: list[str] | None,
    cost_matrix: SankoffCostMatrix | None,
    allow_asymmetric_costs: bool,
    character_weights: ParsimonyCharacterWeights,
) -> tuple[PhyloTree, float, int]:
    equal_best_trees, best_score = select_equal_best_trees(
        candidate_trees,
        matrix,
        method=method,
        state_order=state_order,
        cost_matrix=cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
        character_weights=character_weights,
    )
    return equal_best_trees[0].copy(), best_score, len(equal_best_trees)


def select_equal_best_trees(
    candidate_trees: tuple[PhyloTree, ...],
    matrix: FitchCharacterMatrix,
    *,
    method: str,
    state_order: list[str] | None,
    cost_matrix: SankoffCostMatrix | None,
    allow_asymmetric_costs: bool,
    character_weights: ParsimonyCharacterWeights,
) -> tuple[tuple[PhyloTree, ...], float]:
    best_score: float | None = None
    best_trees_by_newick: dict[str, PhyloTree] = {}
    for tree in candidate_trees:
        score_report = tree_length(
            tree,
            matrix,
            method=method,
            state_order=state_order,
            cost_matrix=cost_matrix,
            allow_asymmetric_costs=allow_asymmetric_costs,
            character_weights=character_weights,
        )
        score = score_report.total_score
        candidate_newick = dumps_newick(tree)
        if (
            best_score is None
            or score < best_score
            and not math.isclose(score, best_score)
        ):
            best_score = score
            best_trees_by_newick = {candidate_newick: tree.copy()}
            continue
        if best_score is not None and math.isclose(score, best_score):
            best_trees_by_newick.setdefault(candidate_newick, tree.copy())
    if best_score is None or not best_trees_by_newick:
        raise AssertionError("exact parsimony search produced no candidate trees")
    return (
        tuple(best_trees_by_newick[key].copy() for key in sorted(best_trees_by_newick)),
        best_score,
    )


def build_resampled_matrix(
    matrix: FitchCharacterMatrix,
    *,
    sampled_character_ids: list[str],
    replicate_index: int,
    workflow_label: str,
) -> FitchCharacterMatrix:
    replicate_character_ids = [
        f"{workflow_label}_replicate_{replicate_index:04d}_draw_{draw_index:04d}"
        for draw_index in range(1, len(sampled_character_ids) + 1)
    ]
    states_by_taxon = {
        taxon: {
            replicate_character_id: states_by_character[source_character_id]
            for replicate_character_id, source_character_id in zip(
                replicate_character_ids,
                sampled_character_ids,
                strict=True,
            )
        }
        for taxon, states_by_character in matrix.states_by_taxon.items()
    }
    return FitchCharacterMatrix(
        matrix_path=matrix.matrix_path,
        taxon_column=matrix.taxon_column,
        character_ids=replicate_character_ids,
        states_by_taxon=states_by_taxon,
    )


def build_resampled_weights(
    weights: ParsimonyCharacterWeights,
    *,
    sampled_character_ids: list[str],
    replicate_index: int,
    workflow_label: str,
) -> ParsimonyCharacterWeights:
    replicate_weights = {
        f"{workflow_label}_replicate_{replicate_index:04d}_draw_{draw_index:04d}": (
            weights.weights_by_character[source_character_id]
        )
        for draw_index, source_character_id in enumerate(sampled_character_ids, start=1)
    }
    return ParsimonyCharacterWeights(
        weights_path=weights.weights_path,
        weights_by_character=replicate_weights,
    )


def build_clade_support_rows(
    reference_tree: PhyloTree,
    replicate_trees: list[PhyloTree],
) -> list[dict[str, object]]:
    shared_taxa = {
        node.name for node in reference_tree.iter_leaves() if node.name is not None
    }
    replicate_clade_sets = [
        set(informative_rooted_clade_nodes(tree, shared_taxa).keys())
        for tree in replicate_trees
    ]
    rows: list[dict[str, object]] = []
    reference_clades = informative_rooted_clade_nodes(reference_tree, shared_taxa)
    for clade_signature, node in sorted(
        reference_clades.items(),
        key=lambda item: split_sort_key(item[0]),
    ):
        supporting_tree_count = sum(
            1
            for replicate_clades in replicate_clade_sets
            if clade_signature in replicate_clades
        )
        clade_frequency = round(supporting_tree_count / len(replicate_trees), 15)
        rows.append(
            {
                "branch_id": canonical_clade_id(clade_signature),
                "node_name": node.name,
                "descendant_taxa": sorted(clade_signature),
                "supporting_tree_count": supporting_tree_count,
                "clade_frequency": clade_frequency,
                "support_percent": round(clade_frequency * 100.0, 15),
            }
        )
    return rows
