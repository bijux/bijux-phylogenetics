from __future__ import annotations

from collections.abc import Mapping
from functools import lru_cache
import itertools
import math
from pathlib import Path
from random import Random

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
    ParsimonyBootstrapCladeSupport,
    ParsimonyBootstrapReplicate,
    ParsimonyBootstrapReport,
    ParsimonyCharacterWeights,
    SankoffCostMatrix,
)
from .tree_length import tree_length
from .weights import resolve_parsimony_character_weights

_SUPPORTED_BOOTSTRAP_METHODS = frozenset(
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
_DEFAULT_MAX_EXACT_TAXA = 6


def bootstrap_parsimony(
    matrix: FitchCharacterMatrix | Path,
    *,
    method: str,
    replicate_count: int,
    random_seed: int,
    taxon_column: str | None = None,
    state_order: list[str] | None = None,
    cost_matrix: SankoffCostMatrix | Path | None = None,
    allow_asymmetric_costs: bool = False,
    character_weights: (
        ParsimonyCharacterWeights | Mapping[str, float] | Path | None
    ) = None,
    max_exact_taxa: int = _DEFAULT_MAX_EXACT_TAXA,
) -> ParsimonyBootstrapReport:
    """Infer bootstrap replicate trees under exact small-taxon parsimony search."""
    if method not in _SUPPORTED_BOOTSTRAP_METHODS:
        raise ParsimonyAnalysisError(
            "parsimony bootstrap requires one supported parsimony method",
            code="parsimony_bootstrap_method_unsupported",
            details={"method": method, "supported_methods": sorted(_SUPPORTED_BOOTSTRAP_METHODS)},
        )
    if replicate_count <= 0:
        raise ParsimonyAnalysisError(
            "parsimony bootstrap requires at least one replicate",
            code="parsimony_bootstrap_invalid_replicate_count",
            details={"replicate_count": replicate_count},
        )
    if max_exact_taxa < 4:
        raise ParsimonyAnalysisError(
            "parsimony bootstrap exact search requires an exact-search taxon limit of at least four",
            code="parsimony_bootstrap_invalid_exact_taxon_limit",
            details={"max_exact_taxa": max_exact_taxa},
        )
    resolved_matrix = (
        matrix
        if isinstance(matrix, FitchCharacterMatrix)
        else load_parsimony_character_matrix(matrix, taxon_column=taxon_column)
    )
    resolved_cost_matrix = _resolve_cost_matrix(
        method=method,
        cost_matrix=cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
    )
    resolved_weights = resolve_parsimony_character_weights(
        resolved_matrix.character_ids,
        character_weights,
    )
    taxa = sorted(resolved_matrix.states_by_taxon)
    if len(taxa) > max_exact_taxa:
        raise ParsimonyAnalysisError(
            "parsimony bootstrap exact search is currently limited to small taxon sets",
            code="parsimony_bootstrap_taxon_limit_exceeded",
            details={
                "taxon_count": len(taxa),
                "max_exact_taxa": max_exact_taxa,
            },
        )
    candidate_trees = _enumerate_rooted_binary_trees(tuple(taxa))
    reference_tree, reference_score, reference_optimal_tree_count = _select_best_tree(
        candidate_trees,
        resolved_matrix,
        method=method,
        state_order=state_order,
        cost_matrix=resolved_cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
        character_weights=resolved_weights,
    )
    rng = Random(random_seed)
    replicate_rows: list[ParsimonyBootstrapReplicate] = []
    replicate_trees: list[PhyloTree] = []
    for replicate_index in range(1, replicate_count + 1):
        sampled_character_ids = [
            rng.choice(resolved_matrix.character_ids)
            for _ in range(resolved_matrix.character_count)
        ]
        replicate_matrix = _build_replicate_matrix(
            resolved_matrix,
            sampled_character_ids=sampled_character_ids,
            replicate_index=replicate_index,
        )
        replicate_weights = _build_replicate_weights(
            resolved_weights,
            sampled_character_ids=sampled_character_ids,
            replicate_index=replicate_index,
        )
        best_tree, best_score, optimal_tree_count = _select_best_tree(
            candidate_trees,
            replicate_matrix,
            method=method,
            state_order=state_order,
            cost_matrix=resolved_cost_matrix,
            allow_asymmetric_costs=allow_asymmetric_costs,
            character_weights=replicate_weights,
        )
        replicate_trees.append(best_tree)
        replicate_rows.append(
            ParsimonyBootstrapReplicate(
                replicate_index=replicate_index,
                sampled_character_ids=sampled_character_ids,
                best_score=best_score,
                optimal_tree_count=optimal_tree_count,
                tree_newick=dumps_newick(best_tree),
            )
        )
    clade_support_rows = _build_clade_support_rows(reference_tree, replicate_trees)
    return ParsimonyBootstrapReport(
        algorithm="parsimony-bootstrap",
        method=method,
        matrix_path=resolved_matrix.matrix_path,
        cost_matrix_path=None
        if resolved_cost_matrix is None
        else resolved_cost_matrix.matrix_path,
        weights_path=resolved_weights.weights_path,
        taxon_column=resolved_matrix.taxon_column,
        taxon_count=resolved_matrix.taxon_count,
        character_count=resolved_matrix.character_count,
        replicate_count=replicate_count,
        random_seed=random_seed,
        candidate_tree_count=len(candidate_trees),
        max_exact_taxa=max_exact_taxa,
        reference_score=reference_score,
        reference_optimal_tree_count=reference_optimal_tree_count,
        reference_tree_newick=dumps_newick(reference_tree),
        replicate_rows=replicate_rows,
        clade_support_rows=clade_support_rows,
    )


def _resolve_cost_matrix(
    *,
    method: str,
    cost_matrix: SankoffCostMatrix | Path | None,
    allow_asymmetric_costs: bool,
) -> SankoffCostMatrix | None:
    if method != "sankoff":
        return None
    if cost_matrix is None:
        raise ParsimonyAnalysisError(
            "parsimony bootstrap requires a Sankoff cost matrix when method='sankoff'",
            code="parsimony_bootstrap_cost_matrix_required",
            details={"method": method},
        )
    if isinstance(cost_matrix, SankoffCostMatrix):
        return cost_matrix
    return load_sankoff_cost_matrix(
        cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
    )


@lru_cache(maxsize=None)
def _enumerate_rooted_binary_trees(taxa: tuple[str, ...]) -> tuple[PhyloTree, ...]:
    if len(taxa) < 2:
        raise ValueError("rooted binary tree enumeration requires at least two taxa")
    return tuple(
        PhyloTree(root=root, source_format="newick", rooted=True)
        for root in _enumerate_rooted_binary_subtrees(taxa)
    )


@lru_cache(maxsize=None)
def _enumerate_rooted_binary_subtrees(taxa: tuple[str, ...]) -> tuple[TreeNode, ...]:
    if len(taxa) == 1:
        return (TreeNode(name=taxa[0]),)
    anchor = taxa[0]
    trailing_taxa = taxa[1:]
    roots_by_newick: dict[str, TreeNode] = {}
    for left_size in range(0, len(trailing_taxa)):
        for left_subset in itertools.combinations(trailing_taxa, left_size):
            left_taxa = tuple(sorted((anchor, *left_subset)))
            right_taxa = tuple(taxon for taxon in taxa if taxon not in set(left_taxa))
            if not right_taxa:
                continue
            for left_root in _enumerate_rooted_binary_subtrees(left_taxa):
                for right_root in _enumerate_rooted_binary_subtrees(right_taxa):
                    children = _order_child_subtrees(left_root.copy(), right_root.copy())
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


def _select_best_tree(
    candidate_trees: tuple[PhyloTree, ...],
    matrix: FitchCharacterMatrix,
    *,
    method: str,
    state_order: list[str] | None,
    cost_matrix: SankoffCostMatrix | None,
    allow_asymmetric_costs: bool,
    character_weights: ParsimonyCharacterWeights,
) -> tuple[PhyloTree, float, int]:
    best_tree: PhyloTree | None = None
    best_newick: str | None = None
    best_score: float | None = None
    optimal_tree_count = 0
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
        if best_score is None or score < best_score and not math.isclose(score, best_score):
            best_tree = tree
            best_newick = candidate_newick
            best_score = score
            optimal_tree_count = 1
            continue
        if best_score is not None and math.isclose(score, best_score):
            optimal_tree_count += 1
            if best_newick is None or candidate_newick < best_newick:
                best_tree = tree
                best_newick = candidate_newick
                best_score = score
    if best_tree is None or best_score is None:
        raise AssertionError("exact parsimony search produced no candidate trees")
    return best_tree.copy(), best_score, optimal_tree_count


def _build_replicate_matrix(
    matrix: FitchCharacterMatrix,
    *,
    sampled_character_ids: list[str],
    replicate_index: int,
) -> FitchCharacterMatrix:
    replicate_character_ids = [
        f"bootstrap_replicate_{replicate_index:04d}_draw_{draw_index:04d}"
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


def _build_replicate_weights(
    weights: ParsimonyCharacterWeights,
    *,
    sampled_character_ids: list[str],
    replicate_index: int,
) -> ParsimonyCharacterWeights:
    replicate_weights = {
        f"bootstrap_replicate_{replicate_index:04d}_draw_{draw_index:04d}": (
            weights.weights_by_character[source_character_id]
        )
        for draw_index, source_character_id in enumerate(sampled_character_ids, start=1)
    }
    return ParsimonyCharacterWeights(
        weights_path=weights.weights_path,
        weights_by_character=replicate_weights,
    )


def _build_clade_support_rows(
    reference_tree: PhyloTree,
    replicate_trees: list[PhyloTree],
) -> list[ParsimonyBootstrapCladeSupport]:
    shared_taxa = set(node.name for node in reference_tree.iter_leaves() if node.name is not None)
    replicate_clade_sets = [
        set(informative_rooted_clade_nodes(tree, shared_taxa).keys())
        for tree in replicate_trees
    ]
    rows: list[ParsimonyBootstrapCladeSupport] = []
    reference_clades = informative_rooted_clade_nodes(reference_tree, shared_taxa)
    for clade_signature, node in sorted(
        reference_clades.items(),
        key=lambda item: split_sort_key(item[0]),
    ):
        supporting_tree_count = sum(
            1 for replicate_clades in replicate_clade_sets if clade_signature in replicate_clades
        )
        clade_frequency = round(supporting_tree_count / len(replicate_trees), 15)
        rows.append(
            ParsimonyBootstrapCladeSupport(
                branch_id=canonical_clade_id(clade_signature),
                node_name=node.name,
                descendant_taxa=sorted(clade_signature),
                supporting_tree_count=supporting_tree_count,
                clade_frequency=clade_frequency,
                support_percent=round(clade_frequency * 100.0, 15),
            )
        )
    return rows
