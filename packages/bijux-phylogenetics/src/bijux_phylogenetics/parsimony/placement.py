from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .fitch import score_fitch
from .matrix import load_fitch_character_matrix
from .models import (
    FitchCharacterMatrix,
    ParsimonyCharacterWeights,
    ParsimonyPlacementAlternativeRow,
    ParsimonyPlacementQuerySummary,
    ParsimonyPlacementReport,
)
from .weights import resolve_parsimony_character_weights

_PLACEMENT_WEIGHT_TOLERANCE = 1e-12


def place_parsimony_queries(
    tree: PhyloTree | Path,
    matrix: FitchCharacterMatrix | Path,
    query_matrix: FitchCharacterMatrix | Path,
    *,
    method: str = "fitch",
    character_weights: (
        ParsimonyCharacterWeights | Mapping[str, float] | Path | None
    ) = None,
) -> ParsimonyPlacementReport:
    """Place one or more query taxa on every edge with unordered Fitch parsimony."""
    normalized_method = method.strip().lower()
    if normalized_method != "fitch":
        raise ValueError("parsimony placement currently supports only the Fitch method")
    resolved_tree, tree_path = _resolve_tree(tree)
    resolved_matrix = (
        matrix
        if isinstance(matrix, FitchCharacterMatrix)
        else load_fitch_character_matrix(matrix)
    )
    resolved_query_matrix = (
        query_matrix
        if isinstance(query_matrix, FitchCharacterMatrix)
        else load_fitch_character_matrix(query_matrix)
    )
    resolved_weights = resolve_parsimony_character_weights(
        resolved_matrix.character_ids,
        character_weights,
    )
    _validate_reference_matrix_against_tree(resolved_tree, resolved_matrix)
    _validate_query_matrix_against_reference(
        resolved_matrix,
        resolved_query_matrix,
    )

    reference_report = score_fitch(
        resolved_tree,
        resolved_matrix,
        character_weights=resolved_weights,
    )
    edge_rows = _iter_candidate_edges(resolved_tree)
    query_summaries: list[ParsimonyPlacementQuerySummary] = []
    alternative_rows: list[ParsimonyPlacementAlternativeRow] = []

    for query_id in sorted(resolved_query_matrix.states_by_taxon):
        placement_rows = _score_query_placements(
            resolved_tree,
            resolved_matrix,
            query_id=query_id,
            query_states=resolved_query_matrix.states_by_taxon[query_id],
            reference_total_steps=reference_report.total_steps,
            reference_total_weighted_score=reference_report.total_weighted_score,
            edge_rows=edge_rows,
            character_weights=resolved_weights,
        )
        best_row = placement_rows[0]
        equally_best_count = sum(1 for row in placement_rows if row.is_equally_best)
        query_summaries.append(
            ParsimonyPlacementQuerySummary(
                query_id=query_id,
                character_count=resolved_matrix.character_count,
                best_edge_id=best_row.edge_id,
                best_child_name=best_row.child_name,
                best_descendant_taxa=list(best_row.descendant_taxa),
                best_total_steps=best_row.total_steps,
                best_additional_steps=best_row.additional_steps,
                best_total_weighted_score=best_row.total_weighted_score,
                best_additional_weighted_score=best_row.additional_weighted_score,
                candidate_placement_count=len(placement_rows),
                equally_best_placement_count=equally_best_count,
                selected_best_tree_newick=best_row.placed_tree_newick,
            )
        )
        alternative_rows.extend(placement_rows)

    return ParsimonyPlacementReport(
        algorithm="parsimony-placement",
        method="unordered-fitch",
        tree_path=tree_path,
        matrix_path=resolved_matrix.matrix_path,
        query_matrix_path=resolved_query_matrix.matrix_path,
        taxon_column=resolved_matrix.taxon_column,
        reference_taxon_count=resolved_matrix.taxon_count,
        character_count=resolved_matrix.character_count,
        edge_count=len(edge_rows),
        query_count=len(resolved_query_matrix.states_by_taxon),
        reference_total_steps=reference_report.total_steps,
        weights_path=resolved_weights.weights_path,
        reference_total_weighted_score=reference_report.total_weighted_score,
        query_summaries=query_summaries,
        alternative_rows=alternative_rows,
    )


def _resolve_tree(tree: PhyloTree | Path) -> tuple[PhyloTree, Path | None]:
    if isinstance(tree, Path):
        return load_tree(tree), tree
    return tree, None


def _validate_reference_matrix_against_tree(
    tree: PhyloTree,
    matrix: FitchCharacterMatrix,
) -> None:
    leaf_taxa = sorted(
        node.name for node in tree.iter_leaves() if node.name is not None
    )
    missing_from_matrix = sorted(set(leaf_taxa) - set(matrix.states_by_taxon))
    extra_in_matrix = sorted(set(matrix.states_by_taxon) - set(leaf_taxa))
    if missing_from_matrix:
        raise ParsimonyAnalysisError(
            "parsimony placement requires every tree taxon to be present in the reference character matrix",
            code="parsimony_placement_reference_missing_taxa",
            details={
                "missing_taxa": missing_from_matrix,
                "matrix_path": None
                if matrix.matrix_path is None
                else str(matrix.matrix_path),
            },
        )
    if extra_in_matrix:
        raise ParsimonyAnalysisError(
            "parsimony placement requires reference matrix taxa to match the tree tips exactly",
            code="parsimony_placement_reference_extra_taxa",
            details={
                "extra_taxa": extra_in_matrix,
                "matrix_path": None
                if matrix.matrix_path is None
                else str(matrix.matrix_path),
            },
        )


def _validate_query_matrix_against_reference(
    reference_matrix: FitchCharacterMatrix,
    query_matrix: FitchCharacterMatrix,
) -> None:
    if reference_matrix.character_ids != query_matrix.character_ids:
        raise ParsimonyAnalysisError(
            "parsimony placement requires query matrix characters to match the reference matrix exactly",
            code="parsimony_placement_query_character_mismatch",
            details={
                "reference_character_ids": reference_matrix.character_ids,
                "query_character_ids": query_matrix.character_ids,
            },
        )
    overlapping_taxa = sorted(
        set(reference_matrix.states_by_taxon) & set(query_matrix.states_by_taxon)
    )
    if overlapping_taxa:
        raise ParsimonyAnalysisError(
            "parsimony placement requires query taxa that are absent from the reference matrix",
            code="parsimony_placement_query_taxon_overlap",
            details={"overlapping_taxa": overlapping_taxa},
        )


def _iter_candidate_edges(
    tree: PhyloTree,
) -> list[tuple[str, str | None, list[str]]]:
    edge_rows: list[tuple[str, str | None, list[str]]] = []
    for _parent, child in tree.iter_edges():
        if child.node_id is None:
            raise ParsimonyAnalysisError(
                "parsimony placement requires stable edge identifiers on the reference tree",
                code="parsimony_placement_missing_edge_id",
            )
        edge_rows.append((child.node_id, child.name, child.descendant_taxa))
    return edge_rows


def _score_query_placements(
    tree: PhyloTree,
    reference_matrix: FitchCharacterMatrix,
    *,
    query_id: str,
    query_states: dict[str, str],
    reference_total_steps: int,
    reference_total_weighted_score: float,
    edge_rows: list[tuple[str, str | None, list[str]]],
    character_weights: ParsimonyCharacterWeights,
) -> list[ParsimonyPlacementAlternativeRow]:
    unresolved_rows: list[ParsimonyPlacementAlternativeRow] = []
    augmented_matrix = _build_augmented_matrix(reference_matrix, query_id, query_states)
    for edge_id, child_name, descendant_taxa in edge_rows:
        placed_tree = _place_query_on_edge(tree, query_id=query_id, edge_id=edge_id)
        score_report = score_fitch(
            placed_tree,
            augmented_matrix,
            character_weights=character_weights,
        )
        unresolved_rows.append(
            ParsimonyPlacementAlternativeRow(
                query_id=query_id,
                placement_rank=0,
                edge_id=edge_id,
                child_name=child_name,
                descendant_taxa=descendant_taxa,
                total_steps=score_report.total_steps,
                additional_steps=score_report.total_steps - reference_total_steps,
                total_weighted_score=score_report.total_weighted_score,
                additional_weighted_score=(
                    score_report.total_weighted_score - reference_total_weighted_score
                ),
                is_equally_best=False,
                placed_tree_newick=dumps_newick(placed_tree),
            )
        )
    return _rank_placement_rows(unresolved_rows)


def _build_augmented_matrix(
    reference_matrix: FitchCharacterMatrix,
    query_id: str,
    query_states: dict[str, str],
) -> FitchCharacterMatrix:
    return FitchCharacterMatrix(
        matrix_path=None,
        taxon_column=reference_matrix.taxon_column,
        character_ids=list(reference_matrix.character_ids),
        states_by_taxon={
            **reference_matrix.states_by_taxon,
            query_id: {
                character_id: query_states[character_id]
                for character_id in reference_matrix.character_ids
            },
        },
    )


def _place_query_on_edge(
    tree: PhyloTree,
    *,
    query_id: str,
    edge_id: str,
) -> PhyloTree:
    placed_tree = tree.copy().refresh()
    child = placed_tree.node_by_id(edge_id)
    parent = child.parent
    if parent is None:
        raise ParsimonyAnalysisError(
            "parsimony placement cannot target the root as a candidate edge",
            code="parsimony_placement_root_edge",
        )
    original_branch_length = child.branch_length
    placement_node = TreeNode(
        branch_length=original_branch_length,
        children=[
            child,
            TreeNode(
                name=query_id,
                branch_length=0.0 if original_branch_length is not None else None,
            ),
        ],
    )
    child.branch_length = 0.0 if original_branch_length is not None else None
    replacement_children = list(parent.children)
    for index, observed_child in enumerate(replacement_children):
        if observed_child.node_id == edge_id:
            replacement_children[index] = placement_node
            break
    parent.replace_children(replacement_children)
    return placed_tree.refresh()


def _rank_placement_rows(
    rows: list[ParsimonyPlacementAlternativeRow],
) -> list[ParsimonyPlacementAlternativeRow]:
    ordered_rows = sorted(
        rows,
        key=lambda row: (
            row.additional_weighted_score,
            row.additional_steps,
            row.edge_id,
        ),
    )
    best_weight = ordered_rows[0].additional_weighted_score
    best_steps = ordered_rows[0].additional_steps
    ranked_rows: list[ParsimonyPlacementAlternativeRow] = []
    for placement_rank, row in enumerate(ordered_rows, start=1):
        ranked_rows.append(
            ParsimonyPlacementAlternativeRow(
                query_id=row.query_id,
                placement_rank=placement_rank,
                edge_id=row.edge_id,
                child_name=row.child_name,
                descendant_taxa=list(row.descendant_taxa),
                total_steps=row.total_steps,
                additional_steps=row.additional_steps,
                total_weighted_score=row.total_weighted_score,
                additional_weighted_score=row.additional_weighted_score,
                is_equally_best=(
                    row.additional_steps == best_steps
                    and abs(row.additional_weighted_score - best_weight)
                    <= _PLACEMENT_WEIGHT_TOLERANCE
                ),
                placed_tree_newick=row.placed_tree_newick,
            )
        )
    return ranked_rows
