from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_clade_id,
    informative_rooted_clade_nodes,
    split_sort_key,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .exact_resampling import (
    DEFAULT_MAX_EXACT_TAXA,
    resolve_resampled_parsimony_context,
    resolve_resampled_parsimony_matrix,
)
from .models import (
    FitchCharacterMatrix,
    ParsimonyBremerSupportReport,
    ParsimonyBremerSupportRow,
    ParsimonyCharacterWeights,
    SankoffCostMatrix,
)
from .topology_search import (
    resolve_topology_search_tree,
    score_topology_search_tree,
    validate_topology_search_tree,
)


@dataclass(slots=True)
class _ShortestLackingTree:
    score: float | None = None
    tree_count: int = 0
    tree_newick: str | None = None


def compute_parsimony_bremer_support(
    tree: PhyloTree | Path,
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
) -> ParsimonyBremerSupportReport:
    """Compute exact Bremer support for every informative clade on one reference tree."""
    resolved_tree, resolved_tree_path = resolve_topology_search_tree(tree)
    validate_topology_search_tree(
        resolved_tree,
        workflow_name="parsimony Bremer support",
    )
    resolved_matrix = resolve_resampled_parsimony_matrix(
        matrix,
        taxon_column=taxon_column,
    )
    reference_taxa = set(resolved_tree.tip_names)
    matrix_taxa = set(resolved_matrix.states_by_taxon)
    if reference_taxa != matrix_taxa:
        raise ParsimonyAnalysisError(
            "parsimony Bremer support requires matrix taxa to match the tree tips exactly",
            code="parsimony_bremer_support_taxa_mismatch",
            details={
                "missing_from_matrix": sorted(reference_taxa - matrix_taxa),
                "extra_in_matrix": sorted(matrix_taxa - reference_taxa),
            },
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
        workflow_name="parsimony Bremer support",
        unsupported_method_code="parsimony_bremer_support_method_unsupported",
        invalid_taxon_limit_code="parsimony_bremer_support_invalid_exact_taxon_limit",
        taxon_limit_exceeded_code="parsimony_bremer_support_taxon_limit_exceeded",
        cost_matrix_required_code="parsimony_bremer_support_cost_matrix_required",
    )
    reference_tree = resolved_tree.copy().refresh()
    reference_tree_newick = dumps_newick(reference_tree)
    reference_tree_score = score_topology_search_tree(
        reference_tree,
        resolved_matrix,
        method=method,
        state_order=state_order,
        cost_matrix=resolved_cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
        character_weights=resolved_weights,
    )
    reference_clades = informative_rooted_clade_nodes(reference_tree, matrix_taxa)
    shortest_lacking_by_clade = {
        clade_signature: _ShortestLackingTree() for clade_signature in reference_clades
    }
    optimal_score: float | None = None
    optimal_tree_count = 0
    optimal_tree_newick: str | None = None
    for candidate_tree in candidate_trees:
        candidate_score = score_topology_search_tree(
            candidate_tree,
            resolved_matrix,
            method=method,
            state_order=state_order,
            cost_matrix=resolved_cost_matrix,
            allow_asymmetric_costs=allow_asymmetric_costs,
            character_weights=resolved_weights,
        )
        candidate_newick = dumps_newick(candidate_tree)
        if optimal_score is None or (
            candidate_score < optimal_score
            and not math.isclose(candidate_score, optimal_score)
        ):
            optimal_score = candidate_score
            optimal_tree_count = 1
            optimal_tree_newick = candidate_newick
        elif optimal_score is not None and math.isclose(candidate_score, optimal_score):
            optimal_tree_count += 1
            if optimal_tree_newick is None or candidate_newick < optimal_tree_newick:
                optimal_tree_newick = candidate_newick
        candidate_clades = set(
            informative_rooted_clade_nodes(candidate_tree, matrix_taxa)
        )
        for clade_signature, best_lacking in shortest_lacking_by_clade.items():
            if clade_signature in candidate_clades:
                continue
            if best_lacking.score is None or (
                candidate_score < best_lacking.score
                and not math.isclose(candidate_score, best_lacking.score)
            ):
                best_lacking.score = candidate_score
                best_lacking.tree_count = 1
                best_lacking.tree_newick = candidate_newick
                continue
            if best_lacking.score is not None and math.isclose(
                candidate_score,
                best_lacking.score,
            ):
                best_lacking.tree_count += 1
                if (
                    best_lacking.tree_newick is None
                    or candidate_newick < best_lacking.tree_newick
                ):
                    best_lacking.tree_newick = candidate_newick
    if optimal_score is None or optimal_tree_newick is None:
        raise AssertionError("exact Bremer support search produced no candidate trees")
    bremer_rows: list[ParsimonyBremerSupportRow] = []
    for clade_signature, node in sorted(
        reference_clades.items(),
        key=lambda item: split_sort_key(item[0]),
    ):
        best_lacking = shortest_lacking_by_clade[clade_signature]
        if best_lacking.score is None or best_lacking.tree_newick is None:
            raise AssertionError(
                "every informative reference-tree clade should be absent from at least one candidate tree"
            )
        bremer_rows.append(
            ParsimonyBremerSupportRow(
                branch_id=canonical_clade_id(clade_signature),
                node_name=node.name,
                descendant_taxa=sorted(clade_signature),
                shortest_lacking_score=best_lacking.score,
                decay_index=best_lacking.score - optimal_score,
                shortest_lacking_tree_count=best_lacking.tree_count,
                shortest_lacking_tree_newick=best_lacking.tree_newick,
            )
        )
    reference_tree_score_delta_from_optimal = reference_tree_score - optimal_score
    return ParsimonyBremerSupportReport(
        algorithm="parsimony-bremer-support",
        method=method,
        tree_path=resolved_tree_path,
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
        reference_tree_newick=reference_tree_newick,
        reference_tree_score=reference_tree_score,
        optimal_score=optimal_score,
        optimal_tree_count=optimal_tree_count,
        optimal_tree_newick=optimal_tree_newick,
        reference_tree_score_delta_from_optimal=reference_tree_score_delta_from_optimal,
        reference_tree_is_optimal=math.isclose(
            reference_tree_score_delta_from_optimal,
            0.0,
        ),
        bremer_rows=bremer_rows,
    )
