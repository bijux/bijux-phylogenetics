from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from random import Random

from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .models import (
    FitchCharacterMatrix,
    ParsimonyCharacterWeights,
    ParsimonyRatchetBestTreeHistory,
    ParsimonyRatchetCycle,
    ParsimonyRatchetReport,
    SankoffCostMatrix,
)
from .spr import search_parsimony_spr
from .topology_search import (
    resolve_topology_search_cost_matrix,
    resolve_topology_search_matrix,
    resolve_topology_search_method,
    resolve_topology_search_tree,
    resolve_topology_search_weights,
    topology_search_prefer_score,
    validate_topology_search_tree,
)


def run_parsimony_ratchet(
    tree: PhyloTree | Path,
    matrix: FitchCharacterMatrix | Path,
    *,
    method: str,
    cycle_count: int,
    random_seed: int,
    perturbed_character_count: int,
    perturbation_factor: float,
    taxon_column: str | None = None,
    state_order: list[str] | None = None,
    cost_matrix: SankoffCostMatrix | Path | None = None,
    allow_asymmetric_costs: bool = False,
    character_weights: ParsimonyCharacterWeights
    | Mapping[str, float]
    | Path
    | None = None,
) -> ParsimonyRatchetReport:
    """Run one deterministic parsimony ratchet over rooted SPR search cycles."""
    if cycle_count <= 0:
        raise ParsimonyAnalysisError(
            "parsimony ratchet requires at least one cycle",
            code="parsimony_ratchet_cycle_count_invalid",
            details={"cycle_count": cycle_count},
        )
    if perturbed_character_count <= 0:
        raise ParsimonyAnalysisError(
            "parsimony ratchet requires at least one perturbed character per cycle",
            code="parsimony_ratchet_perturbed_character_count_invalid",
            details={"perturbed_character_count": perturbed_character_count},
        )
    if perturbation_factor <= 0.0:
        raise ParsimonyAnalysisError(
            "parsimony ratchet requires a positive perturbation factor",
            code="parsimony_ratchet_perturbation_factor_invalid",
            details={"perturbation_factor": perturbation_factor},
        )
    resolved_tree, resolved_tree_path = resolve_topology_search_tree(tree)
    validate_topology_search_tree(
        resolved_tree,
        workflow_name="parsimony ratchet",
    )
    resolved_matrix = resolve_topology_search_matrix(
        matrix,
        taxon_column=taxon_column,
    )
    if perturbed_character_count > resolved_matrix.character_count:
        raise ParsimonyAnalysisError(
            "parsimony ratchet cannot perturb more characters than the matrix contains",
            code="parsimony_ratchet_perturbed_character_count_exceeds_matrix",
            details={
                "perturbed_character_count": perturbed_character_count,
                "character_count": resolved_matrix.character_count,
            },
        )
    resolved_method = resolve_topology_search_method(
        method,
        workflow_name="parsimony ratchet",
    )
    resolved_cost_matrix = resolve_topology_search_cost_matrix(
        method=resolved_method,
        cost_matrix=cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
        workflow_name="parsimony ratchet",
    )
    base_weights = resolve_topology_search_weights(
        resolved_matrix,
        character_weights,
    )
    start_report = search_parsimony_spr(
        resolved_tree,
        resolved_matrix,
        method=resolved_method,
        state_order=state_order,
        cost_matrix=resolved_cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
        character_weights=base_weights,
    )
    current_tree_newick = start_report.final_tree_newick
    current_score = start_report.final_score
    best_tree_newick = current_tree_newick
    best_score = current_score
    best_tree_history_rows = [
        ParsimonyRatchetBestTreeHistory(
            history_index=1,
            cycle_index=0,
            best_score=best_score,
            best_tree_newick=best_tree_newick,
        )
    ]
    cycle_rows: list[ParsimonyRatchetCycle] = []
    rng = Random(random_seed)  # nosec B311
    for cycle_index in range(1, cycle_count + 1):
        cycle_start_tree_newick = current_tree_newick
        cycle_start_score = current_score
        perturbed_character_ids = sorted(
            rng.sample(resolved_matrix.character_ids, perturbed_character_count)
        )
        perturbed_weights = _build_perturbed_weights(
            base_weights,
            perturbed_character_ids=perturbed_character_ids,
            perturbation_factor=perturbation_factor,
        )
        perturbed_report = search_parsimony_spr(
            loads_newick(current_tree_newick),
            resolved_matrix,
            method=resolved_method,
            state_order=state_order,
            cost_matrix=resolved_cost_matrix,
            allow_asymmetric_costs=allow_asymmetric_costs,
            character_weights=perturbed_weights,
        )
        restored_report = search_parsimony_spr(
            loads_newick(perturbed_report.final_tree_newick),
            resolved_matrix,
            method=resolved_method,
            state_order=state_order,
            cost_matrix=resolved_cost_matrix,
            allow_asymmetric_costs=allow_asymmetric_costs,
            character_weights=base_weights,
        )
        current_tree_newick = restored_report.final_tree_newick
        current_score = restored_report.final_score
        best_tree_improved = topology_search_prefer_score(
            current_score,
            current_tree_newick,
            best_score,
            best_tree_newick,
        )
        if best_tree_improved:
            best_score = current_score
            best_tree_newick = current_tree_newick
            best_tree_history_rows.append(
                ParsimonyRatchetBestTreeHistory(
                    history_index=len(best_tree_history_rows) + 1,
                    cycle_index=cycle_index,
                    best_score=best_score,
                    best_tree_newick=best_tree_newick,
                )
            )
        cycle_rows.append(
            ParsimonyRatchetCycle(
                cycle_index=cycle_index,
                start_score=cycle_start_score,
                start_tree_newick=cycle_start_tree_newick,
                perturbed_character_ids=perturbed_character_ids,
                perturbation_factor=perturbation_factor,
                perturbed_score=perturbed_report.final_score,
                perturbed_tree_newick=perturbed_report.final_tree_newick,
                perturbed_accepted_move_count=perturbed_report.accepted_move_count,
                restored_score=restored_report.final_score,
                restored_tree_newick=restored_report.final_tree_newick,
                restored_accepted_move_count=restored_report.accepted_move_count,
                best_score_after_cycle=best_score,
                best_tree_after_cycle=best_tree_newick,
                best_tree_improved=best_tree_improved,
            )
        )
    return ParsimonyRatchetReport(
        algorithm="parsimony-ratchet",
        method=resolved_method,
        tree_path=resolved_tree_path,
        matrix_path=resolved_matrix.matrix_path,
        cost_matrix_path=None
        if resolved_cost_matrix is None
        else resolved_cost_matrix.matrix_path,
        weights_path=base_weights.weights_path,
        taxon_column=resolved_matrix.taxon_column,
        taxon_count=resolved_matrix.taxon_count,
        character_count=resolved_matrix.character_count,
        cycle_count=cycle_count,
        random_seed=random_seed,
        perturbed_character_count=perturbed_character_count,
        perturbation_factor=perturbation_factor,
        start_tree_newick=start_report.final_tree_newick,
        start_score=start_report.final_score,
        final_tree_newick=current_tree_newick,
        final_score=current_score,
        best_tree_newick=best_tree_newick,
        best_score=best_score,
        cycle_rows=cycle_rows,
        best_tree_history_rows=best_tree_history_rows,
    )


def _build_perturbed_weights(
    base_weights: ParsimonyCharacterWeights,
    *,
    perturbed_character_ids: list[str],
    perturbation_factor: float,
) -> ParsimonyCharacterWeights:
    perturbed_character_set = set(perturbed_character_ids)
    return ParsimonyCharacterWeights(
        weights_path=base_weights.weights_path,
        weights_by_character={
            character_id: (
                weight * perturbation_factor
                if character_id in perturbed_character_set
                else weight
            )
            for character_id, weight in base_weights.weights_by_character.items()
        },
    )
