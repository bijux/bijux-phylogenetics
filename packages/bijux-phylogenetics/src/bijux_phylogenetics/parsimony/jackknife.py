from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from random import Random

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .exact_resampling import (
    DEFAULT_MAX_EXACT_TAXA,
    build_clade_support_rows,
    build_resampled_matrix,
    build_resampled_weights,
    resolve_resampled_parsimony_context,
    resolve_resampled_parsimony_matrix,
    select_best_tree,
)
from .models import (
    FitchCharacterMatrix,
    ParsimonyCharacterWeights,
    ParsimonyJackknifeCladeSupport,
    ParsimonyJackknifeReplicate,
    ParsimonyJackknifeReport,
    SankoffCostMatrix,
)


def jackknife_parsimony(
    matrix: FitchCharacterMatrix | Path,
    *,
    method: str,
    replicate_count: int,
    random_seed: int,
    retain_probability: float = 0.75,
    taxon_column: str | None = None,
    state_order: list[str] | None = None,
    cost_matrix: SankoffCostMatrix | Path | None = None,
    allow_asymmetric_costs: bool = False,
    character_weights: (
        ParsimonyCharacterWeights | Mapping[str, float] | Path | None
    ) = None,
    max_exact_taxa: int = DEFAULT_MAX_EXACT_TAXA,
) -> ParsimonyJackknifeReport:
    """Infer jackknife replicate trees under exact small-taxon parsimony search."""
    if replicate_count <= 0:
        raise ParsimonyAnalysisError(
            "parsimony jackknife requires at least one replicate",
            code="parsimony_jackknife_invalid_replicate_count",
            details={"replicate_count": replicate_count},
        )
    if not 0.0 < retain_probability <= 1.0:
        raise ParsimonyAnalysisError(
            "parsimony jackknife retain probability must be greater than 0 and at most 1",
            code="parsimony_jackknife_invalid_retain_probability",
            details={"retain_probability": retain_probability},
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
        workflow_name="parsimony jackknife",
        unsupported_method_code="parsimony_jackknife_method_unsupported",
        invalid_taxon_limit_code="parsimony_jackknife_invalid_exact_taxon_limit",
        taxon_limit_exceeded_code="parsimony_jackknife_taxon_limit_exceeded",
        cost_matrix_required_code="parsimony_jackknife_cost_matrix_required",
    )
    reference_tree, reference_score, reference_optimal_tree_count = select_best_tree(
        candidate_trees,
        resolved_matrix,
        method=method,
        state_order=state_order,
        cost_matrix=resolved_cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
        character_weights=resolved_weights,
    )
    rng = Random(random_seed)  # nosec B311
    replicate_rows: list[ParsimonyJackknifeReplicate] = []
    replicate_trees = []
    for replicate_index in range(1, replicate_count + 1):
        retained_character_ids = _sample_retained_characters(
            resolved_matrix.character_ids,
            retain_probability=retain_probability,
            rng=rng,
        )
        replicate_matrix = build_resampled_matrix(
            resolved_matrix,
            sampled_character_ids=retained_character_ids,
            replicate_index=replicate_index,
            workflow_label="jackknife",
        )
        replicate_weights = build_resampled_weights(
            resolved_weights,
            sampled_character_ids=retained_character_ids,
            replicate_index=replicate_index,
            workflow_label="jackknife",
        )
        best_tree, best_score, optimal_tree_count = select_best_tree(
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
            ParsimonyJackknifeReplicate(
                replicate_index=replicate_index,
                retained_character_count=len(retained_character_ids),
                retained_character_ids=retained_character_ids,
                best_score=best_score,
                optimal_tree_count=optimal_tree_count,
                tree_newick=dumps_newick(best_tree),
            )
        )
    clade_support_rows = [
        ParsimonyJackknifeCladeSupport(**row)
        for row in build_clade_support_rows(reference_tree, replicate_trees)
    ]
    return ParsimonyJackknifeReport(
        algorithm="parsimony-jackknife",
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
        retain_probability=retain_probability,
        candidate_tree_count=len(candidate_trees),
        max_exact_taxa=max_exact_taxa,
        reference_score=reference_score,
        reference_optimal_tree_count=reference_optimal_tree_count,
        reference_tree_newick=dumps_newick(reference_tree),
        replicate_rows=replicate_rows,
        clade_support_rows=clade_support_rows,
    )


def _sample_retained_characters(
    character_ids: list[str],
    *,
    retain_probability: float,
    rng: Random,
) -> list[str]:
    retained = [
        character_id
        for character_id in character_ids
        if rng.random() <= retain_probability
    ]
    if retained:
        return retained
    return [character_ids[rng.randrange(len(character_ids))]]
