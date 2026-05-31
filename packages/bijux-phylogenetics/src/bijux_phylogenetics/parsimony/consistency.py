from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .homoplasy import character_kind, minimum_possible_steps
from .matrix import load_parsimony_character_matrix
from .models import (
    FitchCharacterMatrix,
    ParsimonyConsistencyCharacterIndex,
    ParsimonyConsistencyIndexReport,
)
from .tree_length import tree_length

_CONSISTENCY_INDEX_METHODS = frozenset(
    {
        "fitch",
        "wagner",
        "dollo",
        "camin-sokal",
        "acctran",
        "deltran",
    }
)


def consistency_index(
    tree: Path,
    matrix: FitchCharacterMatrix | Path,
    *,
    method: str,
    state_order: list[str] | None = None,
) -> ParsimonyConsistencyIndexReport:
    """Compute per-character and aggregate consistency index for one matrix."""
    resolved_method = _resolve_consistency_method(method)
    resolved_matrix = (
        matrix
        if isinstance(matrix, FitchCharacterMatrix)
        else load_parsimony_character_matrix(matrix)
    )
    tree_length_report = tree_length(
        tree,
        resolved_matrix,
        method=resolved_method,
        state_order=state_order,
    )

    character_rows: list[ParsimonyConsistencyCharacterIndex] = []
    included_character_count = 0
    excluded_character_count = 0
    minimum_possible_steps_total = 0.0
    observed_steps_total = 0.0
    raw_scores_by_character = {
        row.character_id: row.raw_score for row in tree_length_report.step_rows
    }
    for character_id in resolved_matrix.character_ids:
        observed_states = sorted(
            {
                resolved_matrix.states_by_taxon[taxon][character_id]
                for taxon in resolved_matrix.states_by_taxon
            }
        )
        minimum_steps = minimum_possible_steps(
            resolved_method,
            observed_states,
            character_id=character_id,
            state_order=state_order,
        )
        observed_steps = raw_scores_by_character[character_id]
        row_kind = character_kind(resolved_matrix, character_id, observed_states)
        if minimum_steps == 0.0 and observed_steps == 0.0:
            excluded_character_count += 1
            character_rows.append(
                ParsimonyConsistencyCharacterIndex(
                    character_id=character_id,
                    character_kind=row_kind,
                    observed_states=observed_states,
                    minimum_possible_steps=minimum_steps,
                    observed_steps=observed_steps,
                    consistency_index=None,
                    undefined_reason="constant_character",
                )
            )
            continue
        included_character_count += 1
        minimum_possible_steps_total += minimum_steps
        observed_steps_total += observed_steps
        character_rows.append(
            ParsimonyConsistencyCharacterIndex(
                character_id=character_id,
                character_kind=row_kind,
                observed_states=observed_states,
                minimum_possible_steps=minimum_steps,
                observed_steps=observed_steps,
                consistency_index=minimum_steps / observed_steps,
                undefined_reason=None,
            )
        )

    aggregate_consistency_index = None
    undefined_reason = None
    if observed_steps_total > 0.0:
        aggregate_consistency_index = (
            minimum_possible_steps_total / observed_steps_total
        )
    else:
        undefined_reason = "no_variable_characters"
    return ParsimonyConsistencyIndexReport(
        algorithm="parsimony-consistency-index",
        method=resolved_method,
        tree_path=tree_length_report.tree_path,
        matrix_path=tree_length_report.matrix_path,
        taxon_column=tree_length_report.taxon_column,
        taxon_count=tree_length_report.taxon_count,
        character_count=tree_length_report.character_count,
        included_character_count=included_character_count,
        excluded_character_count=excluded_character_count,
        minimum_possible_steps_total=minimum_possible_steps_total,
        observed_steps_total=observed_steps_total,
        consistency_index=aggregate_consistency_index,
        undefined_reason=undefined_reason,
        character_rows=character_rows,
    )


def _resolve_consistency_method(method: str) -> str:
    resolved_method = method.strip().lower()
    if resolved_method == "sankoff":
        raise ParsimonyAnalysisError(
            "consistency index does not yet support arbitrary Sankoff step matrices",
            code="parsimony_consistency_index_method_unsupported",
            details={
                "method": resolved_method,
                "reason": "minimum possible steps for arbitrary transition-cost matrices are not yet owned separately",
            },
        )
    if resolved_method not in _CONSISTENCY_INDEX_METHODS:
        raise ParsimonyAnalysisError(
            "consistency index requires one supported parsimony method",
            code="parsimony_consistency_index_method_unknown",
            details={
                "method": method,
                "supported_methods": sorted(_CONSISTENCY_INDEX_METHODS),
            },
        )
    return resolved_method
