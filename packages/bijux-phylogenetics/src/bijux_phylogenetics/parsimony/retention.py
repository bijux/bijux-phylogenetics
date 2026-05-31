from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .homoplasy import (
    character_kind,
    minimum_possible_steps,
    unordered_maximum_possible_steps,
)
from .matrix import load_parsimony_character_matrix
from .models import (
    FitchCharacterMatrix,
    ParsimonyRetentionCharacterIndex,
    ParsimonyRetentionIndexReport,
)
from .tree_length import tree_length

_RETENTION_INDEX_METHODS = frozenset({"fitch", "acctran", "deltran"})


def retention_index(
    tree: Path,
    matrix: FitchCharacterMatrix | Path,
    *,
    method: str,
) -> ParsimonyRetentionIndexReport:
    """Compute per-character and aggregate retention index for one matrix."""
    resolved_method = _resolve_retention_method(method)
    resolved_matrix = (
        matrix
        if isinstance(matrix, FitchCharacterMatrix)
        else load_parsimony_character_matrix(matrix)
    )
    tree_length_report = tree_length(
        tree,
        resolved_matrix,
        method=resolved_method,
    )
    raw_scores_by_character = {
        row.character_id: row.raw_score for row in tree_length_report.step_rows
    }

    character_rows: list[ParsimonyRetentionCharacterIndex] = []
    included_character_count = 0
    excluded_character_count = 0
    minimum_possible_steps_total = 0.0
    maximum_possible_steps_total = 0.0
    observed_steps_total = 0.0
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
            state_order=None,
        )
        maximum_steps = unordered_maximum_possible_steps(resolved_matrix, character_id)
        observed_steps = raw_scores_by_character[character_id]
        row_kind = character_kind(resolved_matrix, character_id, observed_states)
        denominator = maximum_steps - minimum_steps
        if denominator == 0.0:
            excluded_character_count += 1
            character_rows.append(
                ParsimonyRetentionCharacterIndex(
                    character_id=character_id,
                    character_kind=row_kind,
                    observed_states=observed_states,
                    minimum_possible_steps=minimum_steps,
                    maximum_possible_steps=maximum_steps,
                    observed_steps=observed_steps,
                    retention_index=None,
                    undefined_reason="zero_range_character",
                )
            )
            continue
        included_character_count += 1
        minimum_possible_steps_total += minimum_steps
        maximum_possible_steps_total += maximum_steps
        observed_steps_total += observed_steps
        character_rows.append(
            ParsimonyRetentionCharacterIndex(
                character_id=character_id,
                character_kind=row_kind,
                observed_states=observed_states,
                minimum_possible_steps=minimum_steps,
                maximum_possible_steps=maximum_steps,
                observed_steps=observed_steps,
                retention_index=(maximum_steps - observed_steps) / denominator,
                undefined_reason=None,
            )
        )

    aggregate_denominator = maximum_possible_steps_total - minimum_possible_steps_total
    aggregate_retention_index = None
    undefined_reason = None
    if aggregate_denominator > 0.0:
        aggregate_retention_index = (
            maximum_possible_steps_total - observed_steps_total
        ) / aggregate_denominator
    else:
        undefined_reason = "no_defined_retention_characters"
    return ParsimonyRetentionIndexReport(
        algorithm="parsimony-retention-index",
        method=resolved_method,
        tree_path=tree_length_report.tree_path,
        matrix_path=tree_length_report.matrix_path,
        taxon_column=tree_length_report.taxon_column,
        taxon_count=tree_length_report.taxon_count,
        character_count=tree_length_report.character_count,
        included_character_count=included_character_count,
        excluded_character_count=excluded_character_count,
        minimum_possible_steps_total=minimum_possible_steps_total,
        maximum_possible_steps_total=maximum_possible_steps_total,
        observed_steps_total=observed_steps_total,
        retention_index=aggregate_retention_index,
        undefined_reason=undefined_reason,
        character_rows=character_rows,
    )


def _resolve_retention_method(method: str) -> str:
    resolved_method = method.strip().lower()
    if resolved_method not in _RETENTION_INDEX_METHODS:
        raise ParsimonyAnalysisError(
            "retention index currently supports only unordered Fitch-style parsimony methods",
            code="parsimony_retention_index_method_unsupported",
            details={
                "method": method,
                "supported_methods": sorted(_RETENTION_INDEX_METHODS),
                "reason": "maximum possible step counts are currently owned only for unordered Fitch-style characters",
            },
        )
    return resolved_method
