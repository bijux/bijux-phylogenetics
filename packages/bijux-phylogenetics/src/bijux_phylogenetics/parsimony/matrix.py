from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .models import FitchCharacterMatrix

_UNKNOWN_STATE_MARKERS = frozenset({"", "?", "-", ".", "na", "n/a", "null"})


def load_parsimony_character_matrix(
    path: Path,
    *,
    taxon_column: str | None = None,
) -> FitchCharacterMatrix:
    """Load a taxon-keyed discrete character matrix for unordered Fitch scoring."""
    table = load_taxon_table(path, taxon_column=taxon_column)
    character_ids = [column for column in table.columns if column != table.taxon_column]
    if not character_ids or table.row_count == 0:
        raise ParsimonyAnalysisError(
            "unordered fitch scoring requires at least one character column and one taxon row",
            code="parsimony_matrix_empty",
            details={
                "matrix_path": str(path),
                "row_count": table.row_count,
                "character_count": len(character_ids),
            },
        )
    states_by_taxon: dict[str, dict[str, str]] = {}
    for row in table.rows:
        taxon = row[table.taxon_column]
        states_by_character: dict[str, str] = {}
        for character_id in character_ids:
            state = row[character_id].strip()
            if state.lower() in _UNKNOWN_STATE_MARKERS:
                raise ParsimonyAnalysisError(
                    "unordered fitch scoring does not accept unknown character states",
                    code="parsimony_matrix_unknown_state",
                    details={
                        "matrix_path": str(path),
                        "taxon": taxon,
                        "character_id": character_id,
                        "value": row[character_id],
                    },
                )
            states_by_character[character_id] = state
        states_by_taxon[taxon] = states_by_character
    return FitchCharacterMatrix(
        matrix_path=path,
        taxon_column=table.taxon_column,
        character_ids=character_ids,
        states_by_taxon=states_by_taxon,
    )


def load_fitch_character_matrix(
    path: Path,
    *,
    taxon_column: str | None = None,
) -> FitchCharacterMatrix:
    """Backward-compatible alias for the parsimony character-matrix loader."""
    return load_parsimony_character_matrix(path, taxon_column=taxon_column)
