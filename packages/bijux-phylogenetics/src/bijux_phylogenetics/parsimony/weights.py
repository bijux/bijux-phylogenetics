from __future__ import annotations

from collections.abc import Mapping
import csv
import math
from pathlib import Path

from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .models import ParsimonyCharacterWeights


def load_parsimony_character_weights(path: Path) -> ParsimonyCharacterWeights:
    """Load one governed per-character weight table."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = [] if reader.fieldnames is None else list(reader.fieldnames)
        missing_columns = [
            column for column in ("character_id", "weight") if column not in fieldnames
        ]
        if missing_columns:
            raise ParsimonyAnalysisError(
                "parsimony character weights require character_id and weight columns",
                code="parsimony_character_weights_invalid_columns",
                details={
                    "weights_path": str(path),
                    "missing_columns": missing_columns,
                    "columns": fieldnames,
                },
            )
        weights_by_character: dict[str, float] = {}
        for line_number, row in enumerate(reader, start=2):
            character_id = (row.get("character_id") or "").strip()
            if not character_id:
                raise ParsimonyAnalysisError(
                    "parsimony character weights require one nonempty character_id per row",
                    code="parsimony_character_weight_missing_character_id",
                    details={
                        "weights_path": str(path),
                        "line_number": line_number,
                    },
                )
            if character_id in weights_by_character:
                raise ParsimonyAnalysisError(
                    "parsimony character weights require one unique row per character",
                    code="parsimony_character_weight_duplicate_character",
                    details={
                        "weights_path": str(path),
                        "character_id": character_id,
                        "line_number": line_number,
                    },
                )
            weights_by_character[character_id] = _parse_weight_value(
                row.get("weight"),
                weights_path=path,
                character_id=character_id,
                line_number=line_number,
            )
    if not weights_by_character:
        raise ParsimonyAnalysisError(
            "parsimony character weights require at least one weight row",
            code="parsimony_character_weights_empty",
            details={"weights_path": str(path)},
        )
    return ParsimonyCharacterWeights(
        weights_path=path,
        weights_by_character=weights_by_character,
    )


def resolve_parsimony_character_weights(
    character_ids: list[str],
    character_weights: ParsimonyCharacterWeights | Mapping[str, float] | Path | None,
) -> ParsimonyCharacterWeights:
    """Resolve explicit character weights for one matrix and validate coverage."""
    if character_weights is None:
        return ParsimonyCharacterWeights(
            weights_path=None,
            weights_by_character=dict.fromkeys(character_ids, 1.0),
        )
    if isinstance(character_weights, Path):
        resolved_weights = load_parsimony_character_weights(character_weights)
    elif isinstance(character_weights, ParsimonyCharacterWeights):
        resolved_weights = character_weights
    else:
        resolved_weights = ParsimonyCharacterWeights(
            weights_path=None,
            weights_by_character={
                character_id: _parse_weight_value(
                    weight,
                    weights_path=None,
                    character_id=character_id,
                    line_number=None,
                )
                for character_id, weight in character_weights.items()
            },
        )
    missing_characters = sorted(
        set(character_ids) - set(resolved_weights.weights_by_character)
    )
    extra_characters = sorted(
        set(resolved_weights.weights_by_character) - set(character_ids)
    )
    if missing_characters:
        raise ParsimonyAnalysisError(
            "parsimony character weights must cover every matrix character",
            code="parsimony_character_weight_missing_character",
            details={
                "weights_path": None
                if resolved_weights.weights_path is None
                else str(resolved_weights.weights_path),
                "missing_characters": missing_characters,
            },
        )
    if extra_characters:
        raise ParsimonyAnalysisError(
            "parsimony character weights cannot reference characters absent from the matrix",
            code="parsimony_character_weight_extra_character",
            details={
                "weights_path": None
                if resolved_weights.weights_path is None
                else str(resolved_weights.weights_path),
                "extra_characters": extra_characters,
            },
        )
    return ParsimonyCharacterWeights(
        weights_path=resolved_weights.weights_path,
        weights_by_character={
            character_id: resolved_weights.weights_by_character[character_id]
            for character_id in character_ids
        },
    )


def _parse_weight_value(
    value: object,
    *,
    weights_path: Path | None,
    character_id: str,
    line_number: int | None,
) -> float:
    try:
        numeric_value = float(value)
    except (TypeError, ValueError) as error:
        raise ParsimonyAnalysisError(
            "parsimony character weights require finite numeric values",
            code="parsimony_character_weight_invalid_value",
            details={
                "weights_path": None if weights_path is None else str(weights_path),
                "character_id": character_id,
                "line_number": line_number,
                "value": value,
            },
        ) from error
    if not math.isfinite(numeric_value):
        raise ParsimonyAnalysisError(
            "parsimony character weights require finite numeric values",
            code="parsimony_character_weight_invalid_value",
            details={
                "weights_path": None if weights_path is None else str(weights_path),
                "character_id": character_id,
                "line_number": line_number,
                "value": value,
            },
        )
    if numeric_value < 0:
        raise ParsimonyAnalysisError(
            "parsimony character weights cannot be negative",
            code="parsimony_character_weight_negative_value",
            details={
                "weights_path": None if weights_path is None else str(weights_path),
                "character_id": character_id,
                "line_number": line_number,
                "value": numeric_value,
            },
        )
    return numeric_value
