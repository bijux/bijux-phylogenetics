from __future__ import annotations

import csv
import math
from collections.abc import Mapping
from pathlib import Path

from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .camin_sokal import score_camin_sokal
from .cost_matrix import load_sankoff_cost_matrix
from .dollo import score_dollo
from .fitch import score_fitch
from .matrix import load_parsimony_character_matrix
from .models import (
    FitchCharacterMatrix,
    ParsimonyCharacterWeights,
    ParsimonyTreeLengthCharacterScore,
    ParsimonyTreeLengthReport,
    SankoffCostMatrix,
)
from .reconstruction import reconstruct_acctran, reconstruct_deltran
from .sankoff import score_sankoff
from .wagner import score_wagner

_TREE_LENGTH_METHODS = frozenset(
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


def tree_length(
    tree: Path,
    matrix: FitchCharacterMatrix | Path,
    *,
    method: str,
    state_order: list[str] | None = None,
    cost_matrix: SankoffCostMatrix | Path | None = None,
    character_weights: ParsimonyCharacterWeights | Mapping[str, float] | Path | None = None,
) -> ParsimonyTreeLengthReport:
    """Score one tree length surface with optional explicit per-character weights."""
    resolved_method = _resolve_method(method)
    resolved_matrix = (
        matrix
        if isinstance(matrix, FitchCharacterMatrix)
        else load_parsimony_character_matrix(matrix)
    )
    resolved_weights = _resolve_character_weights(
        resolved_matrix.character_ids,
        character_weights,
    )
    if resolved_method == "fitch":
        score_report = score_fitch(tree, resolved_matrix)
        raw_scores = {
            row.character_id: float(row.step_count) for row in score_report.step_rows
        }
        cost_matrix_path = None
    elif resolved_method == "wagner":
        score_report = score_wagner(tree, resolved_matrix, state_order=state_order)
        raw_scores = {
            row.character_id: float(row.weighted_step_count)
            for row in score_report.step_rows
        }
        cost_matrix_path = None
    elif resolved_method == "sankoff":
        if cost_matrix is None:
            raise ParsimonyAnalysisError(
                "sankoff tree-length scoring requires one cost matrix",
                code="parsimony_tree_length_cost_matrix_required",
                details={"method": resolved_method},
            )
        resolved_cost_matrix = (
            cost_matrix
            if isinstance(cost_matrix, SankoffCostMatrix)
            else load_sankoff_cost_matrix(cost_matrix)
        )
        score_report = score_sankoff(tree, resolved_matrix, resolved_cost_matrix)
        raw_scores = {
            row.character_id: float(row.minimum_cost) for row in score_report.step_rows
        }
        cost_matrix_path = score_report.cost_matrix_path
    elif resolved_method == "dollo":
        score_report = score_dollo(tree, resolved_matrix)
        raw_scores = {
            row.character_id: float((1 if row.derived_taxon_count > 0 else 0) + row.total_losses)
            for row in score_report.step_rows
        }
        cost_matrix_path = None
    elif resolved_method == "camin-sokal":
        score_report = score_camin_sokal(tree, resolved_matrix)
        raw_scores = {
            row.character_id: float(row.gain_count) for row in score_report.step_rows
        }
        cost_matrix_path = None
    elif resolved_method == "acctran":
        score_report = reconstruct_acctran(tree, resolved_matrix)
        raw_scores = {
            row.character_id: float(row.step_count) for row in score_report.step_rows
        }
        cost_matrix_path = None
    else:
        score_report = reconstruct_deltran(tree, resolved_matrix)
        raw_scores = {
            row.character_id: float(row.step_count) for row in score_report.step_rows
        }
        cost_matrix_path = None

    step_rows: list[ParsimonyTreeLengthCharacterScore] = []
    raw_total_score = 0.0
    total_score = 0.0
    for character_id in resolved_matrix.character_ids:
        raw_score = raw_scores[character_id]
        weight = resolved_weights.weights_by_character[character_id]
        weighted_score = raw_score * weight
        raw_total_score += raw_score
        total_score += weighted_score
        step_rows.append(
            ParsimonyTreeLengthCharacterScore(
                character_id=character_id,
                raw_score=raw_score,
                character_weight=weight,
                weighted_score=weighted_score,
            )
        )
    return ParsimonyTreeLengthReport(
        algorithm="parsimony-tree-length",
        method=resolved_method,
        tree_path=score_report.tree_path,
        matrix_path=score_report.matrix_path,
        cost_matrix_path=cost_matrix_path,
        weights_path=resolved_weights.weights_path,
        taxon_column=score_report.taxon_column,
        taxon_count=score_report.taxon_count,
        character_count=score_report.character_count,
        raw_total_score=raw_total_score,
        total_score=total_score,
        step_rows=step_rows,
    )


def _resolve_method(method: str) -> str:
    resolved_method = method.strip().lower()
    if resolved_method not in _TREE_LENGTH_METHODS:
        raise ParsimonyAnalysisError(
            "tree-length scoring requires one supported parsimony method",
            code="parsimony_tree_length_method_unknown",
            details={
                "method": method,
                "supported_methods": sorted(_TREE_LENGTH_METHODS),
            },
        )
    return resolved_method


def _resolve_character_weights(
    character_ids: list[str],
    character_weights: ParsimonyCharacterWeights | Mapping[str, float] | Path | None,
) -> ParsimonyCharacterWeights:
    if character_weights is None:
        return ParsimonyCharacterWeights(
            weights_path=None,
            weights_by_character={character_id: 1.0 for character_id in character_ids},
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
