from __future__ import annotations

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
from .weights import (
    load_parsimony_character_weights as _load_parsimony_character_weights,
)
from .weights import (
    resolve_parsimony_character_weights,
)

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
    return _load_parsimony_character_weights(path)


def tree_length(
    tree: Path,
    matrix: FitchCharacterMatrix | Path,
    *,
    method: str,
    state_order: list[str] | None = None,
    cost_matrix: SankoffCostMatrix | Path | None = None,
    allow_asymmetric_costs: bool = False,
    character_weights: ParsimonyCharacterWeights | Path | None = None,
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
        score_report = score_fitch(
            tree,
            resolved_matrix,
            character_weights=resolved_weights,
        )
        raw_scores = {
            row.character_id: float(row.step_count) for row in score_report.step_rows
        }
        cost_matrix_path = None
    elif resolved_method == "wagner":
        score_report = score_wagner(
            tree,
            resolved_matrix,
            state_order=state_order,
            character_weights=resolved_weights,
        )
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
            else load_sankoff_cost_matrix(
                cost_matrix,
                allow_asymmetric_costs=allow_asymmetric_costs,
            )
        )
        score_report = score_sankoff(
            tree,
            resolved_matrix,
            resolved_cost_matrix,
            allow_asymmetric_costs=allow_asymmetric_costs,
            character_weights=resolved_weights,
        )
        raw_scores = {
            row.character_id: float(row.minimum_cost) for row in score_report.step_rows
        }
        cost_matrix_path = score_report.cost_matrix_path
    elif resolved_method == "dollo":
        score_report = score_dollo(
            tree,
            resolved_matrix,
            character_weights=resolved_weights,
        )
        raw_scores = {
            row.character_id: float(row.step_count) for row in score_report.step_rows
        }
        cost_matrix_path = None
    elif resolved_method == "camin-sokal":
        score_report = score_camin_sokal(
            tree,
            resolved_matrix,
            character_weights=resolved_weights,
        )
        raw_scores = {
            row.character_id: float(row.gain_count) for row in score_report.step_rows
        }
        cost_matrix_path = None
    elif resolved_method == "acctran":
        score_report = reconstruct_acctran(
            tree,
            resolved_matrix,
            character_weights=resolved_weights,
        )
        raw_scores = {
            row.character_id: float(row.step_count) for row in score_report.step_rows
        }
        cost_matrix_path = None
    else:
        score_report = reconstruct_deltran(
            tree,
            resolved_matrix,
            character_weights=resolved_weights,
        )
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
    character_weights: ParsimonyCharacterWeights | Path | None,
) -> ParsimonyCharacterWeights:
    return resolve_parsimony_character_weights(
        character_ids,
        character_weights,
    )
