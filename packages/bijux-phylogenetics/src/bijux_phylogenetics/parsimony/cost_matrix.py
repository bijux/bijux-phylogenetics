from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .models import SankoffCostMatrix, SankoffCostMatrixWarning


def load_sankoff_cost_matrix(
    path: Path,
    *,
    observed_states: list[str] | None = None,
    allow_asymmetric_costs: bool = False,
) -> SankoffCostMatrix:
    """Load and validate one square Sankoff transition-cost matrix."""
    return validate_sankoff_cost_matrix(
        _load_sankoff_cost_matrix_rows(path),
        observed_states=observed_states,
        allow_asymmetric_costs=allow_asymmetric_costs,
    )


def validate_sankoff_cost_matrix(
    cost_matrix: SankoffCostMatrix | Path,
    *,
    observed_states: list[str] | None = None,
    allow_asymmetric_costs: bool = False,
) -> SankoffCostMatrix:
    """Validate one Sankoff cost matrix against symmetry policy and observed states."""
    resolved_cost_matrix = (
        _load_sankoff_cost_matrix_rows(cost_matrix)
        if isinstance(cost_matrix, Path)
        else cost_matrix
    )
    warnings: list[SankoffCostMatrixWarning] = []
    states = resolved_cost_matrix.states
    costs = resolved_cost_matrix.costs

    if not allow_asymmetric_costs:
        asymmetric_pairs: list[dict[str, object]] = []
        for index, row_state in enumerate(states):
            for column_state in states[index + 1 :]:
                forward_cost = costs[row_state][column_state]
                reverse_cost = costs[column_state][row_state]
                if not math.isclose(forward_cost, reverse_cost):
                    asymmetric_pairs.append(
                        {
                            "from_state": row_state,
                            "to_state": column_state,
                            "forward_cost": forward_cost,
                            "reverse_cost": reverse_cost,
                        }
                    )
        if asymmetric_pairs:
            raise ParsimonyAnalysisError(
                "sankoff scoring requires a symmetric cost matrix unless asymmetric costs are allowed explicitly",
                code="parsimony_cost_matrix_asymmetric",
                details={
                    "cost_matrix_path": None
                    if resolved_cost_matrix.matrix_path is None
                    else str(resolved_cost_matrix.matrix_path),
                    "asymmetric_pairs": asymmetric_pairs,
                },
            )

    diagonal_nonzero_states = sorted(
        state for state in states if not math.isclose(costs[state][state], 0.0)
    )
    if diagonal_nonzero_states:
        warnings.append(
            SankoffCostMatrixWarning(
                code="parsimony_cost_matrix_diagonal_nonzero",
                message="sankoff cost matrix diagonal costs are nonzero for at least one state",
                details={
                    "cost_matrix_path": None
                    if resolved_cost_matrix.matrix_path is None
                    else str(resolved_cost_matrix.matrix_path),
                    "states": diagonal_nonzero_states,
                    "diagonal_costs": {
                        state: costs[state][state] for state in diagonal_nonzero_states
                    },
                },
            )
        )

    if observed_states is not None:
        missing_states = sorted(set(observed_states) - set(states))
        if missing_states:
            raise ParsimonyAnalysisError(
                "sankoff scoring requires every observed state to exist in the cost matrix",
                code="parsimony_cost_matrix_missing_states",
                details={
                    "cost_matrix_path": None
                    if resolved_cost_matrix.matrix_path is None
                    else str(resolved_cost_matrix.matrix_path),
                    "missing_states": missing_states,
                    "cost_matrix_states": states,
                },
            )
        unused_states = sorted(set(states) - set(observed_states))
        if unused_states:
            warnings.append(
                SankoffCostMatrixWarning(
                    code="parsimony_cost_matrix_unused_states",
                    message="sankoff cost matrix includes states that are not observed in the current character matrix",
                    details={
                        "cost_matrix_path": None
                        if resolved_cost_matrix.matrix_path is None
                        else str(resolved_cost_matrix.matrix_path),
                        "unused_states": unused_states,
                        "observed_states": observed_states,
                    },
                )
            )

    return SankoffCostMatrix(
        matrix_path=resolved_cost_matrix.matrix_path,
        states=states,
        costs=costs,
        validation_warnings=warnings,
    )


def _load_sankoff_cost_matrix_rows(path: Path) -> SankoffCostMatrix:
    table = load_taxon_table(path, taxon_column="state")
    states = [column for column in table.columns if column != table.taxon_column]
    if not states:
        raise ParsimonyAnalysisError(
            "sankoff scoring requires at least one state in the cost matrix",
            code="parsimony_cost_matrix_not_square",
            details={"cost_matrix_path": str(path), "state_count": 0},
        )
    if table.row_count != len(states):
        raise ParsimonyAnalysisError(
            "sankoff scoring requires a square cost matrix",
            code="parsimony_cost_matrix_not_square",
            details={
                "cost_matrix_path": str(path),
                "row_count": table.row_count,
                "column_count": len(states),
            },
        )

    row_labels = [row[table.taxon_column] for row in table.rows]
    if set(row_labels) != set(states):
        raise ParsimonyAnalysisError(
            "sankoff scoring requires matching row and column state labels",
            code="parsimony_cost_matrix_inconsistent_labels",
            details={
                "cost_matrix_path": str(path),
                "row_labels": row_labels,
                "column_labels": states,
            },
        )

    costs: dict[str, dict[str, float]] = {}
    for row in table.rows:
        row_state = row[table.taxon_column]
        costs[row_state] = {}
        for column_state in states:
            raw_value = row[column_state]
            try:
                cost = float(raw_value)
            except ValueError as error:
                raise ParsimonyAnalysisError(
                    "sankoff scoring requires numeric transition costs",
                    code="parsimony_cost_matrix_invalid_value",
                    details={
                        "cost_matrix_path": str(path),
                        "row_state": row_state,
                        "column_state": column_state,
                        "value": raw_value,
                    },
                ) from error
            if not math.isfinite(cost):
                raise ParsimonyAnalysisError(
                    "sankoff scoring requires finite transition costs",
                    code="parsimony_cost_matrix_invalid_value",
                    details={
                        "cost_matrix_path": str(path),
                        "row_state": row_state,
                        "column_state": column_state,
                        "value": raw_value,
                    },
                )
            if cost < 0:
                raise ParsimonyAnalysisError(
                    "sankoff scoring does not accept negative transition costs",
                    code="parsimony_cost_matrix_negative_cost",
                    details={
                        "cost_matrix_path": str(path),
                        "row_state": row_state,
                        "column_state": column_state,
                        "cost": cost,
                    },
                )
            costs[row_state][column_state] = cost

    return SankoffCostMatrix(
        matrix_path=path,
        states=states,
        costs=costs,
        validation_warnings=[],
    )
