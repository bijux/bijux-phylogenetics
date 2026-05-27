from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .models import SankoffCostMatrix


def load_sankoff_cost_matrix(path: Path) -> SankoffCostMatrix:
    """Load and validate one square Sankoff transition-cost matrix."""
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
    )
