from __future__ import annotations

import math

from bijux_phylogenetics.phylo.likelihood.models import (
    BranchLengthOptimizationRow,
    LikelihoodOptimizationBoundaryWarning,
    SubstitutionParameterOptimizationRow,
)
from bijux_phylogenetics.phylo.topology.clades import canonical_clade_id

_BOUNDARY_TOLERANCE = 1e-9
_NEAR_ZERO_BRANCH_LENGTH = 1e-8


def build_substitution_parameter_boundary_warnings(
    parameter_rows: list[SubstitutionParameterOptimizationRow],
) -> list[LikelihoodOptimizationBoundaryWarning]:
    """Render explicit boundary warnings from optimized substitution parameters."""
    warnings: list[LikelihoodOptimizationBoundaryWarning] = []
    for row in parameter_rows:
        if row.hit_lower_bound:
            warnings.append(
                LikelihoodOptimizationBoundaryWarning(
                    warning_kind="parameter-boundary",
                    affected_parameter=row.parameter_name,
                    boundary_side="lower",
                    observed_value=row.optimized_value,
                    lower_bound=row.lower_bound,
                    upper_bound=row.upper_bound,
                    affected_branch_id=None,
                    affected_branch_clade_id=None,
                    message=f"{row.parameter_name} hit lower search boundary",
                )
            )
        if row.hit_upper_bound:
            warnings.append(
                LikelihoodOptimizationBoundaryWarning(
                    warning_kind="parameter-boundary",
                    affected_parameter=row.parameter_name,
                    boundary_side="upper",
                    observed_value=row.optimized_value,
                    lower_bound=row.lower_bound,
                    upper_bound=row.upper_bound,
                    affected_branch_id=None,
                    affected_branch_clade_id=None,
                    message=f"{row.parameter_name} hit upper search boundary",
                )
            )
    return warnings


def build_base_frequency_boundary_warnings(
    *,
    base_frequency_source: str | None,
    base_frequency_a: float | None,
    base_frequency_c: float | None,
    base_frequency_g: float | None,
    base_frequency_t: float | None,
) -> list[LikelihoodOptimizationBoundaryWarning]:
    """Render explicit lower-bound frequency warnings for estimated DNA frequencies."""
    if base_frequency_source != "estimated":
        return []
    warnings: list[LikelihoodOptimizationBoundaryWarning] = []
    for parameter_name, value in (
        ("base_frequency_a", base_frequency_a),
        ("base_frequency_c", base_frequency_c),
        ("base_frequency_g", base_frequency_g),
        ("base_frequency_t", base_frequency_t),
    ):
        if value is None:
            continue
        if math.isclose(value, 0.0, rel_tol=0.0, abs_tol=_BOUNDARY_TOLERANCE):
            warnings.append(
                LikelihoodOptimizationBoundaryWarning(
                    warning_kind="frequency-boundary",
                    affected_parameter=parameter_name,
                    boundary_side="lower",
                    observed_value=value,
                    lower_bound=0.0,
                    upper_bound=1.0,
                    affected_branch_id=None,
                    affected_branch_clade_id=None,
                    message=f"{parameter_name} hit lower frequency boundary",
                )
            )
    return warnings


def build_branch_length_boundary_warnings(
    branch_rows: list[BranchLengthOptimizationRow],
    *,
    lower_branch_length_bound: float,
    upper_branch_length_bound: float,
) -> list[LikelihoodOptimizationBoundaryWarning]:
    """Render explicit lower-, upper-, and near-zero branch boundary warnings."""
    warnings: list[LikelihoodOptimizationBoundaryWarning] = []
    for row in branch_rows:
        clade_id = canonical_clade_id(frozenset(row.descendant_taxa))
        if math.isclose(
            row.optimized_branch_length,
            lower_branch_length_bound,
            rel_tol=0.0,
            abs_tol=_BOUNDARY_TOLERANCE,
        ):
            warnings.append(
                LikelihoodOptimizationBoundaryWarning(
                    warning_kind="branch-boundary",
                    affected_parameter="branch_length",
                    boundary_side="lower",
                    observed_value=row.optimized_branch_length,
                    lower_bound=lower_branch_length_bound,
                    upper_bound=upper_branch_length_bound,
                    affected_branch_id=row.branch_id,
                    affected_branch_clade_id=clade_id,
                    message=f"branch {clade_id} hit lower search boundary",
                )
            )
            continue
        if math.isclose(
            row.optimized_branch_length,
            upper_branch_length_bound,
            rel_tol=0.0,
            abs_tol=_BOUNDARY_TOLERANCE,
        ):
            warnings.append(
                LikelihoodOptimizationBoundaryWarning(
                    warning_kind="branch-boundary",
                    affected_parameter="branch_length",
                    boundary_side="upper",
                    observed_value=row.optimized_branch_length,
                    lower_bound=lower_branch_length_bound,
                    upper_bound=upper_branch_length_bound,
                    affected_branch_id=row.branch_id,
                    affected_branch_clade_id=clade_id,
                    message=f"branch {clade_id} hit upper search boundary",
                )
            )
            continue
        if row.optimized_branch_length < _NEAR_ZERO_BRANCH_LENGTH:
            warnings.append(
                LikelihoodOptimizationBoundaryWarning(
                    warning_kind="near-zero-branch",
                    affected_parameter="branch_length",
                    boundary_side="near-zero",
                    observed_value=row.optimized_branch_length,
                    lower_bound=lower_branch_length_bound,
                    upper_bound=upper_branch_length_bound,
                    affected_branch_id=row.branch_id,
                    affected_branch_clade_id=clade_id,
                    message=f"branch {clade_id} is near zero length",
                )
            )
    return warnings


def boundary_warning_messages(
    warnings: list[LikelihoodOptimizationBoundaryWarning],
) -> list[str]:
    """Project explicit boundary warnings to a deterministic message list."""
    return [warning.message for warning in warnings]


__all__ = [
    "LikelihoodOptimizationBoundaryWarning",
    "boundary_warning_messages",
    "build_base_frequency_boundary_warnings",
    "build_branch_length_boundary_warnings",
    "build_substitution_parameter_boundary_warnings",
]
