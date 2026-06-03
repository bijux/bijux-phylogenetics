from __future__ import annotations

import math

import numpy

from bijux_phylogenetics.phylo.likelihood.protein import (
    PROTEIN_GAP_CHARACTER,
    PROTEIN_MISSING_CHARACTER,
    PROTEIN_STATE_INDEX,
    PROTEIN_STATE_ORDER,
    validate_protein_gap_policy,
    validate_protein_missing_policy,
    validate_protein_root_prior,
)


def validate_invariant_proportion(
    invariant_proportion: float,
    *,
    model_name: str,
) -> float:
    if not math.isfinite(invariant_proportion):
        raise ValueError(f"{model_name} invariant proportion must be finite")
    if invariant_proportion < 0.0 or invariant_proportion >= 1.0:
        raise ValueError(
            f"{model_name} invariant proportion must be in the half-open interval [0, 1)"
        )
    return float(invariant_proportion)


def validate_invariant_proportion_bounds(
    *,
    initial_invariant_proportion: float,
    lower_invariant_proportion_bound: float,
    upper_invariant_proportion_bound: float,
    model_name: str,
) -> tuple[float, float, float]:
    validated_lower_bound = validate_invariant_proportion(
        lower_invariant_proportion_bound,
        model_name=model_name,
    )
    validated_upper_bound = validate_invariant_proportion(
        upper_invariant_proportion_bound,
        model_name=model_name,
    )
    if validated_upper_bound <= validated_lower_bound:
        raise ValueError(
            f"{model_name} invariant proportion bounds must be strictly increasing"
        )
    validated_initial_proportion = validate_invariant_proportion(
        initial_invariant_proportion,
        model_name=model_name,
    )
    if not (
        validated_lower_bound <= validated_initial_proportion <= validated_upper_bound
    ):
        raise ValueError(
            f"{model_name} initial invariant proportion must lie within the requested bounds"
        )
    return (
        validated_initial_proportion,
        validated_lower_bound,
        validated_upper_bound,
    )


def invariant_component_site_likelihood(
    states: tuple[str, ...],
    *,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...],
    model_name: str,
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> float:
    """Return the invariant-site component likelihood for one protein site."""
    validated_root_prior = validate_protein_root_prior(
        root_prior,
        model_name=model_name,
    )
    validated_gap_policy = validate_protein_gap_policy(
        gap_policy,
        model_name=model_name,
    )
    validated_missing_policy = validate_protein_missing_policy(
        missing_policy,
        model_name=model_name,
    )
    compatible_state_indices: set[int] = set(range(len(PROTEIN_STATE_ORDER)))
    for state in states:
        if state in PROTEIN_STATE_INDEX:
            compatible_state_indices &= {PROTEIN_STATE_INDEX[state]}
            continue
        if state == PROTEIN_GAP_CHARACTER:
            if validated_gap_policy == "reject":
                raise ValueError(
                    f"{model_name} invariant component gap policy rejects '-'"
                )
            continue
        if state == PROTEIN_MISSING_CHARACTER:
            if validated_missing_policy == "reject":
                raise ValueError(
                    f"{model_name} invariant component missing-state policy rejects '?'"
                )
            continue
        raise ValueError(
            f"{model_name} invariant component encountered unsupported protein state '{state}'"
        )
    if not compatible_state_indices:
        return 0.0
    return float(
        sum(
            validated_root_prior[state_index]
            for state_index in compatible_state_indices
        )
    )


def invariant_mixture_site_likelihood(
    *,
    invariant_proportion: float,
    invariant_component_likelihood: float,
    variable_component_likelihood: float,
    model_name: str,
) -> float:
    validated_invariant_proportion = validate_invariant_proportion(
        invariant_proportion,
        model_name=model_name,
    )
    return validated_invariant_proportion * invariant_component_likelihood + (
        (1.0 - validated_invariant_proportion) * variable_component_likelihood
    )


def invariant_proportion_boundary_warnings(
    *,
    invariant_proportion: float,
    lower_invariant_proportion_bound: float,
    upper_invariant_proportion_bound: float,
    boundary_tolerance: float = 1e-9,
) -> tuple[bool, bool, list[str]]:
    hit_lower_boundary = math.isclose(
        invariant_proportion,
        lower_invariant_proportion_bound,
        rel_tol=0.0,
        abs_tol=boundary_tolerance,
    )
    hit_upper_boundary = math.isclose(
        invariant_proportion,
        upper_invariant_proportion_bound,
        rel_tol=0.0,
        abs_tol=boundary_tolerance,
    )
    warnings: list[str] = []
    if hit_lower_boundary:
        warnings.append("invariant proportion hit lower search boundary")
    if hit_upper_boundary:
        warnings.append("invariant proportion hit upper search boundary")
    return hit_lower_boundary, hit_upper_boundary, warnings
