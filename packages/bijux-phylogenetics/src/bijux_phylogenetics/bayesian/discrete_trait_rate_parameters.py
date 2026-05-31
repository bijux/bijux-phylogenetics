from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math

from bijux_phylogenetics.ancestral.discrete import DiscreteTransitionRateRow
from bijux_phylogenetics.ancestral.discrete.policy import resolve_discrete_model_name
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

DISCRETE_TRAIT_RATE_PARAMETER_MODELS = (
    "equal-rates",
    "symmetric",
    "all-rates-different",
)
_RATE_GROUP_TOLERANCE = 1e-12


@dataclass(frozen=True, slots=True)
class DiscreteTraitRateParameterGroup:
    """One governed discrete-trait rate parameter grouped over transitions."""

    parameter_name: str
    transition_pairs: list[tuple[str, str]]
    rate_value: float


@dataclass(frozen=True, slots=True)
class DiscreteTraitRateParameterization:
    """One grouped Mk transition-rate parameterization."""

    model: str
    parameter_count: int
    parameter_values: dict[str, float]
    groups: list[DiscreteTraitRateParameterGroup]


def parameterize_discrete_trait_rate_rows(
    *,
    model: str,
    transition_rate_rows: Sequence[DiscreteTransitionRateRow],
) -> DiscreteTraitRateParameterization:
    """Group governed discrete-trait transition rows into named Mk parameters."""
    resolved_model = _resolve_supported_discrete_trait_rate_model(model)
    groups = [
        DiscreteTraitRateParameterGroup(
            parameter_name=parameter_name,
            transition_pairs=transition_pairs,
            rate_value=float(format(rate_value, ".15g")),
        )
        for parameter_name, transition_pairs, rate_value in _resolve_rate_parameter_groups(
            model=resolved_model,
            transition_rate_rows=transition_rate_rows,
        )
    ]
    return DiscreteTraitRateParameterization(
        model=resolved_model,
        parameter_count=len(groups),
        parameter_values={group.parameter_name: group.rate_value for group in groups},
        groups=groups,
    )


def resolve_discrete_trait_rate_rows(
    *,
    model: str,
    transition_rate_rows: Sequence[DiscreteTransitionRateRow],
    parameter_values: Mapping[str, float],
) -> list[DiscreteTransitionRateRow]:
    """Resolve one Mk parameterization back onto one transition-row surface."""
    parameterization = parameterize_discrete_trait_rate_rows(
        model=model,
        transition_rate_rows=transition_rate_rows,
    )
    expected_parameter_names = set(parameterization.parameter_values)
    provided_parameter_names = {
        parameter_name.strip() for parameter_name in parameter_values
    }
    if provided_parameter_names != expected_parameter_names:
        raise PhylogeneticsError(
            "discrete-trait rate parameter resolution requires exactly the grouped Mk parameter names",
            code="discrete_trait_rate_parameter_names_mismatch",
            details={
                "expected_parameter_names": sorted(expected_parameter_names),
                "provided_parameter_names": sorted(provided_parameter_names),
            },
        )
    validated_parameter_values = {
        parameter_name: _validate_positive_finite_value(
            parameter_name=parameter_name,
            value=parameter_value,
            owner_name="discrete-trait rate parameter resolution",
        )
        for parameter_name, parameter_value in parameter_values.items()
    }
    resolved_rate_by_transition_pair = {
        transition_pair: validated_parameter_values[group.parameter_name]
        for group in parameterization.groups
        for transition_pair in group.transition_pairs
    }
    resolved_rows: list[DiscreteTransitionRateRow] = []
    for row in transition_rate_rows:
        transition_pair = (row.source_state, row.target_state)
        resolved_rows.append(
            DiscreteTransitionRateRow(
                source_state=row.source_state,
                target_state=row.target_state,
                transition_allowed=row.transition_allowed,
                step_distance=row.step_distance,
                rate=(
                    float(
                        format(
                            resolved_rate_by_transition_pair[transition_pair], ".15g"
                        )
                    )
                    if row.transition_allowed
                    else row.rate
                ),
            )
        )
    return resolved_rows


def _resolve_supported_discrete_trait_rate_model(model: str) -> str:
    resolved_model = resolve_discrete_model_name(model)
    if resolved_model not in DISCRETE_TRAIT_RATE_PARAMETER_MODELS:
        raise PhylogeneticsError(
            "discrete-trait rate parameters support only ER, SYM, and ARD transition-rate surfaces",
            code="discrete_trait_rate_parameter_model_unsupported",
            details={
                "model": model,
                "resolved_model": resolved_model,
                "allowed_models": list(DISCRETE_TRAIT_RATE_PARAMETER_MODELS),
            },
        )
    return resolved_model


def _resolve_rate_parameter_groups(
    *,
    model: str,
    transition_rate_rows: Sequence[DiscreteTransitionRateRow],
) -> list[tuple[str, list[tuple[str, str]], float]]:
    allowed_rows = [row for row in transition_rate_rows if row.transition_allowed]
    if not allowed_rows:
        raise PhylogeneticsError(
            "discrete-trait rate parameterization requires at least one allowed transition rate",
            code="discrete_trait_rate_parameter_no_allowed_rates",
        )
    validated_rows = [
        (
            row.source_state,
            row.target_state,
            _validate_positive_finite_value(
                parameter_name=f"{row.source_state}->{row.target_state} rate",
                value=row.rate,
                owner_name="discrete-trait rate parameterization",
            ),
        )
        for row in allowed_rows
    ]
    if model == "equal-rates":
        shared_rate = validated_rows[0][2]
        if any(
            not _rates_close(rate_value, shared_rate)
            for _source_state, _target_state, rate_value in validated_rows[1:]
        ):
            raise PhylogeneticsError(
                "equal-rates parameterization requires one shared transition rate",
                code="discrete_trait_rate_parameter_equal_rates_inconsistent",
            )
        return [
            (
                "shared-rate",
                sorted(
                    (source_state, target_state)
                    for source_state, target_state, _rate_value in validated_rows
                ),
                shared_rate,
            )
        ]
    if model == "symmetric":
        pair_groups: dict[tuple[str, str], list[tuple[str, str, float]]] = {}
        for source_state, target_state, rate_value in validated_rows:
            pair_groups.setdefault(
                tuple(sorted((source_state, target_state))),
                [],
            ).append((source_state, target_state, rate_value))
        resolved_groups: list[tuple[str, list[tuple[str, str]], float]] = []
        for pair_key in sorted(pair_groups):
            grouped_rows = pair_groups[pair_key]
            transition_pairs = sorted(
                (source_state, target_state)
                for source_state, target_state, _rate_value in grouped_rows
            )
            if len(grouped_rows) != 2:
                raise PhylogeneticsError(
                    "symmetric parameterization requires one forward and one reverse transition per rate parameter",
                    code="discrete_trait_rate_parameter_symmetric_pair_incomplete",
                    details={"pair": list(pair_key)},
                )
            first_rate = grouped_rows[0][2]
            if any(
                not _rates_close(rate_value, first_rate)
                for _source_state, _target_state, rate_value in grouped_rows[1:]
            ):
                raise PhylogeneticsError(
                    "symmetric parameterization requires matched forward and reverse transition rates",
                    code="discrete_trait_rate_parameter_symmetric_pair_mismatched",
                    details={"pair": list(pair_key)},
                )
            resolved_groups.append(
                (
                    f"{pair_key[0]}<->{pair_key[1]}",
                    transition_pairs,
                    first_rate,
                )
            )
        return resolved_groups
    return [
        (
            f"{source_state}->{target_state}",
            [(source_state, target_state)],
            rate_value,
        )
        for source_state, target_state, rate_value in sorted(validated_rows)
    ]


def _rates_close(left: float, right: float) -> bool:
    return math.isclose(
        left,
        right,
        rel_tol=_RATE_GROUP_TOLERANCE,
        abs_tol=_RATE_GROUP_TOLERANCE,
    )


def _validate_positive_finite_value(
    *,
    parameter_name: str,
    value: float,
    owner_name: str,
) -> float:
    validated_value = float(value)
    if not math.isfinite(validated_value):
        raise PhylogeneticsError(
            f"{owner_name} requires finite {parameter_name}",
            code="discrete_trait_rate_parameter_nonfinite",
            details={
                "parameter_name": parameter_name,
                "parameter_value": value,
                "owner_name": owner_name,
            },
        )
    if validated_value <= 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires {parameter_name} > 0",
            code="discrete_trait_rate_parameter_nonpositive",
            details={
                "parameter_name": parameter_name,
                "parameter_value": value,
                "owner_name": owner_name,
            },
        )
    return validated_value
