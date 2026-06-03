from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy

from bijux_phylogenetics.phylo.likelihood.dna import (
    DNA_EXCHANGEABILITY_ORDER,
    DNA_STATE_ORDER,
    validate_dna_exchangeabilities,
)
from bijux_phylogenetics.phylo.likelihood.simplex_coordinates import (
    SimplexCoordinateParameterization,
    parameterize_named_simplex,
    resolve_named_simplex_from_unconstrained,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

DNA_EXCHANGEABILITY_LABELS = tuple("".join(pair) for pair in DNA_EXCHANGEABILITY_ORDER)


def parameterize_dna_base_frequency_simplex(
    base_frequencies: Mapping[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...],
    *,
    normalization_tolerance: float = 1e-9,
) -> SimplexCoordinateParameterization:
    if isinstance(base_frequencies, Mapping):
        values_by_name = {
            str(state): float(probability)
            for state, probability in base_frequencies.items()
        }
    else:
        vector = numpy.asarray(base_frequencies, dtype=float)
        if vector.shape != (len(DNA_STATE_ORDER),):
            raise PhylogeneticsError(
                "DNA base-frequency simplex parameterization requires exactly four values in A/C/G/T order",
                code="dna_base_frequency_simplex_dimension_invalid",
            )
        values_by_name = dict(zip(DNA_STATE_ORDER, vector.tolist(), strict=True))
    return parameterize_named_simplex(
        values_by_name,
        expected_component_names=DNA_STATE_ORDER,
        owner_name="DNA base-frequency simplex parameterization",
        normalization_tolerance=normalization_tolerance,
    )


def resolve_dna_base_frequency_simplex_from_unconstrained(
    unconstrained_values: Sequence[float],
    *,
    normalization_tolerance: float = 1e-9,
) -> SimplexCoordinateParameterization:
    return resolve_named_simplex_from_unconstrained(
        unconstrained_values,
        component_names=DNA_STATE_ORDER,
        owner_name="DNA base-frequency simplex parameterization",
        normalization_tolerance=normalization_tolerance,
    )


def resolve_dna_base_frequencies_from_unconstrained(
    unconstrained_values: Sequence[float],
    *,
    normalization_tolerance: float = 1e-9,
) -> numpy.ndarray:
    parameterization = resolve_dna_base_frequency_simplex_from_unconstrained(
        unconstrained_values,
        normalization_tolerance=normalization_tolerance,
    )
    return numpy.array(parameterization.constrained_values, dtype=float)


def parameterize_dna_exchangeability_simplex(
    exchangeabilities: (
        Mapping[tuple[str, str], float]
        | Mapping[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
    ),
    *,
    normalization_tolerance: float = 1e-9,
) -> SimplexCoordinateParameterization:
    vector = validate_dna_exchangeabilities(
        exchangeabilities,
        model_name="DNA exchangeability simplex parameterization",
    )
    total = float(vector.sum())
    values_by_name = dict(
        zip(DNA_EXCHANGEABILITY_LABELS, (vector / total).tolist(), strict=True)
    )
    return parameterize_named_simplex(
        values_by_name,
        expected_component_names=DNA_EXCHANGEABILITY_LABELS,
        owner_name="DNA exchangeability simplex parameterization",
        normalization_tolerance=normalization_tolerance,
    )


def resolve_dna_exchangeability_simplex_from_unconstrained(
    unconstrained_values: Sequence[float],
    *,
    normalization_tolerance: float = 1e-9,
) -> SimplexCoordinateParameterization:
    return resolve_named_simplex_from_unconstrained(
        unconstrained_values,
        component_names=DNA_EXCHANGEABILITY_LABELS,
        owner_name="DNA exchangeability simplex parameterization",
        normalization_tolerance=normalization_tolerance,
    )


def resolve_anchor_normalized_dna_exchangeabilities_from_unconstrained(
    unconstrained_values: Sequence[float],
    *,
    normalization_tolerance: float = 1e-9,
    anchor_pair: tuple[str, str] = ("A", "C"),
) -> numpy.ndarray:
    parameterization = resolve_dna_exchangeability_simplex_from_unconstrained(
        unconstrained_values,
        normalization_tolerance=normalization_tolerance,
    )
    try:
        anchor_index = DNA_EXCHANGEABILITY_LABELS.index("".join(anchor_pair))
    except ValueError as error:
        raise ValueError(
            f"unsupported DNA exchangeability anchor {anchor_pair}"
        ) from error
    constrained_values = numpy.array(parameterization.constrained_values, dtype=float)
    anchor_value = float(constrained_values[anchor_index])
    if anchor_value <= 0.0:
        raise PhylogeneticsError(
            "DNA exchangeability simplex reconstruction requires a positive anchor component",
            code="dna_exchangeability_simplex_anchor_not_positive",
            details={
                "anchor_pair": "".join(anchor_pair),
                "anchor_value": anchor_value,
            },
        )
    return constrained_values / anchor_value


__all__ = [
    "DNA_EXCHANGEABILITY_LABELS",
    "parameterize_dna_base_frequency_simplex",
    "parameterize_dna_exchangeability_simplex",
    "resolve_anchor_normalized_dna_exchangeabilities_from_unconstrained",
    "resolve_dna_base_frequencies_from_unconstrained",
    "resolve_dna_base_frequency_simplex_from_unconstrained",
    "resolve_dna_exchangeability_simplex_from_unconstrained",
]
