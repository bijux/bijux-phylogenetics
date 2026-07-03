from __future__ import annotations

from functools import cache
from pathlib import Path

import numpy

from bijux_phylogenetics.phylo.likelihood.ctmc import (
    compute_ctmc_expected_substitution_rate,
    normalize_ctmc_rate_matrix_by_expected_substitution_rate,
)
from bijux_phylogenetics.phylo.likelihood.protein import PROTEIN_STATE_ORDER

_RESOURCE_ROOT = (
    Path(__file__).resolve().parent / "resources" / "protein_empirical_models"
)
_EMPIRICAL_MODEL_RESOURCE_BY_NAME = {
    "JTT": "jtt.paml.txt",
    "WAG": "wag.paml.txt",
    "LG": "lg.paml.txt",
}
_EMPIRICAL_MODEL_SOURCE_URL_BY_NAME = {
    "JTT": "https://raw.githubusercontent.com/abacus-gene/paml/master/dat/jones.dat",
    "WAG": "https://www.ebi.ac.uk/goldman-srv/WAG/wag.dat",
    "LG": "http://www.atgc-montpellier.fr/download/datasets/models/lg_LG.PAML.txt",
}
_EXPECTED_TRIANGLE_VALUE_COUNT = (
    len(PROTEIN_STATE_ORDER) * (len(PROTEIN_STATE_ORDER) - 1)
) // 2
_EXPECTED_PAML_VALUE_COUNT = _EXPECTED_TRIANGLE_VALUE_COUNT + len(PROTEIN_STATE_ORDER)


def list_supported_empirical_protein_model_names() -> tuple[str, ...]:
    """Return the built-in empirical protein models owned by the native runtime."""
    return tuple(_EMPIRICAL_MODEL_RESOURCE_BY_NAME)


def validate_empirical_protein_model_name(model_name: str) -> str:
    """Normalize one built-in empirical protein model name."""
    normalized_model_name = model_name.strip().upper()
    if normalized_model_name not in _EMPIRICAL_MODEL_RESOURCE_BY_NAME:
        raise ValueError(
            "empirical protein model name must be one of "
            f"{sorted(_EMPIRICAL_MODEL_RESOURCE_BY_NAME)}"
        )
    return normalized_model_name


def empirical_protein_model_source_url(model_name: str) -> str:
    """Return the source URL used to govern one built-in empirical protein model."""
    normalized_model_name = validate_empirical_protein_model_name(model_name)
    return _EMPIRICAL_MODEL_SOURCE_URL_BY_NAME[normalized_model_name]


def empirical_protein_exchangeability_matrix(model_name: str) -> numpy.ndarray:
    """Return one symmetric empirical protein exchangeability matrix."""
    normalized_model_name = validate_empirical_protein_model_name(model_name)
    return numpy.array(
        _load_empirical_protein_model_payload(normalized_model_name)[
            "exchangeability_matrix"
        ],
        copy=True,
    )


def empirical_protein_stationary_frequencies(model_name: str) -> numpy.ndarray:
    """Return one empirical protein stationary-frequency vector."""
    normalized_model_name = validate_empirical_protein_model_name(model_name)
    return numpy.array(
        _load_empirical_protein_model_payload(normalized_model_name)[
            "stationary_frequencies"
        ],
        copy=True,
    )


def empirical_protein_rate_matrix(model_name: str) -> numpy.ndarray:
    """Return one normalized empirical protein rate matrix."""
    normalized_model_name = validate_empirical_protein_model_name(model_name)
    return numpy.array(
        _load_empirical_protein_model_payload(normalized_model_name)["rate_matrix"],
        copy=True,
    )


@cache
def _load_empirical_protein_model_payload(model_name: str) -> dict[str, numpy.ndarray]:
    numeric_values = _parse_paml_numeric_values(
        _RESOURCE_ROOT / _EMPIRICAL_MODEL_RESOURCE_BY_NAME[model_name]
    )
    exchangeability_values = numeric_values[:_EXPECTED_TRIANGLE_VALUE_COUNT]
    stationary_frequencies = numpy.asarray(
        numeric_values[_EXPECTED_TRIANGLE_VALUE_COUNT:],
        dtype=float,
    )
    stationary_frequencies /= float(stationary_frequencies.sum())
    exchangeability_matrix = _build_symmetric_exchangeability_matrix(
        exchangeability_values
    )
    rate_matrix = _build_empirical_rate_matrix(
        exchangeability_matrix=exchangeability_matrix,
        stationary_frequencies=stationary_frequencies,
    )
    exchangeability_matrix.setflags(write=False)
    stationary_frequencies.setflags(write=False)
    rate_matrix.setflags(write=False)
    return {
        "exchangeability_matrix": exchangeability_matrix,
        "stationary_frequencies": stationary_frequencies,
        "rate_matrix": rate_matrix,
    }


def _parse_paml_numeric_values(path: Path) -> tuple[float, ...]:
    values: list[float] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        for token in line.split():
            values.append(float(token))
            if len(values) == _EXPECTED_PAML_VALUE_COUNT:
                return tuple(values)
    raise ValueError(
        f"empirical protein model resource '{path.name}' is missing required PAML values"
    )


def _build_symmetric_exchangeability_matrix(
    lower_triangle_values: tuple[float, ...],
) -> numpy.ndarray:
    state_count = len(PROTEIN_STATE_ORDER)
    exchangeability_matrix = numpy.zeros((state_count, state_count), dtype=float)
    value_index = 0
    for row_index in range(1, state_count):
        for column_index in range(row_index):
            value = lower_triangle_values[value_index]
            exchangeability_matrix[row_index, column_index] = value
            exchangeability_matrix[column_index, row_index] = value
            value_index += 1
    return exchangeability_matrix


def _build_empirical_rate_matrix(
    *,
    exchangeability_matrix: numpy.ndarray,
    stationary_frequencies: numpy.ndarray,
) -> numpy.ndarray:
    rate_matrix = numpy.zeros_like(exchangeability_matrix, dtype=float)
    state_count = len(PROTEIN_STATE_ORDER)
    for row_index in range(state_count):
        for column_index in range(state_count):
            if row_index == column_index:
                continue
            rate_matrix[row_index, column_index] = (
                exchangeability_matrix[row_index, column_index]
                * stationary_frequencies[column_index]
            )
    for row_index in range(state_count):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    normalized_rate_matrix = normalize_ctmc_rate_matrix_by_expected_substitution_rate(
        rate_matrix,
        stationary_frequencies,
        state_labels=PROTEIN_STATE_ORDER,
    )
    expected_rate = compute_ctmc_expected_substitution_rate(
        normalized_rate_matrix,
        stationary_frequencies,
        state_labels=PROTEIN_STATE_ORDER,
    )
    if not numpy.isclose(expected_rate, 1.0, rtol=0.0, atol=1e-12):
        raise ValueError(
            "empirical protein model normalization failed to reach rate one"
        )
    return normalized_rate_matrix
