from __future__ import annotations

import math

import numpy

from bijux_phylogenetics.phylo.likelihood.codon_states import (
    build_equal_rate_codon_ctmc_rate_matrix,
    resolve_codon_state_space,
    validate_codon_frequency_vector,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError


def test_standard_codon_state_space_excludes_stop_codons() -> None:
    state_space = resolve_codon_state_space()

    assert state_space.genetic_code_id == 1
    assert len(state_space.state_order) == 61
    assert "TAA" not in state_space.state_order
    assert "TAG" not in state_space.state_order
    assert "TGA" not in state_space.state_order
    assert state_space.state_order[0] == "AAA"
    assert state_space.state_order[-1] == "TTT"


def test_validate_codon_frequency_vector_normalizes_mapping() -> None:
    state_space = resolve_codon_state_space()
    frequencies, source = validate_codon_frequency_vector(
        {codon: 2.0 if codon == "AAA" else 1.0 for codon in state_space.state_order},
        state_space=state_space,
        owner_name="test codon likelihood",
    )

    assert source == "provided"
    assert math.isclose(float(frequencies.sum()), 1.0, rel_tol=0.0, abs_tol=1e-12)
    assert (
        frequencies[state_space.state_index["AAA"]]
        > frequencies[state_space.state_index["AAC"]]
    )


def test_validate_codon_frequency_vector_rejects_stop_codon_mapping() -> None:
    state_space = resolve_codon_state_space()
    with_stop = dict.fromkeys(state_space.state_order, 1.0)
    with_stop["TAA"] = 1.0

    try:
        validate_codon_frequency_vector(
            with_stop,
            state_space=state_space,
            owner_name="test codon likelihood",
        )
    except InvalidAlignmentError as error:
        assert "exactly the resolved sense-codon state order" in str(error)
    else:
        raise AssertionError("stop-codon mapping should be rejected")


def test_equal_rate_codon_ctmc_rate_matrix_is_normalized() -> None:
    state_space, rate_matrix, frequencies, source = (
        build_equal_rate_codon_ctmc_rate_matrix()
    )

    assert source == "uniform"
    assert rate_matrix.shape == (61, 61)
    assert numpy.allclose(
        rate_matrix.sum(axis=1),
        numpy.zeros(61, dtype=float),
        rtol=0.0,
        atol=1e-12,
    )
    expected_rate = -float(numpy.sum(frequencies * numpy.diag(rate_matrix)))
    assert state_space.genetic_code_name == "Standard"
    assert math.isclose(expected_rate, 1.0, rel_tol=0.0, abs_tol=1e-12)
