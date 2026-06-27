from __future__ import annotations

import numpy
import pytest

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    solve_ctmc_stationary_distribution,
)
from bijux_phylogenetics.phylo.likelihood.f81 import f81_rate_matrix
from bijux_phylogenetics.phylo.likelihood.gtr import gtr_rate_matrix
from bijux_phylogenetics.phylo.likelihood.hky85 import hky85_rate_matrix
from bijux_phylogenetics.phylo.likelihood.jc69 import (
    jc69_rate_matrix,
    jc69_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.k80 import k80_rate_matrix
from bijux_phylogenetics.phylo.likelihood.nucleotide_finite_state import (
    build_selected_nucleotide_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_models import (
    resolve_selected_nucleotide_likelihood_specification,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    transition_probability_matrix,
)


def _nucleotide_records() -> list[AlignmentRecord]:
    return [
        AlignmentRecord(identifier="taxon-a", sequence="ACGTACGT"),
        AlignmentRecord(identifier="taxon-b", sequence="TGCATGCA"),
    ]


@pytest.mark.parametrize(
    ("model_name", "model_kwargs", "rate_matrix"),
    [
        ("jc69", {}, jc69_rate_matrix()),
        ("k80", {"kappa": 4.0}, k80_rate_matrix(4.0)),
        (
            "f81",
            {"base_frequencies": {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}},
            f81_rate_matrix(numpy.array([0.4, 0.1, 0.2, 0.3], dtype=float)),
        ),
        (
            "hky85",
            {
                "kappa": 4.0,
                "base_frequencies": {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
            },
            hky85_rate_matrix(
                numpy.array([0.4, 0.1, 0.2, 0.3], dtype=float),
                kappa=4.0,
            ),
        ),
        (
            "gtr",
            {
                "base_frequencies": {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
                "exchangeabilities": {
                    "AC": 1.0,
                    "AG": 4.5,
                    "AT": 0.8,
                    "CG": 1.6,
                    "CT": 2.4,
                    "GT": 3.1,
                },
            },
            gtr_rate_matrix(
                {
                    "AC": 1.0,
                    "AG": 4.5,
                    "AT": 0.8,
                    "CG": 1.6,
                    "CT": 2.4,
                    "GT": 3.1,
                },
                base_frequencies=numpy.array([0.4, 0.1, 0.2, 0.3], dtype=float),
            ),
        ),
    ],
)
def test_nucleotide_transition_matrices_stay_valid_from_tiny_to_long_branches(
    model_name: str,
    model_kwargs: dict[str, object],
    rate_matrix: numpy.ndarray,
) -> None:
    tiny_branch_length = 1e-9
    long_branch_length = 25.0
    stationary = numpy.array(
        solve_ctmc_stationary_distribution(
            rate_matrix,
            state_labels=("A", "C", "G", "T"),
        ).probabilities,
        dtype=float,
    )

    tiny_transition = transition_probability_matrix(rate_matrix, tiny_branch_length)
    long_transition = transition_probability_matrix(rate_matrix, long_branch_length)

    assert numpy.allclose(
        tiny_transition,
        numpy.eye(4, dtype=float) + (rate_matrix * tiny_branch_length),
        rtol=0.0,
        atol=1e-9,
    ), model_name
    assert numpy.all(tiny_transition >= 0.0), model_name
    assert numpy.all(long_transition >= 0.0), model_name
    assert numpy.allclose(
        tiny_transition.sum(axis=1),
        numpy.ones(4, dtype=float),
        rtol=0.0,
        atol=1e-12,
    ), model_name
    assert numpy.allclose(
        long_transition.sum(axis=1),
        numpy.ones(4, dtype=float),
        rtol=0.0,
        atol=1e-12,
    ), model_name
    assert numpy.allclose(
        long_transition,
        numpy.repeat(stationary[numpy.newaxis, :], 4, axis=0),
        rtol=0.0,
        atol=2e-8,
    ), model_name


@pytest.mark.parametrize(
    ("model_name", "model_kwargs"),
    [
        ("jc69", {}),
        ("k80", {"kappa": 4.0}),
        (
            "f81",
            {"base_frequencies": {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}},
        ),
        (
            "hky85",
            {
                "kappa": 4.0,
                "base_frequencies": {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
            },
        ),
        (
            "gtr",
            {
                "base_frequencies": {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
                "exchangeabilities": {
                    "AC": 1.0,
                    "AG": 4.5,
                    "AT": 0.8,
                    "CG": 1.6,
                    "CT": 2.4,
                    "GT": 3.1,
                },
            },
        ),
    ],
)
def test_nucleotide_finite_state_surface_uses_generic_transition_engine(
    model_name: str,
    model_kwargs: dict[str, object],
) -> None:
    records = _nucleotide_records()
    specification = resolve_selected_nucleotide_likelihood_specification(
        records,
        model_name=model_name,
        owner_name="nucleotide transition fixture",
        **model_kwargs,
    )
    surface = build_selected_nucleotide_finite_state_surface(
        records,
        model_name=model_name,
        owner_name="nucleotide transition fixture",
        **model_kwargs,
    )
    branch_length = 0.75

    assert numpy.allclose(
        surface.transition_matrix_for_branch_length(branch_length),
        transition_probability_matrix(specification.rate_matrix, branch_length),
        rtol=0.0,
        atol=1e-12,
    )


@pytest.mark.parametrize("branch_length", [0.0, 1e-9, 0.1, 0.75, 5.0, 25.0])
def test_jc69_closed_form_transition_matches_normalized_ctmc_engine(
    branch_length: float,
) -> None:
    closed_form = jc69_transition_probability_matrix(branch_length)
    matrix_exponential = transition_probability_matrix(
        jc69_rate_matrix(),
        branch_length,
    )

    assert numpy.allclose(
        closed_form,
        matrix_exponential,
        rtol=0.0,
        atol=1e-12,
    )
