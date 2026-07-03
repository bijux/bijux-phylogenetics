from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    FiniteStateFrequencyOptimizationReport,
    StateFrequencyOptimizationRow,
    build_protein_equal_input_finite_state_surface,
    optimize_codon_ctmc_state_frequencies,
    optimize_codon_ctmc_state_frequencies_from_alignment,
    optimize_discrete_equal_input_state_frequencies,
    optimize_protein_equal_input_state_frequencies,
)


def test_public_likelihood_exports_finite_state_frequency_optimization_surface() -> (
    None
):
    assert (
        likelihood_api.FiniteStateFrequencyOptimizationReport
        is FiniteStateFrequencyOptimizationReport
    )
    assert likelihood_api.StateFrequencyOptimizationRow is StateFrequencyOptimizationRow
    assert (
        likelihood_api.build_protein_equal_input_finite_state_surface
        is build_protein_equal_input_finite_state_surface
    )
    assert (
        likelihood_api.optimize_codon_ctmc_state_frequencies
        is optimize_codon_ctmc_state_frequencies
    )
    assert (
        likelihood_api.optimize_codon_ctmc_state_frequencies_from_alignment
        is optimize_codon_ctmc_state_frequencies_from_alignment
    )
    assert (
        likelihood_api.optimize_discrete_equal_input_state_frequencies
        is optimize_discrete_equal_input_state_frequencies
    )
    assert (
        likelihood_api.optimize_protein_equal_input_state_frequencies
        is optimize_protein_equal_input_state_frequencies
    )
