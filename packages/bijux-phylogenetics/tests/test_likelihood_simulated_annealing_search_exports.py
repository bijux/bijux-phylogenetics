from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    NucleotideLikelihoodSimulatedAnnealingSearchReport,
    NucleotideLikelihoodSimulatedAnnealingTraceRow,
    search_nucleotide_likelihood_simulated_annealing,
    search_nucleotide_likelihood_simulated_annealing_from_alignment,
    validate_nucleotide_likelihood_annealing_branch_reoptimization_policy,
    validate_nucleotide_likelihood_annealing_cooling_rate,
    validate_nucleotide_likelihood_annealing_initial_temperature,
    validate_nucleotide_likelihood_annealing_iteration_count,
    validate_nucleotide_likelihood_annealing_move_family,
    write_nucleotide_likelihood_simulated_annealing_artifacts,
    write_nucleotide_likelihood_simulated_annealing_run_json,
    write_nucleotide_likelihood_simulated_annealing_trace_table,
)


def test_public_likelihood_exports_simulated_annealing_surface() -> None:
    assert (
        likelihood_api.NucleotideLikelihoodSimulatedAnnealingSearchReport
        is NucleotideLikelihoodSimulatedAnnealingSearchReport
    )
    assert (
        likelihood_api.NucleotideLikelihoodSimulatedAnnealingTraceRow
        is NucleotideLikelihoodSimulatedAnnealingTraceRow
    )
    assert (
        likelihood_api.search_nucleotide_likelihood_simulated_annealing
        is search_nucleotide_likelihood_simulated_annealing
    )
    assert (
        likelihood_api.search_nucleotide_likelihood_simulated_annealing_from_alignment
        is search_nucleotide_likelihood_simulated_annealing_from_alignment
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_annealing_branch_reoptimization_policy
        is validate_nucleotide_likelihood_annealing_branch_reoptimization_policy
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_annealing_cooling_rate
        is validate_nucleotide_likelihood_annealing_cooling_rate
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_annealing_initial_temperature
        is validate_nucleotide_likelihood_annealing_initial_temperature
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_annealing_iteration_count
        is validate_nucleotide_likelihood_annealing_iteration_count
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_annealing_move_family
        is validate_nucleotide_likelihood_annealing_move_family
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_simulated_annealing_artifacts
        is write_nucleotide_likelihood_simulated_annealing_artifacts
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_simulated_annealing_run_json
        is write_nucleotide_likelihood_simulated_annealing_run_json
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_simulated_annealing_trace_table
        is write_nucleotide_likelihood_simulated_annealing_trace_table
    )
