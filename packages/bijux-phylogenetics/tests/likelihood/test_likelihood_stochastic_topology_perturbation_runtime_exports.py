from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    NucleotideLikelihoodStochasticTopologyPerturbationSearchReport,
    NucleotideLikelihoodTopologyPerturbationStep,
    search_nucleotide_likelihood_stochastic_topology_perturbation,
    search_nucleotide_likelihood_stochastic_topology_perturbation_from_alignment,
    validate_nucleotide_likelihood_local_search_method,
    validate_nucleotide_likelihood_perturbation_branch_reoptimization_policy,
    validate_nucleotide_likelihood_perturbation_move_count,
    validate_nucleotide_likelihood_perturbation_move_family,
    write_nucleotide_likelihood_stochastic_topology_perturbation_artifacts,
    write_nucleotide_likelihood_stochastic_topology_perturbation_run_json,
    write_nucleotide_likelihood_topology_perturbation_trace_table,
)


def test_public_runtime_exports_stochastic_topology_perturbation_surface() -> None:
    assert (
        likelihood_api.NucleotideLikelihoodStochasticTopologyPerturbationSearchReport
        is NucleotideLikelihoodStochasticTopologyPerturbationSearchReport
    )
    assert (
        likelihood_api.NucleotideLikelihoodTopologyPerturbationStep
        is NucleotideLikelihoodTopologyPerturbationStep
    )
    assert (
        likelihood_api.search_nucleotide_likelihood_stochastic_topology_perturbation
        is search_nucleotide_likelihood_stochastic_topology_perturbation
    )
    assert (
        likelihood_api.search_nucleotide_likelihood_stochastic_topology_perturbation_from_alignment
        is search_nucleotide_likelihood_stochastic_topology_perturbation_from_alignment
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_local_search_method
        is validate_nucleotide_likelihood_local_search_method
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_perturbation_branch_reoptimization_policy
        is validate_nucleotide_likelihood_perturbation_branch_reoptimization_policy
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_perturbation_move_count
        is validate_nucleotide_likelihood_perturbation_move_count
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_perturbation_move_family
        is validate_nucleotide_likelihood_perturbation_move_family
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_stochastic_topology_perturbation_artifacts
        is write_nucleotide_likelihood_stochastic_topology_perturbation_artifacts
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_stochastic_topology_perturbation_run_json
        is write_nucleotide_likelihood_stochastic_topology_perturbation_run_json
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_topology_perturbation_trace_table
        is write_nucleotide_likelihood_topology_perturbation_trace_table
    )
