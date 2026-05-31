from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    NucleotideLikelihoodRatchetBestTreeHistory,
    NucleotideLikelihoodRatchetCycle,
    NucleotideLikelihoodRatchetReport,
    search_nucleotide_likelihood_ratchet,
    search_nucleotide_likelihood_ratchet_from_alignment,
    validate_nucleotide_likelihood_ratchet_branch_reoptimization_policy,
    validate_nucleotide_likelihood_ratchet_cycle_count,
    validate_nucleotide_likelihood_ratchet_evaluation_budget,
    validate_nucleotide_likelihood_ratchet_local_search_method,
    validate_nucleotide_likelihood_ratchet_perturbation_factor,
    validate_nucleotide_likelihood_ratchet_perturbed_site_count,
    write_nucleotide_likelihood_ratchet_artifacts,
    write_nucleotide_likelihood_ratchet_best_tree_history_table,
    write_nucleotide_likelihood_ratchet_cycle_table,
    write_nucleotide_likelihood_ratchet_run_json,
)


def test_public_likelihood_exports_ratchet_search_surface() -> None:
    assert (
        likelihood_api.NucleotideLikelihoodRatchetBestTreeHistory
        is NucleotideLikelihoodRatchetBestTreeHistory
    )
    assert likelihood_api.NucleotideLikelihoodRatchetCycle is NucleotideLikelihoodRatchetCycle
    assert likelihood_api.NucleotideLikelihoodRatchetReport is NucleotideLikelihoodRatchetReport
    assert likelihood_api.search_nucleotide_likelihood_ratchet is search_nucleotide_likelihood_ratchet
    assert (
        likelihood_api.search_nucleotide_likelihood_ratchet_from_alignment
        is search_nucleotide_likelihood_ratchet_from_alignment
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_ratchet_branch_reoptimization_policy
        is validate_nucleotide_likelihood_ratchet_branch_reoptimization_policy
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_ratchet_cycle_count
        is validate_nucleotide_likelihood_ratchet_cycle_count
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_ratchet_evaluation_budget
        is validate_nucleotide_likelihood_ratchet_evaluation_budget
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_ratchet_local_search_method
        is validate_nucleotide_likelihood_ratchet_local_search_method
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_ratchet_perturbation_factor
        is validate_nucleotide_likelihood_ratchet_perturbation_factor
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_ratchet_perturbed_site_count
        is validate_nucleotide_likelihood_ratchet_perturbed_site_count
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_ratchet_artifacts
        is write_nucleotide_likelihood_ratchet_artifacts
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_ratchet_best_tree_history_table
        is write_nucleotide_likelihood_ratchet_best_tree_history_table
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_ratchet_cycle_table
        is write_nucleotide_likelihood_ratchet_cycle_table
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_ratchet_run_json
        is write_nucleotide_likelihood_ratchet_run_json
    )
