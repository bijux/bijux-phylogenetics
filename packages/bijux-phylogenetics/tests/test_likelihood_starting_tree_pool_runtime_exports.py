from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    NucleotideLikelihoodStartingTreePoolReport,
    NucleotideLikelihoodStartingTreeSummary,
    build_nucleotide_likelihood_starting_tree_pool,
    build_nucleotide_likelihood_starting_tree_pool_from_alignment,
    validate_nucleotide_likelihood_random_start_tree_count,
    validate_nucleotide_likelihood_starting_tree_pool_model,
    write_nucleotide_likelihood_starting_tree_pool_artifacts,
    write_nucleotide_likelihood_starting_tree_pool_run_json,
    write_nucleotide_likelihood_starting_tree_score_table,
)


def test_public_runtime_exports_likelihood_starting_tree_pool_surface() -> None:
    assert (
        likelihood_api.NucleotideLikelihoodStartingTreePoolReport
        is NucleotideLikelihoodStartingTreePoolReport
    )
    assert (
        likelihood_api.NucleotideLikelihoodStartingTreeSummary
        is NucleotideLikelihoodStartingTreeSummary
    )
    assert (
        likelihood_api.build_nucleotide_likelihood_starting_tree_pool
        is build_nucleotide_likelihood_starting_tree_pool
    )
    assert (
        likelihood_api.build_nucleotide_likelihood_starting_tree_pool_from_alignment
        is build_nucleotide_likelihood_starting_tree_pool_from_alignment
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_random_start_tree_count
        is validate_nucleotide_likelihood_random_start_tree_count
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_starting_tree_pool_model
        is validate_nucleotide_likelihood_starting_tree_pool_model
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_starting_tree_score_table
        is write_nucleotide_likelihood_starting_tree_score_table
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_starting_tree_pool_run_json
        is write_nucleotide_likelihood_starting_tree_pool_run_json
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_starting_tree_pool_artifacts
        is write_nucleotide_likelihood_starting_tree_pool_artifacts
    )
