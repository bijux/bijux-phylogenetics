from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    NucleotideLikelihoodBootstrapCladeSupportRow,
    NucleotideLikelihoodBootstrapReplicateRow,
    NucleotideLikelihoodBootstrapTreeInferenceReport,
    bootstrap_nucleotide_likelihood_tree_inference,
    bootstrap_nucleotide_likelihood_tree_inference_from_alignment,
    validate_nucleotide_likelihood_bootstrap_replicate_count,
    write_nucleotide_likelihood_bootstrap_artifacts,
    write_nucleotide_likelihood_bootstrap_clade_support_table,
    write_nucleotide_likelihood_bootstrap_replicate_draws_table,
    write_nucleotide_likelihood_bootstrap_run_json,
)


def test_public_likelihood_exports_bootstrap_inference_surface() -> None:
    assert (
        likelihood_api.NucleotideLikelihoodBootstrapReplicateRow
        is NucleotideLikelihoodBootstrapReplicateRow
    )
    assert (
        likelihood_api.NucleotideLikelihoodBootstrapCladeSupportRow
        is NucleotideLikelihoodBootstrapCladeSupportRow
    )
    assert (
        likelihood_api.NucleotideLikelihoodBootstrapTreeInferenceReport
        is NucleotideLikelihoodBootstrapTreeInferenceReport
    )
    assert (
        likelihood_api.bootstrap_nucleotide_likelihood_tree_inference
        is bootstrap_nucleotide_likelihood_tree_inference
    )
    assert (
        likelihood_api.bootstrap_nucleotide_likelihood_tree_inference_from_alignment
        is bootstrap_nucleotide_likelihood_tree_inference_from_alignment
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_bootstrap_replicate_count
        is validate_nucleotide_likelihood_bootstrap_replicate_count
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_bootstrap_artifacts
        is write_nucleotide_likelihood_bootstrap_artifacts
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_bootstrap_clade_support_table
        is write_nucleotide_likelihood_bootstrap_clade_support_table
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_bootstrap_replicate_draws_table
        is write_nucleotide_likelihood_bootstrap_replicate_draws_table
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_bootstrap_run_json
        is write_nucleotide_likelihood_bootstrap_run_json
    )
