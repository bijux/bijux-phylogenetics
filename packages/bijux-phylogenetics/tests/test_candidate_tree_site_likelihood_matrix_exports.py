from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    CandidateTreeSiteLikelihoodMatrixReport,
    CandidateTreeSiteLikelihoodRow,
    CandidateTreeSiteLikelihoodSummary,
    evaluate_nucleotide_candidate_tree_site_likelihood_matrix,
    evaluate_nucleotide_candidate_tree_site_likelihood_matrix_from_alignment,
    resolve_candidate_tree_alignment_records,
    resolve_candidate_tree_records,
    write_candidate_tree_likelihood_summary_table,
    write_candidate_tree_site_likelihood_matrix_artifacts,
    write_candidate_tree_site_likelihood_matrix_run_json,
    write_candidate_tree_site_likelihood_matrix_table,
)


def test_public_likelihood_exports_candidate_tree_comparison_surface() -> None:
    assert (
        likelihood_api.CandidateTreeSiteLikelihoodSummary
        is CandidateTreeSiteLikelihoodSummary
    )
    assert (
        likelihood_api.CandidateTreeSiteLikelihoodRow is CandidateTreeSiteLikelihoodRow
    )
    assert (
        likelihood_api.CandidateTreeSiteLikelihoodMatrixReport
        is CandidateTreeSiteLikelihoodMatrixReport
    )
    assert (
        likelihood_api.evaluate_nucleotide_candidate_tree_site_likelihood_matrix
        is evaluate_nucleotide_candidate_tree_site_likelihood_matrix
    )
    assert (
        likelihood_api.evaluate_nucleotide_candidate_tree_site_likelihood_matrix_from_alignment
        is evaluate_nucleotide_candidate_tree_site_likelihood_matrix_from_alignment
    )
    assert (
        likelihood_api.resolve_candidate_tree_alignment_records
        is resolve_candidate_tree_alignment_records
    )
    assert (
        likelihood_api.resolve_candidate_tree_records is resolve_candidate_tree_records
    )
    assert (
        likelihood_api.write_candidate_tree_likelihood_summary_table
        is write_candidate_tree_likelihood_summary_table
    )
    assert (
        likelihood_api.write_candidate_tree_site_likelihood_matrix_artifacts
        is write_candidate_tree_site_likelihood_matrix_artifacts
    )
    assert (
        likelihood_api.write_candidate_tree_site_likelihood_matrix_run_json
        is write_candidate_tree_site_likelihood_matrix_run_json
    )
    assert (
        likelihood_api.write_candidate_tree_site_likelihood_matrix_table
        is write_candidate_tree_site_likelihood_matrix_table
    )
