from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    NucleotideShLikeBranchSupportLocalTopologyRow,
    NucleotideShLikeBranchSupportReport,
    NucleotideShLikeBranchSupportResamplingRow,
    NucleotideShLikeBranchSupportRow,
    evaluate_nucleotide_sh_like_branch_support,
    evaluate_nucleotide_sh_like_branch_support_from_alignment,
    validate_nucleotide_sh_like_branch_support_replicate_count,
    write_nucleotide_sh_like_branch_support_artifacts,
    write_nucleotide_sh_like_branch_support_local_topology_table,
    write_nucleotide_sh_like_branch_support_resampling_table,
    write_nucleotide_sh_like_branch_support_run_json,
    write_nucleotide_sh_like_branch_support_table,
)


def test_public_likelihood_exports_sh_like_branch_support_surface() -> None:
    assert likelihood_api.NucleotideShLikeBranchSupportRow is NucleotideShLikeBranchSupportRow
    assert (
        likelihood_api.NucleotideShLikeBranchSupportLocalTopologyRow
        is NucleotideShLikeBranchSupportLocalTopologyRow
    )
    assert (
        likelihood_api.NucleotideShLikeBranchSupportResamplingRow
        is NucleotideShLikeBranchSupportResamplingRow
    )
    assert (
        likelihood_api.NucleotideShLikeBranchSupportReport
        is NucleotideShLikeBranchSupportReport
    )
    assert (
        likelihood_api.evaluate_nucleotide_sh_like_branch_support
        is evaluate_nucleotide_sh_like_branch_support
    )
    assert (
        likelihood_api.evaluate_nucleotide_sh_like_branch_support_from_alignment
        is evaluate_nucleotide_sh_like_branch_support_from_alignment
    )
    assert (
        likelihood_api.validate_nucleotide_sh_like_branch_support_replicate_count
        is validate_nucleotide_sh_like_branch_support_replicate_count
    )
    assert (
        likelihood_api.write_nucleotide_sh_like_branch_support_artifacts
        is write_nucleotide_sh_like_branch_support_artifacts
    )
    assert (
        likelihood_api.write_nucleotide_sh_like_branch_support_local_topology_table
        is write_nucleotide_sh_like_branch_support_local_topology_table
    )
    assert (
        likelihood_api.write_nucleotide_sh_like_branch_support_resampling_table
        is write_nucleotide_sh_like_branch_support_resampling_table
    )
    assert (
        likelihood_api.write_nucleotide_sh_like_branch_support_run_json
        is write_nucleotide_sh_like_branch_support_run_json
    )
    assert (
        likelihood_api.write_nucleotide_sh_like_branch_support_table
        is write_nucleotide_sh_like_branch_support_table
    )
