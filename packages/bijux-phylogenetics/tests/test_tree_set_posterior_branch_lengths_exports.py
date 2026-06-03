from __future__ import annotations

import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    PosteriorBranchLengthSummaryReport,
    PosteriorBranchLengthSummaryRow,
    summarize_posterior_branch_lengths,
    write_posterior_branch_length_summary_table,
)


def test_package_tree_gateway_exports_posterior_branch_length_surface() -> None:
    assert trees_api.PosteriorBranchLengthSummaryRow is PosteriorBranchLengthSummaryRow
    assert (
        trees_api.PosteriorBranchLengthSummaryReport
        is PosteriorBranchLengthSummaryReport
    )
    assert (
        trees_api.summarize_posterior_branch_lengths
        is summarize_posterior_branch_lengths
    )
    assert (
        trees_api.write_posterior_branch_length_summary_table
        is write_posterior_branch_length_summary_table
    )
