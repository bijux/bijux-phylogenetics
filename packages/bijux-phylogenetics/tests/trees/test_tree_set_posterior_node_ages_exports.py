from __future__ import annotations

import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    PosteriorNodeAgeSummaryReport,
    PosteriorNodeAgeSummaryRow,
    summarize_posterior_node_ages,
    write_posterior_node_age_summary_table,
)


def test_package_tree_gateway_exports_posterior_node_age_surface() -> None:
    assert trees_api.PosteriorNodeAgeSummaryRow is PosteriorNodeAgeSummaryRow
    assert trees_api.PosteriorNodeAgeSummaryReport is PosteriorNodeAgeSummaryReport
    assert trees_api.summarize_posterior_node_ages is summarize_posterior_node_ages
    assert (
        trees_api.write_posterior_node_age_summary_table
        is write_posterior_node_age_summary_table
    )
