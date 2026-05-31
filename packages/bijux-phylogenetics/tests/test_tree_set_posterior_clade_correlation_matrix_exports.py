from __future__ import annotations

import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    PosteriorCladeCorrelationReport,
    PosteriorCladeCorrelationRow,
    compute_posterior_clade_correlation_matrix,
    write_posterior_clade_correlation_artifacts,
    write_posterior_clade_correlation_matrix_table,
    write_posterior_clade_correlation_pair_table,
)


def test_package_tree_gateway_exports_posterior_clade_correlation_surface() -> None:
    assert trees_api.PosteriorCladeCorrelationReport is PosteriorCladeCorrelationReport
    assert trees_api.PosteriorCladeCorrelationRow is PosteriorCladeCorrelationRow
    assert (
        trees_api.compute_posterior_clade_correlation_matrix
        is compute_posterior_clade_correlation_matrix
    )
    assert (
        trees_api.write_posterior_clade_correlation_artifacts
        is write_posterior_clade_correlation_artifacts
    )
    assert (
        trees_api.write_posterior_clade_correlation_matrix_table
        is write_posterior_clade_correlation_matrix_table
    )
    assert (
        trees_api.write_posterior_clade_correlation_pair_table
        is write_posterior_clade_correlation_pair_table
    )
