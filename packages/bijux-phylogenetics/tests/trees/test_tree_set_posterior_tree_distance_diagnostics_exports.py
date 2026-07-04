from __future__ import annotations

import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    PosteriorTreeDistanceDiagnosticRow,
    PosteriorTreeDistanceDiagnosticsReport,
    PosteriorTreeDistanceDistributionRow,
    compute_posterior_tree_distance_diagnostics,
    write_posterior_tree_distance_artifacts,
    write_posterior_tree_distance_diagnostic_table,
    write_posterior_tree_distance_distribution_table,
)


def test_package_tree_gateway_exports_posterior_tree_distance_surface() -> None:
    assert (
        trees_api.PosteriorTreeDistanceDiagnosticRow
        is PosteriorTreeDistanceDiagnosticRow
    )
    assert (
        trees_api.PosteriorTreeDistanceDiagnosticsReport
        is PosteriorTreeDistanceDiagnosticsReport
    )
    assert (
        trees_api.PosteriorTreeDistanceDistributionRow
        is PosteriorTreeDistanceDistributionRow
    )
    assert (
        trees_api.compute_posterior_tree_distance_diagnostics
        is compute_posterior_tree_distance_diagnostics
    )
    assert (
        trees_api.write_posterior_tree_distance_artifacts
        is write_posterior_tree_distance_artifacts
    )
    assert (
        trees_api.write_posterior_tree_distance_diagnostic_table
        is write_posterior_tree_distance_diagnostic_table
    )
    assert (
        trees_api.write_posterior_tree_distance_distribution_table
        is write_posterior_tree_distance_distribution_table
    )
