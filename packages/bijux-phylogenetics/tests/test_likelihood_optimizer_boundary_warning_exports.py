from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    LikelihoodOptimizationBoundaryWarning,
    boundary_warning_messages,
    build_base_frequency_boundary_warnings,
    build_branch_length_boundary_warnings,
    build_substitution_parameter_boundary_warnings,
)


def test_public_likelihood_exports_optimizer_boundary_warning_surface() -> None:
    assert (
        likelihood_api.LikelihoodOptimizationBoundaryWarning
        is LikelihoodOptimizationBoundaryWarning
    )
    assert (
        likelihood_api.boundary_warning_messages
        is boundary_warning_messages
    )
    assert (
        likelihood_api.build_base_frequency_boundary_warnings
        is build_base_frequency_boundary_warnings
    )
    assert (
        likelihood_api.build_branch_length_boundary_warnings
        is build_branch_length_boundary_warnings
    )
    assert (
        likelihood_api.build_substitution_parameter_boundary_warnings
        is build_substitution_parameter_boundary_warnings
    )
