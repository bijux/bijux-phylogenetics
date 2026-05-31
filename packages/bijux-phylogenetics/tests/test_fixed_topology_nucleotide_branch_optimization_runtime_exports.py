from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    FixedTopologyNucleotideBranchLengthOptimizationReport,
    optimize_fixed_topology_nucleotide_branch_lengths,
    optimize_fixed_topology_nucleotide_branch_lengths_from_alignment,
    validate_fixed_topology_nucleotide_branch_length_model,
)


def test_public_runtime_exports_fixed_topology_nucleotide_branch_optimization_surface() -> (
    None
):
    assert (
        likelihood_api.FixedTopologyNucleotideBranchLengthOptimizationReport
        is FixedTopologyNucleotideBranchLengthOptimizationReport
    )
    assert (
        likelihood_api.optimize_fixed_topology_nucleotide_branch_lengths
        is optimize_fixed_topology_nucleotide_branch_lengths
    )
    assert (
        likelihood_api.optimize_fixed_topology_nucleotide_branch_lengths_from_alignment
        is optimize_fixed_topology_nucleotide_branch_lengths_from_alignment
    )
    assert (
        likelihood_api.validate_fixed_topology_nucleotide_branch_length_model
        is validate_fixed_topology_nucleotide_branch_length_model
    )
