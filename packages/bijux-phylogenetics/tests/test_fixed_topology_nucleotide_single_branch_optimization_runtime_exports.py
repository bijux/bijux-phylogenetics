from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    FixedTopologyNucleotideSingleBranchOptimizationReport,
    optimize_fixed_topology_nucleotide_single_branch_length,
    optimize_fixed_topology_nucleotide_single_branch_length_from_alignment,
    validate_fixed_topology_nucleotide_single_branch_id,
)


def test_public_runtime_exports_fixed_topology_nucleotide_single_branch_optimization_surface() -> (
    None
):
    assert (
        likelihood_api.FixedTopologyNucleotideSingleBranchOptimizationReport
        is FixedTopologyNucleotideSingleBranchOptimizationReport
    )
    assert (
        likelihood_api.optimize_fixed_topology_nucleotide_single_branch_length
        is optimize_fixed_topology_nucleotide_single_branch_length
    )
    assert (
        likelihood_api.optimize_fixed_topology_nucleotide_single_branch_length_from_alignment
        is optimize_fixed_topology_nucleotide_single_branch_length_from_alignment
    )
    assert (
        likelihood_api.validate_fixed_topology_nucleotide_single_branch_id
        is validate_fixed_topology_nucleotide_single_branch_id
    )
