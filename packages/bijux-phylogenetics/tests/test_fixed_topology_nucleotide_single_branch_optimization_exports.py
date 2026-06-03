from __future__ import annotations

from bijux_phylogenetics.phylo.likelihood import (
    FixedTopologyNucleotideSingleBranchOptimizationReport,
    optimize_fixed_topology_nucleotide_single_branch_length,
    optimize_fixed_topology_nucleotide_single_branch_length_from_alignment,
    validate_fixed_topology_nucleotide_single_branch_id,
)
from bijux_phylogenetics.phylo.likelihood.fixed_topology_branch_lengths import (
    optimize_fixed_topology_nucleotide_single_branch_length as optimize_fixed_topology_nucleotide_single_branch_length_impl,
)
from bijux_phylogenetics.phylo.likelihood.fixed_topology_branch_lengths import (
    optimize_fixed_topology_nucleotide_single_branch_length_from_alignment as optimize_fixed_topology_nucleotide_single_branch_length_from_alignment_impl,
)
from bijux_phylogenetics.phylo.likelihood.fixed_topology_branch_lengths import (
    validate_fixed_topology_nucleotide_single_branch_id as validate_fixed_topology_nucleotide_single_branch_id_impl,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    FixedTopologyNucleotideSingleBranchOptimizationReport as FixedTopologyNucleotideSingleBranchOptimizationReportImpl,
)


def test_likelihood_exports_fixed_topology_nucleotide_single_branch_optimization_surface() -> (
    None
):
    assert (
        FixedTopologyNucleotideSingleBranchOptimizationReport
        is FixedTopologyNucleotideSingleBranchOptimizationReportImpl
    )
    assert (
        optimize_fixed_topology_nucleotide_single_branch_length
        is optimize_fixed_topology_nucleotide_single_branch_length_impl
    )
    assert (
        optimize_fixed_topology_nucleotide_single_branch_length_from_alignment
        is optimize_fixed_topology_nucleotide_single_branch_length_from_alignment_impl
    )
    assert (
        validate_fixed_topology_nucleotide_single_branch_id
        is validate_fixed_topology_nucleotide_single_branch_id_impl
    )
