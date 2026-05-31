from __future__ import annotations

from bijux_phylogenetics.phylo.likelihood import (
    FixedTopologyNucleotideBranchLengthOptimizationReport,
    optimize_fixed_topology_nucleotide_branch_lengths,
    optimize_fixed_topology_nucleotide_branch_lengths_from_alignment,
    validate_fixed_topology_nucleotide_branch_length_model,
)
from bijux_phylogenetics.phylo.likelihood.fixed_topology_branch_lengths import (
    optimize_fixed_topology_nucleotide_branch_lengths as optimize_fixed_topology_nucleotide_branch_lengths_impl,
)
from bijux_phylogenetics.phylo.likelihood.fixed_topology_branch_lengths import (
    optimize_fixed_topology_nucleotide_branch_lengths_from_alignment as optimize_fixed_topology_nucleotide_branch_lengths_from_alignment_impl,
)
from bijux_phylogenetics.phylo.likelihood.fixed_topology_branch_lengths import (
    validate_fixed_topology_nucleotide_branch_length_model as validate_fixed_topology_nucleotide_branch_length_model_impl,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    FixedTopologyNucleotideBranchLengthOptimizationReport as FixedTopologyNucleotideBranchLengthOptimizationReportImpl,
)


def test_likelihood_exports_fixed_topology_nucleotide_branch_optimization_surface() -> None:
    assert (
        FixedTopologyNucleotideBranchLengthOptimizationReport
        is FixedTopologyNucleotideBranchLengthOptimizationReportImpl
    )
    assert (
        optimize_fixed_topology_nucleotide_branch_lengths
        is optimize_fixed_topology_nucleotide_branch_lengths_impl
    )
    assert (
        optimize_fixed_topology_nucleotide_branch_lengths_from_alignment
        is optimize_fixed_topology_nucleotide_branch_lengths_from_alignment_impl
    )
    assert (
        validate_fixed_topology_nucleotide_branch_length_model
        is validate_fixed_topology_nucleotide_branch_length_model_impl
    )
