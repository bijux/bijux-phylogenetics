from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    FixedTopologyNucleotideJointOptimizationReport,
    JointNucleotideOptimizationUpdateRow,
    optimize_fixed_topology_nucleotide_branches_and_model,
    optimize_fixed_topology_nucleotide_branches_and_model_from_alignment,
    validate_fixed_topology_nucleotide_joint_optimization_model,
)


def test_public_likelihood_exports_fixed_topology_nucleotide_joint_optimization_surface() -> (
    None
):
    assert (
        likelihood_api.FixedTopologyNucleotideJointOptimizationReport
        is FixedTopologyNucleotideJointOptimizationReport
    )
    assert (
        likelihood_api.JointNucleotideOptimizationUpdateRow
        is JointNucleotideOptimizationUpdateRow
    )
    assert (
        likelihood_api.optimize_fixed_topology_nucleotide_branches_and_model
        is optimize_fixed_topology_nucleotide_branches_and_model
    )
    assert (
        likelihood_api.optimize_fixed_topology_nucleotide_branches_and_model_from_alignment
        is optimize_fixed_topology_nucleotide_branches_and_model_from_alignment
    )
    assert (
        likelihood_api.validate_fixed_topology_nucleotide_joint_optimization_model
        is validate_fixed_topology_nucleotide_joint_optimization_model
    )
