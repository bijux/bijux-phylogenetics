from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    FixedTopologyNucleotideJointOptimizationRestartReport,
    JointNucleotideOptimizationRestartAttemptRow,
    optimize_fixed_topology_nucleotide_branches_and_model_with_restarts,
    optimize_fixed_topology_nucleotide_branches_and_model_with_restarts_from_alignment,
    validate_fixed_topology_nucleotide_joint_restart_policy,
)


def test_public_likelihood_exports_fixed_topology_nucleotide_joint_restart_surface() -> (
    None
):
    assert (
        likelihood_api.FixedTopologyNucleotideJointOptimizationRestartReport
        is FixedTopologyNucleotideJointOptimizationRestartReport
    )
    assert (
        likelihood_api.JointNucleotideOptimizationRestartAttemptRow
        is JointNucleotideOptimizationRestartAttemptRow
    )
    assert (
        likelihood_api.optimize_fixed_topology_nucleotide_branches_and_model_with_restarts
        is optimize_fixed_topology_nucleotide_branches_and_model_with_restarts
    )
    assert (
        likelihood_api.optimize_fixed_topology_nucleotide_branches_and_model_with_restarts_from_alignment
        is optimize_fixed_topology_nucleotide_branches_and_model_with_restarts_from_alignment
    )
    assert (
        likelihood_api.validate_fixed_topology_nucleotide_joint_restart_policy
        is validate_fixed_topology_nucleotide_joint_restart_policy
    )
