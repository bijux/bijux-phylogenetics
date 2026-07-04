from __future__ import annotations

import math

import pytest

from bijux_phylogenetics.bayesian.time_tree_priors import (
    SkylineCoalescentEpoch,
    build_constant_population_coalescent_tree_prior,
    build_crown_conditioned_birth_death_tree_prior,
    build_crown_conditioned_yule_tree_prior,
    build_skyline_coalescent_tree_prior,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    evaluate_time_tree_topology_log_prior,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


@pytest.mark.parametrize(
    "prior_model",
    [
        build_crown_conditioned_yule_tree_prior(speciation_rate=0.5),
        build_crown_conditioned_birth_death_tree_prior(
            speciation_rate=1.2,
            extinction_rate=0.4,
            sampling_fraction=0.75,
        ),
        build_constant_population_coalescent_tree_prior(
            effective_population_size=2.0,
        ),
        build_skyline_coalescent_tree_prior(
            [
                SkylineCoalescentEpoch(
                    younger_boundary_age=0.0,
                    older_boundary_age=1.5,
                    effective_population_size=2.0,
                ),
                SkylineCoalescentEpoch(
                    younger_boundary_age=1.5,
                    older_boundary_age=None,
                    effective_population_size=4.0,
                ),
            ]
        ),
    ],
)
def test_time_tree_topology_priors_ignore_child_order(
    prior_model: object,
) -> None:
    left_tree = PhyloTree(
        TreeNode(
            children=[
                TreeNode(
                    branch_length=1.0,
                    children=[
                        TreeNode(
                            branch_length=1.0,
                            children=[
                                TreeNode(name="A", branch_length=1.0),
                                TreeNode(name="B", branch_length=1.0),
                            ],
                        ),
                        TreeNode(name="C", branch_length=2.0),
                    ],
                ),
                TreeNode(name="D", branch_length=3.0),
            ]
        ),
        rooted=True,
    )
    right_tree = PhyloTree(
        TreeNode(
            children=[
                TreeNode(name="D", branch_length=3.0),
                TreeNode(
                    branch_length=1.0,
                    children=[
                        TreeNode(name="C", branch_length=2.0),
                        TreeNode(
                            branch_length=1.0,
                            children=[
                                TreeNode(name="B", branch_length=1.0),
                                TreeNode(name="A", branch_length=1.0),
                            ],
                        ),
                    ],
                ),
            ]
        ),
        rooted=True,
    )

    left_report = evaluate_time_tree_topology_log_prior(left_tree, prior_model)
    right_report = evaluate_time_tree_topology_log_prior(right_tree, prior_model)

    assert left_report.family == right_report.family
    assert math.isclose(
        left_report.log_prior,
        right_report.log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
