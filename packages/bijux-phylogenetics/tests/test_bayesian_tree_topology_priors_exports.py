from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    TREE_TOPOLOGY_PRIOR_FAMILIES,
    TreeTopologyPriorEvaluationReport,
    TreeTopologyPriorModel,
    build_uniform_rooted_tree_topology_prior,
    count_rooted_labeled_bifurcating_topologies,
    evaluate_tree_topology_log_prior,
    validate_tree_topology_prior_taxa,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    TREE_TOPOLOGY_PRIOR_FAMILIES as TREE_TOPOLOGY_PRIOR_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    TreeTopologyPriorEvaluationReport as TreeTopologyPriorEvaluationReportImpl,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    TreeTopologyPriorModel as TreeTopologyPriorModelImpl,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    build_uniform_rooted_tree_topology_prior as build_uniform_rooted_tree_topology_prior_impl,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    count_rooted_labeled_bifurcating_topologies as count_rooted_labeled_bifurcating_topologies_impl,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    evaluate_tree_topology_log_prior as evaluate_tree_topology_log_prior_impl,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    validate_tree_topology_prior_taxa as validate_tree_topology_prior_taxa_impl,
)


def test_bayesian_exports_tree_topology_prior_surface() -> None:
    assert TREE_TOPOLOGY_PRIOR_FAMILIES == TREE_TOPOLOGY_PRIOR_FAMILIES_IMPL
    assert TreeTopologyPriorModel is TreeTopologyPriorModelImpl
    assert TreeTopologyPriorEvaluationReport is TreeTopologyPriorEvaluationReportImpl
    assert (
        build_uniform_rooted_tree_topology_prior
        is build_uniform_rooted_tree_topology_prior_impl
    )
    assert (
        count_rooted_labeled_bifurcating_topologies
        is count_rooted_labeled_bifurcating_topologies_impl
    )
    assert evaluate_tree_topology_log_prior is evaluate_tree_topology_log_prior_impl
    assert validate_tree_topology_prior_taxa is validate_tree_topology_prior_taxa_impl
