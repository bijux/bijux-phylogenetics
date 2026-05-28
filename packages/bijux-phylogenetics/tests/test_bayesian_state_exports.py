from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    BayesianModelParameterState,
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    BayesianStateBranchRow,
    BayesianTreeState,
    build_bayesian_model_parameter_state,
    build_bayesian_phylogenetic_state,
    build_bayesian_phylogenetic_state_from_prior_only_sample,
    build_bayesian_prior_component_state,
    build_bayesian_tree_state,
    deserialize_bayesian_phylogenetic_state,
    deserialize_bayesian_phylogenetic_state_json,
    serialize_bayesian_phylogenetic_state,
    serialize_bayesian_phylogenetic_state_json,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianModelParameterState as BayesianModelParameterStateImpl,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState as BayesianPhylogeneticStateImpl,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPriorComponentState as BayesianPriorComponentStateImpl,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianStateBranchRow as BayesianStateBranchRowImpl,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianTreeState as BayesianTreeStateImpl,
)
from bijux_phylogenetics.bayesian.state import (
    build_bayesian_model_parameter_state as build_bayesian_model_parameter_state_impl,
)
from bijux_phylogenetics.bayesian.state import (
    build_bayesian_phylogenetic_state as build_bayesian_phylogenetic_state_impl,
)
from bijux_phylogenetics.bayesian.state import (
    build_bayesian_phylogenetic_state_from_prior_only_sample as build_bayesian_phylogenetic_state_from_prior_only_sample_impl,
)
from bijux_phylogenetics.bayesian.state import (
    build_bayesian_prior_component_state as build_bayesian_prior_component_state_impl,
)
from bijux_phylogenetics.bayesian.state import (
    build_bayesian_tree_state as build_bayesian_tree_state_impl,
)
from bijux_phylogenetics.bayesian.state import (
    deserialize_bayesian_phylogenetic_state as deserialize_bayesian_phylogenetic_state_impl,
)
from bijux_phylogenetics.bayesian.state import (
    deserialize_bayesian_phylogenetic_state_json as deserialize_bayesian_phylogenetic_state_json_impl,
)
from bijux_phylogenetics.bayesian.state import (
    serialize_bayesian_phylogenetic_state as serialize_bayesian_phylogenetic_state_impl,
)
from bijux_phylogenetics.bayesian.state import (
    serialize_bayesian_phylogenetic_state_json as serialize_bayesian_phylogenetic_state_json_impl,
)


def test_bayesian_exports_phylogenetic_state_surface() -> None:
    assert BayesianStateBranchRow is BayesianStateBranchRowImpl
    assert BayesianTreeState is BayesianTreeStateImpl
    assert BayesianModelParameterState is BayesianModelParameterStateImpl
    assert BayesianPriorComponentState is BayesianPriorComponentStateImpl
    assert BayesianPhylogeneticState is BayesianPhylogeneticStateImpl
    assert build_bayesian_tree_state is build_bayesian_tree_state_impl
    assert (
        build_bayesian_model_parameter_state
        is build_bayesian_model_parameter_state_impl
    )
    assert build_bayesian_prior_component_state is (
        build_bayesian_prior_component_state_impl
    )
    assert build_bayesian_phylogenetic_state is build_bayesian_phylogenetic_state_impl
    assert (
        build_bayesian_phylogenetic_state_from_prior_only_sample
        is build_bayesian_phylogenetic_state_from_prior_only_sample_impl
    )
    assert (
        serialize_bayesian_phylogenetic_state
        is serialize_bayesian_phylogenetic_state_impl
    )
    assert (
        serialize_bayesian_phylogenetic_state_json
        is serialize_bayesian_phylogenetic_state_json_impl
    )
    assert (
        deserialize_bayesian_phylogenetic_state
        is deserialize_bayesian_phylogenetic_state_impl
    )
    assert (
        deserialize_bayesian_phylogenetic_state_json
        is deserialize_bayesian_phylogenetic_state_json_impl
    )
