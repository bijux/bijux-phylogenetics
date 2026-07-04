from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    build_partition_model_parameter_state,
    resolve_partition_parameter_linkage_plan_from_model_parameters,
    resolve_partition_parameter_states_from_model_parameters,
    strip_partition_model_parameter_state,
)
from bijux_phylogenetics.bayesian.partition_model_state import (
    build_partition_model_parameter_state as build_partition_model_parameter_state_impl,
)
from bijux_phylogenetics.bayesian.partition_model_state import (
    resolve_partition_parameter_linkage_plan_from_model_parameters as resolve_partition_parameter_linkage_plan_from_model_parameters_impl,
)
from bijux_phylogenetics.bayesian.partition_model_state import (
    resolve_partition_parameter_states_from_model_parameters as resolve_partition_parameter_states_from_model_parameters_impl,
)
from bijux_phylogenetics.bayesian.partition_model_state import (
    strip_partition_model_parameter_state as strip_partition_model_parameter_state_impl,
)


def test_bayesian_exports_partition_model_state_helpers() -> None:
    assert (
        build_partition_model_parameter_state
        is build_partition_model_parameter_state_impl
    )
    assert (
        resolve_partition_parameter_linkage_plan_from_model_parameters
        is resolve_partition_parameter_linkage_plan_from_model_parameters_impl
    )
    assert (
        resolve_partition_parameter_states_from_model_parameters
        is resolve_partition_parameter_states_from_model_parameters_impl
    )
    assert (
        strip_partition_model_parameter_state
        is strip_partition_model_parameter_state_impl
    )
