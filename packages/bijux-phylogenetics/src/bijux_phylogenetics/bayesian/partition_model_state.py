from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
import math

from bijux_phylogenetics.bayesian.partition_model_priors import (
    PARTITION_MODEL_PRIOR_TARGETS,
    PartitionParameterLinkagePlan,
    PartitionSubstitutionModelDefinition,
    PartitionSubstitutionParameterState,
    build_partition_parameter_linkage_plan,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianModelParameterState,
    build_bayesian_model_parameter_state,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

_PARTITION_LINKAGE_PREFIX = "partition-linkage:"
_PARTITION_PARAMETER_PREFIX = "partition-parameter:"
_PARTITION_STATE_VALUE_TOLERANCE = 1e-12
_SCALAR_PARTITION_TARGETS = frozenset(
    {
        "kappa",
        "gamma-alpha",
        "invariant-proportion",
    }
)
_VECTOR_PARTITION_TARGETS = frozenset(
    {
        "exchangeabilities",
        "base-frequencies",
    }
)
_DNA_BASE_FREQUENCY_COMPONENT_NAMES = ("A", "C", "G", "T")
_DNA_EXCHANGEABILITY_COMPONENT_NAMES = ("AC", "AG", "AT", "CG", "CT", "GT")


def build_partition_model_parameter_state(
    *,
    partition_models: Sequence[PartitionSubstitutionModelDefinition],
    linkage_plan: PartitionParameterLinkagePlan,
    partition_parameter_states: Sequence[PartitionSubstitutionParameterState],
    preserved_categorical_parameters: Mapping[str, str] | None = None,
    preserved_scalar_parameters: Mapping[str, float] | None = None,
    preserved_vector_parameters: Mapping[str, Mapping[str, float]] | None = None,
) -> BayesianModelParameterState:
    """Encode one partition-model parameter surface into Bayesian sampler state."""
    validated_partition_models = tuple(partition_models)
    expected_partition_names = _validate_unique_partition_names(
        model.partition_name for model in validated_partition_models
    )
    if expected_partition_names != linkage_plan.partition_names:
        raise PhylogeneticsError(
            "partition model state encoding requires linkage coverage for the exact partition set",
            code="partition_model_state_linkage_partition_mismatch",
            details={
                "model_partition_names": list(expected_partition_names),
                "linkage_partition_names": list(linkage_plan.partition_names),
            },
        )
    state_by_partition_name = _build_partition_state_lookup(
        partition_parameter_states=partition_parameter_states,
        expected_partition_names=expected_partition_names,
    )
    if set(state_by_partition_name) != set(expected_partition_names):
        raise PhylogeneticsError(
            "partition model state encoding requires states for the exact partition set",
            code="partition_model_state_partition_mismatch",
            details={
                "expected_partition_names": list(expected_partition_names),
                "observed_partition_names": sorted(state_by_partition_name),
            },
        )
    categorical_parameters = dict(preserved_categorical_parameters or {})
    scalar_parameters = dict(preserved_scalar_parameters or {})
    vector_parameters = {
        parameter_name: dict(component_values)
        for parameter_name, component_values in (
            preserved_vector_parameters or {}
        ).items()
    }
    _validate_reserved_parameter_collisions(
        categorical_parameters=categorical_parameters,
        scalar_parameters=scalar_parameters,
        vector_parameters=vector_parameters,
    )
    for target_name in PARTITION_MODEL_PRIOR_TARGETS:
        target_groups = linkage_plan.groups_for_target(target_name)
        for partition_name, group_name in target_groups.items():
            categorical_parameters[
                _partition_linkage_key(
                    target_name=target_name,
                    partition_name=partition_name,
                )
            ] = group_name
        required_partition_names = tuple(
            model.partition_name
            for model in partition_models
            if target_name in model.required_targets()
        )
        if not required_partition_names:
            continue
        grouped_partition_names: dict[str, list[str]] = defaultdict(list)
        for partition_name in required_partition_names:
            grouped_partition_names[target_groups[partition_name]].append(
                partition_name
            )
        for group_name, grouped_names in grouped_partition_names.items():
            representative_value = _extract_partition_target_value(
                target_name=target_name,
                state=state_by_partition_name[grouped_names[0]],
            )
            for partition_name in grouped_names[1:]:
                compared_value = _extract_partition_target_value(
                    target_name=target_name,
                    state=state_by_partition_name[partition_name],
                )
                if not _partition_target_values_equal(
                    representative_value,
                    compared_value,
                ):
                    raise PhylogeneticsError(
                        "partition model state encoding requires linked groups to share one realized value",
                        code="partition_model_state_linked_values_mismatched",
                        details={
                            "target_name": target_name,
                            "group_name": group_name,
                            "representative_partition_name": grouped_names[0],
                            "compared_partition_name": partition_name,
                        },
                    )
            parameter_key = _partition_parameter_key(
                target_name=target_name,
                group_name=group_name,
            )
            if target_name in _SCALAR_PARTITION_TARGETS:
                scalar_parameters[parameter_key] = float(representative_value)
                continue
            vector_parameters[parameter_key] = {
                component_name: float(component_value)
                for component_name, component_value in representative_value.items()
            }
    return build_bayesian_model_parameter_state(
        categorical_parameters=categorical_parameters,
        scalar_parameters=scalar_parameters,
        vector_parameters=vector_parameters,
    )


def strip_partition_model_parameter_state(
    model_parameters: BayesianModelParameterState,
) -> BayesianModelParameterState:
    """Return one model-parameter state without reserved partition-model keys."""
    return build_bayesian_model_parameter_state(
        categorical_parameters={
            parameter_name: parameter_value
            for parameter_name, parameter_value in model_parameters.categorical_parameters.items()
            if not parameter_name.startswith(_PARTITION_LINKAGE_PREFIX)
        },
        scalar_parameters={
            parameter_name: parameter_value
            for parameter_name, parameter_value in model_parameters.scalar_parameters.items()
            if not parameter_name.startswith(_PARTITION_PARAMETER_PREFIX)
        },
        vector_parameters={
            parameter_name: component_values
            for parameter_name, component_values in model_parameters.vector_parameters.items()
            if not parameter_name.startswith(_PARTITION_PARAMETER_PREFIX)
        },
    )


def resolve_partition_parameter_linkage_plan_from_model_parameters(
    *,
    model_parameters: BayesianModelParameterState,
    partition_names: Sequence[str],
) -> PartitionParameterLinkagePlan:
    """Decode one explicit partition linkage plan from Bayesian sampler state."""
    validated_partition_names = _validate_unique_partition_names(partition_names)
    group_assignments: dict[str, dict[str, str]] = {}
    for target_name in PARTITION_MODEL_PRIOR_TARGETS:
        target_assignments: dict[str, str] = {}
        for partition_name in validated_partition_names:
            linkage_key = _partition_linkage_key(
                target_name=target_name,
                partition_name=partition_name,
            )
            group_name = model_parameters.categorical_parameters.get(linkage_key)
            if group_name is None:
                raise PhylogeneticsError(
                    "partition model state decoding requires explicit linkage assignment for every target and partition",
                    code="partition_model_state_linkage_assignment_missing",
                    details={
                        "target_name": target_name,
                        "partition_name": partition_name,
                        "linkage_key": linkage_key,
                    },
                )
            target_assignments[partition_name] = group_name
        group_assignments[target_name] = target_assignments
    return build_partition_parameter_linkage_plan(
        partition_names=validated_partition_names,
        group_assignments=group_assignments,
    )


def resolve_partition_parameter_states_from_model_parameters(
    *,
    model_parameters: BayesianModelParameterState,
    partition_models: Sequence[PartitionSubstitutionModelDefinition],
    linkage_plan: PartitionParameterLinkagePlan,
) -> tuple[PartitionSubstitutionParameterState, ...]:
    """Decode one partition-parameter realization from Bayesian sampler state."""
    validated_partition_models = tuple(partition_models)
    _validate_unique_partition_names(
        partition_model.partition_name for partition_model in validated_partition_models
    )
    resolved_states: list[PartitionSubstitutionParameterState] = []
    for partition_model in validated_partition_models:
        target_groups = {
            target_name: linkage_plan.groups_for_target(target_name)
            for target_name in PARTITION_MODEL_PRIOR_TARGETS
        }
        partition_name = partition_model.partition_name
        state_kwargs: dict[str, object] = {"partition_name": partition_name}
        if "kappa" in partition_model.required_targets():
            state_kwargs["kappa"] = _require_partition_scalar_value(
                model_parameters=model_parameters,
                target_name="kappa",
                group_name=target_groups["kappa"][partition_name],
            )
        if "exchangeabilities" in partition_model.required_targets():
            state_kwargs["exchangeabilities"] = _require_partition_vector_value(
                model_parameters=model_parameters,
                target_name="exchangeabilities",
                group_name=target_groups["exchangeabilities"][partition_name],
            )
        if "base-frequencies" in partition_model.required_targets():
            state_kwargs["base_frequencies"] = _require_partition_vector_value(
                model_parameters=model_parameters,
                target_name="base-frequencies",
                group_name=target_groups["base-frequencies"][partition_name],
            )
        if "gamma-alpha" in partition_model.required_targets():
            state_kwargs["gamma_alpha"] = _require_partition_scalar_value(
                model_parameters=model_parameters,
                target_name="gamma-alpha",
                group_name=target_groups["gamma-alpha"][partition_name],
            )
        if "invariant-proportion" in partition_model.required_targets():
            state_kwargs["invariant_proportion"] = _require_partition_scalar_value(
                model_parameters=model_parameters,
                target_name="invariant-proportion",
                group_name=target_groups["invariant-proportion"][partition_name],
            )
        resolved_states.append(PartitionSubstitutionParameterState(**state_kwargs))
    return tuple(resolved_states)


def _partition_linkage_key(*, target_name: str, partition_name: str) -> str:
    return f"{_PARTITION_LINKAGE_PREFIX}{target_name}:{partition_name}"


def _partition_parameter_key(*, target_name: str, group_name: str) -> str:
    return f"{_PARTITION_PARAMETER_PREFIX}{target_name}:{group_name}"


def _extract_partition_target_value(
    *,
    target_name: str,
    state: PartitionSubstitutionParameterState,
) -> float | Mapping[str, float]:
    if target_name == "kappa":
        if state.kappa is None:
            raise PhylogeneticsError(
                "partition model state encoding requires kappa on every partition that uses it",
                code="partition_model_state_kappa_missing",
                details={"partition_name": state.partition_name},
            )
        return state.kappa
    if target_name == "exchangeabilities":
        if state.exchangeabilities is None:
            raise PhylogeneticsError(
                "partition model state encoding requires exchangeabilities on every partition that uses them",
                code="partition_model_state_exchangeabilities_missing",
                details={"partition_name": state.partition_name},
            )
        return _normalize_vector_value(
            target_name=target_name,
            raw_value=state.exchangeabilities,
        )
    if target_name == "base-frequencies":
        if state.base_frequencies is None:
            raise PhylogeneticsError(
                "partition model state encoding requires base frequencies on every partition that uses them",
                code="partition_model_state_base_frequencies_missing",
                details={"partition_name": state.partition_name},
            )
        return _normalize_vector_value(
            target_name=target_name,
            raw_value=state.base_frequencies,
        )
    if target_name == "gamma-alpha":
        if state.gamma_alpha is None:
            raise PhylogeneticsError(
                "partition model state encoding requires gamma alpha on every partition that uses it",
                code="partition_model_state_gamma_alpha_missing",
                details={"partition_name": state.partition_name},
            )
        return state.gamma_alpha
    if target_name == "invariant-proportion":
        if state.invariant_proportion is None:
            raise PhylogeneticsError(
                "partition model state encoding requires invariant proportion on every partition that uses it",
                code="partition_model_state_invariant_proportion_missing",
                details={"partition_name": state.partition_name},
            )
        return state.invariant_proportion
    raise AssertionError(f"unexpected partition target: {target_name}")


def _partition_target_values_equal(
    left: float | Mapping[str, float],
    right: float | Mapping[str, float],
) -> bool:
    if isinstance(left, Mapping):
        if not isinstance(right, Mapping) or set(left) != set(right):
            return False
        return all(
            math.isclose(
                float(left[component_name]),
                float(right[component_name]),
                rel_tol=_PARTITION_STATE_VALUE_TOLERANCE,
                abs_tol=_PARTITION_STATE_VALUE_TOLERANCE,
            )
            for component_name in left
        )
    return math.isclose(
        float(left),
        float(right),
        rel_tol=_PARTITION_STATE_VALUE_TOLERANCE,
        abs_tol=_PARTITION_STATE_VALUE_TOLERANCE,
    )


def _normalize_vector_value(
    *,
    target_name: str,
    raw_value: Mapping[tuple[str, str], float] | Mapping[str, float] | Sequence[float],
) -> dict[str, float]:
    if isinstance(raw_value, Mapping):
        normalized = {
            (
                f"{component_name[0]}{component_name[1]}"
                if isinstance(component_name, tuple)
                else str(component_name)
            ): float(component_value)
            for component_name, component_value in raw_value.items()
        }
        return dict(sorted(normalized.items(), key=lambda item: item[0]))
    expected_component_names = (
        _DNA_EXCHANGEABILITY_COMPONENT_NAMES
        if target_name == "exchangeabilities"
        else _DNA_BASE_FREQUENCY_COMPONENT_NAMES
    )
    if len(raw_value) != len(expected_component_names):
        raise PhylogeneticsError(
            "partition model state encoding received one vector with an unexpected component count",
            code="partition_model_state_vector_component_count_invalid",
            details={
                "target_name": target_name,
                "expected_component_count": len(expected_component_names),
                "observed_component_count": len(raw_value),
            },
        )
    return {
        component_name: float(component_value)
        for component_name, component_value in zip(
            expected_component_names,
            raw_value,
            strict=True,
        )
    }


def _require_partition_scalar_value(
    *,
    model_parameters: BayesianModelParameterState,
    target_name: str,
    group_name: str,
) -> float:
    parameter_key = _partition_parameter_key(
        target_name=target_name,
        group_name=group_name,
    )
    parameter_value = model_parameters.scalar_parameters.get(parameter_key)
    if parameter_value is None:
        raise PhylogeneticsError(
            "partition model state decoding requires one realized scalar value for every active linkage group",
            code="partition_model_state_scalar_value_missing",
            details={
                "target_name": target_name,
                "group_name": group_name,
                "parameter_key": parameter_key,
            },
        )
    return parameter_value


def _require_partition_vector_value(
    *,
    model_parameters: BayesianModelParameterState,
    target_name: str,
    group_name: str,
) -> dict[str, float]:
    parameter_key = _partition_parameter_key(
        target_name=target_name,
        group_name=group_name,
    )
    parameter_value = model_parameters.vector_parameters.get(parameter_key)
    if parameter_value is None:
        raise PhylogeneticsError(
            "partition model state decoding requires one realized vector value for every active linkage group",
            code="partition_model_state_vector_value_missing",
            details={
                "target_name": target_name,
                "group_name": group_name,
                "parameter_key": parameter_key,
            },
        )
    return dict(parameter_value)


def _validate_reserved_parameter_collisions(
    *,
    categorical_parameters: Mapping[str, str],
    scalar_parameters: Mapping[str, float],
    vector_parameters: Mapping[str, Mapping[str, float]],
) -> None:
    reserved_keys = sorted(
        [
            *[
                parameter_name
                for parameter_name in categorical_parameters
                if parameter_name.startswith(_PARTITION_LINKAGE_PREFIX)
            ],
            *[
                parameter_name
                for parameter_name in scalar_parameters
                if parameter_name.startswith(_PARTITION_PARAMETER_PREFIX)
            ],
            *[
                parameter_name
                for parameter_name in vector_parameters
                if parameter_name.startswith(_PARTITION_PARAMETER_PREFIX)
            ],
        ]
    )
    if reserved_keys:
        raise PhylogeneticsError(
            "partition model state encoding does not allow preserved parameters to reuse reserved partition-model keys",
            code="partition_model_state_reserved_key_collision",
            details={"reserved_keys": reserved_keys},
        )


def _validate_unique_partition_names(
    partition_names: Sequence[str],
) -> tuple[str, ...]:
    normalized_partition_names = tuple(
        _validate_partition_name(partition_name) for partition_name in partition_names
    )
    duplicate_partition_names = sorted(
        {
            partition_name
            for partition_name in normalized_partition_names
            if normalized_partition_names.count(partition_name) > 1
        }
    )
    if duplicate_partition_names:
        raise PhylogeneticsError(
            "partition model state encoding requires unique partition names",
            code="partition_model_state_partition_names_duplicated",
            details={"duplicate_partition_names": duplicate_partition_names},
        )
    return normalized_partition_names


def _build_partition_state_lookup(
    *,
    partition_parameter_states: Sequence[PartitionSubstitutionParameterState],
    expected_partition_names: tuple[str, ...],
) -> dict[str, PartitionSubstitutionParameterState]:
    state_by_partition_name: dict[str, PartitionSubstitutionParameterState] = {}
    for state in partition_parameter_states:
        partition_name = _validate_partition_name(state.partition_name)
        if partition_name in state_by_partition_name:
            raise PhylogeneticsError(
                "partition model state encoding received more than one state for the same partition",
                code="partition_model_state_partition_duplicated",
                details={"partition_name": partition_name},
            )
        state_by_partition_name[partition_name] = state
    if set(state_by_partition_name) != set(expected_partition_names):
        raise PhylogeneticsError(
            "partition model state encoding requires states for the exact partition set",
            code="partition_model_state_partition_mismatch",
            details={
                "expected_partition_names": list(expected_partition_names),
                "observed_partition_names": sorted(state_by_partition_name),
            },
        )
    return state_by_partition_name


def _validate_partition_name(partition_name: str) -> str:
    normalized_partition_name = partition_name.strip()
    if not normalized_partition_name:
        raise PhylogeneticsError(
            "partition model state partition name must not be empty",
            code="partition_model_state_partition_name_empty",
        )
    return normalized_partition_name
