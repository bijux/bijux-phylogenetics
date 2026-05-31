from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
import json
import math

from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    PhylogeneticsError,
)


@dataclass(frozen=True, slots=True)
class BayesianStateBranchRow:
    """One explicit branch snapshot inside one Bayesian tree state."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    branch_length: float


@dataclass(frozen=True, slots=True)
class BayesianTreeState:
    """One serialized tree snapshot for Bayesian posterior state tracking."""

    tree_newick: str
    rooted: bool | None
    topology_id: str
    tip_names: list[str]
    branch_rows: list[BayesianStateBranchRow]

    def to_tree(self) -> PhyloTree:
        tree = PhyloTree.from_newick(self.tree_newick)
        tree.rooted = self.rooted
        return tree


@dataclass(frozen=True, slots=True)
class BayesianModelParameterState:
    """One model-parameter snapshot for a Bayesian state."""

    categorical_parameters: dict[str, str] = field(default_factory=dict)
    scalar_parameters: dict[str, float] = field(default_factory=dict)
    vector_parameters: dict[str, dict[str, float]] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class BayesianPriorComponentState:
    """One prior contribution tracked inside a Bayesian state."""

    component_name: str
    family: str | None
    log_prior: float
    parameter_values: dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class BayesianPhylogeneticState:
    """One coherent Bayesian sampler state over tree, parameters, priors, and scores."""

    tree: BayesianTreeState
    model_parameters: BayesianModelParameterState
    prior_components: list[BayesianPriorComponentState]
    total_log_prior: float
    log_likelihood: float
    posterior_log_score: float


def build_bayesian_tree_state(tree: PhyloTree) -> BayesianTreeState:
    """Build one serializable Bayesian tree state from one explicit tree."""
    tree.refresh()
    branch_rows: list[BayesianStateBranchRow] = []
    for _parent, child in tree.iter_edges():
        if child.node_id is None:
            raise InvalidBranchLengthError(
                "bayesian state tree snapshots require stable branch identifiers"
            )
        if child.branch_length is None:
            raise InvalidBranchLengthError(
                "bayesian state tree snapshots require explicit branch lengths on every edge"
            )
        branch_rows.append(
            BayesianStateBranchRow(
                branch_id=child.node_id,
                child_name=child.name,
                descendant_taxa=list(child.descendant_taxa),
                branch_length=_round_float(child.branch_length),
            )
        )
    return BayesianTreeState(
        tree_newick=tree.to_newick(),
        rooted=tree.rooted,
        topology_id=_build_tree_topology_id(tree),
        tip_names=sorted(tree.tip_names),
        branch_rows=branch_rows,
    )


def build_bayesian_model_parameter_state(
    *,
    categorical_parameters: Mapping[str, str] | None = None,
    scalar_parameters: Mapping[str, float] | None = None,
    vector_parameters: Mapping[str, Mapping[str, float]] | None = None,
) -> BayesianModelParameterState:
    """Build one validated Bayesian model-parameter snapshot."""
    normalized_categorical_parameters = {
        parameter_name: _validate_nonblank_name(
            field_name=parameter_name,
            value=parameter_value,
            owner_name="bayesian categorical parameter state",
        )
        for parameter_name, parameter_value in sorted(
            (categorical_parameters or {}).items(),
            key=lambda item: item[0],
        )
    }
    normalized_scalar_parameters = {
        parameter_name: _validate_finite_float(
            parameter_name=parameter_name,
            value=parameter_value,
            owner_name="bayesian scalar parameter state",
        )
        for parameter_name, parameter_value in sorted(
            (scalar_parameters or {}).items(),
            key=lambda item: item[0],
        )
    }
    normalized_vector_parameters = {
        parameter_name: {
            component_name: _validate_finite_float(
                parameter_name=f"{parameter_name}.{component_name}",
                value=component_value,
                owner_name="bayesian vector parameter state",
            )
            for component_name, component_value in sorted(
                parameter_components.items(),
                key=lambda item: item[0],
            )
        }
        for parameter_name, parameter_components in sorted(
            (vector_parameters or {}).items(),
            key=lambda item: item[0],
        )
    }
    _require_nonblank_parameter_names(
        list(normalized_categorical_parameters),
        owner_name="bayesian categorical parameter state",
    )
    _require_nonblank_parameter_names(
        list(normalized_scalar_parameters),
        owner_name="bayesian scalar parameter state",
    )
    _require_nonblank_parameter_names(
        list(normalized_vector_parameters),
        owner_name="bayesian vector parameter state",
    )
    for parameter_name, component_values in normalized_vector_parameters.items():
        if not component_values:
            raise PhylogeneticsError(
                "bayesian vector parameter state requires at least one component per named parameter",
                code="bayesian_vector_parameter_state_empty",
                details={"parameter_name": parameter_name},
            )
        _require_nonblank_parameter_names(
            list(component_values),
            owner_name=f"bayesian vector parameter state for '{parameter_name}'",
        )
    return BayesianModelParameterState(
        categorical_parameters=normalized_categorical_parameters,
        scalar_parameters=normalized_scalar_parameters,
        vector_parameters=normalized_vector_parameters,
    )


def build_bayesian_prior_component_state(
    *,
    component_name: str,
    family: str | None,
    log_prior: float,
    parameter_values: Mapping[str, float] | None = None,
) -> BayesianPriorComponentState:
    """Build one validated prior contribution for Bayesian state tracking."""
    validated_component_name = _validate_nonblank_name(
        field_name="component_name",
        value=component_name,
        owner_name="bayesian prior component state",
    )
    validated_family = (
        _validate_nonblank_name(
            field_name="family",
            value=family,
            owner_name="bayesian prior component state",
        )
        if family is not None
        else None
    )
    normalized_parameter_values = {
        parameter_name: _validate_finite_float(
            parameter_name=parameter_name,
            value=parameter_value,
            owner_name="bayesian prior component parameter values",
        )
        for parameter_name, parameter_value in sorted(
            (parameter_values or {}).items(),
            key=lambda item: item[0],
        )
    }
    _require_nonblank_parameter_names(
        list(normalized_parameter_values),
        owner_name="bayesian prior component parameter values",
    )
    return BayesianPriorComponentState(
        component_name=validated_component_name,
        family=validated_family,
        log_prior=_validate_finite_float(
            parameter_name="log_prior",
            value=log_prior,
            owner_name="bayesian prior component state",
        ),
        parameter_values=normalized_parameter_values,
    )


def build_bayesian_phylogenetic_state(
    *,
    tree: PhyloTree | BayesianTreeState,
    model_parameters: BayesianModelParameterState,
    prior_components: list[BayesianPriorComponentState],
    log_likelihood: float,
) -> BayesianPhylogeneticState:
    """Build one coherent scored Bayesian state."""
    if not prior_components:
        raise PhylogeneticsError(
            "bayesian phylogenetic state requires at least one prior component",
            code="bayesian_phylogenetic_state_prior_components_empty",
        )
    tree_state = (
        tree if isinstance(tree, BayesianTreeState) else build_bayesian_tree_state(tree)
    )
    total_log_prior = math.fsum(component.log_prior for component in prior_components)
    validated_log_likelihood = _validate_finite_float(
        parameter_name="log_likelihood",
        value=log_likelihood,
        owner_name="bayesian phylogenetic state",
    )
    posterior_log_score = total_log_prior + validated_log_likelihood
    return BayesianPhylogeneticState(
        tree=tree_state,
        model_parameters=model_parameters,
        prior_components=list(prior_components),
        total_log_prior=total_log_prior,
        log_likelihood=validated_log_likelihood,
        posterior_log_score=posterior_log_score,
    )


def build_bayesian_phylogenetic_state_from_prior_only_sample(
    sample,
    *,
    log_likelihood: float = 0.0,
) -> BayesianPhylogeneticState:
    """Build one coherent Bayesian state from one prior-only simulation sample."""
    from bijux_phylogenetics.bayesian.prior_sampling import PriorOnlyPhylogeneticSample

    if not isinstance(sample, PriorOnlyPhylogeneticSample):
        raise PhylogeneticsError(
            "bayesian prior-only state builder requires one PriorOnlyPhylogeneticSample",
            code="bayesian_state_prior_only_sample_type_invalid",
        )
    scalar_parameters: dict[str, float] = {}
    if sample.substitution_parameter_state.kappa is not None:
        scalar_parameters["kappa"] = sample.substitution_parameter_state.kappa
    if sample.substitution_parameter_state.gamma_alpha is not None:
        scalar_parameters["gamma-alpha"] = (
            sample.substitution_parameter_state.gamma_alpha
        )
    if sample.substitution_parameter_state.invariant_proportion is not None:
        scalar_parameters["invariant-proportion"] = (
            sample.substitution_parameter_state.invariant_proportion
        )
    vector_parameters: dict[str, dict[str, float]] = {}
    if sample.substitution_parameter_state.exchangeabilities is not None:
        vector_parameters["exchangeabilities"] = (
            sample.substitution_parameter_state.exchangeabilities
        )
    if sample.substitution_parameter_state.base_frequencies is not None:
        vector_parameters["base-frequencies"] = (
            sample.substitution_parameter_state.base_frequencies
        )
    prior_components = [
        build_bayesian_prior_component_state(
            component_name="tree-topology",
            family=sample.tree_topology_prior_family,
            log_prior=sample.topology_log_prior,
        ),
        build_bayesian_prior_component_state(
            component_name="branch-lengths",
            family=sample.branch_length_prior_family,
            log_prior=sample.branch_length_log_prior,
        ),
        *[
            build_bayesian_prior_component_state(
                component_name=f"substitution:{row.target_name}",
                family=row.family,
                log_prior=row.log_prior_contribution,
                parameter_values=row.hyperparameter_values,
            )
            for row in sample.substitution_prior_rows
        ],
    ]
    return build_bayesian_phylogenetic_state(
        tree=PhyloTree.from_newick(sample.tree_newick),
        model_parameters=build_bayesian_model_parameter_state(
            scalar_parameters=scalar_parameters,
            vector_parameters=vector_parameters,
        ),
        prior_components=prior_components,
        log_likelihood=log_likelihood,
    )


def serialize_bayesian_phylogenetic_state(
    state: BayesianPhylogeneticState,
) -> dict[str, object]:
    """Serialize one Bayesian state into one JSON-safe payload."""
    return asdict(state)


def deserialize_bayesian_phylogenetic_state(
    payload: Mapping[str, object],
) -> BayesianPhylogeneticState:
    """Deserialize one Bayesian state payload with score validation."""
    tree_payload = _require_mapping(payload, key="tree")
    model_parameter_payload = _require_mapping(payload, key="model_parameters")
    prior_component_payloads = _require_list(payload, key="prior_components")
    tree_state = BayesianTreeState(
        tree_newick=_require_string(tree_payload, key="tree_newick"),
        rooted=tree_payload.get("rooted"),
        topology_id=_require_string(tree_payload, key="topology_id"),
        tip_names=_require_string_list(tree_payload, key="tip_names"),
        branch_rows=[
            BayesianStateBranchRow(
                branch_id=_require_string(branch_row, key="branch_id"),
                child_name=_optional_string(branch_row.get("child_name")),
                descendant_taxa=_require_string_list(branch_row, key="descendant_taxa"),
                branch_length=_require_float(branch_row, key="branch_length"),
            )
            for branch_row in (
                _require_mapping(item, owner_name="bayesian tree state branch row")
                for item in _require_list(tree_payload, key="branch_rows")
            )
        ],
    )
    model_parameter_state = build_bayesian_model_parameter_state(
        categorical_parameters=_require_string_mapping(
            model_parameter_payload,
            key="categorical_parameters",
        ),
        scalar_parameters=_require_float_mapping(
            model_parameter_payload, key="scalar_parameters"
        ),
        vector_parameters=_require_nested_float_mapping(
            model_parameter_payload,
            key="vector_parameters",
        ),
    )
    prior_components = [
        build_bayesian_prior_component_state(
            component_name=_require_string(
                prior_component_payload, key="component_name"
            ),
            family=_optional_string(prior_component_payload.get("family")),
            log_prior=_require_float(prior_component_payload, key="log_prior"),
            parameter_values=_require_float_mapping(
                prior_component_payload,
                key="parameter_values",
            ),
        )
        for prior_component_payload in (
            _require_mapping(item, owner_name="bayesian prior component payload")
            for item in prior_component_payloads
        )
    ]
    deserialized_state = build_bayesian_phylogenetic_state(
        tree=tree_state,
        model_parameters=model_parameter_state,
        prior_components=prior_components,
        log_likelihood=_require_float(payload, key="log_likelihood"),
    )
    serialized_total_log_prior = _require_float(payload, key="total_log_prior")
    serialized_posterior_log_score = _require_float(payload, key="posterior_log_score")
    if deserialized_state.total_log_prior != serialized_total_log_prior:
        raise PhylogeneticsError(
            "deserialized bayesian state total_log_prior does not match serialized payload",
            code="bayesian_state_total_log_prior_mismatch",
            details={
                "serialized_total_log_prior": serialized_total_log_prior,
                "recomputed_total_log_prior": deserialized_state.total_log_prior,
            },
        )
    if deserialized_state.posterior_log_score != serialized_posterior_log_score:
        raise PhylogeneticsError(
            "deserialized bayesian state posterior_log_score does not match serialized payload",
            code="bayesian_state_posterior_log_score_mismatch",
            details={
                "serialized_posterior_log_score": serialized_posterior_log_score,
                "recomputed_posterior_log_score": deserialized_state.posterior_log_score,
            },
        )
    return deserialized_state


def serialize_bayesian_phylogenetic_state_json(
    state: BayesianPhylogeneticState,
) -> str:
    """Serialize one Bayesian state into canonical JSON text."""
    return json.dumps(
        serialize_bayesian_phylogenetic_state(state),
        indent=2,
        sort_keys=True,
    )


def deserialize_bayesian_phylogenetic_state_json(
    payload: str,
) -> BayesianPhylogeneticState:
    """Deserialize one Bayesian state from JSON text."""
    raw_payload = json.loads(payload)
    if not isinstance(raw_payload, dict):
        raise PhylogeneticsError(
            "bayesian state json payload must decode to one mapping",
            code="bayesian_state_json_payload_type_invalid",
        )
    return deserialize_bayesian_phylogenetic_state(raw_payload)


def _build_tree_topology_id(tree: PhyloTree) -> str:
    internal_clades = sorted(
        {
            tuple(node.descendant_taxa)
            for node in tree.iter_internal_nodes(order="preorder")
            if node is not tree.root and node.descendant_taxa
        }
    )
    rooted_token = (
        "rooted"
        if tree.rooted is True
        else "unrooted"
        if tree.rooted is False
        else "unspecified"
    )
    taxon_token = ",".join(sorted(tree.tip_names))
    clade_token = ";".join(",".join(clade) for clade in internal_clades)
    return f"{rooted_token}:{taxon_token}:{clade_token}"


def _require_mapping(
    payload: Mapping[str, object],
    *,
    key: str | None = None,
    owner_name: str | None = None,
) -> Mapping[str, object]:
    raw_value: object = payload if key is None else payload.get(key)
    if not isinstance(raw_value, Mapping):
        subject = owner_name or f"'{key}'"
        raise PhylogeneticsError(
            f"bayesian state deserialization requires {subject} to be one mapping",
            code="bayesian_state_mapping_required",
        )
    return raw_value


def _require_list(payload: Mapping[str, object], *, key: str) -> list[object]:
    raw_value = payload.get(key)
    if not isinstance(raw_value, list):
        raise PhylogeneticsError(
            f"bayesian state deserialization requires '{key}' to be one list",
            code="bayesian_state_list_required",
        )
    return raw_value


def _require_string(payload: Mapping[str, object], *, key: str) -> str:
    raw_value = payload.get(key)
    if not isinstance(raw_value, str):
        raise PhylogeneticsError(
            f"bayesian state deserialization requires '{key}' to be one string",
            code="bayesian_state_string_required",
        )
    return raw_value


def _optional_string(raw_value: object) -> str | None:
    if raw_value is None:
        return None
    if not isinstance(raw_value, str):
        raise PhylogeneticsError(
            "bayesian state deserialization requires optional string fields to be strings when present",
            code="bayesian_state_optional_string_required",
        )
    return raw_value


def _require_string_list(payload: Mapping[str, object], *, key: str) -> list[str]:
    raw_values = _require_list(payload, key=key)
    string_values: list[str] = []
    for raw_value in raw_values:
        if not isinstance(raw_value, str):
            raise PhylogeneticsError(
                f"bayesian state deserialization requires every '{key}' item to be one string",
                code="bayesian_state_string_list_required",
            )
        string_values.append(raw_value)
    return string_values


def _require_float(payload: Mapping[str, object], *, key: str) -> float:
    raw_value = payload.get(key)
    if isinstance(raw_value, bool) or not isinstance(raw_value, (int, float)):
        raise PhylogeneticsError(
            f"bayesian state deserialization requires '{key}' to be one number",
            code="bayesian_state_float_required",
        )
    return float(raw_value)


def _require_float_mapping(
    payload: Mapping[str, object],
    *,
    key: str,
) -> dict[str, float]:
    raw_mapping = _require_mapping(payload, key=key)
    normalized_mapping: dict[str, float] = {}
    for item_key, item_value in sorted(
        raw_mapping.items(), key=lambda item: str(item[0])
    ):
        if not isinstance(item_key, str):
            raise PhylogeneticsError(
                f"bayesian state deserialization requires '{key}' mapping keys to be strings",
                code="bayesian_state_float_mapping_key_invalid",
            )
        if isinstance(item_value, bool) or not isinstance(item_value, (int, float)):
            raise PhylogeneticsError(
                f"bayesian state deserialization requires '{key}' mapping values to be numeric",
                code="bayesian_state_float_mapping_value_invalid",
            )
        normalized_mapping[item_key] = float(item_value)
    return normalized_mapping


def _require_string_mapping(
    payload: Mapping[str, object],
    *,
    key: str,
) -> dict[str, str]:
    raw_mapping = _require_mapping(payload, key=key)
    normalized_mapping: dict[str, str] = {}
    for item_key, item_value in sorted(
        raw_mapping.items(), key=lambda item: str(item[0])
    ):
        if not isinstance(item_key, str):
            raise PhylogeneticsError(
                f"bayesian state deserialization requires '{key}' mapping keys to be strings",
                code="bayesian_state_string_mapping_key_invalid",
            )
        if not isinstance(item_value, str):
            raise PhylogeneticsError(
                f"bayesian state deserialization requires '{key}' mapping values to be strings",
                code="bayesian_state_string_mapping_value_invalid",
            )
        normalized_mapping[item_key] = item_value
    return normalized_mapping


def _require_nested_float_mapping(
    payload: Mapping[str, object],
    *,
    key: str,
) -> dict[str, dict[str, float]]:
    raw_mapping = _require_mapping(payload, key=key)
    normalized_mapping: dict[str, dict[str, float]] = {}
    for parameter_name, component_values in sorted(
        raw_mapping.items(),
        key=lambda item: str(item[0]),
    ):
        if not isinstance(parameter_name, str):
            raise PhylogeneticsError(
                f"bayesian state deserialization requires '{key}' mapping keys to be strings",
                code="bayesian_state_nested_mapping_key_invalid",
            )
        if not isinstance(component_values, Mapping):
            raise PhylogeneticsError(
                f"bayesian state deserialization requires '{key}' mapping values to be nested mappings",
                code="bayesian_state_nested_mapping_value_invalid",
            )
        normalized_mapping[parameter_name] = _require_float_mapping(
            {"value": component_values},
            key="value",
        )
    return normalized_mapping


def _validate_nonblank_name(
    *,
    field_name: str,
    value: str,
    owner_name: str,
) -> str:
    if not value.strip():
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be non-blank",
            code="bayesian_state_nonblank_name_required",
            details={"field_name": field_name},
        )
    return value


def _require_nonblank_parameter_names(
    parameter_names: list[str],
    *,
    owner_name: str,
) -> None:
    blank_parameter_names = sorted(
        parameter_name
        for parameter_name in parameter_names
        if not parameter_name.strip()
    )
    if blank_parameter_names:
        raise PhylogeneticsError(
            f"{owner_name} does not allow blank parameter names",
            code="bayesian_state_blank_parameter_name",
            details={"blank_parameter_names": blank_parameter_names},
        )


def _validate_finite_float(
    *,
    parameter_name: str,
    value: float,
    owner_name: str,
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise PhylogeneticsError(
            f"{owner_name} requires '{parameter_name}' to be numeric",
            code="bayesian_state_numeric_value_required",
            details={"parameter_name": parameter_name},
        )
    normalized_value = float(value)
    if not math.isfinite(normalized_value):
        raise PhylogeneticsError(
            f"{owner_name} requires '{parameter_name}' to be finite",
            code="bayesian_state_finite_value_required",
            details={"parameter_name": parameter_name},
        )
    return normalized_value


def _round_float(value: float) -> float:
    return float(format(value, ".15g"))
