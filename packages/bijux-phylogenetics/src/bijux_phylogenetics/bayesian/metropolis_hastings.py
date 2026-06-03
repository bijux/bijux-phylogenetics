from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass
import json
import math
import random

from bijux_phylogenetics.bayesian.partition_model_priors import (
    PARTITION_MODEL_PRIOR_TARGETS,
    PartitionSubstitutionModelDefinition,
    build_partition_parameter_linkage_plan,
)
from bijux_phylogenetics.bayesian.partition_model_state import (
    build_partition_model_parameter_state,
    resolve_partition_parameter_linkage_plan_from_model_parameters,
    resolve_partition_parameter_states_from_model_parameters,
    strip_partition_model_parameter_state,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianModelParameterState,
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
    build_bayesian_phylogenetic_state,
    deserialize_bayesian_phylogenetic_state,
    serialize_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    APE_ULTRAMETRIC_TOLERANCE,
)
from bijux_phylogenetics.phylo.likelihood.dna import validate_positive_kappa
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    parameterize_dna_base_frequency_simplex,
    parameterize_dna_exchangeability_simplex,
    resolve_dna_base_frequency_simplex_from_unconstrained,
    resolve_dna_exchangeability_simplex_from_unconstrained,
)
from bijux_phylogenetics.phylo.likelihood.gamma import (
    validate_discrete_gamma_alpha,
)
from bijux_phylogenetics.phylo.likelihood.invariant import (
    validate_invariant_proportion,
)
from bijux_phylogenetics.phylo.topology.clades import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.rooted_nni import (
    RootedNniMoveCandidate,
    apply_rooted_nni_move,
    iter_rooted_nni_move_candidates,
    validate_rooted_nni_tree,
)
from bijux_phylogenetics.phylo.topology.rooted_spr import (
    enumerate_rooted_spr_neighbors,
    validate_rooted_spr_tree,
)
from bijux_phylogenetics.phylo.topology.rooted_tbr import (
    enumerate_rooted_tbr_neighbors,
    validate_rooted_tbr_tree,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import (
    InvalidAlignmentError,
    PhylogeneticsError,
)

BayesianPriorUpdate = Callable[
    [BayesianPhylogeneticState],
    list[BayesianPriorComponentState],
]
BayesianLikelihoodUpdate = Callable[[BayesianPhylogeneticState], float]
BayesianStateProposal = Callable[
    [BayesianPhylogeneticState, random.Random],
    "MetropolisHastingsProposal",
]


@dataclass(frozen=True, slots=True)
class MetropolisHastingsProposal:
    """One validated proposal result for a Metropolis-Hastings iteration."""

    changed_fields: tuple[str, ...]
    log_forward_density: float
    log_reverse_density: float
    is_valid: bool
    invalid_reason: str | None
    proposed_tree: object | None = None
    proposed_model_parameters: BayesianModelParameterState | None = None

    @property
    def log_hastings_ratio(self) -> float:
        return self.log_reverse_density - self.log_forward_density


@dataclass(frozen=True, slots=True)
class MetropolisHastingsStepRow:
    """One iteration-level Metropolis-Hastings trace row."""

    iteration_index: int
    proposal_changed_fields: tuple[str, ...]
    proposal_valid: bool
    proposal_invalid_reason: str | None
    log_forward_density: float
    log_reverse_density: float
    accepted: bool
    log_hastings_ratio: float
    current_posterior_log_score: float
    proposed_posterior_log_score: float | None
    log_acceptance_ratio: float | None
    recorded_posterior_log_score: float


@dataclass(frozen=True, slots=True)
class MetropolisHastingsRunReport:
    """One completed Metropolis-Hastings chain report."""

    iteration_count: int
    sample_every: int
    seed: int
    accepted_count: int
    rejected_count: int
    acceptance_rate: float
    initial_state: BayesianPhylogeneticState
    final_state: BayesianPhylogeneticState
    sampled_states: list[BayesianPhylogeneticState]
    step_rows: list[MetropolisHastingsStepRow]


@dataclass(frozen=True, slots=True)
class MetropolisHastingsRandomState:
    """One serialized Python random-generator state for exact chain resumption."""

    version: int
    internal_state: list[int]
    gaussian_spare: float | None


@dataclass(frozen=True, slots=True)
class MetropolisHastingsCheckpoint:
    """One resumable Metropolis-Hastings sampler checkpoint."""

    iteration_count: int
    sample_every: int
    seed: int
    completed_iteration_count: int
    accepted_count: int
    initial_state: BayesianPhylogeneticState
    current_state: BayesianPhylogeneticState
    sampled_states: list[BayesianPhylogeneticState]
    step_rows: list[MetropolisHastingsStepRow]
    random_state: MetropolisHastingsRandomState


@dataclass(frozen=True, slots=True)
class CheckpointedMetropolisHastingsRunReport:
    """One checkpoint-aware Metropolis-Hastings execution report."""

    iteration_count: int
    completed_iteration_count: int
    sample_every: int
    seed: int
    resumed: bool
    completed: bool
    accepted_count: int
    rejected_count: int
    acceptance_rate: float
    initial_state: BayesianPhylogeneticState
    current_state: BayesianPhylogeneticState
    sampled_states: list[BayesianPhylogeneticState]
    step_rows: list[MetropolisHastingsStepRow]
    checkpoints: list[MetropolisHastingsCheckpoint]

    def to_chain_report(self) -> MetropolisHastingsRunReport:
        """Convert one completed checkpoint-aware execution into a chain report."""
        if self.completed is False:
            raise PhylogeneticsError(
                "checkpointed metropolis-hastings execution can convert to one chain report only after all iterations complete",
                code="metropolis_hastings_checkpointed_execution_incomplete",
                details={
                    "iteration_count": self.iteration_count,
                    "completed_iteration_count": self.completed_iteration_count,
                },
            )
        return MetropolisHastingsRunReport(
            iteration_count=self.iteration_count,
            sample_every=self.sample_every,
            seed=self.seed,
            accepted_count=self.accepted_count,
            rejected_count=self.rejected_count,
            acceptance_rate=self.acceptance_rate,
            initial_state=self.initial_state,
            final_state=self.current_state,
            sampled_states=list(self.sampled_states),
            step_rows=list(self.step_rows),
        )


@dataclass(frozen=True, slots=True)
class _NodeHeightSlideCandidate:
    node_id: str
    current_height: float
    lower_height_bound: float
    upper_height_bound: float


@dataclass(frozen=True, slots=True)
class _ReversibleRootedTbrNeighbor:
    neighbor_row: object
    proposed_tree: PhyloTree
    reverse_neighbor_row: object


_REVERSIBLE_JUMP_MODEL_SWITCH_FAMILIES = ("nucleotide-substitution-model",)


def build_metropolis_hastings_proposal(
    *,
    changed_fields: list[str] | tuple[str, ...],
    log_forward_density: float,
    log_reverse_density: float,
    is_valid: bool,
    invalid_reason: str | None = None,
    proposed_tree: object | None = None,
    proposed_model_parameters: BayesianModelParameterState | None = None,
) -> MetropolisHastingsProposal:
    """Build one validated Metropolis-Hastings proposal result."""
    return _validate_metropolis_hastings_proposal(
        MetropolisHastingsProposal(
            changed_fields=tuple(changed_fields),
            log_forward_density=log_forward_density,
            log_reverse_density=log_reverse_density,
            is_valid=is_valid,
            invalid_reason=invalid_reason,
            proposed_tree=proposed_tree,
            proposed_model_parameters=proposed_model_parameters,
        )
    )


def build_metropolis_hastings_random_state(
    *,
    version: int,
    internal_state: Sequence[int],
    gaussian_spare: float | None,
) -> MetropolisHastingsRandomState:
    """Build one validated serialized random-generator state."""
    validated_version = _validate_integer_field(
        value=version,
        field_name="version",
        owner_name="metropolis-hastings random state",
    )
    validated_internal_state = _validate_integer_sequence(
        values=internal_state,
        field_name="internal_state",
        owner_name="metropolis-hastings random state",
        minimum_length=1,
    )
    validated_gaussian_spare = (
        _validate_finite_float(
            value=gaussian_spare,
            field_name="gaussian_spare",
            owner_name="metropolis-hastings random state",
        )
        if gaussian_spare is not None
        else None
    )
    return MetropolisHastingsRandomState(
        version=validated_version,
        internal_state=validated_internal_state,
        gaussian_spare=validated_gaussian_spare,
    )


def build_metropolis_hastings_checkpoint(
    *,
    iteration_count: int,
    sample_every: int,
    seed: int,
    completed_iteration_count: int,
    accepted_count: int,
    initial_state: BayesianPhylogeneticState,
    current_state: BayesianPhylogeneticState,
    sampled_states: Sequence[BayesianPhylogeneticState],
    step_rows: Sequence[MetropolisHastingsStepRow],
    random_state: MetropolisHastingsRandomState,
) -> MetropolisHastingsCheckpoint:
    """Build one validated sampler checkpoint for exact chain resumption."""
    validated_iteration_count = _validate_positive_integer(
        value=iteration_count,
        field_name="iteration_count",
        owner_name="metropolis-hastings checkpoint",
    )
    validated_sample_every = _validate_positive_integer(
        value=sample_every,
        field_name="sample_every",
        owner_name="metropolis-hastings checkpoint",
    )
    validated_seed = _validate_integer_seed(seed)
    validated_completed_iteration_count = _validate_nonnegative_integer(
        value=completed_iteration_count,
        field_name="completed_iteration_count",
        owner_name="metropolis-hastings checkpoint",
    )
    if validated_completed_iteration_count > validated_iteration_count:
        raise PhylogeneticsError(
            "metropolis-hastings checkpoint requires 'completed_iteration_count' to be less than or equal to 'iteration_count'",
            code="metropolis_hastings_checkpoint_iteration_range_invalid",
            details={
                "iteration_count": validated_iteration_count,
                "completed_iteration_count": validated_completed_iteration_count,
            },
        )
    validated_accepted_count = _validate_nonnegative_integer(
        value=accepted_count,
        field_name="accepted_count",
        owner_name="metropolis-hastings checkpoint",
    )
    if validated_accepted_count > validated_completed_iteration_count:
        raise PhylogeneticsError(
            "metropolis-hastings checkpoint requires 'accepted_count' to be less than or equal to 'completed_iteration_count'",
            code="metropolis_hastings_checkpoint_acceptance_count_invalid",
            details={
                "accepted_count": validated_accepted_count,
                "completed_iteration_count": validated_completed_iteration_count,
            },
        )
    if not isinstance(initial_state, BayesianPhylogeneticState):
        raise PhylogeneticsError(
            "metropolis-hastings checkpoint requires one BayesianPhylogeneticState initial_state",
            code="metropolis_hastings_checkpoint_initial_state_type_invalid",
        )
    if not isinstance(current_state, BayesianPhylogeneticState):
        raise PhylogeneticsError(
            "metropolis-hastings checkpoint requires one BayesianPhylogeneticState current_state",
            code="metropolis_hastings_checkpoint_current_state_type_invalid",
        )
    validated_sampled_states = _validate_sampled_states(
        sampled_states=sampled_states,
        initial_state=initial_state,
        current_state=current_state,
        completed_iteration_count=validated_completed_iteration_count,
        sample_every=validated_sample_every,
        owner_name="metropolis-hastings checkpoint",
    )
    validated_step_rows = _validate_step_rows(
        step_rows=step_rows,
        completed_iteration_count=validated_completed_iteration_count,
        current_state=current_state,
        owner_name="metropolis-hastings checkpoint",
    )
    if validated_accepted_count != sum(
        1 for step_row in validated_step_rows if step_row.accepted
    ):
        raise PhylogeneticsError(
            "metropolis-hastings checkpoint requires 'accepted_count' to match accepted step rows",
            code="metropolis_hastings_checkpoint_accepted_count_mismatch",
            details={
                "accepted_count": validated_accepted_count,
                "accepted_step_row_count": sum(
                    1 for step_row in validated_step_rows if step_row.accepted
                ),
            },
        )
    if not isinstance(random_state, MetropolisHastingsRandomState):
        raise PhylogeneticsError(
            "metropolis-hastings checkpoint requires one MetropolisHastingsRandomState",
            code="metropolis_hastings_checkpoint_random_state_type_invalid",
        )
    if validated_completed_iteration_count == 0 and current_state != initial_state:
        raise PhylogeneticsError(
            "metropolis-hastings checkpoint requires zero-iteration checkpoints to preserve the initial state as the current state",
            code="metropolis_hastings_checkpoint_zero_iteration_state_invalid",
        )
    return MetropolisHastingsCheckpoint(
        iteration_count=validated_iteration_count,
        sample_every=validated_sample_every,
        seed=validated_seed,
        completed_iteration_count=validated_completed_iteration_count,
        accepted_count=validated_accepted_count,
        initial_state=initial_state,
        current_state=current_state,
        sampled_states=validated_sampled_states,
        step_rows=validated_step_rows,
        random_state=random_state,
    )


def serialize_metropolis_hastings_checkpoint(
    checkpoint: MetropolisHastingsCheckpoint,
) -> dict[str, object]:
    """Serialize one Metropolis-Hastings checkpoint into one JSON-safe payload."""
    if not isinstance(checkpoint, MetropolisHastingsCheckpoint):
        raise PhylogeneticsError(
            "metropolis-hastings checkpoint serialization requires one MetropolisHastingsCheckpoint",
            code="metropolis_hastings_checkpoint_serialization_type_invalid",
        )
    return {
        "iteration_count": checkpoint.iteration_count,
        "sample_every": checkpoint.sample_every,
        "seed": checkpoint.seed,
        "completed_iteration_count": checkpoint.completed_iteration_count,
        "accepted_count": checkpoint.accepted_count,
        "initial_state": serialize_bayesian_phylogenetic_state(
            checkpoint.initial_state
        ),
        "current_state": serialize_bayesian_phylogenetic_state(
            checkpoint.current_state
        ),
        "sampled_states": [
            serialize_bayesian_phylogenetic_state(sampled_state)
            for sampled_state in checkpoint.sampled_states
        ],
        "step_rows": [
            {
                **asdict(step_row),
                "proposal_changed_fields": list(step_row.proposal_changed_fields),
            }
            for step_row in checkpoint.step_rows
        ],
        "random_state": asdict(checkpoint.random_state),
    }


def deserialize_metropolis_hastings_checkpoint(
    payload: Mapping[str, object],
) -> MetropolisHastingsCheckpoint:
    """Deserialize one Metropolis-Hastings checkpoint payload with consistency checks."""
    initial_state_payload = _require_mapping_payload(
        payload,
        key="initial_state",
        owner_name="metropolis-hastings checkpoint deserialization",
    )
    current_state_payload = _require_mapping_payload(
        payload,
        key="current_state",
        owner_name="metropolis-hastings checkpoint deserialization",
    )
    sampled_state_payloads = _require_list_payload(
        payload,
        key="sampled_states",
        owner_name="metropolis-hastings checkpoint deserialization",
    )
    step_row_payloads = _require_list_payload(
        payload,
        key="step_rows",
        owner_name="metropolis-hastings checkpoint deserialization",
    )
    random_state_payload = _require_mapping_payload(
        payload,
        key="random_state",
        owner_name="metropolis-hastings checkpoint deserialization",
    )
    return build_metropolis_hastings_checkpoint(
        iteration_count=_require_integer_payload(
            payload,
            key="iteration_count",
            owner_name="metropolis-hastings checkpoint deserialization",
        ),
        sample_every=_require_integer_payload(
            payload,
            key="sample_every",
            owner_name="metropolis-hastings checkpoint deserialization",
        ),
        seed=_require_integer_payload(
            payload,
            key="seed",
            owner_name="metropolis-hastings checkpoint deserialization",
        ),
        completed_iteration_count=_require_integer_payload(
            payload,
            key="completed_iteration_count",
            owner_name="metropolis-hastings checkpoint deserialization",
        ),
        accepted_count=_require_integer_payload(
            payload,
            key="accepted_count",
            owner_name="metropolis-hastings checkpoint deserialization",
        ),
        initial_state=deserialize_bayesian_phylogenetic_state(initial_state_payload),
        current_state=deserialize_bayesian_phylogenetic_state(current_state_payload),
        sampled_states=[
            deserialize_bayesian_phylogenetic_state(
                _require_mapping_payload(
                    sampled_state_payload,
                    owner_name="metropolis-hastings sampled state payload",
                )
            )
            for sampled_state_payload in sampled_state_payloads
        ],
        step_rows=[
            _deserialize_metropolis_hastings_step_row(
                _require_mapping_payload(
                    step_row_payload,
                    owner_name="metropolis-hastings step row payload",
                )
            )
            for step_row_payload in step_row_payloads
        ],
        random_state=build_metropolis_hastings_random_state(
            version=_require_integer_payload(
                random_state_payload,
                key="version",
                owner_name="metropolis-hastings random state deserialization",
            ),
            internal_state=_require_integer_list_payload(
                random_state_payload,
                key="internal_state",
                owner_name="metropolis-hastings random state deserialization",
            ),
            gaussian_spare=_optional_float_payload(
                random_state_payload,
                key="gaussian_spare",
                owner_name="metropolis-hastings random state deserialization",
            ),
        ),
    )


def serialize_metropolis_hastings_checkpoint_json(
    checkpoint: MetropolisHastingsCheckpoint,
) -> str:
    """Serialize one Metropolis-Hastings checkpoint into canonical JSON text."""
    return json.dumps(
        serialize_metropolis_hastings_checkpoint(checkpoint),
        indent=2,
        sort_keys=True,
    )


def deserialize_metropolis_hastings_checkpoint_json(
    payload: str,
) -> MetropolisHastingsCheckpoint:
    """Deserialize one Metropolis-Hastings checkpoint from JSON text."""
    raw_payload = json.loads(payload)
    if not isinstance(raw_payload, dict):
        raise PhylogeneticsError(
            "metropolis-hastings checkpoint json payload must decode to one mapping",
            code="metropolis_hastings_checkpoint_json_payload_type_invalid",
        )
    return deserialize_metropolis_hastings_checkpoint(raw_payload)


def list_reversible_jump_model_switch_families() -> tuple[str, ...]:
    """List declared reversible-jump model-switch proposal families."""
    return _REVERSIBLE_JUMP_MODEL_SWITCH_FAMILIES


def validate_reversible_jump_model_switch_family(model_family: str) -> str:
    """Validate one declared reversible-jump model-switch proposal family."""
    validated_model_family = _validate_parameter_name(
        value=model_family,
        field_name="model_family",
        owner_name="reversible-jump model-switch proposal",
    )
    if validated_model_family not in _REVERSIBLE_JUMP_MODEL_SWITCH_FAMILIES:
        raise PhylogeneticsError(
            "reversible-jump model-switch proposal requires one declared model family",
            code="reversible_jump_model_switch_family_invalid",
            details={
                "model_family": model_family,
                "allowed_model_families": list(_REVERSIBLE_JUMP_MODEL_SWITCH_FAMILIES),
            },
        )
    return validated_model_family


def propose_reversible_jump_model_switch_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    *,
    model_family: str = "nucleotide-substitution-model",
    log_kappa_standard_deviation: float = 0.5,
) -> MetropolisHastingsProposal:
    """Propose one reversible-jump model switch within one declared family."""
    validated_model_family = validate_reversible_jump_model_switch_family(model_family)
    if validated_model_family != "nucleotide-substitution-model":
        raise AssertionError(
            "reversible-jump model-switch proposal reached one unsupported family handler"
        )
    validated_log_kappa_standard_deviation = _validate_positive_finite_float(
        value=log_kappa_standard_deviation,
        field_name="log_kappa_standard_deviation",
        owner_name="reversible-jump model-switch proposal",
    )
    current_model_name = current_state.model_parameters.categorical_parameters.get(
        "substitution-model"
    )
    changed_fields = (
        "categorical_parameters.substitution-model",
        "scalar_parameters.kappa",
    )
    if current_model_name is None:
        return build_metropolis_hastings_proposal(
            changed_fields=changed_fields,
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "reversible-jump model-switch proposal requires one "
                "'substitution-model' categorical parameter"
            ),
        )
    normalized_current_model_name = current_model_name.strip().upper()
    if normalized_current_model_name not in {"JC69", "K80"}:
        return build_metropolis_hastings_proposal(
            changed_fields=changed_fields,
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "reversible-jump model-switch proposal supports only JC69 and K80 "
                "within the nucleotide-substitution-model family"
            ),
        )
    current_tree = current_state.tree.to_tree()
    current_tree.rooted = current_state.tree.rooted
    proposed_categorical_parameters = dict(
        current_state.model_parameters.categorical_parameters
    )
    proposed_scalar_parameters = dict(current_state.model_parameters.scalar_parameters)
    proposed_vector_parameters = {
        parameter_name: dict(component_values)
        for parameter_name, component_values in current_state.model_parameters.vector_parameters.items()
    }
    if normalized_current_model_name == "JC69":
        if "kappa" in proposed_scalar_parameters:
            return build_metropolis_hastings_proposal(
                changed_fields=changed_fields,
                log_forward_density=0.0,
                log_reverse_density=0.0,
                is_valid=False,
                invalid_reason=(
                    "reversible-jump model-switch proposal requires JC69 states to omit "
                    "the standalone 'kappa' scalar parameter"
                ),
            )
        try:
            proposed_kappa = math.exp(
                rng.gauss(0.0, validated_log_kappa_standard_deviation)
            )
        except OverflowError:
            return build_metropolis_hastings_proposal(
                changed_fields=changed_fields,
                log_forward_density=0.0,
                log_reverse_density=0.0,
                is_valid=False,
                invalid_reason="reversible-jump model-switch proposal overflowed while drawing K80 kappa",
            )
        try:
            validated_proposed_kappa = float(
                format(
                    validate_positive_kappa(
                        proposed_kappa,
                        model_name="reversible-jump model-switch proposal",
                    ),
                    ".15g",
                )
            )
        except ValueError as error:
            return build_metropolis_hastings_proposal(
                changed_fields=changed_fields,
                log_forward_density=0.0,
                log_reverse_density=0.0,
                is_valid=False,
                invalid_reason=str(error),
            )
        proposed_categorical_parameters["substitution-model"] = "K80"
        proposed_scalar_parameters["kappa"] = validated_proposed_kappa
        log_forward_density = _lognormal_positive_draw_density(
            proposed_value=validated_proposed_kappa,
            log_standard_deviation=validated_log_kappa_standard_deviation,
        )
        log_reverse_density = 0.0
    else:
        current_kappa = proposed_scalar_parameters.get("kappa")
        if current_kappa is None:
            return build_metropolis_hastings_proposal(
                changed_fields=changed_fields,
                log_forward_density=0.0,
                log_reverse_density=0.0,
                is_valid=False,
                invalid_reason=(
                    "reversible-jump model-switch proposal requires K80 states to include "
                    "one positive 'kappa' scalar parameter"
                ),
            )
        try:
            validated_current_kappa = float(
                format(
                    validate_positive_kappa(
                        current_kappa,
                        model_name="reversible-jump model-switch proposal",
                    ),
                    ".15g",
                )
            )
        except ValueError as error:
            return build_metropolis_hastings_proposal(
                changed_fields=changed_fields,
                log_forward_density=0.0,
                log_reverse_density=0.0,
                is_valid=False,
                invalid_reason=str(error),
            )
        proposed_categorical_parameters["substitution-model"] = "JC69"
        proposed_scalar_parameters.pop("kappa", None)
        log_forward_density = 0.0
        log_reverse_density = _lognormal_positive_draw_density(
            proposed_value=validated_current_kappa,
            log_standard_deviation=validated_log_kappa_standard_deviation,
        )
    return build_metropolis_hastings_proposal(
        changed_fields=changed_fields,
        log_forward_density=log_forward_density,
        log_reverse_density=log_reverse_density,
        is_valid=True,
        proposed_tree=current_tree,
        proposed_model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters=proposed_categorical_parameters,
            scalar_parameters=proposed_scalar_parameters,
            vector_parameters=proposed_vector_parameters,
        ),
    )


def propose_partition_linking_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    *,
    partition_models: Sequence[PartitionSubstitutionModelDefinition],
    target_name: str | None = None,
) -> MetropolisHastingsProposal:
    """Propose one linkage toggle between linked and unlinked partition states."""
    validated_partition_models = _validate_partition_linking_partition_models(
        partition_models
    )
    validated_target_name = (
        _validate_partition_linking_target_name(target_name)
        if target_name is not None
        else None
    )
    eligible_target_names = _eligible_partition_linking_targets(
        validated_partition_models
    )
    if not eligible_target_names:
        return build_metropolis_hastings_proposal(
            changed_fields=("categorical_parameters.partition-linkage",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "partition-linking proposal requires at least one target used by "
                "more than one partition"
            ),
        )
    if validated_target_name is None:
        selected_target_name = eligible_target_names[
            rng.randrange(len(eligible_target_names))
        ]
        selection_log_density = -math.log(len(eligible_target_names))
    else:
        selected_target_name = validated_target_name
        if selected_target_name not in eligible_target_names:
            return build_metropolis_hastings_proposal(
                changed_fields=(
                    f"categorical_parameters.partition-linkage:{selected_target_name}",
                ),
                log_forward_density=0.0,
                log_reverse_density=0.0,
                is_valid=False,
                invalid_reason=(
                    "partition-linking proposal can switch only targets used by "
                    "more than one partition"
                ),
            )
        selection_log_density = 0.0
    current_tree = current_state.tree.to_tree()
    current_tree.rooted = current_state.tree.rooted
    partition_names = tuple(
        partition_model.partition_name for partition_model in validated_partition_models
    )
    try:
        current_linkage_plan = (
            resolve_partition_parameter_linkage_plan_from_model_parameters(
                model_parameters=current_state.model_parameters,
                partition_names=partition_names,
            )
        )
        current_partition_states = (
            resolve_partition_parameter_states_from_model_parameters(
                model_parameters=current_state.model_parameters,
                partition_models=validated_partition_models,
                linkage_plan=current_linkage_plan,
            )
        )
    except PhylogeneticsError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=(
                f"categorical_parameters.partition-linkage:{selected_target_name}",
            ),
            log_forward_density=selection_log_density,
            log_reverse_density=selection_log_density,
            is_valid=False,
            invalid_reason=str(error),
        )
    required_partition_names = tuple(
        partition_model.partition_name
        for partition_model in validated_partition_models
        if selected_target_name in partition_model.required_targets()
    )
    current_target_groups = current_linkage_plan.groups_for_target(selected_target_name)
    target_linkage_state = _classify_partition_linkage_state(
        group_by_partition_name={
            partition_name: current_target_groups[partition_name]
            for partition_name in required_partition_names
        }
    )
    if target_linkage_state is None:
        return build_metropolis_hastings_proposal(
            changed_fields=(
                f"categorical_parameters.partition-linkage:{selected_target_name}",
            ),
            log_forward_density=selection_log_density,
            log_reverse_density=selection_log_density,
            is_valid=False,
            invalid_reason=(
                "partition-linking proposal currently supports only fully linked or "
                "fully unlinked states for one target"
            ),
        )
    proposed_group_assignments = {
        candidate_target_name: current_linkage_plan.groups_for_target(
            candidate_target_name
        )
        for candidate_target_name in PARTITION_MODEL_PRIOR_TARGETS
    }
    if target_linkage_state == "linked":
        for partition_name in required_partition_names:
            proposed_group_assignments[selected_target_name][partition_name] = (
                partition_name
            )
    else:
        for partition_name in required_partition_names:
            proposed_group_assignments[selected_target_name][partition_name] = (
                f"{selected_target_name}-shared"
            )
    stripped_model_parameters = strip_partition_model_parameter_state(
        current_state.model_parameters
    )
    proposed_linkage_plan = build_partition_parameter_linkage_plan(
        partition_names=partition_names,
        group_assignments=proposed_group_assignments,
    )
    try:
        proposed_model_parameters = build_partition_model_parameter_state(
            partition_models=validated_partition_models,
            linkage_plan=proposed_linkage_plan,
            partition_parameter_states=current_partition_states,
            preserved_categorical_parameters=(
                stripped_model_parameters.categorical_parameters
            ),
            preserved_scalar_parameters=stripped_model_parameters.scalar_parameters,
            preserved_vector_parameters=stripped_model_parameters.vector_parameters,
        )
    except PhylogeneticsError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=(
                f"categorical_parameters.partition-linkage:{selected_target_name}",
            ),
            log_forward_density=selection_log_density,
            log_reverse_density=selection_log_density,
            is_valid=False,
            invalid_reason=str(error),
        )
    changed_fields = _diff_model_parameter_changed_fields(
        current_state.model_parameters,
        proposed_model_parameters,
    )
    return build_metropolis_hastings_proposal(
        changed_fields=changed_fields,
        log_forward_density=selection_log_density,
        log_reverse_density=selection_log_density,
        is_valid=True,
        proposed_tree=current_tree,
        proposed_model_parameters=proposed_model_parameters,
    )


def propose_branch_length_scaling_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    *,
    log_scale_standard_deviation: float,
) -> MetropolisHastingsProposal:
    """Propose one multiplicative change to one positive branch length."""
    validated_log_scale_standard_deviation = _validate_positive_finite_float(
        value=log_scale_standard_deviation,
        field_name="log_scale_standard_deviation",
        owner_name="branch-length scaling proposal",
    )
    current_tree = current_state.tree.to_tree()
    current_tree.rooted = current_state.tree.rooted
    candidate_branch_rows = [
        (child.node_id, float(child.branch_length))
        for _parent, child in current_tree.iter_edges()
        if child.node_id is not None
        and child.branch_length is not None
        and math.isfinite(child.branch_length)
        and child.branch_length > 0.0
    ]
    if not candidate_branch_rows:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.branch_lengths",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="tree has no strictly positive finite branch lengths",
        )
    selected_branch_id, current_branch_length = candidate_branch_rows[
        rng.randrange(len(candidate_branch_rows))
    ]
    try:
        scale_factor = math.exp(rng.gauss(0.0, validated_log_scale_standard_deviation))
    except OverflowError:
        return build_metropolis_hastings_proposal(
            changed_fields=(f"tree.branch_length:{selected_branch_id}",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="branch-length scaling factor overflowed",
        )
    proposed_branch_length = current_branch_length * scale_factor
    if not math.isfinite(proposed_branch_length) or proposed_branch_length <= 0.0:
        return build_metropolis_hastings_proposal(
            changed_fields=(f"tree.branch_length:{selected_branch_id}",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="branch-length scaling produced a non-positive finite branch length",
        )
    proposed_tree = _copy_tree_with_scaled_branch_length(
        current_tree=current_tree,
        branch_id=selected_branch_id,
        proposed_branch_length=proposed_branch_length,
    )
    branch_selection_log_probability = -math.log(len(candidate_branch_rows))
    log_forward_density = branch_selection_log_probability + _lognormal_scaling_density(
        current_branch_length=current_branch_length,
        proposed_branch_length=proposed_branch_length,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    log_reverse_density = branch_selection_log_probability + _lognormal_scaling_density(
        current_branch_length=proposed_branch_length,
        proposed_branch_length=current_branch_length,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    return build_metropolis_hastings_proposal(
        changed_fields=(f"tree.branch_length:{selected_branch_id}",),
        log_forward_density=log_forward_density,
        log_reverse_density=log_reverse_density,
        is_valid=True,
        proposed_tree=proposed_tree,
        proposed_model_parameters=current_state.model_parameters,
    )


def propose_global_tree_height_scaling_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    *,
    log_scale_standard_deviation: float,
) -> MetropolisHastingsProposal:
    """Propose one coherent multiplicative scale change across every branch."""
    validated_log_scale_standard_deviation = _validate_positive_finite_float(
        value=log_scale_standard_deviation,
        field_name="log_scale_standard_deviation",
        owner_name="global tree-height scaling proposal",
    )
    current_tree = current_state.tree.to_tree()
    current_tree.rooted = current_state.tree.rooted
    current_branch_lengths = [
        float(child.branch_length)
        for _parent, child in current_tree.iter_edges()
        if child.branch_length is not None
    ]
    if not current_branch_lengths:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.branch_lengths",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="global tree-height scaling requires explicit branch lengths on every edge",
        )
    if len(current_branch_lengths) != sum(1 for _ in current_tree.iter_edges()):
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.branch_lengths",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="global tree-height scaling requires explicit branch lengths on every edge",
        )
    if any(
        not math.isfinite(branch_length) or branch_length <= 0.0
        for branch_length in current_branch_lengths
    ):
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.branch_lengths",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="global tree-height scaling requires strictly positive finite branch lengths on every edge",
        )
    try:
        scale_factor = math.exp(rng.gauss(0.0, validated_log_scale_standard_deviation))
    except OverflowError:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.branch_lengths",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="global tree-height scaling factor overflowed",
        )
    if not math.isfinite(scale_factor) or scale_factor <= 0.0:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.branch_lengths",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="global tree-height scaling factor was not strictly positive and finite",
        )
    proposed_tree = _copy_tree_with_globally_scaled_branch_lengths(
        current_tree=current_tree,
        scale_factor=scale_factor,
    )
    proposed_branch_lengths = [
        float(child.branch_length)
        for _parent, child in proposed_tree.iter_edges()
        if child.branch_length is not None
    ]
    log_forward_density = _global_tree_scaling_density(
        current_branch_lengths=current_branch_lengths,
        proposed_branch_lengths=proposed_branch_lengths,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    log_reverse_density = _global_tree_scaling_density(
        current_branch_lengths=proposed_branch_lengths,
        proposed_branch_lengths=current_branch_lengths,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    return build_metropolis_hastings_proposal(
        changed_fields=("tree.branch_lengths",),
        log_forward_density=log_forward_density,
        log_reverse_density=log_reverse_density,
        is_valid=True,
        proposed_tree=proposed_tree,
        proposed_model_parameters=current_state.model_parameters,
    )


def propose_gtr_exchangeability_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    *,
    unconstrained_coordinate_standard_deviation: float,
) -> MetropolisHastingsProposal:
    """Propose one simplex-preserving GTR exchangeability parameter move."""
    validated_coordinate_standard_deviation = _validate_positive_finite_float(
        value=unconstrained_coordinate_standard_deviation,
        field_name="unconstrained_coordinate_standard_deviation",
        owner_name="GTR exchangeability proposal",
    )
    current_exchangeabilities = current_state.model_parameters.vector_parameters.get(
        "exchangeabilities"
    )
    if current_exchangeabilities is None:
        return build_metropolis_hastings_proposal(
            changed_fields=("vector_parameters.exchangeabilities",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="GTR exchangeability proposal requires one 'exchangeabilities' vector parameter",
        )
    try:
        current_parameterization = parameterize_dna_exchangeability_simplex(
            current_exchangeabilities
        )
    except (InvalidAlignmentError, PhylogeneticsError) as error:
        return build_metropolis_hastings_proposal(
            changed_fields=("vector_parameters.exchangeabilities",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    current_unconstrained_values = list(current_parameterization.unconstrained_values)
    coordinate_component_names = _simplex_coordinate_component_names(
        current_parameterization.component_names,
        reference_component_name=current_parameterization.reference_component_name,
    )
    selected_coordinate_index = rng.randrange(len(current_unconstrained_values))
    current_coordinate_value = current_unconstrained_values[selected_coordinate_index]
    proposed_coordinate_value = current_coordinate_value + rng.gauss(
        0.0,
        validated_coordinate_standard_deviation,
    )
    changed_field = (
        "vector_parameters.exchangeabilities."
        f"{coordinate_component_names[selected_coordinate_index]}"
    )
    if not math.isfinite(proposed_coordinate_value):
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="GTR exchangeability proposal produced one non-finite simplex coordinate",
        )
    current_unconstrained_values[selected_coordinate_index] = proposed_coordinate_value
    proposed_parameterization = resolve_dna_exchangeability_simplex_from_unconstrained(
        current_unconstrained_values
    )
    if (
        proposed_parameterization.constrained_values
        == current_parameterization.constrained_values
    ):
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="GTR exchangeability proposal must change normalized exchangeabilities",
        )
    proposed_vector_parameters = {
        parameter_name: dict(component_values)
        for parameter_name, component_values in current_state.model_parameters.vector_parameters.items()
    }
    proposed_vector_parameters["exchangeabilities"] = (
        proposed_parameterization.constrained_mapping()
    )
    proposal_density = -math.log(
        len(current_unconstrained_values)
    ) + _gaussian_random_walk_density(
        coordinate_change=proposed_coordinate_value - current_coordinate_value,
        standard_deviation=validated_coordinate_standard_deviation,
    )
    return build_metropolis_hastings_proposal(
        changed_fields=(changed_field,),
        log_forward_density=proposal_density,
        log_reverse_density=proposal_density,
        is_valid=True,
        proposed_tree=current_state.tree.to_tree(),
        proposed_model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters=current_state.model_parameters.categorical_parameters,
            scalar_parameters=current_state.model_parameters.scalar_parameters,
            vector_parameters=proposed_vector_parameters,
        ),
    )


def propose_base_frequency_simplex_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    *,
    unconstrained_coordinate_standard_deviation: float,
) -> MetropolisHastingsProposal:
    """Propose one simplex-preserving DNA base-frequency parameter move."""
    validated_coordinate_standard_deviation = _validate_positive_finite_float(
        value=unconstrained_coordinate_standard_deviation,
        field_name="unconstrained_coordinate_standard_deviation",
        owner_name="base-frequency simplex proposal",
    )
    current_base_frequencies = current_state.model_parameters.vector_parameters.get(
        "base-frequencies"
    )
    if current_base_frequencies is None:
        return build_metropolis_hastings_proposal(
            changed_fields=("vector_parameters.base-frequencies",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "base-frequency simplex proposal requires one 'base-frequencies' vector parameter"
            ),
        )
    try:
        current_parameterization = parameterize_dna_base_frequency_simplex(
            current_base_frequencies
        )
    except (InvalidAlignmentError, PhylogeneticsError) as error:
        return build_metropolis_hastings_proposal(
            changed_fields=("vector_parameters.base-frequencies",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    current_unconstrained_values = list(current_parameterization.unconstrained_values)
    coordinate_component_names = _simplex_coordinate_component_names(
        current_parameterization.component_names,
        reference_component_name=current_parameterization.reference_component_name,
    )
    selected_coordinate_index = rng.randrange(len(current_unconstrained_values))
    current_coordinate_value = current_unconstrained_values[selected_coordinate_index]
    proposed_coordinate_value = current_coordinate_value + rng.gauss(
        0.0,
        validated_coordinate_standard_deviation,
    )
    changed_field = (
        "vector_parameters.base-frequencies."
        f"{coordinate_component_names[selected_coordinate_index]}"
    )
    if not math.isfinite(proposed_coordinate_value):
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "base-frequency simplex proposal produced one non-finite simplex coordinate"
            ),
        )
    current_unconstrained_values[selected_coordinate_index] = proposed_coordinate_value
    try:
        proposed_parameterization = (
            resolve_dna_base_frequency_simplex_from_unconstrained(
                current_unconstrained_values
            )
        )
    except (InvalidAlignmentError, PhylogeneticsError) as error:
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    if (
        proposed_parameterization.constrained_values
        == current_parameterization.constrained_values
    ):
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "base-frequency simplex proposal must change normalized base frequencies"
            ),
        )
    proposed_vector_parameters = {
        parameter_name: dict(component_values)
        for parameter_name, component_values in current_state.model_parameters.vector_parameters.items()
    }
    proposed_vector_parameters["base-frequencies"] = (
        proposed_parameterization.constrained_mapping()
    )
    proposal_density = -math.log(
        len(current_unconstrained_values)
    ) + _gaussian_random_walk_density(
        coordinate_change=proposed_coordinate_value - current_coordinate_value,
        standard_deviation=validated_coordinate_standard_deviation,
    )
    return build_metropolis_hastings_proposal(
        changed_fields=(changed_field,),
        log_forward_density=proposal_density,
        log_reverse_density=proposal_density,
        is_valid=True,
        proposed_tree=current_state.tree.to_tree(),
        proposed_model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters=current_state.model_parameters.categorical_parameters,
            scalar_parameters=current_state.model_parameters.scalar_parameters,
            vector_parameters=proposed_vector_parameters,
        ),
    )


def propose_gamma_alpha_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    *,
    log_scale_standard_deviation: float,
) -> MetropolisHastingsProposal:
    """Propose one multiplicative change to one positive discrete-gamma alpha."""
    validated_log_scale_standard_deviation = _validate_positive_finite_float(
        value=log_scale_standard_deviation,
        field_name="log_scale_standard_deviation",
        owner_name="gamma-alpha proposal",
    )
    current_gamma_alpha = current_state.model_parameters.scalar_parameters.get(
        "gamma-alpha"
    )
    if current_gamma_alpha is None:
        return build_metropolis_hastings_proposal(
            changed_fields=("scalar_parameters.gamma-alpha",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "gamma-alpha proposal requires one 'gamma-alpha' scalar parameter"
            ),
        )
    try:
        validated_current_gamma_alpha = validate_discrete_gamma_alpha(
            current_gamma_alpha
        )
    except ValueError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=("scalar_parameters.gamma-alpha",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    try:
        scale_factor = math.exp(rng.gauss(0.0, validated_log_scale_standard_deviation))
    except OverflowError:
        return build_metropolis_hastings_proposal(
            changed_fields=("scalar_parameters.gamma-alpha",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="gamma-alpha scaling factor overflowed",
        )
    proposed_gamma_alpha = validated_current_gamma_alpha * scale_factor
    try:
        validated_proposed_gamma_alpha = validate_discrete_gamma_alpha(
            proposed_gamma_alpha
        )
    except ValueError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=("scalar_parameters.gamma-alpha",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    proposed_scalar_parameters = dict(current_state.model_parameters.scalar_parameters)
    proposed_scalar_parameters["gamma-alpha"] = validated_proposed_gamma_alpha
    log_forward_density = _lognormal_scaling_density(
        current_branch_length=validated_current_gamma_alpha,
        proposed_branch_length=validated_proposed_gamma_alpha,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    log_reverse_density = _lognormal_scaling_density(
        current_branch_length=validated_proposed_gamma_alpha,
        proposed_branch_length=validated_current_gamma_alpha,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    current_tree = current_state.tree.to_tree()
    current_tree.rooted = current_state.tree.rooted
    return build_metropolis_hastings_proposal(
        changed_fields=("scalar_parameters.gamma-alpha",),
        log_forward_density=log_forward_density,
        log_reverse_density=log_reverse_density,
        is_valid=True,
        proposed_tree=current_tree,
        proposed_model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters=current_state.model_parameters.categorical_parameters,
            scalar_parameters=proposed_scalar_parameters,
            vector_parameters=current_state.model_parameters.vector_parameters,
        ),
    )


def propose_kappa_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    *,
    log_scale_standard_deviation: float,
) -> MetropolisHastingsProposal:
    """Propose one multiplicative change to one positive nucleotide kappa parameter."""
    validated_log_scale_standard_deviation = _validate_positive_finite_float(
        value=log_scale_standard_deviation,
        field_name="log_scale_standard_deviation",
        owner_name="kappa proposal",
    )
    current_kappa = current_state.model_parameters.scalar_parameters.get("kappa")
    if current_kappa is None:
        return build_metropolis_hastings_proposal(
            changed_fields=("scalar_parameters.kappa",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="kappa proposal requires one 'kappa' scalar parameter",
        )
    try:
        validated_current_kappa = validate_positive_kappa(
            current_kappa,
            model_name="kappa proposal",
        )
    except ValueError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=("scalar_parameters.kappa",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    try:
        scale_factor = math.exp(rng.gauss(0.0, validated_log_scale_standard_deviation))
    except OverflowError:
        return build_metropolis_hastings_proposal(
            changed_fields=("scalar_parameters.kappa",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="kappa scaling factor overflowed",
        )
    proposed_kappa = validated_current_kappa * scale_factor
    try:
        validated_proposed_kappa = validate_positive_kappa(
            proposed_kappa,
            model_name="kappa proposal",
        )
    except ValueError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=("scalar_parameters.kappa",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    proposed_scalar_parameters = dict(current_state.model_parameters.scalar_parameters)
    proposed_scalar_parameters["kappa"] = validated_proposed_kappa
    log_forward_density = _lognormal_scaling_density(
        current_branch_length=validated_current_kappa,
        proposed_branch_length=validated_proposed_kappa,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    log_reverse_density = _lognormal_scaling_density(
        current_branch_length=validated_proposed_kappa,
        proposed_branch_length=validated_current_kappa,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    current_tree = current_state.tree.to_tree()
    current_tree.rooted = current_state.tree.rooted
    return build_metropolis_hastings_proposal(
        changed_fields=("scalar_parameters.kappa",),
        log_forward_density=log_forward_density,
        log_reverse_density=log_reverse_density,
        is_valid=True,
        proposed_tree=current_tree,
        proposed_model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters=current_state.model_parameters.categorical_parameters,
            scalar_parameters=proposed_scalar_parameters,
            vector_parameters=current_state.model_parameters.vector_parameters,
        ),
    )


def propose_clock_rate_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    *,
    log_scale_standard_deviation: float,
    parameter_name: str = "clock-rate",
) -> MetropolisHastingsProposal:
    """Propose one multiplicative change to one positive scalar clock-rate parameter."""
    validated_log_scale_standard_deviation = _validate_positive_finite_float(
        value=log_scale_standard_deviation,
        field_name="log_scale_standard_deviation",
        owner_name="clock-rate proposal",
    )
    validated_parameter_name = _validate_parameter_name(
        value=parameter_name,
        field_name="parameter_name",
        owner_name="clock-rate proposal",
    )
    changed_field = f"scalar_parameters.{validated_parameter_name}"
    current_clock_rate = current_state.model_parameters.scalar_parameters.get(
        validated_parameter_name
    )
    if current_clock_rate is None:
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "clock-rate proposal requires one "
                f"'{validated_parameter_name}' scalar parameter"
            ),
        )
    try:
        validated_current_clock_rate = _validate_positive_finite_float(
            value=current_clock_rate,
            field_name=validated_parameter_name,
            owner_name="clock-rate proposal",
        )
    except PhylogeneticsError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    try:
        scale_factor = math.exp(rng.gauss(0.0, validated_log_scale_standard_deviation))
    except OverflowError:
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="clock-rate scaling factor overflowed",
        )
    proposed_clock_rate = validated_current_clock_rate * scale_factor
    try:
        validated_proposed_clock_rate = _validate_positive_finite_float(
            value=proposed_clock_rate,
            field_name=validated_parameter_name,
            owner_name="clock-rate proposal",
        )
    except PhylogeneticsError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    proposed_scalar_parameters = dict(current_state.model_parameters.scalar_parameters)
    proposed_scalar_parameters[validated_parameter_name] = validated_proposed_clock_rate
    log_forward_density = _lognormal_scaling_density(
        current_branch_length=validated_current_clock_rate,
        proposed_branch_length=validated_proposed_clock_rate,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    log_reverse_density = _lognormal_scaling_density(
        current_branch_length=validated_proposed_clock_rate,
        proposed_branch_length=validated_current_clock_rate,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    current_tree = current_state.tree.to_tree()
    current_tree.rooted = current_state.tree.rooted
    return build_metropolis_hastings_proposal(
        changed_fields=(changed_field,),
        log_forward_density=log_forward_density,
        log_reverse_density=log_reverse_density,
        is_valid=True,
        proposed_tree=current_tree,
        proposed_model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters=current_state.model_parameters.categorical_parameters,
            scalar_parameters=proposed_scalar_parameters,
            vector_parameters=current_state.model_parameters.vector_parameters,
        ),
    )


def propose_continuous_trait_location_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    *,
    standard_deviation: float,
    parameter_name: str = "root-state",
) -> MetropolisHastingsProposal:
    """Propose one additive Gaussian move on one real-valued continuous-trait location parameter."""
    validated_standard_deviation = _validate_positive_finite_float(
        value=standard_deviation,
        field_name="standard_deviation",
        owner_name="continuous-trait location proposal",
    )
    validated_parameter_name = _validate_parameter_name(
        value=parameter_name,
        field_name="parameter_name",
        owner_name="continuous-trait location proposal",
    )
    changed_field = f"scalar_parameters.{validated_parameter_name}"
    current_location = current_state.model_parameters.scalar_parameters.get(
        validated_parameter_name
    )
    if current_location is None:
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "continuous-trait location proposal requires one "
                f"'{validated_parameter_name}' scalar parameter"
            ),
        )
    try:
        validated_current_location = _validate_finite_float(
            value=current_location,
            field_name=validated_parameter_name,
            owner_name="continuous-trait location proposal",
        )
    except PhylogeneticsError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    proposed_location = validated_current_location + rng.gauss(
        0.0,
        validated_standard_deviation,
    )
    try:
        validated_proposed_location = _validate_finite_float(
            value=proposed_location,
            field_name=validated_parameter_name,
            owner_name="continuous-trait location proposal",
        )
    except PhylogeneticsError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    proposed_scalar_parameters = dict(current_state.model_parameters.scalar_parameters)
    proposed_scalar_parameters[validated_parameter_name] = validated_proposed_location
    current_tree = current_state.tree.to_tree()
    current_tree.rooted = current_state.tree.rooted
    return build_metropolis_hastings_proposal(
        changed_fields=(changed_field,),
        log_forward_density=0.0,
        log_reverse_density=0.0,
        is_valid=True,
        proposed_tree=current_tree,
        proposed_model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters=current_state.model_parameters.categorical_parameters,
            scalar_parameters=proposed_scalar_parameters,
            vector_parameters=current_state.model_parameters.vector_parameters,
        ),
    )


def propose_discrete_trait_rate_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    *,
    log_scale_standard_deviation: float,
    parameter_name: str = "discrete-trait-rates",
) -> MetropolisHastingsProposal:
    """Propose one multiplicative change to one positive Mk transition-rate parameter."""
    validated_log_scale_standard_deviation = _validate_positive_finite_float(
        value=log_scale_standard_deviation,
        field_name="log_scale_standard_deviation",
        owner_name="discrete-trait rate proposal",
    )
    validated_parameter_name = _validate_parameter_name(
        value=parameter_name,
        field_name="parameter_name",
        owner_name="discrete-trait rate proposal",
    )
    current_rate_parameters = current_state.model_parameters.vector_parameters.get(
        validated_parameter_name
    )
    if current_rate_parameters is None:
        return build_metropolis_hastings_proposal(
            changed_fields=(f"vector_parameters.{validated_parameter_name}",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "discrete-trait rate proposal requires one "
                f"'{validated_parameter_name}' vector parameter"
            ),
        )
    if not current_rate_parameters:
        return build_metropolis_hastings_proposal(
            changed_fields=(f"vector_parameters.{validated_parameter_name}",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "discrete-trait rate proposal requires at least one named transition-rate parameter"
            ),
        )
    validated_current_rate_parameters = {
        component_name: _validate_positive_finite_float(
            value=component_value,
            field_name=f"{validated_parameter_name}.{component_name}",
            owner_name="discrete-trait rate proposal",
        )
        for component_name, component_value in current_rate_parameters.items()
    }
    component_names = sorted(validated_current_rate_parameters)
    selected_component_name = component_names[rng.randrange(len(component_names))]
    current_rate_value = validated_current_rate_parameters[selected_component_name]
    try:
        scale_factor = math.exp(rng.gauss(0.0, validated_log_scale_standard_deviation))
    except OverflowError:
        return build_metropolis_hastings_proposal(
            changed_fields=(
                f"vector_parameters.{validated_parameter_name}.{selected_component_name}",
            ),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="discrete-trait rate scaling factor overflowed",
        )
    proposed_rate_value = current_rate_value * scale_factor
    try:
        validated_proposed_rate_value = _validate_positive_finite_float(
            value=proposed_rate_value,
            field_name=f"{validated_parameter_name}.{selected_component_name}",
            owner_name="discrete-trait rate proposal",
        )
    except PhylogeneticsError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=(
                f"vector_parameters.{validated_parameter_name}.{selected_component_name}",
            ),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    proposed_vector_parameters = {
        parameter_label: dict(component_values)
        for parameter_label, component_values in current_state.model_parameters.vector_parameters.items()
    }
    proposed_vector_parameters[validated_parameter_name][selected_component_name] = (
        validated_proposed_rate_value
    )
    log_forward_density = -math.log(len(component_names)) + _lognormal_scaling_density(
        current_branch_length=current_rate_value,
        proposed_branch_length=validated_proposed_rate_value,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    log_reverse_density = -math.log(len(component_names)) + _lognormal_scaling_density(
        current_branch_length=validated_proposed_rate_value,
        proposed_branch_length=current_rate_value,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    current_tree = current_state.tree.to_tree()
    current_tree.rooted = current_state.tree.rooted
    return build_metropolis_hastings_proposal(
        changed_fields=(
            f"vector_parameters.{validated_parameter_name}.{selected_component_name}",
        ),
        log_forward_density=log_forward_density,
        log_reverse_density=log_reverse_density,
        is_valid=True,
        proposed_tree=current_tree,
        proposed_model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters=current_state.model_parameters.categorical_parameters,
            scalar_parameters=current_state.model_parameters.scalar_parameters,
            vector_parameters=proposed_vector_parameters,
        ),
    )


def propose_invariant_proportion_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    *,
    logit_standard_deviation: float,
) -> MetropolisHastingsProposal:
    """Propose one bounded change to one invariant-site mixture proportion."""
    validated_logit_standard_deviation = _validate_positive_finite_float(
        value=logit_standard_deviation,
        field_name="logit_standard_deviation",
        owner_name="invariant-proportion proposal",
    )
    current_invariant_proportion = current_state.model_parameters.scalar_parameters.get(
        "invariant-proportion"
    )
    if current_invariant_proportion is None:
        return build_metropolis_hastings_proposal(
            changed_fields=("scalar_parameters.invariant-proportion",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "invariant-proportion proposal requires one 'invariant-proportion' scalar parameter"
            ),
        )
    try:
        validated_current_invariant_proportion = validate_invariant_proportion(
            current_invariant_proportion,
            model_name="invariant-proportion proposal",
        )
    except ValueError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=("scalar_parameters.invariant-proportion",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    if not 0.0 < validated_current_invariant_proportion < 1.0:
        return build_metropolis_hastings_proposal(
            changed_fields=("scalar_parameters.invariant-proportion",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "invariant-proportion proposal requires one interior invariant proportion in (0, 1)"
            ),
        )
    current_logit = _logit_probability(validated_current_invariant_proportion)
    proposed_logit = current_logit + rng.gauss(0.0, validated_logit_standard_deviation)
    if not math.isfinite(proposed_logit):
        return build_metropolis_hastings_proposal(
            changed_fields=("scalar_parameters.invariant-proportion",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "invariant-proportion proposal produced one non-finite logit coordinate"
            ),
        )
    proposed_invariant_proportion = _inverse_logit_probability(proposed_logit)
    try:
        validated_proposed_invariant_proportion = validate_invariant_proportion(
            proposed_invariant_proportion,
            model_name="invariant-proportion proposal",
        )
    except ValueError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=("scalar_parameters.invariant-proportion",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    proposed_scalar_parameters = dict(current_state.model_parameters.scalar_parameters)
    proposed_scalar_parameters["invariant-proportion"] = (
        validated_proposed_invariant_proportion
    )
    coordinate_change = proposed_logit - current_logit
    gaussian_density = _gaussian_random_walk_density(
        coordinate_change=coordinate_change,
        standard_deviation=validated_logit_standard_deviation,
    )
    log_forward_density = gaussian_density + _log_probability_logit_jacobian(
        validated_proposed_invariant_proportion
    )
    log_reverse_density = gaussian_density + _log_probability_logit_jacobian(
        validated_current_invariant_proportion
    )
    current_tree = current_state.tree.to_tree()
    current_tree.rooted = current_state.tree.rooted
    return build_metropolis_hastings_proposal(
        changed_fields=("scalar_parameters.invariant-proportion",),
        log_forward_density=log_forward_density,
        log_reverse_density=log_reverse_density,
        is_valid=True,
        proposed_tree=current_tree,
        proposed_model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters=current_state.model_parameters.categorical_parameters,
            scalar_parameters=proposed_scalar_parameters,
            vector_parameters=current_state.model_parameters.vector_parameters,
        ),
    )


def propose_node_height_sliding_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    *,
    height_slide_standard_deviation: float,
    ultrametric_tolerance: float = APE_ULTRAMETRIC_TOLERANCE,
) -> MetropolisHastingsProposal:
    """Propose one additive slide on one internal node height."""
    validated_height_slide_standard_deviation = _validate_positive_finite_float(
        value=height_slide_standard_deviation,
        field_name="height_slide_standard_deviation",
        owner_name="node-height sliding proposal",
    )
    validated_ultrametric_tolerance = _validate_positive_finite_float(
        value=ultrametric_tolerance,
        field_name="ultrametric_tolerance",
        owner_name="node-height sliding proposal",
    )
    current_tree = current_state.tree.to_tree()
    current_tree.rooted = current_state.tree.rooted
    try:
        current_node_heights = _compute_rooted_ultrametric_node_heights(
            current_tree=current_tree,
            ultrametric_tolerance=validated_ultrametric_tolerance,
        )
    except PhylogeneticsError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.node_heights",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    candidate_nodes = _collect_node_height_slide_candidates(
        current_tree=current_tree,
        node_heights=current_node_heights,
    )
    if not candidate_nodes:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.node_heights",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "node-height sliding requires one non-root internal node with a positive age interval"
            ),
        )
    selected_candidate = candidate_nodes[rng.randrange(len(candidate_nodes))]
    proposed_height = selected_candidate.current_height + rng.gauss(
        0.0,
        validated_height_slide_standard_deviation,
    )
    changed_field = f"tree.node_height:{selected_candidate.node_id}"
    if not math.isfinite(proposed_height):
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="node-height sliding produced a non-finite proposed height",
        )
    if not (
        selected_candidate.lower_height_bound
        < proposed_height
        < selected_candidate.upper_height_bound
    ):
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "node-height sliding would violate ancestor-descendant age order"
            ),
        )
    proposed_tree = _copy_tree_with_slid_node_height(
        current_tree=current_tree,
        node_id=selected_candidate.node_id,
        proposed_height=proposed_height,
        node_heights=current_node_heights,
    )
    proposed_node_heights = _compute_rooted_ultrametric_node_heights(
        current_tree=proposed_tree,
        ultrametric_tolerance=validated_ultrametric_tolerance,
    )
    proposed_candidate_nodes = _collect_node_height_slide_candidates(
        current_tree=proposed_tree,
        node_heights=proposed_node_heights,
    )
    log_forward_density = -math.log(
        len(candidate_nodes)
    ) + _normal_node_height_slide_density(
        current_height=selected_candidate.current_height,
        proposed_height=proposed_height,
        height_slide_standard_deviation=validated_height_slide_standard_deviation,
    )
    log_reverse_density = -math.log(
        len(proposed_candidate_nodes)
    ) + _normal_node_height_slide_density(
        current_height=proposed_height,
        proposed_height=selected_candidate.current_height,
        height_slide_standard_deviation=validated_height_slide_standard_deviation,
    )
    return build_metropolis_hastings_proposal(
        changed_fields=(changed_field,),
        log_forward_density=log_forward_density,
        log_reverse_density=log_reverse_density,
        is_valid=True,
        proposed_tree=proposed_tree,
        proposed_model_parameters=current_state.model_parameters,
    )


def propose_nni_topology_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
) -> MetropolisHastingsProposal:
    """Propose one rooted nearest-neighbor interchange topology move."""
    current_tree = current_state.tree.to_tree()
    current_tree.rooted = current_state.tree.rooted
    try:
        validate_rooted_nni_tree(current_tree)
    except ValueError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.topology",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    current_neighbors = _enumerate_rooted_nni_neighbors(current_tree)
    if not current_neighbors:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.topology",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="NNI topology proposal requires at least one legal rooted NNI move",
        )
    current_topology_fingerprint = rooted_topology_fingerprint(current_tree)
    selected_candidate, proposed_tree, proposed_topology_fingerprint = (
        current_neighbors[rng.randrange(len(current_neighbors))]
    )
    if proposed_topology_fingerprint == current_topology_fingerprint:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.topology",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "NNI topology proposal must change rooted topology rather than child order"
            ),
        )
    proposed_neighbors = _enumerate_rooted_nni_neighbors(proposed_tree)
    reverse_match_count = sum(
        1
        for _candidate, _neighbor_tree, neighbor_topology_fingerprint in proposed_neighbors
        if neighbor_topology_fingerprint == current_topology_fingerprint
    )
    if reverse_match_count <= 0:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.topology",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="NNI topology proposal could not find one reverse rooted interchange",
        )
    log_forward_density = -math.log(len(current_neighbors))
    log_reverse_density = math.log(reverse_match_count) - math.log(
        len(proposed_neighbors)
    )
    changed_field = (
        "tree.topology:nni:"
        f"{selected_candidate.parent_node_id}:{selected_candidate.child_node_id}:"
        f"{selected_candidate.exchanged_child_node_id}"
    )
    return build_metropolis_hastings_proposal(
        changed_fields=(changed_field,),
        log_forward_density=log_forward_density,
        log_reverse_density=log_reverse_density,
        is_valid=True,
        proposed_tree=proposed_tree,
        proposed_model_parameters=current_state.model_parameters,
    )


def propose_spr_topology_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
) -> MetropolisHastingsProposal:
    """Propose one rooted subtree-prune-regraft topology move."""
    current_tree = current_state.tree.to_tree()
    current_tree.rooted = current_state.tree.rooted
    try:
        validate_rooted_spr_tree(current_tree)
    except ValueError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.topology",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    current_report = enumerate_rooted_spr_neighbors(current_tree)
    if not current_report.neighbor_rows:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.topology",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="SPR topology proposal requires at least one legal rooted SPR move",
        )
    current_legal_move_count = _rooted_spr_legal_move_count(current_report)
    selected_neighbor_row = current_report.neighbor_rows[
        rng.randrange(len(current_report.neighbor_rows))
    ]
    proposed_tree = PhyloTree.from_newick(selected_neighbor_row.neighbor_tree_newick)
    proposed_tree.rooted = current_state.tree.rooted
    current_topology_fingerprint = rooted_topology_fingerprint(current_tree)
    if (
        selected_neighbor_row.neighbor_topology_fingerprint
        == current_topology_fingerprint
    ):
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.topology",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "SPR topology proposal must change rooted topology rather than child order"
            ),
        )
    proposed_report = enumerate_rooted_spr_neighbors(proposed_tree)
    reverse_neighbor_row = next(
        (
            row
            for row in proposed_report.neighbor_rows
            if row.neighbor_topology_fingerprint == current_topology_fingerprint
        ),
        None,
    )
    if reverse_neighbor_row is None:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.topology",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="SPR topology proposal could not find one reverse rooted prune-regraft move",
        )
    proposed_legal_move_count = _rooted_spr_legal_move_count(proposed_report)
    log_forward_density = math.log(
        selected_neighbor_row.supporting_move_count
    ) - math.log(current_legal_move_count)
    log_reverse_density = math.log(
        reverse_neighbor_row.supporting_move_count
    ) - math.log(proposed_legal_move_count)
    changed_field = (
        "tree.topology:spr:"
        f"{selected_neighbor_row.representative_pruned_node_id}:"
        f"{selected_neighbor_row.representative_regraft_target_branch_id}"
    )
    return build_metropolis_hastings_proposal(
        changed_fields=(changed_field,),
        log_forward_density=log_forward_density,
        log_reverse_density=log_reverse_density,
        is_valid=True,
        proposed_tree=proposed_tree,
        proposed_model_parameters=current_state.model_parameters,
    )


def propose_tbr_topology_move(
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
) -> MetropolisHastingsProposal:
    """Propose one rooted tree-bisection-reconnection topology move."""
    current_tree = current_state.tree.to_tree()
    current_tree.rooted = current_state.tree.rooted
    try:
        validate_rooted_tbr_tree(current_tree)
    except ValueError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.topology",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    reversible_neighbors = _enumerate_reversible_rooted_tbr_neighbors(current_tree)
    if not reversible_neighbors:
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.topology",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="TBR topology proposal requires at least one legal rooted TBR move",
        )
    current_legal_reconnection_count = _reversible_rooted_tbr_move_count(
        reversible_neighbors
    )
    selected_neighbor = reversible_neighbors[rng.randrange(len(reversible_neighbors))]
    selected_neighbor_row = selected_neighbor.neighbor_row
    proposed_tree = selected_neighbor.proposed_tree
    current_topology_fingerprint = rooted_topology_fingerprint(current_tree)
    if (
        selected_neighbor_row.neighbor_topology_fingerprint
        == current_topology_fingerprint
    ):
        return build_metropolis_hastings_proposal(
            changed_fields=("tree.topology",),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "TBR topology proposal must change rooted topology rather than child order"
            ),
        )
    reverse_neighbor_row = selected_neighbor.reverse_neighbor_row
    proposed_legal_reconnection_count = _reversible_rooted_tbr_move_count(
        _enumerate_reversible_rooted_tbr_neighbors(proposed_tree)
    )
    log_forward_density = math.log(
        selected_neighbor_row.supporting_reconnection_count
    ) - math.log(current_legal_reconnection_count)
    log_reverse_density = math.log(
        reverse_neighbor_row.supporting_reconnection_count
    ) - math.log(proposed_legal_reconnection_count)
    changed_field = (
        "tree.topology:tbr:"
        f"{selected_neighbor_row.representative_cut_edge_id}:"
        f"{selected_neighbor_row.representative_left_attachment_branch_id}:"
        f"{selected_neighbor_row.representative_right_attachment_branch_id}"
    )
    return build_metropolis_hastings_proposal(
        changed_fields=(changed_field,),
        log_forward_density=log_forward_density,
        log_reverse_density=log_reverse_density,
        is_valid=True,
        proposed_tree=proposed_tree,
        proposed_model_parameters=current_state.model_parameters,
    )


def score_bayesian_phylogenetic_state(
    *,
    tree,
    model_parameters: BayesianModelParameterState,
    update_prior_components: BayesianPriorUpdate,
    update_log_likelihood: BayesianLikelihoodUpdate,
) -> BayesianPhylogeneticState:
    """Score one Bayesian state from owned tree and model-parameter surfaces."""
    provisional_state = build_bayesian_phylogenetic_state(
        tree=tree,
        model_parameters=model_parameters,
        prior_components=[
            BayesianPriorComponentState(
                component_name="provisional",
                family=None,
                log_prior=0.0,
            )
        ],
        log_likelihood=0.0,
    )
    updated_prior_components = update_prior_components(provisional_state)
    if not updated_prior_components:
        raise PhylogeneticsError(
            "metropolis-hastings state scoring requires at least one prior component",
            code="metropolis_hastings_prior_components_empty",
        )
    return build_bayesian_phylogenetic_state(
        tree=provisional_state.tree.to_tree(),
        model_parameters=provisional_state.model_parameters,
        prior_components=updated_prior_components,
        log_likelihood=update_log_likelihood(provisional_state),
    )


def run_metropolis_hastings_sampler(
    *,
    initial_state: BayesianPhylogeneticState,
    propose_state: BayesianStateProposal,
    update_prior_components: BayesianPriorUpdate,
    update_log_likelihood: BayesianLikelihoodUpdate,
    iteration_count: int,
    sample_every: int = 1,
    seed: int = 0,
) -> MetropolisHastingsRunReport:
    """Run one native Metropolis-Hastings chain over Bayesian phylogenetic states."""
    return run_checkpointed_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=propose_state,
        update_prior_components=update_prior_components,
        update_log_likelihood=update_log_likelihood,
        iteration_count=iteration_count,
        sample_every=sample_every,
        seed=seed,
    ).to_chain_report()


def run_checkpointed_metropolis_hastings_sampler(
    *,
    initial_state: BayesianPhylogeneticState,
    propose_state: BayesianStateProposal,
    update_prior_components: BayesianPriorUpdate,
    update_log_likelihood: BayesianLikelihoodUpdate,
    iteration_count: int,
    sample_every: int = 1,
    seed: int = 0,
    checkpoint_iteration_indexes: Sequence[int] = (),
    stop_after_iteration_index: int | None = None,
) -> CheckpointedMetropolisHastingsRunReport:
    """Run one checkpoint-aware Metropolis-Hastings execution."""
    validated_iteration_count = _validate_positive_integer(
        value=iteration_count,
        field_name="iteration_count",
        owner_name="checkpointed metropolis-hastings sampler",
    )
    validated_sample_every = _validate_positive_integer(
        value=sample_every,
        field_name="sample_every",
        owner_name="checkpointed metropolis-hastings sampler",
    )
    validated_seed = _validate_integer_seed(seed)
    validated_checkpoint_iteration_indexes = _validate_checkpoint_iteration_indexes(
        checkpoint_iteration_indexes=checkpoint_iteration_indexes,
        iteration_count=validated_iteration_count,
    )
    validated_stop_after_iteration_index = _validate_stop_after_iteration_index(
        stop_after_iteration_index=stop_after_iteration_index,
        checkpoint_iteration_indexes=validated_checkpoint_iteration_indexes,
        completed_iteration_count=0,
        iteration_count=validated_iteration_count,
        owner_name="checkpointed metropolis-hastings sampler",
    )
    rng = random.Random(validated_seed)  # nosec B311
    return _run_metropolis_hastings_execution(
        initial_state=initial_state,
        current_state=initial_state,
        propose_state=propose_state,
        update_prior_components=update_prior_components,
        update_log_likelihood=update_log_likelihood,
        iteration_count=validated_iteration_count,
        sample_every=validated_sample_every,
        seed=validated_seed,
        resumed=False,
        completed_iteration_count=0,
        accepted_count=0,
        sampled_states=[initial_state],
        step_rows=[],
        rng=rng,
        checkpoint_iteration_indexes=validated_checkpoint_iteration_indexes,
        stop_after_iteration_index=validated_stop_after_iteration_index,
    )


def resume_metropolis_hastings_sampler(
    *,
    checkpoint: MetropolisHastingsCheckpoint,
    propose_state: BayesianStateProposal,
    update_prior_components: BayesianPriorUpdate,
    update_log_likelihood: BayesianLikelihoodUpdate,
    checkpoint_iteration_indexes: Sequence[int] = (),
    stop_after_iteration_index: int | None = None,
) -> CheckpointedMetropolisHastingsRunReport:
    """Resume one checkpointed Metropolis-Hastings execution without restarting the chain."""
    validated_checkpoint = build_metropolis_hastings_checkpoint(
        iteration_count=checkpoint.iteration_count,
        sample_every=checkpoint.sample_every,
        seed=checkpoint.seed,
        completed_iteration_count=checkpoint.completed_iteration_count,
        accepted_count=checkpoint.accepted_count,
        initial_state=checkpoint.initial_state,
        current_state=checkpoint.current_state,
        sampled_states=checkpoint.sampled_states,
        step_rows=checkpoint.step_rows,
        random_state=checkpoint.random_state,
    )
    validated_checkpoint_iteration_indexes = _validate_checkpoint_iteration_indexes(
        checkpoint_iteration_indexes=checkpoint_iteration_indexes,
        iteration_count=validated_checkpoint.iteration_count,
    )
    validated_stop_after_iteration_index = _validate_stop_after_iteration_index(
        stop_after_iteration_index=stop_after_iteration_index,
        checkpoint_iteration_indexes=validated_checkpoint_iteration_indexes,
        completed_iteration_count=validated_checkpoint.completed_iteration_count,
        iteration_count=validated_checkpoint.iteration_count,
        owner_name="metropolis-hastings sampler resume",
    )
    rng = random.Random()  # nosec B311
    rng.setstate(_restore_rng_state_tuple(validated_checkpoint.random_state))
    return _run_metropolis_hastings_execution(
        initial_state=validated_checkpoint.initial_state,
        current_state=validated_checkpoint.current_state,
        propose_state=propose_state,
        update_prior_components=update_prior_components,
        update_log_likelihood=update_log_likelihood,
        iteration_count=validated_checkpoint.iteration_count,
        sample_every=validated_checkpoint.sample_every,
        seed=validated_checkpoint.seed,
        resumed=True,
        completed_iteration_count=validated_checkpoint.completed_iteration_count,
        accepted_count=validated_checkpoint.accepted_count,
        sampled_states=list(validated_checkpoint.sampled_states),
        step_rows=list(validated_checkpoint.step_rows),
        rng=rng,
        checkpoint_iteration_indexes=validated_checkpoint_iteration_indexes,
        stop_after_iteration_index=validated_stop_after_iteration_index,
    )


def _run_metropolis_hastings_execution(
    *,
    initial_state: BayesianPhylogeneticState,
    current_state: BayesianPhylogeneticState,
    propose_state: BayesianStateProposal,
    update_prior_components: BayesianPriorUpdate,
    update_log_likelihood: BayesianLikelihoodUpdate,
    iteration_count: int,
    sample_every: int,
    seed: int,
    resumed: bool,
    completed_iteration_count: int,
    accepted_count: int,
    sampled_states: list[BayesianPhylogeneticState],
    step_rows: list[MetropolisHastingsStepRow],
    rng: random.Random,
    checkpoint_iteration_indexes: tuple[int, ...],
    stop_after_iteration_index: int | None,
) -> CheckpointedMetropolisHastingsRunReport:
    final_iteration_index = stop_after_iteration_index or iteration_count
    checkpoint_iteration_index_set = set(checkpoint_iteration_indexes)
    emitted_checkpoints: list[MetropolisHastingsCheckpoint] = []

    for iteration_index in range(
        completed_iteration_count + 1, final_iteration_index + 1
    ):
        previous_state = current_state
        proposal = _validate_metropolis_hastings_proposal(
            propose_state(current_state, rng)
        )
        proposed_state: BayesianPhylogeneticState | None = None
        log_acceptance_ratio: float | None = None
        accepted = False
        if proposal.is_valid:
            proposed_state = score_bayesian_phylogenetic_state(
                tree=proposal.proposed_tree,
                model_parameters=proposal.proposed_model_parameters,
                update_prior_components=update_prior_components,
                update_log_likelihood=update_log_likelihood,
            )
            log_acceptance_ratio = (
                proposed_state.posterior_log_score
                - current_state.posterior_log_score
                + proposal.log_hastings_ratio
            )
            accepted = _accept_metropolis_hastings_proposal(
                log_acceptance_ratio=log_acceptance_ratio,
                rng=rng,
            )
            if accepted:
                current_state = proposed_state
                accepted_count += 1
        if iteration_index % sample_every == 0:
            sampled_states.append(current_state)
        step_rows.append(
            MetropolisHastingsStepRow(
                iteration_index=iteration_index,
                proposal_changed_fields=proposal.changed_fields,
                proposal_valid=proposal.is_valid,
                proposal_invalid_reason=proposal.invalid_reason,
                log_forward_density=proposal.log_forward_density,
                log_reverse_density=proposal.log_reverse_density,
                accepted=accepted,
                log_hastings_ratio=proposal.log_hastings_ratio,
                current_posterior_log_score=previous_state.posterior_log_score,
                proposed_posterior_log_score=(
                    proposed_state.posterior_log_score
                    if proposed_state is not None
                    else None
                ),
                log_acceptance_ratio=log_acceptance_ratio,
                recorded_posterior_log_score=current_state.posterior_log_score,
            )
        )
        if iteration_index in checkpoint_iteration_index_set:
            emitted_checkpoints.append(
                build_metropolis_hastings_checkpoint(
                    iteration_count=iteration_count,
                    sample_every=sample_every,
                    seed=seed,
                    completed_iteration_count=iteration_index,
                    accepted_count=accepted_count,
                    initial_state=initial_state,
                    current_state=current_state,
                    sampled_states=sampled_states,
                    step_rows=step_rows,
                    random_state=_capture_rng_state(rng),
                )
            )

    realized_completed_iteration_count = final_iteration_index
    rejected_count = realized_completed_iteration_count - accepted_count
    return CheckpointedMetropolisHastingsRunReport(
        iteration_count=iteration_count,
        completed_iteration_count=realized_completed_iteration_count,
        sample_every=sample_every,
        seed=seed,
        resumed=resumed,
        completed=realized_completed_iteration_count == iteration_count,
        accepted_count=accepted_count,
        rejected_count=rejected_count,
        acceptance_rate=(
            accepted_count / realized_completed_iteration_count
            if realized_completed_iteration_count > 0
            else 0.0
        ),
        initial_state=initial_state,
        current_state=current_state,
        sampled_states=list(sampled_states),
        step_rows=list(step_rows),
        checkpoints=emitted_checkpoints,
    )


def _accept_metropolis_hastings_proposal(
    *,
    log_acceptance_ratio: float,
    rng: random.Random,
) -> bool:
    if log_acceptance_ratio >= 0.0:
        return True
    return math.log(rng.random()) <= log_acceptance_ratio


def _capture_rng_state(rng: random.Random) -> MetropolisHastingsRandomState:
    version, internal_state, gaussian_spare = rng.getstate()
    if not isinstance(internal_state, tuple):
        raise PhylogeneticsError(
            "metropolis-hastings sampler expected one tuple-backed random internal state",
            code="metropolis_hastings_random_state_internal_type_invalid",
        )
    return build_metropolis_hastings_random_state(
        version=version,
        internal_state=internal_state,
        gaussian_spare=gaussian_spare,
    )


def _restore_rng_state_tuple(
    random_state: MetropolisHastingsRandomState,
) -> tuple[int, tuple[int, ...], float | None]:
    return (
        random_state.version,
        tuple(random_state.internal_state),
        random_state.gaussian_spare,
    )


def _deserialize_metropolis_hastings_step_row(
    payload: Mapping[str, object],
) -> MetropolisHastingsStepRow:
    iteration_index = _require_integer_payload(
        payload,
        key="iteration_index",
        owner_name="metropolis-hastings step row deserialization",
    )
    proposal_changed_fields = _require_string_list_payload(
        payload,
        key="proposal_changed_fields",
        owner_name="metropolis-hastings step row deserialization",
    )
    proposal_valid = _require_boolean_payload(
        payload,
        key="proposal_valid",
        owner_name="metropolis-hastings step row deserialization",
    )
    proposal_invalid_reason = _optional_string_payload(
        payload,
        key="proposal_invalid_reason",
        owner_name="metropolis-hastings step row deserialization",
    )
    if proposal_valid and proposal_invalid_reason is not None:
        raise PhylogeneticsError(
            "metropolis-hastings step row deserialization requires valid rows to omit proposal_invalid_reason",
            code="metropolis_hastings_step_row_invalid_reason_unexpected",
        )
    if proposal_valid is False and proposal_invalid_reason is None:
        raise PhylogeneticsError(
            "metropolis-hastings step row deserialization requires invalid rows to include proposal_invalid_reason",
            code="metropolis_hastings_step_row_invalid_reason_missing",
        )
    proposed_posterior_log_score = _optional_float_payload(
        payload,
        key="proposed_posterior_log_score",
        owner_name="metropolis-hastings step row deserialization",
    )
    log_acceptance_ratio = _optional_float_payload(
        payload,
        key="log_acceptance_ratio",
        owner_name="metropolis-hastings step row deserialization",
    )
    return MetropolisHastingsStepRow(
        iteration_index=iteration_index,
        proposal_changed_fields=tuple(proposal_changed_fields),
        proposal_valid=proposal_valid,
        proposal_invalid_reason=proposal_invalid_reason,
        log_forward_density=_require_float_payload(
            payload,
            key="log_forward_density",
            owner_name="metropolis-hastings step row deserialization",
        ),
        log_reverse_density=_require_float_payload(
            payload,
            key="log_reverse_density",
            owner_name="metropolis-hastings step row deserialization",
        ),
        accepted=_require_boolean_payload(
            payload,
            key="accepted",
            owner_name="metropolis-hastings step row deserialization",
        ),
        log_hastings_ratio=_require_float_payload(
            payload,
            key="log_hastings_ratio",
            owner_name="metropolis-hastings step row deserialization",
        ),
        current_posterior_log_score=_require_float_payload(
            payload,
            key="current_posterior_log_score",
            owner_name="metropolis-hastings step row deserialization",
        ),
        proposed_posterior_log_score=proposed_posterior_log_score,
        log_acceptance_ratio=log_acceptance_ratio,
        recorded_posterior_log_score=_require_float_payload(
            payload,
            key="recorded_posterior_log_score",
            owner_name="metropolis-hastings step row deserialization",
        ),
    )


def _validate_checkpoint_iteration_indexes(
    *,
    checkpoint_iteration_indexes: Sequence[int],
    iteration_count: int,
) -> tuple[int, ...]:
    validated_checkpoint_iteration_indexes = tuple(checkpoint_iteration_indexes)
    if len(validated_checkpoint_iteration_indexes) != len(
        set(validated_checkpoint_iteration_indexes)
    ):
        raise PhylogeneticsError(
            "checkpointed metropolis-hastings sampler requires unique checkpoint iteration indexes",
            code="metropolis_hastings_checkpoint_iteration_duplicate",
        )
    for checkpoint_iteration_index in validated_checkpoint_iteration_indexes:
        _validate_positive_integer(
            value=checkpoint_iteration_index,
            field_name="checkpoint_iteration_index",
            owner_name="checkpointed metropolis-hastings sampler",
        )
        if checkpoint_iteration_index > iteration_count:
            raise PhylogeneticsError(
                "checkpointed metropolis-hastings sampler requires checkpoint iteration indexes to lie within the chain length",
                code="metropolis_hastings_checkpoint_iteration_out_of_range",
                details={
                    "checkpoint_iteration_index": checkpoint_iteration_index,
                    "iteration_count": iteration_count,
                },
            )
    return tuple(sorted(validated_checkpoint_iteration_indexes))


def _validate_stop_after_iteration_index(
    *,
    stop_after_iteration_index: int | None,
    checkpoint_iteration_indexes: tuple[int, ...],
    completed_iteration_count: int,
    iteration_count: int,
    owner_name: str,
) -> int | None:
    if stop_after_iteration_index is None:
        return None
    validated_stop_after_iteration_index = _validate_positive_integer(
        value=stop_after_iteration_index,
        field_name="stop_after_iteration_index",
        owner_name=owner_name,
    )
    if validated_stop_after_iteration_index <= completed_iteration_count:
        raise PhylogeneticsError(
            f"{owner_name} requires 'stop_after_iteration_index' to advance beyond the completed trace position",
            code="metropolis_hastings_stop_after_iteration_not_ahead",
            details={
                "stop_after_iteration_index": validated_stop_after_iteration_index,
                "completed_iteration_count": completed_iteration_count,
            },
        )
    if validated_stop_after_iteration_index > iteration_count:
        raise PhylogeneticsError(
            f"{owner_name} requires 'stop_after_iteration_index' to be less than or equal to 'iteration_count'",
            code="metropolis_hastings_stop_after_iteration_out_of_range",
            details={
                "stop_after_iteration_index": validated_stop_after_iteration_index,
                "iteration_count": iteration_count,
            },
        )
    if (
        validated_stop_after_iteration_index != iteration_count
        and validated_stop_after_iteration_index not in checkpoint_iteration_indexes
    ):
        raise PhylogeneticsError(
            f"{owner_name} requires partial executions to stop only on one declared checkpoint iteration index",
            code="metropolis_hastings_stop_after_iteration_not_checkpointed",
            details={
                "stop_after_iteration_index": validated_stop_after_iteration_index,
                "checkpoint_iteration_indexes": list(checkpoint_iteration_indexes),
            },
        )
    return validated_stop_after_iteration_index


def _validate_sampled_states(
    *,
    sampled_states: Sequence[BayesianPhylogeneticState],
    initial_state: BayesianPhylogeneticState,
    current_state: BayesianPhylogeneticState,
    completed_iteration_count: int,
    sample_every: int,
    owner_name: str,
) -> list[BayesianPhylogeneticState]:
    validated_sampled_states = list(sampled_states)
    if not validated_sampled_states:
        raise PhylogeneticsError(
            f"{owner_name} requires at least one sampled state",
            code="metropolis_hastings_sampled_states_empty",
        )
    if any(
        not isinstance(sampled_state, BayesianPhylogeneticState)
        for sampled_state in validated_sampled_states
    ):
        raise PhylogeneticsError(
            f"{owner_name} requires every sampled state to be one BayesianPhylogeneticState",
            code="metropolis_hastings_sampled_state_type_invalid",
        )
    expected_sampled_state_count = 1 + (completed_iteration_count // sample_every)
    if len(validated_sampled_states) != expected_sampled_state_count:
        raise PhylogeneticsError(
            f"{owner_name} requires sampled state count to match the retained trace positions",
            code="metropolis_hastings_sampled_state_count_invalid",
            details={
                "expected_sampled_state_count": expected_sampled_state_count,
                "actual_sampled_state_count": len(validated_sampled_states),
            },
        )
    if validated_sampled_states[0] != initial_state:
        raise PhylogeneticsError(
            f"{owner_name} requires the first sampled state to equal the initial state",
            code="metropolis_hastings_sampled_state_initial_mismatch",
        )
    if (
        completed_iteration_count % sample_every == 0
        and validated_sampled_states[-1] != current_state
    ):
        raise PhylogeneticsError(
            f"{owner_name} requires checkpoints at sampled iteration boundaries to retain the current state in sampled_states",
            code="metropolis_hastings_sampled_state_current_mismatch",
        )
    return validated_sampled_states


def _validate_step_rows(
    *,
    step_rows: Sequence[MetropolisHastingsStepRow],
    completed_iteration_count: int,
    current_state: BayesianPhylogeneticState,
    owner_name: str,
) -> list[MetropolisHastingsStepRow]:
    validated_step_rows = list(step_rows)
    if len(validated_step_rows) != completed_iteration_count:
        raise PhylogeneticsError(
            f"{owner_name} requires step row count to equal completed_iteration_count",
            code="metropolis_hastings_step_row_count_invalid",
            details={
                "completed_iteration_count": completed_iteration_count,
                "step_row_count": len(validated_step_rows),
            },
        )
    for expected_iteration_index, step_row in enumerate(validated_step_rows, start=1):
        if not isinstance(step_row, MetropolisHastingsStepRow):
            raise PhylogeneticsError(
                f"{owner_name} requires every step row to be one MetropolisHastingsStepRow",
                code="metropolis_hastings_step_row_type_invalid",
            )
        if step_row.iteration_index != expected_iteration_index:
            raise PhylogeneticsError(
                f"{owner_name} requires contiguous step row iteration indexes",
                code="metropolis_hastings_step_row_iteration_gap",
                details={
                    "expected_iteration_index": expected_iteration_index,
                    "actual_iteration_index": step_row.iteration_index,
                },
            )
    if validated_step_rows and (
        validated_step_rows[-1].recorded_posterior_log_score
        != current_state.posterior_log_score
    ):
        raise PhylogeneticsError(
            f"{owner_name} requires the last recorded posterior score to match the current state",
            code="metropolis_hastings_step_row_current_score_mismatch",
            details={
                "last_recorded_posterior_log_score": validated_step_rows[
                    -1
                ].recorded_posterior_log_score,
                "current_posterior_log_score": current_state.posterior_log_score,
            },
        )
    return validated_step_rows


def _require_mapping_payload(
    payload: Mapping[str, object] | object,
    *,
    key: str | None = None,
    owner_name: str,
) -> Mapping[str, object]:
    raw_value = payload if key is None else payload.get(key)
    if not isinstance(raw_value, Mapping):
        subject = owner_name if key is None else f"'{key}'"
        raise PhylogeneticsError(
            f"{owner_name} requires {subject} to be one mapping",
            code="metropolis_hastings_mapping_required",
        )
    return raw_value


def _require_list_payload(
    payload: Mapping[str, object],
    *,
    key: str,
    owner_name: str,
) -> list[object]:
    raw_value = payload.get(key)
    if not isinstance(raw_value, list):
        raise PhylogeneticsError(
            f"{owner_name} requires '{key}' to be one list",
            code="metropolis_hastings_list_required",
        )
    return raw_value


def _require_integer_payload(
    payload: Mapping[str, object],
    *,
    key: str,
    owner_name: str,
) -> int:
    raw_value = payload.get(key)
    if isinstance(raw_value, bool) or not isinstance(raw_value, int):
        raise PhylogeneticsError(
            f"{owner_name} requires '{key}' to be one integer",
            code="metropolis_hastings_integer_payload_required",
        )
    return raw_value


def _require_integer_list_payload(
    payload: Mapping[str, object],
    *,
    key: str,
    owner_name: str,
) -> list[int]:
    raw_values = _require_list_payload(
        payload,
        key=key,
        owner_name=owner_name,
    )
    return _validate_integer_sequence(
        values=raw_values,
        field_name=key,
        owner_name=owner_name,
        minimum_length=1,
    )


def _require_string_list_payload(
    payload: Mapping[str, object],
    *,
    key: str,
    owner_name: str,
) -> list[str]:
    raw_values = _require_list_payload(
        payload,
        key=key,
        owner_name=owner_name,
    )
    string_values: list[str] = []
    for raw_value in raw_values:
        if not isinstance(raw_value, str):
            raise PhylogeneticsError(
                f"{owner_name} requires every '{key}' item to be one string",
                code="metropolis_hastings_string_list_required",
            )
        string_values.append(raw_value)
    return string_values


def _require_float_payload(
    payload: Mapping[str, object],
    *,
    key: str,
    owner_name: str,
) -> float:
    raw_value = payload.get(key)
    if isinstance(raw_value, bool) or not isinstance(raw_value, (int, float)):
        raise PhylogeneticsError(
            f"{owner_name} requires '{key}' to be numeric",
            code="metropolis_hastings_float_payload_required",
        )
    normalized_value = float(raw_value)
    if not math.isfinite(normalized_value):
        raise PhylogeneticsError(
            f"{owner_name} requires '{key}' to be finite",
            code="metropolis_hastings_float_payload_nonfinite",
        )
    return normalized_value


def _optional_float_payload(
    payload: Mapping[str, object],
    *,
    key: str,
    owner_name: str,
) -> float | None:
    raw_value = payload.get(key)
    if raw_value is None:
        return None
    return _require_float_payload(payload, key=key, owner_name=owner_name)


def _require_boolean_payload(
    payload: Mapping[str, object],
    *,
    key: str,
    owner_name: str,
) -> bool:
    raw_value = payload.get(key)
    if not isinstance(raw_value, bool):
        raise PhylogeneticsError(
            f"{owner_name} requires '{key}' to be boolean",
            code="metropolis_hastings_boolean_payload_required",
        )
    return raw_value


def _optional_string_payload(
    payload: Mapping[str, object],
    *,
    key: str,
    owner_name: str,
) -> str | None:
    raw_value = payload.get(key)
    if raw_value is None:
        return None
    if not isinstance(raw_value, str):
        raise PhylogeneticsError(
            f"{owner_name} requires optional '{key}' values to be strings",
            code="metropolis_hastings_optional_string_payload_required",
        )
    return raw_value


def _validate_integer_field(
    *,
    value: int,
    field_name: str,
    owner_name: str,
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one integer",
            code="metropolis_hastings_integer_field_invalid",
            details={"field_name": field_name},
        )
    return value


def _validate_integer_sequence(
    *,
    values: Sequence[object],
    field_name: str,
    owner_name: str,
    minimum_length: int = 0,
) -> list[int]:
    normalized_values: list[int] = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, int):
            raise PhylogeneticsError(
                f"{owner_name} requires every '{field_name}' value to be one integer",
                code="metropolis_hastings_integer_sequence_invalid",
                details={"field_name": field_name},
            )
        normalized_values.append(value)
    if len(normalized_values) < minimum_length:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to contain at least {minimum_length} integers",
            code="metropolis_hastings_integer_sequence_too_short",
            details={
                "field_name": field_name,
                "minimum_length": minimum_length,
            },
        )
    return normalized_values


def _validate_positive_integer(
    *,
    value: int,
    field_name: str,
    owner_name: str,
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one integer",
            code="metropolis_hastings_integer_required",
            details={"field_name": field_name},
        )
    if value <= 0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be positive",
            code="metropolis_hastings_positive_integer_required",
            details={"field_name": field_name},
        )
    return value


def _validate_nonnegative_integer(
    *,
    value: int,
    field_name: str,
    owner_name: str,
) -> int:
    normalized_value = _validate_integer_field(
        value=value,
        field_name=field_name,
        owner_name=owner_name,
    )
    if normalized_value < 0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be nonnegative",
            code="metropolis_hastings_nonnegative_integer_required",
            details={"field_name": field_name},
        )
    return normalized_value


def _lognormal_positive_draw_density(
    *,
    proposed_value: float,
    log_standard_deviation: float,
) -> float:
    return (
        -math.log(proposed_value)
        - math.log(log_standard_deviation)
        - (0.5 * math.log(2.0 * math.pi))
        - (math.log(proposed_value) ** 2 / (2.0 * (log_standard_deviation**2)))
    )


def _validate_metropolis_hastings_proposal(
    proposal: MetropolisHastingsProposal,
) -> MetropolisHastingsProposal:
    if not isinstance(proposal, MetropolisHastingsProposal):
        raise PhylogeneticsError(
            "metropolis-hastings sampler requires proposals to return one MetropolisHastingsProposal",
            code="metropolis_hastings_proposal_type_invalid",
        )
    validated_changed_fields = _validate_changed_fields(proposal.changed_fields)
    validated_log_forward_density = _validate_finite_float(
        value=proposal.log_forward_density,
        field_name="log_forward_density",
        owner_name="metropolis-hastings proposal",
    )
    validated_log_reverse_density = _validate_finite_float(
        value=proposal.log_reverse_density,
        field_name="log_reverse_density",
        owner_name="metropolis-hastings proposal",
    )
    if not isinstance(proposal.is_valid, bool):
        raise PhylogeneticsError(
            "metropolis-hastings proposal requires 'is_valid' to be boolean",
            code="metropolis_hastings_proposal_validity_type_invalid",
        )
    validated_invalid_reason = _validate_invalid_reason(
        invalid_reason=proposal.invalid_reason,
        is_valid=proposal.is_valid,
    )
    if proposal.is_valid and not isinstance(
        proposal.proposed_model_parameters,
        BayesianModelParameterState,
    ):
        raise PhylogeneticsError(
            "metropolis-hastings proposal requires one BayesianModelParameterState",
            code="metropolis_hastings_model_parameter_state_type_invalid",
        )
    return MetropolisHastingsProposal(
        changed_fields=validated_changed_fields,
        log_forward_density=validated_log_forward_density,
        log_reverse_density=validated_log_reverse_density,
        is_valid=proposal.is_valid,
        invalid_reason=validated_invalid_reason,
        proposed_tree=proposal.proposed_tree,
        proposed_model_parameters=proposal.proposed_model_parameters,
    )


def _validate_integer_seed(seed: int) -> int:
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise PhylogeneticsError(
            "metropolis-hastings sampler requires 'seed' to be one integer",
            code="metropolis_hastings_seed_type_invalid",
        )
    return seed


def _validate_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be numeric",
            code="metropolis_hastings_float_required",
            details={"field_name": field_name},
        )
    normalized_value = float(value)
    if not math.isfinite(normalized_value):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be finite",
            code="metropolis_hastings_finite_float_required",
            details={"field_name": field_name},
        )
    return normalized_value


def _validate_positive_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    normalized_value = _validate_finite_float(
        value=value,
        field_name=field_name,
        owner_name=owner_name,
    )
    if normalized_value <= 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be strictly positive",
            code="metropolis_hastings_positive_float_required",
            details={"field_name": field_name},
        )
    return normalized_value


def _validate_changed_fields(
    changed_fields: tuple[str, ...],
) -> tuple[str, ...]:
    if not isinstance(changed_fields, tuple):
        raise PhylogeneticsError(
            "metropolis-hastings proposal requires 'changed_fields' to be one tuple",
            code="metropolis_hastings_changed_fields_type_invalid",
        )
    if not changed_fields:
        raise PhylogeneticsError(
            "metropolis-hastings proposal requires at least one changed field",
            code="metropolis_hastings_changed_fields_empty",
        )
    normalized_changed_fields: list[str] = []
    seen_changed_fields: set[str] = set()
    for changed_field in changed_fields:
        if not isinstance(changed_field, str) or not changed_field.strip():
            raise PhylogeneticsError(
                "metropolis-hastings proposal requires every changed field to be nonblank text",
                code="metropolis_hastings_changed_field_invalid",
            )
        normalized_changed_field = changed_field.strip()
        if normalized_changed_field in seen_changed_fields:
            raise PhylogeneticsError(
                "metropolis-hastings proposal requires changed fields to be unique",
                code="metropolis_hastings_changed_fields_duplicate",
                details={"changed_field": normalized_changed_field},
            )
        seen_changed_fields.add(normalized_changed_field)
        normalized_changed_fields.append(normalized_changed_field)
    return tuple(normalized_changed_fields)


def _validate_invalid_reason(
    *,
    invalid_reason: str | None,
    is_valid: bool,
) -> str | None:
    if is_valid:
        if invalid_reason is not None:
            raise PhylogeneticsError(
                "metropolis-hastings proposal does not allow 'invalid_reason' on valid proposals",
                code="metropolis_hastings_invalid_reason_not_allowed",
            )
        return None
    if not isinstance(invalid_reason, str) or not invalid_reason.strip():
        raise PhylogeneticsError(
            "metropolis-hastings proposal requires one nonblank invalid reason when 'is_valid' is false",
            code="metropolis_hastings_invalid_reason_required",
        )
    return invalid_reason.strip()


def _validate_parameter_name(
    *,
    value: str,
    field_name: str,
    owner_name: str,
) -> str:
    if not isinstance(value, str) or not value.strip():
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one nonblank parameter name",
            code="metropolis_hastings_parameter_name_invalid",
            details={"field_name": field_name},
        )
    return value.strip()


def _validate_partition_linking_partition_models(
    partition_models: Sequence[PartitionSubstitutionModelDefinition],
) -> tuple[PartitionSubstitutionModelDefinition, ...]:
    validated_partition_models = tuple(partition_models)
    if not validated_partition_models:
        raise PhylogeneticsError(
            "partition-linking proposal requires at least one partition model definition",
            code="partition_linking_partition_models_empty",
        )
    partition_names = tuple(
        partition_model.partition_name for partition_model in validated_partition_models
    )
    duplicate_partition_names = sorted(
        {
            partition_name
            for partition_name in partition_names
            if partition_names.count(partition_name) > 1
        }
    )
    if duplicate_partition_names:
        raise PhylogeneticsError(
            "partition-linking proposal requires unique partition names",
            code="partition_linking_partition_names_duplicated",
            details={"duplicate_partition_names": duplicate_partition_names},
        )
    return validated_partition_models


def _validate_partition_linking_target_name(target_name: str) -> str:
    validated_target_name = _validate_parameter_name(
        value=target_name,
        field_name="target_name",
        owner_name="partition-linking proposal",
    )
    if validated_target_name not in PARTITION_MODEL_PRIOR_TARGETS:
        raise PhylogeneticsError(
            "partition-linking proposal requires one supported partition target",
            code="partition_linking_target_invalid",
            details={
                "target_name": target_name,
                "allowed_target_names": list(PARTITION_MODEL_PRIOR_TARGETS),
            },
        )
    return validated_target_name


def _eligible_partition_linking_targets(
    partition_models: Sequence[PartitionSubstitutionModelDefinition],
) -> tuple[str, ...]:
    eligible_target_names: list[str] = []
    for target_name in PARTITION_MODEL_PRIOR_TARGETS:
        required_partition_count = sum(
            target_name in partition_model.required_targets()
            for partition_model in partition_models
        )
        if required_partition_count >= 2:
            eligible_target_names.append(target_name)
    return tuple(eligible_target_names)


def _classify_partition_linkage_state(
    *,
    group_by_partition_name: Mapping[str, str],
) -> str | None:
    group_names = tuple(group_by_partition_name.values())
    if len(set(group_names)) == 1:
        return "linked"
    if len(set(group_names)) == len(group_names):
        return "unlinked"
    return None


def _diff_model_parameter_changed_fields(
    current_model_parameters: BayesianModelParameterState,
    proposed_model_parameters: BayesianModelParameterState,
) -> tuple[str, ...]:
    changed_fields: list[str] = []
    for parameter_name in sorted(
        set(current_model_parameters.categorical_parameters)
        | set(proposed_model_parameters.categorical_parameters)
    ):
        if current_model_parameters.categorical_parameters.get(
            parameter_name
        ) != proposed_model_parameters.categorical_parameters.get(parameter_name):
            changed_fields.append(f"categorical_parameters.{parameter_name}")
    for parameter_name in sorted(
        set(current_model_parameters.scalar_parameters)
        | set(proposed_model_parameters.scalar_parameters)
    ):
        if current_model_parameters.scalar_parameters.get(
            parameter_name
        ) != proposed_model_parameters.scalar_parameters.get(parameter_name):
            changed_fields.append(f"scalar_parameters.{parameter_name}")
    for parameter_name in sorted(
        set(current_model_parameters.vector_parameters)
        | set(proposed_model_parameters.vector_parameters)
    ):
        if current_model_parameters.vector_parameters.get(
            parameter_name
        ) != proposed_model_parameters.vector_parameters.get(parameter_name):
            changed_fields.append(f"vector_parameters.{parameter_name}")
    return tuple(changed_fields)


def _copy_tree_with_scaled_branch_length(
    *,
    current_tree: PhyloTree,
    branch_id: str,
    proposed_branch_length: float,
) -> PhyloTree:
    proposed_tree = current_tree.copy()
    proposed_tree.node_by_id(branch_id).branch_length = proposed_branch_length
    return proposed_tree


def _copy_tree_with_globally_scaled_branch_lengths(
    *,
    current_tree: PhyloTree,
    scale_factor: float,
) -> PhyloTree:
    proposed_tree = current_tree.copy()
    for _parent, child in proposed_tree.iter_edges():
        if child.branch_length is None:
            raise PhylogeneticsError(
                "global tree-height scaling requires explicit branch lengths on every edge",
                code="metropolis_hastings_global_tree_scale_branch_length_missing",
            )
        child.branch_length = child.branch_length * scale_factor
    return proposed_tree


def _copy_tree_with_slid_node_height(
    *,
    current_tree: PhyloTree,
    node_id: str,
    proposed_height: float,
    node_heights: dict[str, float],
) -> PhyloTree:
    proposed_tree = current_tree.copy()
    proposed_node = proposed_tree.node_by_id(node_id)
    if proposed_node.parent is None:
        raise PhylogeneticsError(
            "node-height sliding requires one non-root internal node",
            code="metropolis_hastings_node_height_slide_root_not_allowed",
        )
    parent_height = node_heights[proposed_node.parent.node_id]
    proposed_node.branch_length = parent_height - proposed_height
    for child in proposed_node.children:
        child_height = node_heights[child.node_id]
        child.branch_length = proposed_height - child_height
    return proposed_tree


def _materialize_missing_branch_lengths(
    tree: PhyloTree,
    *,
    default_branch_length: float = 0.0,
) -> PhyloTree:
    proposed_tree = tree.copy()
    for _parent, child in proposed_tree.iter_edges():
        if child.branch_length is None:
            child.branch_length = default_branch_length
    return proposed_tree


def _lognormal_scaling_density(
    *,
    current_branch_length: float,
    proposed_branch_length: float,
    log_scale_standard_deviation: float,
) -> float:
    log_scale_change = math.log(proposed_branch_length / current_branch_length)
    z_score = log_scale_change / log_scale_standard_deviation
    return (
        -math.log(proposed_branch_length)
        - math.log(log_scale_standard_deviation)
        - (math.log(2.0 * math.pi) / 2.0)
        - ((z_score * z_score) / 2.0)
    )


def _gaussian_random_walk_density(
    *,
    coordinate_change: float,
    standard_deviation: float,
) -> float:
    z_score = coordinate_change / standard_deviation
    return (
        -math.log(standard_deviation)
        - (math.log(2.0 * math.pi) / 2.0)
        - ((z_score * z_score) / 2.0)
    )


def _logit_probability(value: float) -> float:
    return math.log(value) - math.log1p(-value)


def _inverse_logit_probability(value: float) -> float:
    if value >= 0.0:
        decay = math.exp(-value)
        return 1.0 / (1.0 + decay)
    growth = math.exp(value)
    return growth / (1.0 + growth)


def _log_probability_logit_jacobian(value: float) -> float:
    return -math.log(value) - math.log1p(-value)


def _global_tree_scaling_density(
    *,
    current_branch_lengths: list[float],
    proposed_branch_lengths: list[float],
    log_scale_standard_deviation: float,
) -> float:
    if len(current_branch_lengths) != len(proposed_branch_lengths):
        raise PhylogeneticsError(
            "global tree-height scaling requires matching branch dimensions",
            code="metropolis_hastings_global_tree_scale_dimension_mismatch",
        )
    log_scale_change = math.log(proposed_branch_lengths[0] / current_branch_lengths[0])
    z_score = log_scale_change / log_scale_standard_deviation
    return (
        -math.fsum(math.log(branch_length) for branch_length in proposed_branch_lengths)
        - math.log(log_scale_standard_deviation)
        - (math.log(2.0 * math.pi) / 2.0)
        - ((z_score * z_score) / 2.0)
    )


def _normal_node_height_slide_density(
    *,
    current_height: float,
    proposed_height: float,
    height_slide_standard_deviation: float,
) -> float:
    height_change = proposed_height - current_height
    z_score = height_change / height_slide_standard_deviation
    return (
        -math.log(height_slide_standard_deviation)
        - (math.log(2.0 * math.pi) / 2.0)
        - ((z_score * z_score) / 2.0)
    )


def _compute_rooted_ultrametric_node_heights(
    *,
    current_tree: PhyloTree,
    ultrametric_tolerance: float,
) -> dict[str, float]:
    if current_tree.rooted is not True:
        raise PhylogeneticsError(
            "node-height sliding requires one rooted ultrametric tree",
            code="metropolis_hastings_node_height_slide_tree_not_rooted",
        )
    node_depths: dict[str, float] = {}
    tip_depths: list[float] = []

    def visit(node, current_depth: float) -> None:
        if node.node_id is None:
            raise PhylogeneticsError(
                "node-height sliding requires stable node identifiers",
                code="metropolis_hastings_node_height_slide_node_id_missing",
            )
        node_depths[node.node_id] = current_depth
        if node.is_leaf():
            tip_depths.append(current_depth)
            return
        for child in node.children:
            branch_length = child.branch_length
            if branch_length is None or not math.isfinite(branch_length):
                raise PhylogeneticsError(
                    "node-height sliding requires explicit finite branch lengths on every edge",
                    code="metropolis_hastings_node_height_slide_branch_length_missing",
                )
            if branch_length < 0.0:
                raise PhylogeneticsError(
                    "node-height sliding requires non-negative branch lengths on every edge",
                    code="metropolis_hastings_node_height_slide_branch_length_negative",
                )
            visit(child, current_depth + branch_length)

    visit(current_tree.root, 0.0)
    if not tip_depths:
        raise PhylogeneticsError(
            "node-height sliding requires at least one tip",
            code="metropolis_hastings_node_height_slide_tip_missing",
        )
    root_age = max(tip_depths)
    if root_age - min(tip_depths) > ultrametric_tolerance:
        raise PhylogeneticsError(
            "node-height sliding requires one rooted ultrametric tree",
            code="metropolis_hastings_node_height_slide_tree_not_ultrametric",
        )
    return {
        node_id: root_age - node_depth for node_id, node_depth in node_depths.items()
    }


def _collect_node_height_slide_candidates(
    *,
    current_tree: PhyloTree,
    node_heights: dict[str, float],
) -> list[_NodeHeightSlideCandidate]:
    candidates: list[_NodeHeightSlideCandidate] = []
    for node in current_tree.iter_internal_nodes(order="preorder"):
        if node.parent is None or node.node_id is None:
            continue
        current_height = node_heights[node.node_id]
        lower_height_bound = max(node_heights[child.node_id] for child in node.children)
        upper_height_bound = node_heights[node.parent.node_id]
        if lower_height_bound < upper_height_bound:
            candidates.append(
                _NodeHeightSlideCandidate(
                    node_id=node.node_id,
                    current_height=current_height,
                    lower_height_bound=lower_height_bound,
                    upper_height_bound=upper_height_bound,
                )
            )
    return candidates


def _enumerate_rooted_nni_neighbors(
    current_tree: PhyloTree,
) -> list[tuple[RootedNniMoveCandidate, PhyloTree, str]]:
    neighbors: list[tuple[RootedNniMoveCandidate, PhyloTree, str]] = []
    for candidate in iter_rooted_nni_move_candidates(current_tree):
        neighbor_tree = apply_rooted_nni_move(current_tree, candidate)
        neighbors.append(
            (
                candidate,
                neighbor_tree,
                rooted_topology_fingerprint(neighbor_tree),
            )
        )
    return neighbors


def _rooted_spr_legal_move_count(report) -> int:
    return (
        report.generated_move_candidate_count
        - report.identity_move_candidate_count
        - report.self_regraft_candidate_count
    )


def _enumerate_reversible_rooted_tbr_neighbors(
    current_tree: PhyloTree,
) -> list[_ReversibleRootedTbrNeighbor]:
    current_topology_fingerprint = rooted_topology_fingerprint(current_tree)
    reversible_neighbors: list[_ReversibleRootedTbrNeighbor] = []
    for neighbor_row in enumerate_rooted_tbr_neighbors(current_tree).neighbor_rows:
        proposed_tree = _materialize_missing_branch_lengths(
            PhyloTree.from_newick(neighbor_row.neighbor_tree_newick)
        )
        proposed_tree.rooted = current_tree.rooted
        reverse_neighbor_row = next(
            (
                row
                for row in enumerate_rooted_tbr_neighbors(proposed_tree).neighbor_rows
                if row.neighbor_topology_fingerprint == current_topology_fingerprint
            ),
            None,
        )
        if reverse_neighbor_row is None:
            continue
        reversible_neighbors.append(
            _ReversibleRootedTbrNeighbor(
                neighbor_row=neighbor_row,
                proposed_tree=proposed_tree,
                reverse_neighbor_row=reverse_neighbor_row,
            )
        )
    return reversible_neighbors


def _reversible_rooted_tbr_move_count(
    neighbors: list[_ReversibleRootedTbrNeighbor],
) -> int:
    return sum(
        neighbor.neighbor_row.supporting_reconnection_count for neighbor in neighbors
    )


def _simplex_coordinate_component_names(
    component_names: tuple[str, ...],
    *,
    reference_component_name: str,
) -> tuple[str, ...]:
    return tuple(
        component_name
        for component_name in component_names
        if component_name != reference_component_name
    )
