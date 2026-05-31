from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import math

from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    FixedTopologyDnaModelDefinition,
    FixedTopologyDnaProposalSchedule,
    _build_fixed_topology_dna_prior_components,
    _build_initial_model_parameters,
    _evaluate_fixed_topology_dna_log_likelihood,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsRunReport,
    propose_base_frequency_simplex_move,
    propose_branch_length_scaling_move,
    propose_gtr_exchangeability_move,
    propose_kappa_move,
    propose_nni_topology_move,
    propose_spr_topology_move,
    propose_tbr_topology_move,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.required_values import require_present
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_prior_component_state,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    TreeTopologyPriorModel,
    evaluate_tree_topology_log_prior,
)
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna_observation_policies import (
    normalize_dna_likelihood_records,
    validate_dna_observation_policy,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

JOINT_TOPOLOGY_DNA_TOPOLOGY_MOVE_KINDS = ("nni", "spr", "tbr")


@dataclass(frozen=True, slots=True)
class JointTopologyDnaModelDefinition:
    """One validated joint topology-and-parameter DNA posterior model."""

    sequence_model_definition: FixedTopologyDnaModelDefinition
    topology_prior: TreeTopologyPriorModel

    @property
    def substitution_model_name(self) -> str:
        return self.sequence_model_definition.substitution_model_name

    @property
    def active_parameter_targets(self) -> tuple[str, ...]:
        return self.sequence_model_definition.active_parameter_targets


@dataclass(frozen=True, slots=True)
class JointTopologyDnaProposalSchedule:
    """One validated proposal schedule for joint DNA topology and parameter sampling."""

    sequence_proposal_schedule: FixedTopologyDnaProposalSchedule
    nni_move_weight: float
    spr_move_weight: float
    tbr_move_weight: float

    @property
    def substitution_model_name(self) -> str:
        return self.sequence_proposal_schedule.substitution_model_name


@dataclass(frozen=True, slots=True)
class JointTopologyDnaPosteriorRow:
    """One sampled posterior row from a joint topology-and-parameter DNA chain."""

    sample_index: int
    iteration_index: int
    topology_id: str
    substitution_model_name: str
    total_log_prior: float
    log_likelihood: float
    posterior_log_score: float
    prior_component_log_priors: dict[str, float]
    branch_lengths: dict[str, float]
    scalar_parameters: dict[str, float]
    vector_parameters: dict[str, dict[str, float]]


@dataclass(frozen=True, slots=True)
class JointTopologyDnaRunReport:
    """One completed joint DNA topology-and-parameter posterior run."""

    model_definition: JointTopologyDnaModelDefinition
    proposal_schedule: JointTopologyDnaProposalSchedule
    observation_policy: str
    chain_report: MetropolisHastingsRunReport
    posterior_rows: list[JointTopologyDnaPosteriorRow]
    distinct_topology_ids: list[str]

    @property
    def distinct_topology_count(self) -> int:
        return len(self.distinct_topology_ids)


def build_joint_topology_dna_model_definition(
    *,
    sequence_model_definition: FixedTopologyDnaModelDefinition,
    topology_prior: TreeTopologyPriorModel,
) -> JointTopologyDnaModelDefinition:
    """Build one validated joint DNA topology-and-parameter posterior model."""
    if not isinstance(sequence_model_definition, FixedTopologyDnaModelDefinition):
        raise PhylogeneticsError(
            "joint topology DNA posterior model requires one FixedTopologyDnaModelDefinition",
            code="joint_topology_dna_sequence_model_definition_type_invalid",
        )
    if not isinstance(topology_prior, TreeTopologyPriorModel):
        raise PhylogeneticsError(
            "joint topology DNA posterior model requires one TreeTopologyPriorModel",
            code="joint_topology_dna_topology_prior_type_invalid",
        )
    return JointTopologyDnaModelDefinition(
        sequence_model_definition=sequence_model_definition,
        topology_prior=topology_prior,
    )


def build_joint_topology_dna_proposal_schedule(
    *,
    sequence_proposal_schedule: FixedTopologyDnaProposalSchedule,
    nni_move_weight: float = 0.0,
    spr_move_weight: float = 0.0,
    tbr_move_weight: float = 0.0,
) -> JointTopologyDnaProposalSchedule:
    """Build one validated joint DNA proposal schedule with explicit topology moves."""
    if not isinstance(sequence_proposal_schedule, FixedTopologyDnaProposalSchedule):
        raise PhylogeneticsError(
            "joint topology DNA proposal schedule requires one FixedTopologyDnaProposalSchedule",
            code="joint_topology_dna_sequence_proposal_schedule_type_invalid",
        )
    validated_nni_move_weight = _validate_nonnegative_finite_float(
        value=nni_move_weight,
        field_name="nni_move_weight",
        owner_name="joint topology DNA proposal schedule",
    )
    validated_spr_move_weight = _validate_nonnegative_finite_float(
        value=spr_move_weight,
        field_name="spr_move_weight",
        owner_name="joint topology DNA proposal schedule",
    )
    validated_tbr_move_weight = _validate_nonnegative_finite_float(
        value=tbr_move_weight,
        field_name="tbr_move_weight",
        owner_name="joint topology DNA proposal schedule",
    )
    topology_weight_total = math.fsum(
        (
            validated_nni_move_weight,
            validated_spr_move_weight,
            validated_tbr_move_weight,
        )
    )
    if topology_weight_total <= 0.0:
        raise PhylogeneticsError(
            "joint topology DNA proposal schedule requires at least one positive topology move weight",
            code="joint_topology_dna_topology_move_weight_missing",
        )
    return JointTopologyDnaProposalSchedule(
        sequence_proposal_schedule=sequence_proposal_schedule,
        nni_move_weight=validated_nni_move_weight,
        spr_move_weight=validated_spr_move_weight,
        tbr_move_weight=validated_tbr_move_weight,
    )


def run_joint_topology_dna_metropolis_hastings(
    *,
    tree: PhyloTree,
    records: Sequence[AlignmentRecord],
    model_definition: JointTopologyDnaModelDefinition,
    proposal_schedule: JointTopologyDnaProposalSchedule,
    iteration_count: int,
    sample_every: int = 1,
    seed: int = 0,
    observation_policy: str = "reject",
) -> JointTopologyDnaRunReport:
    """Run one joint DNA posterior sampler over topology, branch lengths, and model parameters."""
    if not isinstance(tree, PhyloTree):
        raise PhylogeneticsError(
            "joint topology DNA posterior runner requires one PhyloTree",
            code="joint_topology_dna_tree_type_invalid",
        )
    if not isinstance(model_definition, JointTopologyDnaModelDefinition):
        raise PhylogeneticsError(
            "joint topology DNA posterior runner requires one JointTopologyDnaModelDefinition",
            code="joint_topology_dna_model_definition_type_invalid",
        )
    if not isinstance(proposal_schedule, JointTopologyDnaProposalSchedule):
        raise PhylogeneticsError(
            "joint topology DNA posterior runner requires one JointTopologyDnaProposalSchedule",
            code="joint_topology_dna_proposal_schedule_type_invalid",
        )
    if (
        proposal_schedule.substitution_model_name
        != model_definition.substitution_model_name
    ):
        raise PhylogeneticsError(
            "joint topology DNA posterior runner requires the proposal schedule and model definition to use the same substitution model",
            code="joint_topology_dna_model_schedule_mismatch",
            details={
                "model_definition": model_definition.substitution_model_name,
                "proposal_schedule": proposal_schedule.substitution_model_name,
            },
        )
    validated_observation_policy = validate_dna_observation_policy(
        observation_policy,
        owner_name="joint topology DNA posterior runner",
    )
    normalized_records = normalize_dna_likelihood_records(
        list(records),
        model_name=model_definition.substitution_model_name,
        observation_policy=validated_observation_policy,
    )
    working_tree = tree.copy()
    working_tree.rooted = tree.rooted
    initial_model_parameters = _build_initial_model_parameters(
        model_definition=model_definition.sequence_model_definition,
        records=normalized_records,
        observation_policy=validated_observation_policy,
    )
    chain_report = run_metropolis_hastings_sampler(
        initial_state=score_bayesian_phylogenetic_state(
            tree=working_tree,
            model_parameters=initial_model_parameters,
            update_prior_components=lambda state: (
                _build_joint_topology_dna_prior_components(
                    state=state,
                    model_definition=model_definition,
                )
            ),
            update_log_likelihood=lambda state: (
                _evaluate_joint_topology_dna_log_likelihood(
                    state=state,
                    records=normalized_records,
                    model_definition=model_definition,
                    observation_policy=validated_observation_policy,
                )
            ),
        ),
        propose_state=lambda current_state, rng: _propose_joint_topology_dna_state(
            current_state=current_state,
            rng=rng,
            proposal_schedule=proposal_schedule,
        ),
        update_prior_components=lambda state: (
            _build_joint_topology_dna_prior_components(
                state=state,
                model_definition=model_definition,
            )
        ),
        update_log_likelihood=lambda state: _evaluate_joint_topology_dna_log_likelihood(
            state=state,
            records=normalized_records,
            model_definition=model_definition,
            observation_policy=validated_observation_policy,
        ),
        iteration_count=iteration_count,
        sample_every=sample_every,
        seed=seed,
    )
    posterior_rows = _build_joint_topology_dna_posterior_rows(
        chain_report=chain_report,
        substitution_model_name=model_definition.substitution_model_name,
    )
    return JointTopologyDnaRunReport(
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        observation_policy=validated_observation_policy,
        chain_report=chain_report,
        posterior_rows=posterior_rows,
        distinct_topology_ids=sorted({row.topology_id for row in posterior_rows}),
    )


def _build_joint_topology_dna_prior_components(
    *,
    state: BayesianPhylogeneticState,
    model_definition: JointTopologyDnaModelDefinition,
) -> list[BayesianPriorComponentState]:
    topology_prior_report = evaluate_tree_topology_log_prior(
        state.tree.to_tree(),
        model_definition.topology_prior,
    )
    return [
        build_bayesian_prior_component_state(
            component_name="tree-topology",
            family=topology_prior_report.family,
            log_prior=topology_prior_report.log_prior,
            parameter_values={
                "topology_count": float(topology_prior_report.topology_count)
            },
        ),
        *_build_fixed_topology_dna_prior_components(
            state=state,
            model_definition=model_definition.sequence_model_definition,
            fixed_topology_id=None,
        ),
    ]


def _evaluate_joint_topology_dna_log_likelihood(
    *,
    state: BayesianPhylogeneticState,
    records: list[AlignmentRecord],
    model_definition: JointTopologyDnaModelDefinition,
    observation_policy: str,
) -> float:
    return _evaluate_fixed_topology_dna_log_likelihood(
        state=state,
        records=records,
        model_definition=model_definition.sequence_model_definition,
        fixed_topology_id=None,
        observation_policy=observation_policy,
    )


def _propose_joint_topology_dna_state(
    *,
    current_state: BayesianPhylogeneticState,
    rng,
    proposal_schedule: JointTopologyDnaProposalSchedule,
):
    sequence_proposal_schedule = proposal_schedule.sequence_proposal_schedule
    weighted_moves = []
    if proposal_schedule.nni_move_weight > 0.0:
        weighted_moves.append(
            (
                proposal_schedule.nni_move_weight,
                lambda: propose_nni_topology_move(current_state, rng),
            )
        )
    if proposal_schedule.spr_move_weight > 0.0:
        weighted_moves.append(
            (
                proposal_schedule.spr_move_weight,
                lambda: propose_spr_topology_move(current_state, rng),
            )
        )
    if proposal_schedule.tbr_move_weight > 0.0:
        weighted_moves.append(
            (
                proposal_schedule.tbr_move_weight,
                lambda: propose_tbr_topology_move(current_state, rng),
            )
        )
    weighted_moves.append(
        (
            sequence_proposal_schedule.branch_length_move_weight,
            lambda: propose_branch_length_scaling_move(
                current_state,
                rng,
                log_scale_standard_deviation=(
                    sequence_proposal_schedule.branch_length_log_scale_standard_deviation
                ),
            ),
        )
    )
    if sequence_proposal_schedule.kappa_move_weight > 0.0:
        kappa_log_scale_standard_deviation = require_present(
            sequence_proposal_schedule.kappa_log_scale_standard_deviation,
            owner_name="joint-topology DNA sequence proposal schedule",
            field_name="kappa_log_scale_standard_deviation",
        )
        weighted_moves.append(
            (
                sequence_proposal_schedule.kappa_move_weight,
                lambda: propose_kappa_move(
                    current_state,
                    rng,
                    log_scale_standard_deviation=kappa_log_scale_standard_deviation,
                ),
            )
        )
    if sequence_proposal_schedule.base_frequency_move_weight > 0.0:
        base_frequency_coordinate_standard_deviation = require_present(
            sequence_proposal_schedule.base_frequency_coordinate_standard_deviation,
            owner_name="joint-topology DNA sequence proposal schedule",
            field_name="base_frequency_coordinate_standard_deviation",
        )
        weighted_moves.append(
            (
                sequence_proposal_schedule.base_frequency_move_weight,
                lambda: propose_base_frequency_simplex_move(
                    current_state,
                    rng,
                    unconstrained_coordinate_standard_deviation=(
                        base_frequency_coordinate_standard_deviation
                    ),
                ),
            )
        )
    if sequence_proposal_schedule.exchangeability_move_weight > 0.0:
        exchangeability_coordinate_standard_deviation = require_present(
            sequence_proposal_schedule.exchangeability_coordinate_standard_deviation,
            owner_name="joint-topology DNA sequence proposal schedule",
            field_name="exchangeability_coordinate_standard_deviation",
        )
        weighted_moves.append(
            (
                sequence_proposal_schedule.exchangeability_move_weight,
                lambda: propose_gtr_exchangeability_move(
                    current_state,
                    rng,
                    unconstrained_coordinate_standard_deviation=(
                        exchangeability_coordinate_standard_deviation
                    ),
                ),
            )
        )
    total_weight = math.fsum(weight for weight, _move in weighted_moves)
    move_threshold = rng.random() * total_weight
    cumulative_weight = 0.0
    for weight, move in weighted_moves:
        cumulative_weight += weight
        if move_threshold <= cumulative_weight:
            return move()
    return weighted_moves[-1][1]()


def _build_joint_topology_dna_posterior_rows(
    *,
    chain_report: MetropolisHastingsRunReport,
    substitution_model_name: str,
) -> list[JointTopologyDnaPosteriorRow]:
    return [
        JointTopologyDnaPosteriorRow(
            sample_index=sample_index,
            iteration_index=sample_index * chain_report.sample_every,
            topology_id=state.tree.topology_id,
            substitution_model_name=substitution_model_name,
            total_log_prior=state.total_log_prior,
            log_likelihood=state.log_likelihood,
            posterior_log_score=state.posterior_log_score,
            prior_component_log_priors={
                component.component_name: component.log_prior
                for component in state.prior_components
            },
            branch_lengths={
                branch_row.branch_id: branch_row.branch_length
                for branch_row in state.tree.branch_rows
            },
            scalar_parameters=dict(state.model_parameters.scalar_parameters),
            vector_parameters={
                parameter_name: dict(component_values)
                for parameter_name, component_values in state.model_parameters.vector_parameters.items()
            },
        )
        for sample_index, state in enumerate(chain_report.sampled_states)
    ]


def _validate_nonnegative_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    try:
        validated_value = float(value)
    except (TypeError, ValueError) as error:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be a finite float",
            code="joint_topology_dna_float_type_invalid",
            details={"field_name": field_name},
        ) from error
    if not math.isfinite(validated_value) or validated_value < 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be nonnegative and finite",
            code="joint_topology_dna_nonnegative_float_invalid",
            details={"field_name": field_name, "value": value},
        )
    return float(format(validated_value, ".15g"))


__all__ = [
    "JOINT_TOPOLOGY_DNA_TOPOLOGY_MOVE_KINDS",
    "JointTopologyDnaModelDefinition",
    "JointTopologyDnaPosteriorRow",
    "JointTopologyDnaProposalSchedule",
    "JointTopologyDnaRunReport",
    "build_joint_topology_dna_model_definition",
    "build_joint_topology_dna_proposal_schedule",
    "run_joint_topology_dna_metropolis_hastings",
]
