from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
import math
import random

from bijux_phylogenetics.bayesian.branch_length_priors import (
    BranchLengthPriorModel,
    evaluate_tree_branch_length_log_prior,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsProposal,
    MetropolisHastingsRunReport,
    build_metropolis_hastings_proposal,
    propose_branch_length_scaling_move,
    propose_partition_linking_move,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    PARTITION_MODEL_PRIOR_TARGETS,
    PartitionModelPriorBundle,
    PartitionSubstitutionModelDefinition,
    PartitionSubstitutionParameterState,
    evaluate_partition_model_log_prior,
)
from bijux_phylogenetics.bayesian.partition_model_state import (
    build_partition_model_parameter_state,
    resolve_partition_parameter_linkage_plan_from_model_parameters,
    resolve_partition_parameter_states_from_model_parameters,
    strip_partition_model_parameter_state,
)
from bijux_phylogenetics.bayesian.required_values import require_present
from bijux_phylogenetics.bayesian.state import (
    BayesianModelParameterState,
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_prior_component_state,
)
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    normalize_partition_data_type,
    slice_partition_sequence,
    validate_locus_partitions,
)
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_f81_tree_likelihood,
    evaluate_gtr_tree_likelihood,
    evaluate_hky85_tree_likelihood,
    evaluate_jc69_tree_likelihood,
    evaluate_k80_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.dna import validate_positive_kappa
from bijux_phylogenetics.phylo.likelihood.dna_observation_policies import (
    estimate_empirical_dna_base_frequencies_from_records,
    normalize_dna_likelihood_records,
    validate_dna_observation_policy,
)
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    DNA_EXCHANGEABILITY_LABELS,
    parameterize_dna_base_frequency_simplex,
    parameterize_dna_exchangeability_simplex,
    resolve_dna_base_frequency_simplex_from_unconstrained,
    resolve_dna_exchangeability_simplex_from_unconstrained,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXED_TOPOLOGY_PARTITIONED_DNA_SUBSTITUTION_MODELS = (
    "F81",
    "GTR",
    "HKY85",
    "JC69",
    "K80",
)
_MODEL_NAME = "partitioned-dna"
_SUPPORTED_PARTITION_PARAMETER_TARGETS = frozenset(
    {"base-frequencies", "exchangeabilities", "kappa"}
)
_DEFAULT_INITIAL_KAPPA = 2.0
_DEFAULT_INITIAL_EXCHANGEABILITIES = dict.fromkeys(DNA_EXCHANGEABILITY_LABELS, 1.0)
_MINIMUM_SIMPLEX_COMPONENT = 1e-6


@dataclass(frozen=True, slots=True)
class FixedTopologyPartitionedDnaModelDefinition:
    """One validated fixed-topology partitioned DNA posterior model definition."""

    locus_partitions: tuple[LocusPartition, ...]
    partition_prior_bundle: PartitionModelPriorBundle
    branch_length_prior: BranchLengthPriorModel
    initial_partition_parameter_states: (
        tuple[PartitionSubstitutionParameterState, ...] | None
    ) = None

    @property
    def partition_models(self) -> tuple[PartitionSubstitutionModelDefinition, ...]:
        return self.partition_prior_bundle.partition_models

    @property
    def active_parameter_targets(self) -> tuple[str, ...]:
        return tuple(
            target_name
            for target_name in PARTITION_MODEL_PRIOR_TARGETS
            if target_name in _SUPPORTED_PARTITION_PARAMETER_TARGETS
            and any(
                target_name in partition_model.required_targets()
                for partition_model in self.partition_models
            )
        )

    @property
    def linkage_eligible_targets(self) -> tuple[str, ...]:
        return tuple(
            target_name
            for target_name in self.active_parameter_targets
            if sum(
                target_name in partition_model.required_targets()
                for partition_model in self.partition_models
            )
            > 1
        )


@dataclass(frozen=True, slots=True)
class FixedTopologyPartitionedDnaProposalSchedule:
    """One validated proposal schedule for partitioned DNA posterior sampling."""

    branch_length_move_weight: float
    branch_length_log_scale_standard_deviation: float
    kappa_move_weight: float
    kappa_log_scale_standard_deviation: float | None
    base_frequency_move_weight: float
    base_frequency_coordinate_standard_deviation: float | None
    exchangeability_move_weight: float
    exchangeability_coordinate_standard_deviation: float | None
    linkage_move_weight: float


@dataclass(frozen=True, slots=True)
class FixedTopologyPartitionedDnaPartitionRow:
    """One partition-level posterior snapshot inside one sampled state."""

    partition_name: str
    model_name: str
    log_likelihood: float
    linkage_groups: dict[str, str]
    scalar_parameters: dict[str, float]
    vector_parameters: dict[str, dict[str, float]]


@dataclass(frozen=True, slots=True)
class FixedTopologyPartitionedDnaPosteriorRow:
    """One sampled posterior row from a fixed-topology partitioned DNA chain."""

    sample_index: int
    iteration_index: int
    topology_id: str
    total_log_prior: float
    log_likelihood: float
    posterior_log_score: float
    prior_component_log_priors: dict[str, float]
    branch_lengths: dict[str, float]
    partition_rows: list[FixedTopologyPartitionedDnaPartitionRow]


@dataclass(frozen=True, slots=True)
class FixedTopologyPartitionedDnaRunReport:
    """One completed fixed-topology partitioned DNA posterior run."""

    model_definition: FixedTopologyPartitionedDnaModelDefinition
    proposal_schedule: FixedTopologyPartitionedDnaProposalSchedule
    observation_policy: str
    chain_report: MetropolisHastingsRunReport
    posterior_rows: list[FixedTopologyPartitionedDnaPosteriorRow]


def build_fixed_topology_partitioned_dna_model_definition(
    *,
    locus_partitions: Sequence[LocusPartition],
    partition_prior_bundle: PartitionModelPriorBundle,
    branch_length_prior: BranchLengthPriorModel,
    initial_partition_parameter_states: (
        Sequence[PartitionSubstitutionParameterState] | None
    ) = None,
) -> FixedTopologyPartitionedDnaModelDefinition:
    """Build one validated fixed-topology partitioned DNA posterior model."""
    if not isinstance(partition_prior_bundle, PartitionModelPriorBundle):
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior model requires one PartitionModelPriorBundle",
            code="fixed_topology_partitioned_dna_prior_bundle_type_invalid",
        )
    if not isinstance(branch_length_prior, BranchLengthPriorModel):
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior model requires one BranchLengthPriorModel",
            code="fixed_topology_partitioned_dna_branch_length_prior_type_invalid",
        )
    if branch_length_prior.family == "fixed":
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior model does not support fixed branch-length priors because Metropolis-Hastings proposals would leave finite prior support",
            code="fixed_topology_partitioned_dna_branch_length_prior_family_invalid",
        )
    validated_locus_partitions = _validate_locus_partitions(
        locus_partitions=locus_partitions,
        expected_partition_names=tuple(
            partition_model.partition_name
            for partition_model in partition_prior_bundle.partition_models
        ),
    )
    _validate_supported_partition_models(
        partition_models=partition_prior_bundle.partition_models,
    )
    _validate_supported_partition_prior_families(
        partition_prior_bundle=partition_prior_bundle,
    )
    validated_initial_partition_parameter_states = (
        _validate_initial_partition_parameter_states(
            partition_models=partition_prior_bundle.partition_models,
            partition_prior_bundle=partition_prior_bundle,
            initial_partition_parameter_states=initial_partition_parameter_states,
        )
        if initial_partition_parameter_states is not None
        else None
    )
    if not any(
        partition_model.required_targets()
        and any(
            target_name in _SUPPORTED_PARTITION_PARAMETER_TARGETS
            for target_name in partition_model.required_targets()
        )
        for partition_model in partition_prior_bundle.partition_models
    ):
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior model requires at least one partition to sample one supported substitution parameter",
            code="fixed_topology_partitioned_dna_active_targets_missing",
        )
    return FixedTopologyPartitionedDnaModelDefinition(
        locus_partitions=validated_locus_partitions,
        partition_prior_bundle=partition_prior_bundle,
        branch_length_prior=branch_length_prior,
        initial_partition_parameter_states=validated_initial_partition_parameter_states,
    )


def build_fixed_topology_partitioned_dna_proposal_schedule(
    *,
    model_definition: FixedTopologyPartitionedDnaModelDefinition,
    branch_length_move_weight: float,
    branch_length_log_scale_standard_deviation: float,
    kappa_move_weight: float = 0.0,
    kappa_log_scale_standard_deviation: float | None = None,
    base_frequency_move_weight: float = 0.0,
    base_frequency_coordinate_standard_deviation: float | None = None,
    exchangeability_move_weight: float = 0.0,
    exchangeability_coordinate_standard_deviation: float | None = None,
    linkage_move_weight: float = 0.0,
) -> FixedTopologyPartitionedDnaProposalSchedule:
    """Build one validated proposal schedule for partitioned DNA sampling."""
    if not isinstance(model_definition, FixedTopologyPartitionedDnaModelDefinition):
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA proposal schedule requires one FixedTopologyPartitionedDnaModelDefinition",
            code="fixed_topology_partitioned_dna_proposal_schedule_model_definition_type_invalid",
        )
    validated_branch_length_move_weight = _validate_nonnegative_finite_float(
        value=branch_length_move_weight,
        field_name="branch_length_move_weight",
        owner_name="fixed-topology partitioned DNA proposal schedule",
    )
    if validated_branch_length_move_weight <= 0.0:
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA proposal schedule requires 'branch_length_move_weight' to be greater than zero",
            code="fixed_topology_partitioned_dna_branch_length_move_weight_invalid",
        )
    validated_branch_length_log_scale_standard_deviation = (
        _validate_positive_finite_float(
            value=branch_length_log_scale_standard_deviation,
            field_name="branch_length_log_scale_standard_deviation",
            owner_name="fixed-topology partitioned DNA proposal schedule",
        )
    )
    active_targets = set(model_definition.active_parameter_targets)
    validated_kappa_move_weight = _validate_nonnegative_finite_float(
        value=kappa_move_weight,
        field_name="kappa_move_weight",
        owner_name="fixed-topology partitioned DNA proposal schedule",
    )
    validated_base_frequency_move_weight = _validate_nonnegative_finite_float(
        value=base_frequency_move_weight,
        field_name="base_frequency_move_weight",
        owner_name="fixed-topology partitioned DNA proposal schedule",
    )
    validated_exchangeability_move_weight = _validate_nonnegative_finite_float(
        value=exchangeability_move_weight,
        field_name="exchangeability_move_weight",
        owner_name="fixed-topology partitioned DNA proposal schedule",
    )
    validated_linkage_move_weight = _validate_nonnegative_finite_float(
        value=linkage_move_weight,
        field_name="linkage_move_weight",
        owner_name="fixed-topology partitioned DNA proposal schedule",
    )
    validated_kappa_log_scale_standard_deviation = (
        _validate_optional_positive_finite_float(
            value=kappa_log_scale_standard_deviation,
            field_name="kappa_log_scale_standard_deviation",
            owner_name="fixed-topology partitioned DNA proposal schedule",
        )
    )
    validated_base_frequency_coordinate_standard_deviation = (
        _validate_optional_positive_finite_float(
            value=base_frequency_coordinate_standard_deviation,
            field_name="base_frequency_coordinate_standard_deviation",
            owner_name="fixed-topology partitioned DNA proposal schedule",
        )
    )
    validated_exchangeability_coordinate_standard_deviation = (
        _validate_optional_positive_finite_float(
            value=exchangeability_coordinate_standard_deviation,
            field_name="exchangeability_coordinate_standard_deviation",
            owner_name="fixed-topology partitioned DNA proposal schedule",
        )
    )
    _validate_parameter_move_activation(
        target_name="kappa",
        active_targets=active_targets,
        move_weight=validated_kappa_move_weight,
        standard_deviation=validated_kappa_log_scale_standard_deviation,
    )
    _validate_parameter_move_activation(
        target_name="base-frequencies",
        active_targets=active_targets,
        move_weight=validated_base_frequency_move_weight,
        standard_deviation=validated_base_frequency_coordinate_standard_deviation,
    )
    _validate_parameter_move_activation(
        target_name="exchangeabilities",
        active_targets=active_targets,
        move_weight=validated_exchangeability_move_weight,
        standard_deviation=validated_exchangeability_coordinate_standard_deviation,
    )
    if (
        validated_linkage_move_weight > 0.0
        and not model_definition.linkage_eligible_targets
    ):
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA proposal schedule can activate linkage moves only when at least one supported target is shared by more than one partition",
            code="fixed_topology_partitioned_dna_linkage_move_weight_invalid",
        )
    return FixedTopologyPartitionedDnaProposalSchedule(
        branch_length_move_weight=validated_branch_length_move_weight,
        branch_length_log_scale_standard_deviation=(
            validated_branch_length_log_scale_standard_deviation
        ),
        kappa_move_weight=validated_kappa_move_weight,
        kappa_log_scale_standard_deviation=validated_kappa_log_scale_standard_deviation,
        base_frequency_move_weight=validated_base_frequency_move_weight,
        base_frequency_coordinate_standard_deviation=(
            validated_base_frequency_coordinate_standard_deviation
        ),
        exchangeability_move_weight=validated_exchangeability_move_weight,
        exchangeability_coordinate_standard_deviation=(
            validated_exchangeability_coordinate_standard_deviation
        ),
        linkage_move_weight=validated_linkage_move_weight,
    )


def run_fixed_topology_partitioned_dna_metropolis_hastings(
    *,
    tree: PhyloTree,
    records: Sequence[AlignmentRecord],
    model_definition: FixedTopologyPartitionedDnaModelDefinition,
    proposal_schedule: FixedTopologyPartitionedDnaProposalSchedule,
    iteration_count: int,
    sample_every: int = 1,
    seed: int = 0,
    observation_policy: str = "reject",
) -> FixedTopologyPartitionedDnaRunReport:
    """Run one fixed-topology partitioned DNA posterior sampler."""
    if not isinstance(tree, PhyloTree):
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior runner requires one PhyloTree",
            code="fixed_topology_partitioned_dna_tree_type_invalid",
        )
    if not isinstance(model_definition, FixedTopologyPartitionedDnaModelDefinition):
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior runner requires one FixedTopologyPartitionedDnaModelDefinition",
            code="fixed_topology_partitioned_dna_model_definition_type_invalid",
        )
    if not isinstance(proposal_schedule, FixedTopologyPartitionedDnaProposalSchedule):
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior runner requires one FixedTopologyPartitionedDnaProposalSchedule",
            code="fixed_topology_partitioned_dna_proposal_schedule_type_invalid",
        )
    validated_observation_policy = validate_dna_observation_policy(
        observation_policy,
        owner_name="fixed-topology partitioned DNA posterior runner",
    )
    normalized_records = normalize_dna_likelihood_records(
        list(records),
        model_name=_MODEL_NAME,
        observation_policy=validated_observation_policy,
    )
    partition_records_by_name = _build_partition_record_sets(
        records=normalized_records,
        locus_partitions=model_definition.locus_partitions,
    )
    fixed_tree = tree.copy()
    fixed_tree.rooted = tree.rooted
    initial_model_parameters = _build_initial_model_parameters(
        model_definition=model_definition,
        partition_records_by_name=partition_records_by_name,
        observation_policy=validated_observation_policy,
    )
    initial_state = score_bayesian_phylogenetic_state(
        tree=fixed_tree,
        model_parameters=initial_model_parameters,
        update_prior_components=lambda state: (
            _build_fixed_topology_partitioned_dna_prior_components(
                state=state,
                model_definition=model_definition,
                fixed_topology_id=None,
            )
        ),
        update_log_likelihood=lambda state: (
            _evaluate_fixed_topology_partitioned_dna_log_likelihood(
                state=state,
                model_definition=model_definition,
                fixed_topology_id=None,
                partition_records_by_name=partition_records_by_name,
                observation_policy=validated_observation_policy,
            )
        ),
    )
    fixed_topology_id = initial_state.tree.topology_id
    chain_report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=lambda current_state, rng: (
            _propose_fixed_topology_partitioned_dna_state(
                current_state=current_state,
                rng=rng,
                model_definition=model_definition,
                proposal_schedule=proposal_schedule,
            )
        ),
        update_prior_components=lambda state: (
            _build_fixed_topology_partitioned_dna_prior_components(
                state=state,
                model_definition=model_definition,
                fixed_topology_id=fixed_topology_id,
            )
        ),
        update_log_likelihood=lambda state: (
            _evaluate_fixed_topology_partitioned_dna_log_likelihood(
                state=state,
                model_definition=model_definition,
                fixed_topology_id=fixed_topology_id,
                partition_records_by_name=partition_records_by_name,
                observation_policy=validated_observation_policy,
            )
        ),
        iteration_count=iteration_count,
        sample_every=sample_every,
        seed=seed,
    )
    posterior_rows = _build_fixed_topology_partitioned_dna_posterior_rows(
        chain_report=chain_report,
        model_definition=model_definition,
        fixed_topology_id=fixed_topology_id,
        partition_records_by_name=partition_records_by_name,
        observation_policy=validated_observation_policy,
    )
    return FixedTopologyPartitionedDnaRunReport(
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        observation_policy=validated_observation_policy,
        chain_report=chain_report,
        posterior_rows=posterior_rows,
    )


def _build_partition_record_sets(
    *,
    records: Sequence[AlignmentRecord],
    locus_partitions: Sequence[LocusPartition],
) -> dict[str, list[AlignmentRecord]]:
    if not records:
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior runner requires at least one alignment record",
            code="fixed_topology_partitioned_dna_records_empty",
        )
    alignment_length = len(records[0].sequence)
    _assigned_site_count, unassigned_site_count = validate_locus_partitions(
        tuple(locus_partitions),
        alignment_length=alignment_length,
    )
    if unassigned_site_count != 0:
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior runner requires locus partitions to assign every alignment site",
            code="fixed_topology_partitioned_dna_partition_coverage_incomplete",
            details={
                "alignment_length": alignment_length,
                "unassigned_site_count": unassigned_site_count,
            },
        )
    partition_records_by_name: dict[str, list[AlignmentRecord]] = {}
    for locus_partition in locus_partitions:
        partition_records_by_name[locus_partition.name] = [
            AlignmentRecord(
                identifier=record.identifier,
                sequence=slice_partition_sequence(record.sequence, locus_partition),
            )
            for record in records
        ]
    return partition_records_by_name


def _build_initial_model_parameters(
    *,
    model_definition: FixedTopologyPartitionedDnaModelDefinition,
    partition_records_by_name: dict[str, list[AlignmentRecord]],
    observation_policy: str,
) -> BayesianModelParameterState:
    if model_definition.initial_partition_parameter_states is not None:
        partition_parameter_states = model_definition.initial_partition_parameter_states
    else:
        partition_parameter_states = _build_default_partition_parameter_states(
            model_definition=model_definition,
            partition_records_by_name=partition_records_by_name,
            observation_policy=observation_policy,
        )
    return build_partition_model_parameter_state(
        partition_models=model_definition.partition_models,
        linkage_plan=model_definition.partition_prior_bundle.linkage_plan,
        partition_parameter_states=partition_parameter_states,
        preserved_categorical_parameters={"substitution-model": _MODEL_NAME},
    )


def _build_default_partition_parameter_states(
    *,
    model_definition: FixedTopologyPartitionedDnaModelDefinition,
    partition_records_by_name: dict[str, list[AlignmentRecord]],
    observation_policy: str,
) -> tuple[PartitionSubstitutionParameterState, ...]:
    default_states = tuple(
        _build_default_partition_parameter_state(
            partition_model=partition_model,
            partition_records=partition_records_by_name[partition_model.partition_name],
            observation_policy=observation_policy,
        )
        for partition_model in model_definition.partition_models
    )
    return _harmonize_linked_partition_parameter_states(
        partition_models=model_definition.partition_models,
        partition_parameter_states=default_states,
        model_parameters=build_partition_model_parameter_state(
            partition_models=model_definition.partition_models,
            linkage_plan=model_definition.partition_prior_bundle.linkage_plan,
            partition_parameter_states=_coerce_partition_states_to_linkage_representatives(
                partition_models=model_definition.partition_models,
                partition_parameter_states=default_states,
                linkage_plan=model_definition.partition_prior_bundle.linkage_plan,
            ),
            preserved_categorical_parameters={"substitution-model": _MODEL_NAME},
        ),
        linkage_plan=model_definition.partition_prior_bundle.linkage_plan,
    )


def _build_default_partition_parameter_state(
    *,
    partition_model: PartitionSubstitutionModelDefinition,
    partition_records: Sequence[AlignmentRecord],
    observation_policy: str,
) -> PartitionSubstitutionParameterState:
    state_kwargs: dict[str, object] = {"partition_name": partition_model.partition_name}
    if "kappa" in partition_model.required_targets():
        state_kwargs["kappa"] = _DEFAULT_INITIAL_KAPPA
    if "base-frequencies" in partition_model.required_targets():
        state_kwargs["base_frequencies"] = parameterize_dna_base_frequency_simplex(
            _stabilize_positive_simplex_mapping(
                estimate_empirical_dna_base_frequencies_from_records(
                    partition_records,
                    model_name=partition_model.base_model_name,
                    observation_policy=observation_policy,
                )
            )
        ).constrained_mapping()
    if "exchangeabilities" in partition_model.required_targets():
        state_kwargs["exchangeabilities"] = parameterize_dna_exchangeability_simplex(
            _DEFAULT_INITIAL_EXCHANGEABILITIES
        ).constrained_mapping()
    return PartitionSubstitutionParameterState(**state_kwargs)


def _coerce_partition_states_to_linkage_representatives(
    *,
    partition_models: Sequence[PartitionSubstitutionModelDefinition],
    partition_parameter_states: Sequence[PartitionSubstitutionParameterState],
    linkage_plan,
) -> tuple[PartitionSubstitutionParameterState, ...]:
    state_by_partition_name = {
        state.partition_name: state for state in partition_parameter_states
    }
    partition_state_payloads = {
        state.partition_name: _partition_state_payload(state)
        for state in partition_parameter_states
    }
    for target_name in _SUPPORTED_PARTITION_PARAMETER_TARGETS:
        required_partition_names = [
            partition_model.partition_name
            for partition_model in partition_models
            if target_name in partition_model.required_targets()
        ]
        if not required_partition_names:
            continue
        grouped_partition_names: dict[str, list[str]] = defaultdict(list)
        target_groups = linkage_plan.groups_for_target(target_name)
        for partition_name in required_partition_names:
            grouped_partition_names[target_groups[partition_name]].append(
                partition_name
            )
        for grouped_names in grouped_partition_names.values():
            representative_payload = partition_state_payloads[grouped_names[0]]
            representative_value = _partition_target_value_from_payload(
                payload=representative_payload,
                target_name=target_name,
            )
            for partition_name in grouped_names[1:]:
                _assign_partition_target_value_to_payload(
                    payload=partition_state_payloads[partition_name],
                    target_name=target_name,
                    target_value=representative_value,
                )
    return tuple(
        PartitionSubstitutionParameterState(
            **partition_state_payloads[state.partition_name]
        )
        for state in partition_parameter_states
        if state.partition_name in state_by_partition_name
    )


def _harmonize_linked_partition_parameter_states(
    *,
    partition_models: Sequence[PartitionSubstitutionModelDefinition],
    partition_parameter_states: Sequence[PartitionSubstitutionParameterState],
    model_parameters: BayesianModelParameterState,
    linkage_plan,
) -> tuple[PartitionSubstitutionParameterState, ...]:
    return resolve_partition_parameter_states_from_model_parameters(
        model_parameters=model_parameters,
        partition_models=partition_models,
        linkage_plan=linkage_plan,
    )


def _build_fixed_topology_partitioned_dna_prior_components(
    *,
    state: BayesianPhylogeneticState,
    model_definition: FixedTopologyPartitionedDnaModelDefinition,
    fixed_topology_id: str | None,
) -> list[BayesianPriorComponentState]:
    linkage_plan, partition_parameter_states = (
        _require_fixed_topology_partitioned_dna_state_consistency(
            state=state,
            model_definition=model_definition,
            fixed_topology_id=fixed_topology_id,
        )
    )
    branch_length_prior_report = evaluate_tree_branch_length_log_prior(
        state.tree.to_tree(),
        model_definition.branch_length_prior,
    )
    partition_prior_report = evaluate_partition_model_log_prior(
        prior_bundle=PartitionModelPriorBundle(
            partition_models=model_definition.partition_models,
            linkage_plan=linkage_plan,
            substitution_prior_bundle=(
                model_definition.partition_prior_bundle.substitution_prior_bundle
            ),
        ),
        partition_parameter_states=partition_parameter_states,
    )
    prior_components = [
        build_bayesian_prior_component_state(
            component_name="branch-lengths",
            family=branch_length_prior_report.family,
            log_prior=branch_length_prior_report.total_log_prior,
            parameter_values=branch_length_prior_report.parameter_values,
        )
    ]
    prior_components.extend(
        _build_partition_prior_component(row) for row in partition_prior_report.rows
    )
    return prior_components


def _build_partition_prior_component(row) -> BayesianPriorComponentState:
    return build_bayesian_prior_component_state(
        component_name=f"partition-substitution:{row.target_name}:{row.group_name}",
        family=row.family,
        log_prior=row.log_prior_contribution,
        parameter_values=row.hyperparameter_values,
    )


def _evaluate_fixed_topology_partitioned_dna_log_likelihood(
    *,
    state: BayesianPhylogeneticState,
    model_definition: FixedTopologyPartitionedDnaModelDefinition,
    fixed_topology_id: str | None,
    partition_records_by_name: dict[str, list[AlignmentRecord]],
    observation_policy: str,
) -> float:
    _linkage_plan, partition_parameter_states = (
        _require_fixed_topology_partitioned_dna_state_consistency(
            state=state,
            model_definition=model_definition,
            fixed_topology_id=fixed_topology_id,
        )
    )
    partition_state_by_name = {
        partition_state.partition_name: partition_state
        for partition_state in partition_parameter_states
    }
    tree = state.tree.to_tree()
    return float(
        format(
            math.fsum(
                _evaluate_partition_log_likelihood(
                    tree=tree,
                    partition_model=partition_model,
                    partition_state=partition_state_by_name[
                        partition_model.partition_name
                    ],
                    partition_records=partition_records_by_name[
                        partition_model.partition_name
                    ],
                    observation_policy=observation_policy,
                )
                for partition_model in model_definition.partition_models
            ),
            ".15g",
        )
    )


def _evaluate_partition_log_likelihood(
    *,
    tree: PhyloTree,
    partition_model: PartitionSubstitutionModelDefinition,
    partition_state: PartitionSubstitutionParameterState,
    partition_records: Sequence[AlignmentRecord],
    observation_policy: str,
) -> float:
    model_name = partition_model.base_model_name
    if model_name == "JC69":
        return evaluate_jc69_tree_likelihood(
            tree,
            partition_records,
            observation_policy=observation_policy,
        ).log_likelihood
    if model_name == "F81":
        if partition_state.base_frequencies is None:
            raise PhylogeneticsError(
                "fixed-topology partitioned DNA likelihood requires base frequencies for F81 partitions",
                code="fixed_topology_partitioned_dna_f81_base_frequencies_missing",
                details={"partition_name": partition_model.partition_name},
            )
        return evaluate_f81_tree_likelihood(
            tree,
            partition_records,
            base_frequencies=partition_state.base_frequencies,
            observation_policy=observation_policy,
        ).log_likelihood
    if model_name == "K80":
        if partition_state.kappa is None:
            raise PhylogeneticsError(
                "fixed-topology partitioned DNA likelihood requires kappa for K80 partitions",
                code="fixed_topology_partitioned_dna_k80_kappa_missing",
                details={"partition_name": partition_model.partition_name},
            )
        return evaluate_k80_tree_likelihood(
            tree,
            partition_records,
            kappa=partition_state.kappa,
            observation_policy=observation_policy,
        ).log_likelihood
    if model_name == "HKY85":
        if partition_state.kappa is None or partition_state.base_frequencies is None:
            raise PhylogeneticsError(
                "fixed-topology partitioned DNA likelihood requires kappa and base frequencies for HKY85 partitions",
                code="fixed_topology_partitioned_dna_hky85_parameters_missing",
                details={"partition_name": partition_model.partition_name},
            )
        return evaluate_hky85_tree_likelihood(
            tree,
            partition_records,
            kappa=partition_state.kappa,
            base_frequencies=partition_state.base_frequencies,
            observation_policy=observation_policy,
        ).log_likelihood
    if model_name == "GTR":
        if (
            partition_state.exchangeabilities is None
            or partition_state.base_frequencies is None
        ):
            raise PhylogeneticsError(
                "fixed-topology partitioned DNA likelihood requires exchangeabilities and base frequencies for GTR partitions",
                code="fixed_topology_partitioned_dna_gtr_parameters_missing",
                details={"partition_name": partition_model.partition_name},
            )
        return evaluate_gtr_tree_likelihood(
            tree,
            partition_records,
            exchangeabilities=partition_state.exchangeabilities,
            base_frequencies=partition_state.base_frequencies,
            observation_policy=observation_policy,
        ).log_likelihood
    raise AssertionError(
        f"unsupported fixed-topology partitioned DNA model {model_name}"
    )


def _propose_fixed_topology_partitioned_dna_state(
    *,
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    model_definition: FixedTopologyPartitionedDnaModelDefinition,
    proposal_schedule: FixedTopologyPartitionedDnaProposalSchedule,
) -> MetropolisHastingsProposal:
    weighted_moves = [
        (
            proposal_schedule.branch_length_move_weight,
            lambda: propose_branch_length_scaling_move(
                current_state,
                rng,
                log_scale_standard_deviation=(
                    proposal_schedule.branch_length_log_scale_standard_deviation
                ),
            ),
        )
    ]
    if proposal_schedule.kappa_move_weight > 0.0:
        kappa_log_scale_standard_deviation = require_present(
            proposal_schedule.kappa_log_scale_standard_deviation,
            owner_name="fixed-topology partitioned DNA proposal schedule",
            field_name="kappa_log_scale_standard_deviation",
        )
        weighted_moves.append(
            (
                proposal_schedule.kappa_move_weight,
                lambda: _propose_partition_kappa_move(
                    current_state=current_state,
                    rng=rng,
                    partition_models=model_definition.partition_models,
                    log_scale_standard_deviation=kappa_log_scale_standard_deviation,
                ),
            )
        )
    if proposal_schedule.base_frequency_move_weight > 0.0:
        base_frequency_coordinate_standard_deviation = require_present(
            proposal_schedule.base_frequency_coordinate_standard_deviation,
            owner_name="fixed-topology partitioned DNA proposal schedule",
            field_name="base_frequency_coordinate_standard_deviation",
        )
        weighted_moves.append(
            (
                proposal_schedule.base_frequency_move_weight,
                lambda: _propose_partition_base_frequency_move(
                    current_state=current_state,
                    rng=rng,
                    partition_models=model_definition.partition_models,
                    unconstrained_coordinate_standard_deviation=(
                        base_frequency_coordinate_standard_deviation
                    ),
                ),
            )
        )
    if proposal_schedule.exchangeability_move_weight > 0.0:
        exchangeability_coordinate_standard_deviation = require_present(
            proposal_schedule.exchangeability_coordinate_standard_deviation,
            owner_name="fixed-topology partitioned DNA proposal schedule",
            field_name="exchangeability_coordinate_standard_deviation",
        )
        weighted_moves.append(
            (
                proposal_schedule.exchangeability_move_weight,
                lambda: _propose_partition_exchangeability_move(
                    current_state=current_state,
                    rng=rng,
                    partition_models=model_definition.partition_models,
                    unconstrained_coordinate_standard_deviation=(
                        exchangeability_coordinate_standard_deviation
                    ),
                ),
            )
        )
    if proposal_schedule.linkage_move_weight > 0.0:
        weighted_moves.append(
            (
                proposal_schedule.linkage_move_weight,
                lambda: propose_partition_linking_move(
                    current_state,
                    rng,
                    partition_models=model_definition.partition_models,
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


def _propose_partition_kappa_move(
    *,
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    partition_models: Sequence[PartitionSubstitutionModelDefinition],
    log_scale_standard_deviation: float,
) -> MetropolisHastingsProposal:
    validated_log_scale_standard_deviation = _validate_positive_finite_float(
        value=log_scale_standard_deviation,
        field_name="log_scale_standard_deviation",
        owner_name="partitioned DNA kappa proposal",
    )
    selected_group_name, current_value, prepared = _select_partition_target_group(
        current_state=current_state,
        rng=rng,
        partition_models=partition_models,
        target_name="kappa",
    )
    changed_field = f"scalar_parameters.partition-parameter:kappa:{selected_group_name}"
    if selected_group_name is None or current_value is None or prepared is None:
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason="partitioned DNA kappa proposal requires at least one linked or unlinked kappa group",
        )
    try:
        validated_current_kappa = validate_positive_kappa(
            current_value,
            model_name="partitioned DNA kappa proposal",
        )
    except ValueError as error:
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
            invalid_reason="partitioned DNA kappa proposal scaling factor overflowed",
        )
    proposed_kappa = validated_current_kappa * scale_factor
    try:
        validated_proposed_kappa = validate_positive_kappa(
            proposed_kappa,
            model_name="partitioned DNA kappa proposal",
        )
    except ValueError as error:
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    proposed_model_parameters = _rebuild_partition_model_parameters_with_group_value(
        prepared=prepared,
        target_name="kappa",
        group_name=selected_group_name,
        target_value=validated_proposed_kappa,
    )
    log_forward_density = _log_group_selection_density(
        group_count=prepared["group_count"]
    ) + _lognormal_scaling_density(
        current_value=validated_current_kappa,
        proposed_value=validated_proposed_kappa,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    log_reverse_density = _log_group_selection_density(
        group_count=prepared["group_count"]
    ) + _lognormal_scaling_density(
        current_value=validated_proposed_kappa,
        proposed_value=validated_current_kappa,
        log_scale_standard_deviation=validated_log_scale_standard_deviation,
    )
    return build_metropolis_hastings_proposal(
        changed_fields=(changed_field,),
        log_forward_density=log_forward_density,
        log_reverse_density=log_reverse_density,
        is_valid=True,
        proposed_tree=_copy_current_tree(current_state),
        proposed_model_parameters=proposed_model_parameters,
    )


def _propose_partition_base_frequency_move(
    *,
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    partition_models: Sequence[PartitionSubstitutionModelDefinition],
    unconstrained_coordinate_standard_deviation: float,
) -> MetropolisHastingsProposal:
    validated_coordinate_standard_deviation = _validate_positive_finite_float(
        value=unconstrained_coordinate_standard_deviation,
        field_name="unconstrained_coordinate_standard_deviation",
        owner_name="partitioned DNA base-frequency proposal",
    )
    selected_group_name, current_value, prepared = _select_partition_target_group(
        current_state=current_state,
        rng=rng,
        partition_models=partition_models,
        target_name="base-frequencies",
    )
    changed_field = (
        f"vector_parameters.partition-parameter:base-frequencies:{selected_group_name}"
    )
    if selected_group_name is None or current_value is None or prepared is None:
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "partitioned DNA base-frequency proposal requires at least one linked or unlinked base-frequency group"
            ),
        )
    parameterization = parameterize_dna_base_frequency_simplex(current_value)
    current_unconstrained_values = list(parameterization.unconstrained_values)
    coordinate_component_names = _simplex_coordinate_component_names(
        parameterization.component_names,
        reference_component_name=parameterization.reference_component_name,
    )
    selected_coordinate_index = rng.randrange(len(current_unconstrained_values))
    current_coordinate_value = current_unconstrained_values[selected_coordinate_index]
    proposed_coordinate_value = current_coordinate_value + rng.gauss(
        0.0,
        validated_coordinate_standard_deviation,
    )
    changed_field = (
        f"{changed_field}.{coordinate_component_names[selected_coordinate_index]}"
    )
    if not math.isfinite(proposed_coordinate_value):
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "partitioned DNA base-frequency proposal produced one non-finite simplex coordinate"
            ),
        )
    current_unconstrained_values[selected_coordinate_index] = proposed_coordinate_value
    proposed_parameterization = resolve_dna_base_frequency_simplex_from_unconstrained(
        current_unconstrained_values
    )
    proposed_model_parameters = _rebuild_partition_model_parameters_with_group_value(
        prepared=prepared,
        target_name="base-frequencies",
        group_name=selected_group_name,
        target_value=proposed_parameterization.constrained_mapping(),
    )
    proposal_density = (
        _log_group_selection_density(group_count=prepared["group_count"])
        - math.log(len(current_unconstrained_values))
        + _gaussian_random_walk_density(
            coordinate_change=proposed_coordinate_value - current_coordinate_value,
            standard_deviation=validated_coordinate_standard_deviation,
        )
    )
    return build_metropolis_hastings_proposal(
        changed_fields=(changed_field,),
        log_forward_density=proposal_density,
        log_reverse_density=proposal_density,
        is_valid=True,
        proposed_tree=_copy_current_tree(current_state),
        proposed_model_parameters=proposed_model_parameters,
    )


def _propose_partition_exchangeability_move(
    *,
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    partition_models: Sequence[PartitionSubstitutionModelDefinition],
    unconstrained_coordinate_standard_deviation: float,
) -> MetropolisHastingsProposal:
    validated_coordinate_standard_deviation = _validate_positive_finite_float(
        value=unconstrained_coordinate_standard_deviation,
        field_name="unconstrained_coordinate_standard_deviation",
        owner_name="partitioned DNA exchangeability proposal",
    )
    selected_group_name, current_value, prepared = _select_partition_target_group(
        current_state=current_state,
        rng=rng,
        partition_models=partition_models,
        target_name="exchangeabilities",
    )
    changed_field = (
        f"vector_parameters.partition-parameter:exchangeabilities:{selected_group_name}"
    )
    if selected_group_name is None or current_value is None or prepared is None:
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "partitioned DNA exchangeability proposal requires at least one linked or unlinked exchangeability group"
            ),
        )
    parameterization = parameterize_dna_exchangeability_simplex(current_value)
    current_unconstrained_values = list(parameterization.unconstrained_values)
    coordinate_component_names = _simplex_coordinate_component_names(
        parameterization.component_names,
        reference_component_name=parameterization.reference_component_name,
    )
    selected_coordinate_index = rng.randrange(len(current_unconstrained_values))
    current_coordinate_value = current_unconstrained_values[selected_coordinate_index]
    proposed_coordinate_value = current_coordinate_value + rng.gauss(
        0.0,
        validated_coordinate_standard_deviation,
    )
    changed_field = (
        f"{changed_field}.{coordinate_component_names[selected_coordinate_index]}"
    )
    if not math.isfinite(proposed_coordinate_value):
        return build_metropolis_hastings_proposal(
            changed_fields=(changed_field,),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "partitioned DNA exchangeability proposal produced one non-finite simplex coordinate"
            ),
        )
    current_unconstrained_values[selected_coordinate_index] = proposed_coordinate_value
    proposed_parameterization = resolve_dna_exchangeability_simplex_from_unconstrained(
        current_unconstrained_values
    )
    proposed_model_parameters = _rebuild_partition_model_parameters_with_group_value(
        prepared=prepared,
        target_name="exchangeabilities",
        group_name=selected_group_name,
        target_value=proposed_parameterization.constrained_mapping(),
    )
    proposal_density = (
        _log_group_selection_density(group_count=prepared["group_count"])
        - math.log(len(current_unconstrained_values))
        + _gaussian_random_walk_density(
            coordinate_change=proposed_coordinate_value - current_coordinate_value,
            standard_deviation=validated_coordinate_standard_deviation,
        )
    )
    return build_metropolis_hastings_proposal(
        changed_fields=(changed_field,),
        log_forward_density=proposal_density,
        log_reverse_density=proposal_density,
        is_valid=True,
        proposed_tree=_copy_current_tree(current_state),
        proposed_model_parameters=proposed_model_parameters,
    )


def _select_partition_target_group(
    *,
    current_state: BayesianPhylogeneticState,
    rng: random.Random,
    partition_models: Sequence[PartitionSubstitutionModelDefinition],
    target_name: str,
) -> tuple[str | None, object | None, dict[str, object] | None]:
    prepared = _prepare_partition_parameter_surface(
        current_state=current_state,
        partition_models=partition_models,
    )
    grouped_partition_names = _group_partition_names_for_target(
        partition_models=partition_models,
        linkage_plan=prepared["linkage_plan"],
        target_name=target_name,
    )
    if not grouped_partition_names:
        return None, None, None
    group_names = sorted(grouped_partition_names)
    selected_group_name = group_names[rng.randrange(len(group_names))]
    representative_partition_name = grouped_partition_names[selected_group_name][0]
    current_value = _partition_state_value_for_target(
        partition_state=prepared["partition_state_by_name"][
            representative_partition_name
        ],
        target_name=target_name,
    )
    prepared["grouped_partition_names"] = grouped_partition_names
    prepared["group_count"] = len(group_names)
    return selected_group_name, current_value, prepared


def _prepare_partition_parameter_surface(
    *,
    current_state: BayesianPhylogeneticState,
    partition_models: Sequence[PartitionSubstitutionModelDefinition],
) -> dict[str, object]:
    partition_names = tuple(
        partition_model.partition_name for partition_model in partition_models
    )
    linkage_plan = resolve_partition_parameter_linkage_plan_from_model_parameters(
        model_parameters=current_state.model_parameters,
        partition_names=partition_names,
    )
    partition_parameter_states = (
        resolve_partition_parameter_states_from_model_parameters(
            model_parameters=current_state.model_parameters,
            partition_models=partition_models,
            linkage_plan=linkage_plan,
        )
    )
    partition_state_by_name = {
        partition_state.partition_name: partition_state
        for partition_state in partition_parameter_states
    }
    return {
        "partition_models": tuple(partition_models),
        "linkage_plan": linkage_plan,
        "partition_parameter_states": partition_parameter_states,
        "partition_state_by_name": partition_state_by_name,
        "stripped_model_parameters": strip_partition_model_parameter_state(
            current_state.model_parameters
        ),
        "current_state": current_state,
    }


def _group_partition_names_for_target(
    *,
    partition_models: Sequence[PartitionSubstitutionModelDefinition],
    linkage_plan,
    target_name: str,
) -> dict[str, list[str]]:
    target_groups = linkage_plan.groups_for_target(target_name)
    grouped_partition_names: dict[str, list[str]] = defaultdict(list)
    for partition_model in partition_models:
        if target_name not in partition_model.required_targets():
            continue
        grouped_partition_names[target_groups[partition_model.partition_name]].append(
            partition_model.partition_name
        )
    return grouped_partition_names


def _partition_state_value_for_target(
    *,
    partition_state: PartitionSubstitutionParameterState,
    target_name: str,
) -> object:
    if target_name == "kappa":
        return partition_state.kappa
    if target_name == "base-frequencies":
        return partition_state.base_frequencies
    if target_name == "exchangeabilities":
        return partition_state.exchangeabilities
    raise AssertionError(f"unsupported partition target {target_name}")


def _rebuild_partition_model_parameters_with_group_value(
    *,
    prepared: dict[str, object],
    target_name: str,
    group_name: str,
    target_value: object,
) -> BayesianModelParameterState:
    partition_state_payloads = {
        partition_state.partition_name: _partition_state_payload(partition_state)
        for partition_state in prepared["partition_parameter_states"]
    }
    for partition_name in prepared["grouped_partition_names"][group_name]:
        _assign_partition_target_value_to_payload(
            payload=partition_state_payloads[partition_name],
            target_name=target_name,
            target_value=target_value,
        )
    return build_partition_model_parameter_state(
        partition_models=prepared["partition_models"],
        linkage_plan=prepared["linkage_plan"],
        partition_parameter_states=tuple(
            PartitionSubstitutionParameterState(
                **partition_state_payloads[partition_model.partition_name]
            )
            for partition_model in prepared["partition_models"]
        ),
        preserved_categorical_parameters=prepared[
            "stripped_model_parameters"
        ].categorical_parameters,
        preserved_scalar_parameters=prepared[
            "stripped_model_parameters"
        ].scalar_parameters,
        preserved_vector_parameters=prepared[
            "stripped_model_parameters"
        ].vector_parameters,
    )


def _partition_state_payload(
    partition_state: PartitionSubstitutionParameterState,
) -> dict[str, object]:
    payload: dict[str, object] = {"partition_name": partition_state.partition_name}
    if partition_state.kappa is not None:
        payload["kappa"] = partition_state.kappa
    if partition_state.base_frequencies is not None:
        payload["base_frequencies"] = dict(partition_state.base_frequencies)
    if partition_state.exchangeabilities is not None:
        payload["exchangeabilities"] = dict(partition_state.exchangeabilities)
    if partition_state.gamma_alpha is not None:
        payload["gamma_alpha"] = partition_state.gamma_alpha
    if partition_state.invariant_proportion is not None:
        payload["invariant_proportion"] = partition_state.invariant_proportion
    return payload


def _partition_target_value_from_payload(
    *,
    payload: dict[str, object],
    target_name: str,
) -> object:
    if target_name == "kappa":
        return payload.get("kappa")
    if target_name == "base-frequencies":
        return dict(payload["base_frequencies"])
    if target_name == "exchangeabilities":
        return dict(payload["exchangeabilities"])
    raise AssertionError(f"unsupported partition target {target_name}")


def _assign_partition_target_value_to_payload(
    *,
    payload: dict[str, object],
    target_name: str,
    target_value: object,
) -> None:
    if target_name == "kappa":
        payload["kappa"] = float(target_value)
        return
    if target_name == "base-frequencies":
        payload["base_frequencies"] = dict(target_value)
        return
    if target_name == "exchangeabilities":
        payload["exchangeabilities"] = dict(target_value)
        return
    raise AssertionError(f"unsupported partition target {target_name}")


def _build_fixed_topology_partitioned_dna_posterior_rows(
    *,
    chain_report: MetropolisHastingsRunReport,
    model_definition: FixedTopologyPartitionedDnaModelDefinition,
    fixed_topology_id: str,
    partition_records_by_name: dict[str, list[AlignmentRecord]],
    observation_policy: str,
) -> list[FixedTopologyPartitionedDnaPosteriorRow]:
    posterior_rows: list[FixedTopologyPartitionedDnaPosteriorRow] = []
    for sample_index, state in enumerate(chain_report.sampled_states):
        linkage_plan, partition_parameter_states = (
            _require_fixed_topology_partitioned_dna_state_consistency(
                state=state,
                model_definition=model_definition,
                fixed_topology_id=fixed_topology_id,
            )
        )
        partition_rows = _build_partition_posterior_rows(
            state=state,
            model_definition=model_definition,
            partition_records_by_name=partition_records_by_name,
            partition_parameter_states=partition_parameter_states,
            linkage_plan=linkage_plan,
            observation_policy=observation_policy,
        )
        posterior_rows.append(
            FixedTopologyPartitionedDnaPosteriorRow(
                sample_index=sample_index,
                iteration_index=sample_index * chain_report.sample_every,
                topology_id=state.tree.topology_id,
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
                partition_rows=partition_rows,
            )
        )
    return posterior_rows


def _build_partition_posterior_rows(
    *,
    state: BayesianPhylogeneticState,
    model_definition: FixedTopologyPartitionedDnaModelDefinition,
    partition_records_by_name: dict[str, list[AlignmentRecord]],
    partition_parameter_states: Sequence[PartitionSubstitutionParameterState],
    linkage_plan,
    observation_policy: str,
) -> list[FixedTopologyPartitionedDnaPartitionRow]:
    partition_state_by_name = {
        partition_state.partition_name: partition_state
        for partition_state in partition_parameter_states
    }
    tree = state.tree.to_tree()
    partition_rows: list[FixedTopologyPartitionedDnaPartitionRow] = []
    for partition_model in model_definition.partition_models:
        partition_state = partition_state_by_name[partition_model.partition_name]
        partition_rows.append(
            FixedTopologyPartitionedDnaPartitionRow(
                partition_name=partition_model.partition_name,
                model_name=partition_model.model_name,
                log_likelihood=_evaluate_partition_log_likelihood(
                    tree=tree,
                    partition_model=partition_model,
                    partition_state=partition_state,
                    partition_records=partition_records_by_name[
                        partition_model.partition_name
                    ],
                    observation_policy=observation_policy,
                ),
                linkage_groups={
                    target_name: linkage_plan.groups_for_target(target_name)[
                        partition_model.partition_name
                    ]
                    for target_name in partition_model.required_targets()
                    if target_name in _SUPPORTED_PARTITION_PARAMETER_TARGETS
                },
                scalar_parameters=_partition_scalar_parameters(partition_state),
                vector_parameters=_partition_vector_parameters(partition_state),
            )
        )
    return partition_rows


def _partition_scalar_parameters(
    partition_state: PartitionSubstitutionParameterState,
) -> dict[str, float]:
    scalar_parameters: dict[str, float] = {}
    if partition_state.kappa is not None:
        scalar_parameters["kappa"] = partition_state.kappa
    return scalar_parameters


def _partition_vector_parameters(
    partition_state: PartitionSubstitutionParameterState,
) -> dict[str, dict[str, float]]:
    vector_parameters: dict[str, dict[str, float]] = {}
    if partition_state.base_frequencies is not None:
        vector_parameters["base-frequencies"] = dict(partition_state.base_frequencies)
    if partition_state.exchangeabilities is not None:
        vector_parameters["exchangeabilities"] = dict(partition_state.exchangeabilities)
    return vector_parameters


def _require_fixed_topology_partitioned_dna_state_consistency(
    *,
    state: BayesianPhylogeneticState,
    model_definition: FixedTopologyPartitionedDnaModelDefinition,
    fixed_topology_id: str | None,
) -> tuple[object, tuple[PartitionSubstitutionParameterState, ...]]:
    model_name = state.model_parameters.categorical_parameters.get("substitution-model")
    if model_name != _MODEL_NAME:
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior model requires every sampled state to preserve the partitioned substitution-model label",
            code="fixed_topology_partitioned_dna_state_model_label_invalid",
            details={
                "expected_model_name": _MODEL_NAME,
                "observed_model_name": model_name,
            },
        )
    if fixed_topology_id is not None and state.tree.topology_id != fixed_topology_id:
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior model requires topology to remain unchanged across sampled states",
            code="fixed_topology_partitioned_dna_state_topology_changed",
            details={
                "expected_topology_id": fixed_topology_id,
                "observed_topology_id": state.tree.topology_id,
            },
        )
    linkage_plan = resolve_partition_parameter_linkage_plan_from_model_parameters(
        model_parameters=state.model_parameters,
        partition_names=tuple(
            partition_model.partition_name
            for partition_model in model_definition.partition_models
        ),
    )
    partition_parameter_states = (
        resolve_partition_parameter_states_from_model_parameters(
            model_parameters=state.model_parameters,
            partition_models=model_definition.partition_models,
            linkage_plan=linkage_plan,
        )
    )
    return linkage_plan, partition_parameter_states


def _copy_current_tree(current_state: BayesianPhylogeneticState) -> PhyloTree:
    current_tree = current_state.tree.to_tree()
    current_tree.rooted = current_state.tree.rooted
    return current_tree


def _log_group_selection_density(*, group_count: int) -> float:
    return -math.log(group_count)


def _lognormal_scaling_density(
    *,
    current_value: float,
    proposed_value: float,
    log_scale_standard_deviation: float,
) -> float:
    log_scale_change = math.log(proposed_value / current_value)
    z_score = log_scale_change / log_scale_standard_deviation
    return (
        -math.log(proposed_value)
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


def _validate_locus_partitions(
    *,
    locus_partitions: Sequence[LocusPartition],
    expected_partition_names: tuple[str, ...],
) -> tuple[LocusPartition, ...]:
    validated_locus_partitions = tuple(locus_partitions)
    if not validated_locus_partitions:
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior model requires at least one locus partition",
            code="fixed_topology_partitioned_dna_locus_partitions_empty",
        )
    partition_by_name = {}
    duplicate_partition_names: list[str] = []
    for locus_partition in validated_locus_partitions:
        if not isinstance(locus_partition, LocusPartition):
            raise PhylogeneticsError(
                "fixed-topology partitioned DNA posterior model requires every locus partition to be one LocusPartition",
                code="fixed_topology_partitioned_dna_locus_partition_type_invalid",
            )
        if normalize_partition_data_type(locus_partition.data_type) not in {
            None,
            "DNA",
        }:
            raise PhylogeneticsError(
                "fixed-topology partitioned DNA posterior model supports DNA locus partitions only",
                code="fixed_topology_partitioned_dna_locus_partition_data_type_invalid",
                details={
                    "partition_name": locus_partition.name,
                    "data_type": locus_partition.data_type,
                },
            )
        if locus_partition.name in partition_by_name:
            duplicate_partition_names.append(locus_partition.name)
        partition_by_name[locus_partition.name] = locus_partition
    if duplicate_partition_names:
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior model requires unique locus partition names",
            code="fixed_topology_partitioned_dna_locus_partition_names_duplicated",
            details={
                "duplicate_partition_names": sorted(set(duplicate_partition_names))
            },
        )
    observed_partition_names = tuple(partition_by_name)
    if set(observed_partition_names) != set(expected_partition_names):
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior model requires locus partitions and partition priors to cover the exact same partition names",
            code="fixed_topology_partitioned_dna_partition_name_mismatch",
            details={
                "expected_partition_names": list(expected_partition_names),
                "observed_partition_names": sorted(observed_partition_names),
            },
        )
    return tuple(
        partition_by_name[partition_name] for partition_name in expected_partition_names
    )


def _validate_supported_partition_models(
    *,
    partition_models: Sequence[PartitionSubstitutionModelDefinition],
) -> None:
    unsupported_partition_names = [
        partition_model.partition_name
        for partition_model in partition_models
        if partition_model.gamma_enabled or partition_model.invariant_enabled
    ]
    if unsupported_partition_names:
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior model currently supports partition models without gamma-rate or invariant-site modifiers",
            code="fixed_topology_partitioned_dna_partition_model_unsupported",
            details={"partition_names": unsupported_partition_names},
        )


def _validate_supported_partition_prior_families(
    *,
    partition_prior_bundle: PartitionModelPriorBundle,
) -> None:
    prior_models = (
        ("kappa", partition_prior_bundle.substitution_prior_bundle.kappa_prior),
        (
            "exchangeabilities",
            partition_prior_bundle.substitution_prior_bundle.exchangeability_prior,
        ),
        (
            "base-frequencies",
            partition_prior_bundle.substitution_prior_bundle.base_frequency_prior,
        ),
        (
            "gamma-alpha",
            partition_prior_bundle.substitution_prior_bundle.gamma_alpha_prior,
        ),
        (
            "invariant-proportion",
            partition_prior_bundle.substitution_prior_bundle.invariant_proportion_prior,
        ),
    )
    unsupported_prior_targets = [
        target_name
        for target_name, prior_model in prior_models
        if prior_model is not None
        and target_name not in _SUPPORTED_PARTITION_PARAMETER_TARGETS
    ]
    if unsupported_prior_targets:
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior model currently supports priors only for sampled kappa, base-frequency, and exchangeability targets",
            code="fixed_topology_partitioned_dna_prior_targets_unsupported",
            details={"unsupported_prior_targets": sorted(unsupported_prior_targets)},
        )
    fixed_prior_targets = [
        target_name
        for target_name, prior_model in prior_models
        if prior_model is not None and getattr(prior_model, "family", None) == "fixed"
    ]
    if fixed_prior_targets:
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior model does not support fixed substitution priors because Metropolis-Hastings proposals would leave finite prior support",
            code="fixed_topology_partitioned_dna_prior_family_invalid",
            details={"fixed_prior_targets": sorted(fixed_prior_targets)},
        )


def _validate_initial_partition_parameter_states(
    *,
    partition_models: Sequence[PartitionSubstitutionModelDefinition],
    partition_prior_bundle: PartitionModelPriorBundle,
    initial_partition_parameter_states: Sequence[PartitionSubstitutionParameterState],
) -> tuple[PartitionSubstitutionParameterState, ...]:
    validated_initial_partition_parameter_states = tuple(
        initial_partition_parameter_states
    )
    build_partition_model_parameter_state(
        partition_models=partition_models,
        linkage_plan=partition_prior_bundle.linkage_plan,
        partition_parameter_states=validated_initial_partition_parameter_states,
        preserved_categorical_parameters={"substitution-model": _MODEL_NAME},
    )
    return validated_initial_partition_parameter_states


def _validate_parameter_move_activation(
    *,
    target_name: str,
    active_targets: set[str],
    move_weight: float,
    standard_deviation: float | None,
) -> None:
    if target_name in active_targets:
        if move_weight <= 0.0:
            raise PhylogeneticsError(
                "fixed-topology partitioned DNA proposal schedule requires every active parameter target to receive one positive move weight",
                code="fixed_topology_partitioned_dna_parameter_move_weight_invalid",
                details={"target_name": target_name},
            )
        if standard_deviation is None:
            raise PhylogeneticsError(
                "fixed-topology partitioned DNA proposal schedule requires one positive tuning scale for every active parameter target",
                code="fixed_topology_partitioned_dna_parameter_tuning_missing",
                details={"target_name": target_name},
            )
        return
    if move_weight != 0.0 or standard_deviation is not None:
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA proposal schedule received tuning for a parameter target that no partition model uses",
            code="fixed_topology_partitioned_dna_parameter_move_unused",
            details={"target_name": target_name},
        )


def _validate_positive_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    if not math.isfinite(value) or value <= 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be positive and finite",
            code="fixed_topology_partitioned_dna_positive_float_invalid",
            details={field_name: value},
        )
    return float(format(value, ".15g"))


def _validate_nonnegative_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    if not math.isfinite(value) or value < 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be non-negative and finite",
            code="fixed_topology_partitioned_dna_nonnegative_float_invalid",
            details={field_name: value},
        )
    return float(format(value, ".15g"))


def _validate_optional_positive_finite_float(
    *,
    value: float | None,
    field_name: str,
    owner_name: str,
) -> float | None:
    if value is None:
        return None
    return _validate_positive_finite_float(
        value=value,
        field_name=field_name,
        owner_name=owner_name,
    )


def _stabilize_positive_simplex_mapping(
    raw_values: Sequence[float],
) -> dict[str, float]:
    stabilized_values = [
        max(float(component_value), _MINIMUM_SIMPLEX_COMPONENT)
        for component_value in raw_values
    ]
    total = math.fsum(stabilized_values)
    return {
        component_name: float(format(component_value / total, ".15g"))
        for component_name, component_value in zip(
            ("A", "C", "G", "T"),
            stabilized_values,
            strict=True,
        )
    }


__all__ = [
    "FIXED_TOPOLOGY_PARTITIONED_DNA_SUBSTITUTION_MODELS",
    "FixedTopologyPartitionedDnaModelDefinition",
    "FixedTopologyPartitionedDnaPartitionRow",
    "FixedTopologyPartitionedDnaPosteriorRow",
    "FixedTopologyPartitionedDnaProposalSchedule",
    "FixedTopologyPartitionedDnaRunReport",
    "build_fixed_topology_partitioned_dna_model_definition",
    "build_fixed_topology_partitioned_dna_proposal_schedule",
    "run_fixed_topology_partitioned_dna_metropolis_hastings",
]
