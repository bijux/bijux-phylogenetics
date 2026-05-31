from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import math

from bijux_phylogenetics.bayesian.branch_length_priors import BranchLengthPriorModel
from bijux_phylogenetics.bayesian.partition_model_priors import (
    PARTITION_MODEL_PRIOR_TARGETS,
    PartitionModelPriorBundle,
    PartitionSubstitutionModelDefinition,
    PartitionSubstitutionParameterState,
)
from bijux_phylogenetics.bayesian.partition_model_state import (
    build_partition_model_parameter_state,
)
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    normalize_partition_data_type,
)
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
    chain_report: object
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
    validated_kappa_log_scale_standard_deviation = _validate_optional_positive_finite_float(
        value=kappa_log_scale_standard_deviation,
        field_name="kappa_log_scale_standard_deviation",
        owner_name="fixed-topology partitioned DNA proposal schedule",
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
    if validated_linkage_move_weight > 0.0 and not model_definition.linkage_eligible_targets:
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
        if normalize_partition_data_type(locus_partition.data_type) not in {None, "DNA"}:
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
            details={"duplicate_partition_names": sorted(set(duplicate_partition_names))},
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
    return tuple(partition_by_name[partition_name] for partition_name in expected_partition_names)


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
        if prior_model is not None and target_name not in _SUPPORTED_PARTITION_PARAMETER_TARGETS
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
    validated_initial_partition_parameter_states = tuple(initial_partition_parameter_states)
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
