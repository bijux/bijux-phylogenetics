from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math
import random
from statistics import mean, median

import numpy

from bijux_phylogenetics.ancestral.discrete.likelihood.rate_matrix import (
    build_transition_rate_rows,
    rate_matrix_from_log_parameters,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_allowed_transition_pairs,
    resolve_root_prior,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_parameters import (
    parameterize_discrete_trait_rate_rows,
)
from bijux_phylogenetics.bayesian.partition_model_state import (
    resolve_partition_parameter_linkage_plan_from_model_parameters,
    resolve_partition_parameter_states_from_model_parameters,
)
from bijux_phylogenetics.bayesian.state import BayesianPhylogeneticState
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
)
from bijux_phylogenetics.phylo.likelihood.dna import (
    DNA_STATE_ORDER,
    UNIFORM_DNA_ROOT_PRIOR,
)
from bijux_phylogenetics.phylo.likelihood.f81 import f81_transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.gtr import gtr_transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.hky85 import (
    hky85_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.jc69 import jc69_transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.k80 import k80_transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import PhylogeneticsError
from bijux_phylogenetics.simulation.contracts import (
    SimulatedContinuousTrait,
    SimulatedDiscreteTrait,
)

from .brownian_continuous_trait import BrownianContinuousTraitRunReport
from .discrete_trait_mk import DiscreteTraitMkRunReport
from .fixed_topology_dna import FixedTopologyDnaRunReport
from .fixed_topology_partitioned_dna import FixedTopologyPartitionedDnaRunReport
from .joint_topology_dna import JointTopologyDnaRunReport
from .ornstein_uhlenbeck_continuous_trait import (
    OrnsteinUhlenbeckContinuousTraitRunReport,
)

POSTERIOR_PREDICTIVE_SAMPLE_SELECTION_POLICIES = ("uniform-with-replacement",)

_DNA_STATE_INDEX = {state: index for index, state in enumerate(DNA_STATE_ORDER)}
_SEGREGATING_SITE_COUNT = "segregating-site-count"
_GC_FRACTION = "gc-fraction"
_STATE_ENTROPY = "state-entropy"
_MAJORITY_STATE_FREQUENCY = "majority-state-frequency"
_TIP_MEAN = "tip-mean"
_TIP_VARIANCE = "tip-variance"
_TIP_RANGE = "tip-range"
_DISCRETE_TRAIT_STATE_ORDERING = "unordered"


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveSimulationDefinition:
    """Validated posterior predictive simulation configuration."""

    replicate_count: int
    sample_selection_policy: str
    seed: int


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveObservedStatisticRow:
    """One observed-data statistic retained beside posterior predictive replicates."""

    statistic_name: str
    value: float


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveReplicateStatisticRow:
    """One replicate statistic generated from one sampled posterior state."""

    statistic_name: str
    replicate_index: int
    posterior_sample_index: int
    posterior_iteration_index: int
    value: float


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveStatisticSummaryRow:
    """One observed-vs-replicate summary for one posterior predictive statistic."""

    statistic_name: str
    observed_value: float
    replicate_count: int
    replicate_mean: float
    replicate_standard_deviation: float
    replicate_minimum: float
    replicate_median: float
    replicate_maximum: float


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveAlignmentReplicate:
    """One posterior predictive alignment replicate."""

    replicate_index: int
    posterior_sample_index: int
    posterior_iteration_index: int
    records: list[AlignmentRecord]


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveDiscreteTraitReplicate:
    """One posterior predictive discrete-trait replicate."""

    replicate_index: int
    posterior_sample_index: int
    posterior_iteration_index: int
    traits: list[SimulatedDiscreteTrait]


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveContinuousTraitReplicate:
    """One posterior predictive continuous-trait replicate."""

    replicate_index: int
    posterior_sample_index: int
    posterior_iteration_index: int
    traits: list[SimulatedContinuousTrait]


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveAlignmentSimulationReport:
    """Posterior predictive alignment simulation output for one Bayesian run."""

    definition: PosteriorPredictiveSimulationDefinition
    model_name: str
    taxon_count: int
    alignment_length: int
    observed_statistic_rows: list[PosteriorPredictiveObservedStatisticRow]
    replicate_statistic_rows: list[PosteriorPredictiveReplicateStatisticRow]
    statistic_summary_rows: list[PosteriorPredictiveStatisticSummaryRow]
    replicates: list[PosteriorPredictiveAlignmentReplicate]


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveDiscreteTraitSimulationReport:
    """Posterior predictive discrete-trait simulation output for one Bayesian run."""

    definition: PosteriorPredictiveSimulationDefinition
    model_name: str
    taxon_count: int
    state_order: list[str]
    observed_statistic_rows: list[PosteriorPredictiveObservedStatisticRow]
    replicate_statistic_rows: list[PosteriorPredictiveReplicateStatisticRow]
    statistic_summary_rows: list[PosteriorPredictiveStatisticSummaryRow]
    replicates: list[PosteriorPredictiveDiscreteTraitReplicate]


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveContinuousTraitSimulationReport:
    """Posterior predictive continuous-trait simulation output for one Bayesian run."""

    definition: PosteriorPredictiveSimulationDefinition
    model_name: str
    taxon_count: int
    observed_statistic_rows: list[PosteriorPredictiveObservedStatisticRow]
    replicate_statistic_rows: list[PosteriorPredictiveReplicateStatisticRow]
    statistic_summary_rows: list[PosteriorPredictiveStatisticSummaryRow]
    replicates: list[PosteriorPredictiveContinuousTraitReplicate]


def build_posterior_predictive_simulation_definition(
    *,
    replicate_count: int,
    sample_selection_policy: str = "uniform-with-replacement",
    seed: int = 0,
) -> PosteriorPredictiveSimulationDefinition:
    """Build one validated posterior predictive simulation configuration."""
    validated_policy = sample_selection_policy.strip()
    if validated_policy not in POSTERIOR_PREDICTIVE_SAMPLE_SELECTION_POLICIES:
        raise PhylogeneticsError(
            "posterior predictive simulation definition requires one supported posterior-sample selection policy",
            code="posterior_predictive_simulation_policy_invalid",
            details={
                "sample_selection_policy": sample_selection_policy,
                "supported_policies": list(
                    POSTERIOR_PREDICTIVE_SAMPLE_SELECTION_POLICIES
                ),
            },
        )
    if not isinstance(replicate_count, int):
        raise PhylogeneticsError(
            "posterior predictive simulation definition requires replicate_count to be one integer",
            code="posterior_predictive_simulation_replicate_count_type_invalid",
        )
    if replicate_count <= 0:
        raise PhylogeneticsError(
            "posterior predictive simulation definition requires replicate_count to be positive",
            code="posterior_predictive_simulation_replicate_count_invalid",
            details={"replicate_count": replicate_count},
        )
    if not isinstance(seed, int):
        raise PhylogeneticsError(
            "posterior predictive simulation definition requires seed to be one integer",
            code="posterior_predictive_simulation_seed_type_invalid",
        )
    return PosteriorPredictiveSimulationDefinition(
        replicate_count=replicate_count,
        sample_selection_policy=validated_policy,
        seed=seed,
    )


def simulate_fixed_topology_dna_posterior_predictive(
    *,
    run_report: FixedTopologyDnaRunReport,
    records: Sequence[AlignmentRecord],
    definition: PosteriorPredictiveSimulationDefinition,
) -> PosteriorPredictiveAlignmentSimulationReport:
    """Simulate DNA alignments from sampled fixed-topology posterior states."""
    if not isinstance(run_report, FixedTopologyDnaRunReport):
        raise PhylogeneticsError(
            "fixed-topology DNA posterior predictive simulation requires one FixedTopologyDnaRunReport",
            code="posterior_predictive_fixed_topology_dna_run_report_type_invalid",
        )
    _require_reject_only_observation_policy(
        observation_policy=run_report.observation_policy,
        owner_name="fixed-topology DNA posterior predictive simulation",
    )
    validated_definition = _require_posterior_predictive_definition(
        definition,
        owner_name="fixed-topology DNA posterior predictive simulation",
    )
    validated_records = _validate_alignment_records(
        records=records,
        owner_name="fixed-topology DNA posterior predictive simulation",
    )
    sample_tree = run_report.chain_report.sampled_states[0].tree.to_tree()
    record_ids = _validate_alignment_taxa_against_tree(
        records=validated_records,
        tree=sample_tree,
        owner_name="fixed-topology DNA posterior predictive simulation",
    )
    observed_statistic_rows = _build_alignment_observed_statistic_rows(
        validated_records
    )
    rng = random.Random(validated_definition.seed)  # nosec B311
    replicates: list[PosteriorPredictiveAlignmentReplicate] = []
    replicate_statistic_rows: list[PosteriorPredictiveReplicateStatisticRow] = []
    for replicate_index, sample_index in enumerate(
        _choose_posterior_sample_indices(
            sample_count=len(run_report.chain_report.sampled_states),
            definition=validated_definition,
            rng=rng,
        )
    ):
        state = run_report.chain_report.sampled_states[sample_index]
        replicate_records = _simulate_dna_alignment_replicate(
            tree=state.tree.to_tree(),
            record_ids=record_ids,
            alignment_length=len(validated_records[0].sequence),
            model_name=run_report.model_definition.substitution_model_name,
            scalar_parameters=state.model_parameters.scalar_parameters,
            vector_parameters=state.model_parameters.vector_parameters,
            rng=rng,
        )
        posterior_iteration_index = sample_index * run_report.chain_report.sample_every
        replicates.append(
            PosteriorPredictiveAlignmentReplicate(
                replicate_index=replicate_index,
                posterior_sample_index=sample_index,
                posterior_iteration_index=posterior_iteration_index,
                records=replicate_records,
            )
        )
        replicate_statistic_rows.extend(
            _build_alignment_replicate_statistic_rows(
                records=replicate_records,
                replicate_index=replicate_index,
                posterior_sample_index=sample_index,
                posterior_iteration_index=posterior_iteration_index,
            )
        )
    return PosteriorPredictiveAlignmentSimulationReport(
        definition=validated_definition,
        model_name=run_report.model_definition.substitution_model_name,
        taxon_count=len(record_ids),
        alignment_length=len(validated_records[0].sequence),
        observed_statistic_rows=observed_statistic_rows,
        replicate_statistic_rows=replicate_statistic_rows,
        statistic_summary_rows=_build_posterior_predictive_statistic_summary_rows(
            observed_statistic_rows=observed_statistic_rows,
            replicate_statistic_rows=replicate_statistic_rows,
        ),
        replicates=replicates,
    )


def simulate_joint_topology_dna_posterior_predictive(
    *,
    run_report: JointTopologyDnaRunReport,
    records: Sequence[AlignmentRecord],
    definition: PosteriorPredictiveSimulationDefinition,
) -> PosteriorPredictiveAlignmentSimulationReport:
    """Simulate DNA alignments from sampled joint-topology posterior states."""
    if not isinstance(run_report, JointTopologyDnaRunReport):
        raise PhylogeneticsError(
            "joint topology DNA posterior predictive simulation requires one JointTopologyDnaRunReport",
            code="posterior_predictive_joint_topology_dna_run_report_type_invalid",
        )
    _require_reject_only_observation_policy(
        observation_policy=run_report.observation_policy,
        owner_name="joint topology DNA posterior predictive simulation",
    )
    validated_definition = _require_posterior_predictive_definition(
        definition,
        owner_name="joint topology DNA posterior predictive simulation",
    )
    validated_records = _validate_alignment_records(
        records=records,
        owner_name="joint topology DNA posterior predictive simulation",
    )
    sample_tree = run_report.chain_report.sampled_states[0].tree.to_tree()
    record_ids = _validate_alignment_taxa_against_tree(
        records=validated_records,
        tree=sample_tree,
        owner_name="joint topology DNA posterior predictive simulation",
    )
    observed_statistic_rows = _build_alignment_observed_statistic_rows(
        validated_records
    )
    rng = random.Random(validated_definition.seed)  # nosec B311
    replicates: list[PosteriorPredictiveAlignmentReplicate] = []
    replicate_statistic_rows: list[PosteriorPredictiveReplicateStatisticRow] = []
    for replicate_index, sample_index in enumerate(
        _choose_posterior_sample_indices(
            sample_count=len(run_report.chain_report.sampled_states),
            definition=validated_definition,
            rng=rng,
        )
    ):
        state = run_report.chain_report.sampled_states[sample_index]
        replicate_records = _simulate_dna_alignment_replicate(
            tree=state.tree.to_tree(),
            record_ids=record_ids,
            alignment_length=len(validated_records[0].sequence),
            model_name=run_report.model_definition.substitution_model_name,
            scalar_parameters=state.model_parameters.scalar_parameters,
            vector_parameters=state.model_parameters.vector_parameters,
            rng=rng,
        )
        posterior_iteration_index = sample_index * run_report.chain_report.sample_every
        replicates.append(
            PosteriorPredictiveAlignmentReplicate(
                replicate_index=replicate_index,
                posterior_sample_index=sample_index,
                posterior_iteration_index=posterior_iteration_index,
                records=replicate_records,
            )
        )
        replicate_statistic_rows.extend(
            _build_alignment_replicate_statistic_rows(
                records=replicate_records,
                replicate_index=replicate_index,
                posterior_sample_index=sample_index,
                posterior_iteration_index=posterior_iteration_index,
            )
        )
    return PosteriorPredictiveAlignmentSimulationReport(
        definition=validated_definition,
        model_name=run_report.model_definition.substitution_model_name,
        taxon_count=len(record_ids),
        alignment_length=len(validated_records[0].sequence),
        observed_statistic_rows=observed_statistic_rows,
        replicate_statistic_rows=replicate_statistic_rows,
        statistic_summary_rows=_build_posterior_predictive_statistic_summary_rows(
            observed_statistic_rows=observed_statistic_rows,
            replicate_statistic_rows=replicate_statistic_rows,
        ),
        replicates=replicates,
    )


def simulate_fixed_topology_partitioned_dna_posterior_predictive(
    *,
    run_report: FixedTopologyPartitionedDnaRunReport,
    records: Sequence[AlignmentRecord],
    definition: PosteriorPredictiveSimulationDefinition,
) -> PosteriorPredictiveAlignmentSimulationReport:
    """Simulate partitioned DNA alignments from sampled posterior states."""
    if not isinstance(run_report, FixedTopologyPartitionedDnaRunReport):
        raise PhylogeneticsError(
            "fixed-topology partitioned DNA posterior predictive simulation requires one FixedTopologyPartitionedDnaRunReport",
            code="posterior_predictive_partitioned_dna_run_report_type_invalid",
        )
    _require_reject_only_observation_policy(
        observation_policy=run_report.observation_policy,
        owner_name="fixed-topology partitioned DNA posterior predictive simulation",
    )
    validated_definition = _require_posterior_predictive_definition(
        definition,
        owner_name="fixed-topology partitioned DNA posterior predictive simulation",
    )
    validated_records = _validate_alignment_records(
        records=records,
        owner_name="fixed-topology partitioned DNA posterior predictive simulation",
    )
    sample_tree = run_report.chain_report.sampled_states[0].tree.to_tree()
    record_ids = _validate_alignment_taxa_against_tree(
        records=validated_records,
        tree=sample_tree,
        owner_name="fixed-topology partitioned DNA posterior predictive simulation",
    )
    observed_statistic_rows = _build_alignment_observed_statistic_rows(
        validated_records
    )
    partition_length_by_name = {
        partition.name: partition.total_sites
        for partition in run_report.model_definition.locus_partitions
    }
    rng = random.Random(validated_definition.seed)  # nosec B311
    replicates: list[PosteriorPredictiveAlignmentReplicate] = []
    replicate_statistic_rows: list[PosteriorPredictiveReplicateStatisticRow] = []
    for replicate_index, sample_index in enumerate(
        _choose_posterior_sample_indices(
            sample_count=len(run_report.chain_report.sampled_states),
            definition=validated_definition,
            rng=rng,
        )
    ):
        state = run_report.chain_report.sampled_states[sample_index]
        tree = state.tree.to_tree()
        linkage_plan = resolve_partition_parameter_linkage_plan_from_model_parameters(
            model_parameters=state.model_parameters,
            partition_names=tuple(
                partition_model.partition_name
                for partition_model in run_report.model_definition.partition_models
            ),
        )
        partition_parameter_states = (
            resolve_partition_parameter_states_from_model_parameters(
                model_parameters=state.model_parameters,
                partition_models=run_report.model_definition.partition_models,
                linkage_plan=linkage_plan,
            )
        )
        records_by_partition_name = {
            partition_state.partition_name: _simulate_dna_alignment_replicate(
                tree=tree,
                record_ids=record_ids,
                alignment_length=partition_length_by_name[
                    partition_state.partition_name
                ],
                model_name=next(
                    partition_model.model_name
                    for partition_model in run_report.model_definition.partition_models
                    if partition_model.partition_name == partition_state.partition_name
                ),
                scalar_parameters=_partition_scalar_parameters(partition_state),
                vector_parameters=_partition_vector_parameters(partition_state),
                rng=rng,
            )
            for partition_state in partition_parameter_states
        }
        replicate_records = _merge_partition_alignment_records(
            record_ids=record_ids,
            partition_order=run_report.model_definition.locus_partitions,
            records_by_partition_name=records_by_partition_name,
        )
        posterior_iteration_index = sample_index * run_report.chain_report.sample_every
        replicates.append(
            PosteriorPredictiveAlignmentReplicate(
                replicate_index=replicate_index,
                posterior_sample_index=sample_index,
                posterior_iteration_index=posterior_iteration_index,
                records=replicate_records,
            )
        )
        replicate_statistic_rows.extend(
            _build_alignment_replicate_statistic_rows(
                records=replicate_records,
                replicate_index=replicate_index,
                posterior_sample_index=sample_index,
                posterior_iteration_index=posterior_iteration_index,
            )
        )
    return PosteriorPredictiveAlignmentSimulationReport(
        definition=validated_definition,
        model_name="partitioned-dna",
        taxon_count=len(record_ids),
        alignment_length=len(validated_records[0].sequence),
        observed_statistic_rows=observed_statistic_rows,
        replicate_statistic_rows=replicate_statistic_rows,
        statistic_summary_rows=_build_posterior_predictive_statistic_summary_rows(
            observed_statistic_rows=observed_statistic_rows,
            replicate_statistic_rows=replicate_statistic_rows,
        ),
        replicates=replicates,
    )


def simulate_discrete_trait_mk_posterior_predictive(
    *,
    run_report: DiscreteTraitMkRunReport,
    tip_states: Mapping[str, str],
    definition: PosteriorPredictiveSimulationDefinition,
) -> PosteriorPredictiveDiscreteTraitSimulationReport:
    """Simulate discrete traits from sampled Mk posterior states."""
    if not isinstance(run_report, DiscreteTraitMkRunReport):
        raise PhylogeneticsError(
            "discrete-trait Mk posterior predictive simulation requires one DiscreteTraitMkRunReport",
            code="posterior_predictive_discrete_trait_run_report_type_invalid",
        )
    validated_definition = _require_posterior_predictive_definition(
        definition,
        owner_name="discrete-trait Mk posterior predictive simulation",
    )
    normalized_tip_states = _validate_tip_states(
        tip_states=tip_states,
        expected_taxa=run_report.taxa,
        owner_name="discrete-trait Mk posterior predictive simulation",
    )
    observed_statistic_rows = _build_discrete_trait_observed_statistic_rows(
        tip_states=normalized_tip_states,
        state_order=run_report.state_order,
    )
    observed_state_counts = _state_counts(normalized_tip_states.values())
    allowed_transition_pairs = resolve_allowed_transition_pairs(
        run_report.state_order,
        model=run_report.model_definition.transition_model_name,
        state_ordering=_DISCRETE_TRAIT_STATE_ORDERING,
        allowed_transition_pairs=None,
    )
    rng = random.Random(validated_definition.seed)  # nosec B311
    replicates: list[PosteriorPredictiveDiscreteTraitReplicate] = []
    replicate_statistic_rows: list[PosteriorPredictiveReplicateStatisticRow] = []
    for replicate_index, sample_index in enumerate(
        _choose_posterior_sample_indices(
            sample_count=len(run_report.chain_report.sampled_states),
            definition=validated_definition,
            rng=rng,
        )
    ):
        state = run_report.chain_report.sampled_states[sample_index]
        tree = state.tree.to_tree()
        rate_matrix = _build_discrete_trait_rate_matrix(
            state=state,
            state_order=run_report.state_order,
            model_name=run_report.model_definition.transition_model_name,
            allowed_transition_pairs=allowed_transition_pairs,
        )
        root_prior = resolve_root_prior(
            run_report.state_order,
            state_counts=observed_state_counts,
            mode=run_report.model_definition.root_prior_mode,
            fixed_root_state=run_report.model_definition.fixed_root_state,
        )
        simulated_traits = _simulate_discrete_trait_replicate(
            tree=tree,
            taxa=run_report.taxa,
            state_order=run_report.state_order,
            root_prior=root_prior,
            rate_matrix=rate_matrix,
            rng=rng,
        )
        posterior_iteration_index = sample_index * run_report.chain_report.sample_every
        replicates.append(
            PosteriorPredictiveDiscreteTraitReplicate(
                replicate_index=replicate_index,
                posterior_sample_index=sample_index,
                posterior_iteration_index=posterior_iteration_index,
                traits=simulated_traits,
            )
        )
        replicate_statistic_rows.extend(
            _build_discrete_trait_replicate_statistic_rows(
                traits=simulated_traits,
                state_order=run_report.state_order,
                replicate_index=replicate_index,
                posterior_sample_index=sample_index,
                posterior_iteration_index=posterior_iteration_index,
            )
        )
    return PosteriorPredictiveDiscreteTraitSimulationReport(
        definition=validated_definition,
        model_name=run_report.model_definition.transition_model_name,
        taxon_count=len(run_report.taxa),
        state_order=list(run_report.state_order),
        observed_statistic_rows=observed_statistic_rows,
        replicate_statistic_rows=replicate_statistic_rows,
        statistic_summary_rows=_build_posterior_predictive_statistic_summary_rows(
            observed_statistic_rows=observed_statistic_rows,
            replicate_statistic_rows=replicate_statistic_rows,
        ),
        replicates=replicates,
    )


def simulate_brownian_continuous_trait_posterior_predictive(
    *,
    run_report: BrownianContinuousTraitRunReport,
    definition: PosteriorPredictiveSimulationDefinition,
) -> PosteriorPredictiveContinuousTraitSimulationReport:
    """Simulate continuous traits from sampled Brownian posterior states."""
    if not isinstance(run_report, BrownianContinuousTraitRunReport):
        raise PhylogeneticsError(
            "Brownian posterior predictive simulation requires one BrownianContinuousTraitRunReport",
            code="posterior_predictive_brownian_run_report_type_invalid",
        )
    validated_definition = _require_posterior_predictive_definition(
        definition,
        owner_name="Brownian posterior predictive simulation",
    )
    observed_statistic_rows = _build_continuous_trait_observed_statistic_rows(
        tip_values=run_report.tip_values,
        taxa=run_report.taxa,
    )
    rng = random.Random(validated_definition.seed)  # nosec B311
    replicates: list[PosteriorPredictiveContinuousTraitReplicate] = []
    replicate_statistic_rows: list[PosteriorPredictiveReplicateStatisticRow] = []
    for replicate_index, sample_index in enumerate(
        _choose_posterior_sample_indices(
            sample_count=len(run_report.chain_report.sampled_states),
            definition=validated_definition,
            rng=rng,
        )
    ):
        state = run_report.chain_report.sampled_states[sample_index]
        simulated_traits = _simulate_brownian_continuous_trait_replicate(
            tree=state.tree.to_tree(),
            taxa=run_report.taxa,
            root_state=state.model_parameters.scalar_parameters["root-state"],
            sigma_squared=state.model_parameters.scalar_parameters["sigma-squared"],
            rng=rng,
        )
        posterior_iteration_index = sample_index * run_report.chain_report.sample_every
        replicates.append(
            PosteriorPredictiveContinuousTraitReplicate(
                replicate_index=replicate_index,
                posterior_sample_index=sample_index,
                posterior_iteration_index=posterior_iteration_index,
                traits=simulated_traits,
            )
        )
        replicate_statistic_rows.extend(
            _build_continuous_trait_replicate_statistic_rows(
                traits=simulated_traits,
                replicate_index=replicate_index,
                posterior_sample_index=sample_index,
                posterior_iteration_index=posterior_iteration_index,
            )
        )
    return PosteriorPredictiveContinuousTraitSimulationReport(
        definition=validated_definition,
        model_name="brownian",
        taxon_count=len(run_report.taxa),
        observed_statistic_rows=observed_statistic_rows,
        replicate_statistic_rows=replicate_statistic_rows,
        statistic_summary_rows=_build_posterior_predictive_statistic_summary_rows(
            observed_statistic_rows=observed_statistic_rows,
            replicate_statistic_rows=replicate_statistic_rows,
        ),
        replicates=replicates,
    )


def simulate_ornstein_uhlenbeck_continuous_trait_posterior_predictive(
    *,
    run_report: OrnsteinUhlenbeckContinuousTraitRunReport,
    definition: PosteriorPredictiveSimulationDefinition,
) -> PosteriorPredictiveContinuousTraitSimulationReport:
    """Simulate continuous traits from sampled OU posterior states."""
    if not isinstance(run_report, OrnsteinUhlenbeckContinuousTraitRunReport):
        raise PhylogeneticsError(
            "OU posterior predictive simulation requires one OrnsteinUhlenbeckContinuousTraitRunReport",
            code="posterior_predictive_ou_run_report_type_invalid",
        )
    validated_definition = _require_posterior_predictive_definition(
        definition,
        owner_name="OU posterior predictive simulation",
    )
    observed_statistic_rows = _build_continuous_trait_observed_statistic_rows(
        tip_values=run_report.tip_values,
        taxa=run_report.taxa,
    )
    rng = random.Random(validated_definition.seed)  # nosec B311
    replicates: list[PosteriorPredictiveContinuousTraitReplicate] = []
    replicate_statistic_rows: list[PosteriorPredictiveReplicateStatisticRow] = []
    for replicate_index, sample_index in enumerate(
        _choose_posterior_sample_indices(
            sample_count=len(run_report.chain_report.sampled_states),
            definition=validated_definition,
            rng=rng,
        )
    ):
        state = run_report.chain_report.sampled_states[sample_index]
        simulated_traits = _simulate_ou_continuous_trait_replicate(
            tree=state.tree.to_tree(),
            taxa=run_report.taxa,
            alpha=state.model_parameters.scalar_parameters["alpha"],
            optimum=state.model_parameters.scalar_parameters["optimum"],
            sigma_squared=state.model_parameters.scalar_parameters["sigma-squared"],
            rng=rng,
        )
        posterior_iteration_index = sample_index * run_report.chain_report.sample_every
        replicates.append(
            PosteriorPredictiveContinuousTraitReplicate(
                replicate_index=replicate_index,
                posterior_sample_index=sample_index,
                posterior_iteration_index=posterior_iteration_index,
                traits=simulated_traits,
            )
        )
        replicate_statistic_rows.extend(
            _build_continuous_trait_replicate_statistic_rows(
                traits=simulated_traits,
                replicate_index=replicate_index,
                posterior_sample_index=sample_index,
                posterior_iteration_index=posterior_iteration_index,
            )
        )
    return PosteriorPredictiveContinuousTraitSimulationReport(
        definition=validated_definition,
        model_name="ornstein-uhlenbeck",
        taxon_count=len(run_report.taxa),
        observed_statistic_rows=observed_statistic_rows,
        replicate_statistic_rows=replicate_statistic_rows,
        statistic_summary_rows=_build_posterior_predictive_statistic_summary_rows(
            observed_statistic_rows=observed_statistic_rows,
            replicate_statistic_rows=replicate_statistic_rows,
        ),
        replicates=replicates,
    )


def _build_posterior_predictive_statistic_summary_rows(
    *,
    observed_statistic_rows: list[PosteriorPredictiveObservedStatisticRow],
    replicate_statistic_rows: list[PosteriorPredictiveReplicateStatisticRow],
) -> list[PosteriorPredictiveStatisticSummaryRow]:
    replicate_values_by_statistic: dict[str, list[float]] = {}
    for row in replicate_statistic_rows:
        replicate_values_by_statistic.setdefault(row.statistic_name, []).append(
            row.value
        )
    summary_rows: list[PosteriorPredictiveStatisticSummaryRow] = []
    for observed_row in observed_statistic_rows:
        replicate_values = replicate_values_by_statistic.get(
            observed_row.statistic_name
        )
        if not replicate_values:
            raise PhylogeneticsError(
                "posterior predictive statistic summaries require at least one replicate value for every observed statistic",
                code="posterior_predictive_statistic_replicates_missing",
                details={"statistic_name": observed_row.statistic_name},
            )
        summary_rows.append(
            PosteriorPredictiveStatisticSummaryRow(
                statistic_name=observed_row.statistic_name,
                observed_value=observed_row.value,
                replicate_count=len(replicate_values),
                replicate_mean=_round_float(mean(replicate_values)),
                replicate_standard_deviation=_round_float(
                    _sample_standard_deviation(replicate_values)
                ),
                replicate_minimum=_round_float(min(replicate_values)),
                replicate_median=_round_float(median(replicate_values)),
                replicate_maximum=_round_float(max(replicate_values)),
            )
        )
    return summary_rows


def _require_posterior_predictive_definition(
    definition: PosteriorPredictiveSimulationDefinition,
    *,
    owner_name: str,
) -> PosteriorPredictiveSimulationDefinition:
    if not isinstance(definition, PosteriorPredictiveSimulationDefinition):
        raise PhylogeneticsError(
            f"{owner_name} requires one PosteriorPredictiveSimulationDefinition",
            code="posterior_predictive_definition_type_invalid",
        )
    return definition


def _require_reject_only_observation_policy(
    *,
    observation_policy: str,
    owner_name: str,
) -> None:
    if observation_policy != "reject":
        raise PhylogeneticsError(
            f"{owner_name} currently supports observation_policy='reject' only",
            code="posterior_predictive_observation_policy_unsupported",
            details={"observation_policy": observation_policy},
        )


def _choose_posterior_sample_indices(
    *,
    sample_count: int,
    definition: PosteriorPredictiveSimulationDefinition,
    rng: random.Random,
) -> list[int]:
    if sample_count <= 0:
        raise PhylogeneticsError(
            "posterior predictive simulation requires at least one posterior sample",
            code="posterior_predictive_sample_count_invalid",
            details={"sample_count": sample_count},
        )
    if definition.sample_selection_policy != "uniform-with-replacement":
        raise PhylogeneticsError(
            "posterior predictive simulation encountered one unsupported posterior-sample selection policy",
            code="posterior_predictive_simulation_policy_unsupported",
            details={"sample_selection_policy": definition.sample_selection_policy},
        )
    return [rng.randrange(sample_count) for _ in range(definition.replicate_count)]


def _validate_alignment_records(
    *,
    records: Sequence[AlignmentRecord],
    owner_name: str,
) -> list[AlignmentRecord]:
    validated_records = list(records)
    if not validated_records:
        raise PhylogeneticsError(
            f"{owner_name} requires at least one alignment record",
            code="posterior_predictive_alignment_records_empty",
        )
    if any(not isinstance(record, AlignmentRecord) for record in validated_records):
        raise PhylogeneticsError(
            f"{owner_name} requires every record to be one AlignmentRecord",
            code="posterior_predictive_alignment_record_type_invalid",
        )
    identifiers = [record.identifier for record in validated_records]
    if len(set(identifiers)) != len(identifiers):
        raise PhylogeneticsError(
            f"{owner_name} requires unique alignment record identifiers",
            code="posterior_predictive_alignment_record_identifier_duplicated",
        )
    sequence_lengths = {len(record.sequence) for record in validated_records}
    if len(sequence_lengths) != 1:
        raise PhylogeneticsError(
            f"{owner_name} requires one aligned sequence matrix with equal sequence lengths",
            code="posterior_predictive_alignment_length_invalid",
        )
    return validated_records


def _validate_alignment_taxa_against_tree(
    *,
    records: Sequence[AlignmentRecord],
    tree: PhyloTree,
    owner_name: str,
) -> list[str]:
    record_ids = [record.identifier for record in records]
    if sorted(record_ids) != sorted(tree.tip_names):
        raise PhylogeneticsError(
            f"{owner_name} requires alignment identifiers to match the sampled tree tip set exactly",
            code="posterior_predictive_alignment_taxa_mismatch",
            details={
                "record_identifiers": sorted(record_ids),
                "tree_tip_names": sorted(tree.tip_names),
            },
        )
    return record_ids


def _validate_tip_states(
    *,
    tip_states: Mapping[str, str],
    expected_taxa: Sequence[str],
    owner_name: str,
) -> dict[str, str]:
    validated_tip_states = {
        str(taxon): state.strip() for taxon, state in tip_states.items()
    }
    if sorted(validated_tip_states) != sorted(expected_taxa):
        raise PhylogeneticsError(
            f"{owner_name} requires tip states to match the posterior tip set exactly",
            code="posterior_predictive_tip_states_taxa_mismatch",
            details={
                "provided_taxa": sorted(validated_tip_states),
                "expected_taxa": sorted(expected_taxa),
            },
        )
    if any(not state for state in validated_tip_states.values()):
        raise PhylogeneticsError(
            f"{owner_name} requires every tip state to be nonblank",
            code="posterior_predictive_tip_state_blank",
        )
    return validated_tip_states


def _build_alignment_observed_statistic_rows(
    records: Sequence[AlignmentRecord],
) -> list[PosteriorPredictiveObservedStatisticRow]:
    return [
        PosteriorPredictiveObservedStatisticRow(
            statistic_name=_SEGREGATING_SITE_COUNT,
            value=_round_float(_alignment_segregating_site_count(records)),
        ),
        PosteriorPredictiveObservedStatisticRow(
            statistic_name=_GC_FRACTION,
            value=_round_float(_alignment_gc_fraction(records)),
        ),
    ]


def _build_alignment_replicate_statistic_rows(
    *,
    records: Sequence[AlignmentRecord],
    replicate_index: int,
    posterior_sample_index: int,
    posterior_iteration_index: int,
) -> list[PosteriorPredictiveReplicateStatisticRow]:
    return [
        PosteriorPredictiveReplicateStatisticRow(
            statistic_name=_SEGREGATING_SITE_COUNT,
            replicate_index=replicate_index,
            posterior_sample_index=posterior_sample_index,
            posterior_iteration_index=posterior_iteration_index,
            value=_round_float(_alignment_segregating_site_count(records)),
        ),
        PosteriorPredictiveReplicateStatisticRow(
            statistic_name=_GC_FRACTION,
            replicate_index=replicate_index,
            posterior_sample_index=posterior_sample_index,
            posterior_iteration_index=posterior_iteration_index,
            value=_round_float(_alignment_gc_fraction(records)),
        ),
    ]


def _build_discrete_trait_observed_statistic_rows(
    *,
    tip_states: Mapping[str, str],
    state_order: Sequence[str],
) -> list[PosteriorPredictiveObservedStatisticRow]:
    counts = _state_counts(tip_states.values())
    return [
        PosteriorPredictiveObservedStatisticRow(
            statistic_name=_STATE_ENTROPY,
            value=_round_float(
                _categorical_state_entropy(
                    counts=counts,
                    state_order=state_order,
                    total_count=len(tip_states),
                )
            ),
        ),
        PosteriorPredictiveObservedStatisticRow(
            statistic_name=_MAJORITY_STATE_FREQUENCY,
            value=_round_float(_majority_state_frequency(counts, len(tip_states))),
        ),
    ]


def _build_discrete_trait_replicate_statistic_rows(
    *,
    traits: Sequence[SimulatedDiscreteTrait],
    state_order: Sequence[str],
    replicate_index: int,
    posterior_sample_index: int,
    posterior_iteration_index: int,
) -> list[PosteriorPredictiveReplicateStatisticRow]:
    counts = _state_counts(trait.state for trait in traits)
    total_count = len(traits)
    return [
        PosteriorPredictiveReplicateStatisticRow(
            statistic_name=_STATE_ENTROPY,
            replicate_index=replicate_index,
            posterior_sample_index=posterior_sample_index,
            posterior_iteration_index=posterior_iteration_index,
            value=_round_float(
                _categorical_state_entropy(
                    counts=counts,
                    state_order=state_order,
                    total_count=total_count,
                )
            ),
        ),
        PosteriorPredictiveReplicateStatisticRow(
            statistic_name=_MAJORITY_STATE_FREQUENCY,
            replicate_index=replicate_index,
            posterior_sample_index=posterior_sample_index,
            posterior_iteration_index=posterior_iteration_index,
            value=_round_float(_majority_state_frequency(counts, total_count)),
        ),
    ]


def _build_continuous_trait_observed_statistic_rows(
    *,
    tip_values: Mapping[str, float],
    taxa: Sequence[str],
) -> list[PosteriorPredictiveObservedStatisticRow]:
    values = [float(tip_values[taxon]) for taxon in taxa]
    return [
        PosteriorPredictiveObservedStatisticRow(
            statistic_name=_TIP_MEAN,
            value=_round_float(mean(values)),
        ),
        PosteriorPredictiveObservedStatisticRow(
            statistic_name=_TIP_VARIANCE,
            value=_round_float(_sample_variance(values)),
        ),
        PosteriorPredictiveObservedStatisticRow(
            statistic_name=_TIP_RANGE,
            value=_round_float(max(values) - min(values)),
        ),
    ]


def _build_continuous_trait_replicate_statistic_rows(
    *,
    traits: Sequence[SimulatedContinuousTrait],
    replicate_index: int,
    posterior_sample_index: int,
    posterior_iteration_index: int,
) -> list[PosteriorPredictiveReplicateStatisticRow]:
    values = [trait.value for trait in traits]
    return [
        PosteriorPredictiveReplicateStatisticRow(
            statistic_name=_TIP_MEAN,
            replicate_index=replicate_index,
            posterior_sample_index=posterior_sample_index,
            posterior_iteration_index=posterior_iteration_index,
            value=_round_float(mean(values)),
        ),
        PosteriorPredictiveReplicateStatisticRow(
            statistic_name=_TIP_VARIANCE,
            replicate_index=replicate_index,
            posterior_sample_index=posterior_sample_index,
            posterior_iteration_index=posterior_iteration_index,
            value=_round_float(_sample_variance(values)),
        ),
        PosteriorPredictiveReplicateStatisticRow(
            statistic_name=_TIP_RANGE,
            replicate_index=replicate_index,
            posterior_sample_index=posterior_sample_index,
            posterior_iteration_index=posterior_iteration_index,
            value=_round_float(max(values) - min(values)),
        ),
    ]


def _simulate_dna_alignment_replicate(
    *,
    tree: PhyloTree,
    record_ids: Sequence[str],
    alignment_length: int,
    model_name: str,
    scalar_parameters: Mapping[str, float],
    vector_parameters: Mapping[str, Mapping[str, float]],
    rng: random.Random,
) -> list[AlignmentRecord]:
    root_probabilities, transition_by_child_id = _build_dna_transition_surface(
        tree=tree,
        model_name=model_name,
        scalar_parameters=scalar_parameters,
        vector_parameters=vector_parameters,
    )
    sequences_by_taxon = {taxon: [] for taxon in record_ids}
    for _site_index in range(alignment_length):
        tip_states = _simulate_dna_site(
            tree=tree,
            root_probabilities=root_probabilities,
            transition_by_child_id=transition_by_child_id,
            rng=rng,
        )
        for taxon in record_ids:
            sequences_by_taxon[taxon].append(tip_states[taxon])
    return [
        AlignmentRecord(identifier=taxon, sequence="".join(sequences_by_taxon[taxon]))
        for taxon in record_ids
    ]


def _build_dna_transition_surface(
    *,
    tree: PhyloTree,
    model_name: str,
    scalar_parameters: Mapping[str, float],
    vector_parameters: Mapping[str, Mapping[str, float]],
) -> tuple[numpy.ndarray, dict[str, numpy.ndarray]]:
    if model_name == "JC69":
        root_probabilities = numpy.array(UNIFORM_DNA_ROOT_PRIOR, dtype=float)
        transition_by_child_id = {
            _require_node_id(child): jc69_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0)
            )
            for _parent, child in tree.iter_edges()
        }
        return root_probabilities, transition_by_child_id
    if model_name == "K80":
        root_probabilities = numpy.array(UNIFORM_DNA_ROOT_PRIOR, dtype=float)
        kappa = scalar_parameters["kappa"]
        transition_by_child_id = {
            _require_node_id(child): k80_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0),
                kappa=kappa,
            )
            for _parent, child in tree.iter_edges()
        }
        return root_probabilities, transition_by_child_id
    if model_name == "F81":
        base_frequencies = vector_parameters["base-frequencies"]
        root_probabilities = _dna_probability_vector(base_frequencies)
        transition_by_child_id = {
            _require_node_id(child): f81_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0),
                base_frequencies=base_frequencies,
            )
            for _parent, child in tree.iter_edges()
        }
        return root_probabilities, transition_by_child_id
    if model_name == "HKY85":
        base_frequencies = vector_parameters["base-frequencies"]
        kappa = scalar_parameters["kappa"]
        root_probabilities = _dna_probability_vector(base_frequencies)
        transition_by_child_id = {
            _require_node_id(child): hky85_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0),
                base_frequencies=base_frequencies,
                kappa=kappa,
            )
            for _parent, child in tree.iter_edges()
        }
        return root_probabilities, transition_by_child_id
    if model_name == "GTR":
        base_frequencies = vector_parameters["base-frequencies"]
        exchangeabilities = vector_parameters["exchangeabilities"]
        root_probabilities = _dna_probability_vector(base_frequencies)
        transition_by_child_id = {
            _require_node_id(child): gtr_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0),
                exchangeabilities=exchangeabilities,
                base_frequencies=base_frequencies,
            )
            for _parent, child in tree.iter_edges()
        }
        return root_probabilities, transition_by_child_id
    raise PhylogeneticsError(
        "posterior predictive DNA simulation encountered one unsupported substitution model",
        code="posterior_predictive_dna_model_unsupported",
        details={"model_name": model_name},
    )


def _simulate_dna_site(
    *,
    tree: PhyloTree,
    root_probabilities: numpy.ndarray,
    transition_by_child_id: Mapping[str, numpy.ndarray],
    rng: random.Random,
) -> dict[str, str]:
    root_state = DNA_STATE_ORDER[
        _sample_index_from_probabilities(root_probabilities, rng=rng)
    ]
    tip_states: dict[str, str] = {}
    _propagate_dna_state(
        node=tree.root,
        state=root_state,
        transition_by_child_id=transition_by_child_id,
        tip_states=tip_states,
        rng=rng,
    )
    return tip_states


def _propagate_dna_state(
    *,
    node: TreeNode,
    state: str,
    transition_by_child_id: Mapping[str, numpy.ndarray],
    tip_states: dict[str, str],
    rng: random.Random,
) -> None:
    if node.is_leaf():
        if node.name is None:
            raise PhylogeneticsError(
                "posterior predictive DNA simulation requires every sampled tip to have one name",
                code="posterior_predictive_dna_tip_name_missing",
            )
        tip_states[node.name] = state
        return
    parent_state_index = _DNA_STATE_INDEX[state]
    for child in node.children:
        child_node_id = _require_node_id(child)
        child_state_index = _sample_index_from_probabilities(
            transition_by_child_id[child_node_id][parent_state_index],
            rng=rng,
        )
        _propagate_dna_state(
            node=child,
            state=DNA_STATE_ORDER[child_state_index],
            transition_by_child_id=transition_by_child_id,
            tip_states=tip_states,
            rng=rng,
        )


def _simulate_discrete_trait_replicate(
    *,
    tree: PhyloTree,
    taxa: Sequence[str],
    state_order: Sequence[str],
    root_prior: numpy.ndarray,
    rate_matrix: numpy.ndarray,
    rng: random.Random,
) -> list[SimulatedDiscreteTrait]:
    transition_by_child_id = {
        _require_node_id(child): transition_probability_matrix(
            rate_matrix,
            max(float(child.branch_length or 0.0), 0.0),
        )
        for _parent, child in tree.iter_edges()
    }
    root_state = state_order[_sample_index_from_probabilities(root_prior, rng=rng)]
    tip_states: dict[str, str] = {}
    _propagate_discrete_trait_state(
        node=tree.root,
        state=root_state,
        state_order=state_order,
        transition_by_child_id=transition_by_child_id,
        tip_states=tip_states,
        rng=rng,
    )
    return [
        SimulatedDiscreteTrait(taxon=taxon, state=tip_states[taxon]) for taxon in taxa
    ]


def _propagate_discrete_trait_state(
    *,
    node: TreeNode,
    state: str,
    state_order: Sequence[str],
    transition_by_child_id: Mapping[str, numpy.ndarray],
    tip_states: dict[str, str],
    rng: random.Random,
) -> None:
    if node.is_leaf():
        if node.name is None:
            raise PhylogeneticsError(
                "posterior predictive discrete-trait simulation requires every sampled tip to have one name",
                code="posterior_predictive_discrete_tip_name_missing",
            )
        tip_states[node.name] = state
        return
    parent_state_index = state_order.index(state)
    for child in node.children:
        child_node_id = _require_node_id(child)
        child_state_index = _sample_index_from_probabilities(
            transition_by_child_id[child_node_id][parent_state_index],
            rng=rng,
        )
        _propagate_discrete_trait_state(
            node=child,
            state=state_order[child_state_index],
            state_order=state_order,
            transition_by_child_id=transition_by_child_id,
            tip_states=tip_states,
            rng=rng,
        )


def _simulate_brownian_continuous_trait_replicate(
    *,
    tree: PhyloTree,
    taxa: Sequence[str],
    root_state: float,
    sigma_squared: float,
    rng: random.Random,
) -> list[SimulatedContinuousTrait]:
    tip_values: dict[str, float] = {}
    _propagate_brownian_trait_value(
        node=tree.root,
        current_value=root_state,
        sigma_squared=sigma_squared,
        tip_values=tip_values,
        rng=rng,
    )
    return [
        SimulatedContinuousTrait(taxon=taxon, value=_round_float(tip_values[taxon]))
        for taxon in taxa
    ]


def _propagate_brownian_trait_value(
    *,
    node: TreeNode,
    current_value: float,
    sigma_squared: float,
    tip_values: dict[str, float],
    rng: random.Random,
) -> None:
    if node.is_leaf():
        if node.name is None:
            raise PhylogeneticsError(
                "posterior predictive Brownian simulation requires every sampled tip to have one name",
                code="posterior_predictive_brownian_tip_name_missing",
            )
        tip_values[node.name] = current_value
        return
    for child in node.children:
        branch_length = max(float(child.branch_length or 0.0), 0.0)
        child_value = current_value + rng.gauss(
            0.0,
            math.sqrt(sigma_squared * branch_length),
        )
        _propagate_brownian_trait_value(
            node=child,
            current_value=child_value,
            sigma_squared=sigma_squared,
            tip_values=tip_values,
            rng=rng,
        )


def _simulate_ou_continuous_trait_replicate(
    *,
    tree: PhyloTree,
    taxa: Sequence[str],
    alpha: float,
    optimum: float,
    sigma_squared: float,
    rng: random.Random,
) -> list[SimulatedContinuousTrait]:
    root_standard_deviation = math.sqrt(sigma_squared / (2.0 * alpha))
    root_value = rng.gauss(optimum, root_standard_deviation)
    tip_values: dict[str, float] = {}
    _propagate_ou_trait_value(
        node=tree.root,
        current_value=root_value,
        alpha=alpha,
        optimum=optimum,
        sigma_squared=sigma_squared,
        tip_values=tip_values,
        rng=rng,
    )
    return [
        SimulatedContinuousTrait(taxon=taxon, value=_round_float(tip_values[taxon]))
        for taxon in taxa
    ]


def _propagate_ou_trait_value(
    *,
    node: TreeNode,
    current_value: float,
    alpha: float,
    optimum: float,
    sigma_squared: float,
    tip_values: dict[str, float],
    rng: random.Random,
) -> None:
    if node.is_leaf():
        if node.name is None:
            raise PhylogeneticsError(
                "posterior predictive OU simulation requires every sampled tip to have one name",
                code="posterior_predictive_ou_tip_name_missing",
            )
        tip_values[node.name] = current_value
        return
    for child in node.children:
        branch_length = max(float(child.branch_length or 0.0), 0.0)
        decay = math.exp(-alpha * branch_length)
        conditional_mean = optimum + ((current_value - optimum) * decay)
        conditional_variance = (
            sigma_squared * (1.0 - math.exp(-2.0 * alpha * branch_length))
        ) / (2.0 * alpha)
        child_value = rng.gauss(
            conditional_mean,
            math.sqrt(max(conditional_variance, 0.0)),
        )
        _propagate_ou_trait_value(
            node=child,
            current_value=child_value,
            alpha=alpha,
            optimum=optimum,
            sigma_squared=sigma_squared,
            tip_values=tip_values,
            rng=rng,
        )


def _build_discrete_trait_rate_matrix(
    *,
    state: BayesianPhylogeneticState,
    state_order: Sequence[str],
    model_name: str,
    allowed_transition_pairs: set[tuple[int, int]],
) -> numpy.ndarray:
    parameter_values = state.model_parameters.vector_parameters["discrete-trait-rates"]
    template_log_parameters = numpy.zeros(len(parameter_values), dtype=float)
    template_rate_matrix = rate_matrix_from_log_parameters(
        template_log_parameters,
        state_order=list(state_order),
        model=model_name,
        state_ordering=_DISCRETE_TRAIT_STATE_ORDERING,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    template_rows = build_transition_rate_rows(
        state_order=list(state_order),
        state_ordering=_DISCRETE_TRAIT_STATE_ORDERING,
        rate_matrix=template_rate_matrix,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    parameterization = parameterize_discrete_trait_rate_rows(
        model=model_name,
        transition_rate_rows=template_rows,
    )
    ordered_log_parameters = numpy.array(
        [
            math.log(parameter_values[group.parameter_name])
            for group in parameterization.groups
        ],
        dtype=float,
    )
    return rate_matrix_from_log_parameters(
        ordered_log_parameters,
        state_order=list(state_order),
        model=model_name,
        state_ordering=_DISCRETE_TRAIT_STATE_ORDERING,
        allowed_transition_pairs=allowed_transition_pairs,
    )


def _merge_partition_alignment_records(
    *,
    record_ids: Sequence[str],
    partition_order: Sequence[LocusPartition],
    records_by_partition_name: Mapping[str, Sequence[AlignmentRecord]],
) -> list[AlignmentRecord]:
    merged_sequences_by_taxon = {taxon: [] for taxon in record_ids}
    for partition in partition_order:
        partition_records = records_by_partition_name[partition.name]
        partition_sequence_by_taxon = {
            record.identifier: record.sequence for record in partition_records
        }
        for taxon in record_ids:
            merged_sequences_by_taxon[taxon].append(partition_sequence_by_taxon[taxon])
    return [
        AlignmentRecord(
            identifier=taxon,
            sequence="".join(merged_sequences_by_taxon[taxon]),
        )
        for taxon in record_ids
    ]


def _partition_scalar_parameters(partition_state) -> dict[str, float]:
    scalar_parameters: dict[str, float] = {}
    if partition_state.kappa is not None:
        scalar_parameters["kappa"] = partition_state.kappa
    return scalar_parameters


def _partition_vector_parameters(partition_state) -> dict[str, dict[str, float]]:
    vector_parameters: dict[str, dict[str, float]] = {}
    if partition_state.base_frequencies is not None:
        vector_parameters["base-frequencies"] = dict(partition_state.base_frequencies)
    if partition_state.exchangeabilities is not None:
        vector_parameters["exchangeabilities"] = dict(partition_state.exchangeabilities)
    return vector_parameters


def _alignment_segregating_site_count(records: Sequence[AlignmentRecord]) -> float:
    sequences = [record.sequence.upper() for record in records]
    return float(
        sum(
            1
            for site_index in range(len(sequences[0]))
            if len({sequence[site_index] for sequence in sequences}) > 1
        )
    )


def _alignment_gc_fraction(records: Sequence[AlignmentRecord]) -> float:
    sequence_text = "".join(record.sequence.upper() for record in records)
    if not sequence_text:
        return 0.0
    gc_count = sum(state in {"G", "C"} for state in sequence_text)
    return gc_count / len(sequence_text)


def _categorical_state_entropy(
    *,
    counts: Mapping[str, int],
    state_order: Sequence[str],
    total_count: int,
) -> float:
    if total_count <= 0:
        return 0.0
    probabilities = [
        counts.get(state, 0) / total_count
        for state in state_order
        if counts.get(state, 0) > 0
    ]
    return -math.fsum(
        probability * math.log(probability) for probability in probabilities
    )


def _majority_state_frequency(counts: Mapping[str, int], total_count: int) -> float:
    if total_count <= 0:
        return 0.0
    return max(counts.values(), default=0) / total_count


def _state_counts(states: Sequence[str] | object) -> dict[str, int]:
    counts: dict[str, int] = {}
    for state in states:
        counts[state] = counts.get(state, 0) + 1
    return counts


def _dna_probability_vector(base_frequencies: Mapping[str, float]) -> numpy.ndarray:
    return numpy.array(
        [float(base_frequencies[state]) for state in DNA_STATE_ORDER],
        dtype=float,
    )


def _sample_index_from_probabilities(
    probabilities: numpy.ndarray | Sequence[float],
    *,
    rng: random.Random,
) -> int:
    threshold = rng.random()
    cumulative_probability = 0.0
    for index, probability in enumerate(probabilities):
        cumulative_probability += float(probability)
        if threshold <= cumulative_probability:
            return index
    return len(probabilities) - 1


def _sample_variance(values: Sequence[float]) -> float:
    if len(values) <= 1:
        return 0.0
    center = mean(values)
    return math.fsum((value - center) ** 2 for value in values) / (len(values) - 1)


def _require_node_id(node: TreeNode) -> str:
    if node.node_id is None:
        raise PhylogeneticsError(
            "posterior predictive simulation requires stable node identifiers on sampled trees",
            code="posterior_predictive_tree_node_id_missing",
        )
    return node.node_id


def _sample_standard_deviation(values: Sequence[float]) -> float:
    if len(values) <= 1:
        return 0.0
    center = mean(values)
    return math.sqrt(
        math.fsum((value - center) ** 2 for value in values) / (len(values) - 1)
    )


def _round_float(value: float) -> float:
    return float(format(value, ".15g"))


__all__ = [
    "POSTERIOR_PREDICTIVE_SAMPLE_SELECTION_POLICIES",
    "PosteriorPredictiveAlignmentReplicate",
    "PosteriorPredictiveAlignmentSimulationReport",
    "PosteriorPredictiveContinuousTraitReplicate",
    "PosteriorPredictiveContinuousTraitSimulationReport",
    "PosteriorPredictiveDiscreteTraitReplicate",
    "PosteriorPredictiveDiscreteTraitSimulationReport",
    "PosteriorPredictiveObservedStatisticRow",
    "PosteriorPredictiveReplicateStatisticRow",
    "PosteriorPredictiveSimulationDefinition",
    "PosteriorPredictiveStatisticSummaryRow",
    "build_posterior_predictive_simulation_definition",
    "simulate_brownian_continuous_trait_posterior_predictive",
    "simulate_discrete_trait_mk_posterior_predictive",
    "simulate_fixed_topology_dna_posterior_predictive",
    "simulate_fixed_topology_partitioned_dna_posterior_predictive",
    "simulate_joint_topology_dna_posterior_predictive",
    "simulate_ornstein_uhlenbeck_continuous_trait_posterior_predictive",
]
