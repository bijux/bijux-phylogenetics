from __future__ import annotations

import math
from math import exp, sqrt
from pathlib import Path
import random

from bijux_phylogenetics.core.tree import PhyloTree
from .models import (
    AlignmentSimulationReport,
    ContinuousTraitSimulationCollectionReport,
    ContinuousTraitSimulationReport,
    ContinuousTraitSimulationSummaryRow,
    CorrelatedContinuousTraitSimulationCollectionReport,
    CorrelatedContinuousTraitSimulationReport,
    DiscreteHistoryRateRow,
    DiscreteHistorySimulationCollectionReport,
    DiscreteHistorySummaryRow,
    DiscreteTraitSimulationReport,
    SimulatedContinuousNode,
    SimulatedContinuousTrait,
    SimulatedCorrelatedContinuousTrait,
    SimulatedDiscreteBranchHistory,
    SimulatedDiscreteNode,
    SimulatedDiscreteStateSegment,
    SimulatedDiscreteTrait,
    SimulatedDiscreteTransitionEvent,
    SimulatedTreeRecord,
    TreeSimulationEnvelopeMetric,
    TreeSimulationReport,
)
from .propagation import (
    _iter_node_trait_values,
    _iter_tip_trait_values,
    _resolve_brownian_sigma_parameters,
    _simulate_brownian_node_values,
    _tip_values_from_node_map,
)
from .statistics import (
    _mean,
    _median,
    _population_standard_deviation,
    _round_float,
    _sample_correlation,
    _sample_covariance,
    _sample_standard_deviation,
)


def simulate_birth_death_trees(
    *,
    tree_count: int,
    tip_count: int,
    birth_rate: float = 1.0,
    death_rate: float = 0.25,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
) -> tuple[list[PhyloTree], TreeSimulationReport]:
    from .trees import simulate_birth_death_trees as simulate_birth_death_trees_impl

    return simulate_birth_death_trees_impl(
        tree_count=tree_count,
        tip_count=tip_count,
        birth_rate=birth_rate,
        death_rate=death_rate,
        seed=seed,
        taxon_prefix=taxon_prefix,
    )


def simulate_random_trees(
    *,
    tree_count: int,
    tip_count: int,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
    branch_length_model: str = "uniform",
) -> tuple[list[PhyloTree], TreeSimulationReport]:
    from .trees import simulate_random_trees as simulate_random_trees_impl

    return simulate_random_trees_impl(
        tree_count=tree_count,
        tip_count=tip_count,
        seed=seed,
        taxon_prefix=taxon_prefix,
        branch_length_model=branch_length_model,
    )


def simulate_random_tree(
    *,
    tip_count: int,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
    branch_length_model: str = "uniform",
) -> tuple[PhyloTree, TreeSimulationReport]:
    from .trees import simulate_random_tree as simulate_random_tree_impl

    return simulate_random_tree_impl(
        tip_count=tip_count,
        seed=seed,
        taxon_prefix=taxon_prefix,
        branch_length_model=branch_length_model,
)


def _poisson_count(expected_changes: float, rng: random.Random) -> int:
    if expected_changes <= 0.0:
        return 0
    threshold = exp(-expected_changes)
    product = 1.0
    changes = 0
    while product > threshold:
        changes += 1
        product *= rng.random()
    return changes - 1


def _simulate_symmetric_state_trajectory(
    state: str,
    *,
    branch_length: float,
    rate: float,
    states: tuple[str, ...],
    rng: random.Random,
) -> tuple[
    str,
    list[SimulatedDiscreteTransitionEvent],
    list[SimulatedDiscreteStateSegment],
]:
    if rate < 0.0:
        raise ValueError(f"rate must be nonnegative, got {rate}")
    next_state = state
    event_states: list[tuple[str, str]] = []
    segment_boundaries = [0.0]
    for _ in range(_poisson_count(rate * branch_length, rng)):
        alternatives = [candidate for candidate in states if candidate != next_state]
        candidate = rng.choice(alternatives)
        event_states.append((next_state, candidate))
        next_state = candidate
    if event_states:
        event_distances = [
            branch_length * index / (len(event_states) + 1)
            for index in range(1, len(event_states) + 1)
        ]
        segment_boundaries.extend(event_distances)
    segment_boundaries.append(branch_length)
    events: list[SimulatedDiscreteTransitionEvent] = []
    segments: list[SimulatedDiscreteStateSegment] = []
    current_state = state
    parent_node = ""
    child_node = ""
    for index, (start_distance, end_distance) in enumerate(
        zip(segment_boundaries[:-1], segment_boundaries[1:], strict=True),
        start=1,
    ):
        segments.append(
            SimulatedDiscreteStateSegment(
                parent_node=parent_node,
                child_node=child_node,
                state=current_state,
                start_distance=_round_float(start_distance),
                end_distance=_round_float(end_distance),
                duration=_round_float(max(end_distance - start_distance, 0.0)),
            )
        )
        if index <= len(event_states):
            source_state, target_state = event_states[index - 1]
            events.append(
                SimulatedDiscreteTransitionEvent(
                    parent_node=parent_node,
                    child_node=child_node,
                    source_state=source_state,
                    target_state=target_state,
                    event_index=index,
                    branch_distance=_round_float(end_distance),
                )
            )
            current_state = target_state
    return next_state, events, segments


def _normalize_discrete_states(states: list[str]) -> tuple[str, ...]:
    unique_states = tuple(dict.fromkeys(state for state in states if state))
    if len(unique_states) < 2:
        raise ValueError("states must contain at least two distinct non-empty states")
    return unique_states


def _normalize_root_state_probabilities(
    *,
    states: tuple[str, ...],
    root_state: str | None,
    root_state_probabilities: dict[str, float] | None,
) -> dict[str, float]:
    if root_state is not None and root_state_probabilities is not None:
        raise ValueError(
            "root_state and root_state_probabilities cannot be supplied together"
        )
    state_set = set(states)
    if root_state is not None:
        if root_state not in state_set:
            raise ValueError(f"root_state '{root_state}' is not present in states")
        return {state: 1.0 if state == root_state else 0.0 for state in states}
    if root_state_probabilities is None:
        probability = 1.0 / len(states)
        return {state: _round_float(probability) for state in states}
    unknown_states = set(root_state_probabilities).difference(state_set)
    if unknown_states:
        unknown_state = sorted(unknown_states)[0]
        raise ValueError(
            f"root_state_probabilities contains unknown state '{unknown_state}'"
        )
    probabilities = {
        state: float(root_state_probabilities.get(state, 0.0)) for state in states
    }
    if any(value < 0.0 for value in probabilities.values()):
        raise ValueError("root_state_probabilities cannot contain negative values")
    total = sum(probabilities.values())
    if total <= 0.0:
        raise ValueError("root_state_probabilities must sum to a positive value")
    return {
        state: _round_float(value / total) for state, value in probabilities.items()
    }


def simulate_discrete_histories(
    tree_path: Path,
    *,
    states: list[str],
    rate_rows: list[DiscreteHistoryRateRow],
    root_state: str | None = None,
    root_state_probabilities: dict[str, float] | None = None,
    transform: str | None = None,
    transform_parameter_value: float | None = None,
    replicates: int = 1,
    seed: int = 1,
) -> DiscreteHistorySimulationCollectionReport:
    from .discrete_histories import (
        simulate_discrete_histories as simulate_discrete_histories_impl,
    )

    return simulate_discrete_histories_impl(
        tree_path,
        states=states,
        rate_rows=rate_rows,
        root_state=root_state,
        root_state_probabilities=root_state_probabilities,
        transform=transform,
        transform_parameter_value=transform_parameter_value,
        replicates=replicates,
        seed=seed,
    )


def simulate_coalescent_trees(
    *,
    tree_count: int,
    tip_count: int,
    population_size: float = 1.0,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
) -> tuple[list[PhyloTree], TreeSimulationReport]:
    from .trees import simulate_coalescent_trees as simulate_coalescent_trees_impl

    return simulate_coalescent_trees_impl(
        tree_count=tree_count,
        tip_count=tip_count,
        population_size=population_size,
        seed=seed,
        taxon_prefix=taxon_prefix,
    )


def simulate_coalescent_tree(
    *,
    tip_count: int,
    population_size: float = 1.0,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
) -> tuple[PhyloTree, TreeSimulationReport]:
    from .trees import simulate_coalescent_tree as simulate_coalescent_tree_impl

    return simulate_coalescent_tree_impl(
        tip_count=tip_count,
        population_size=population_size,
        seed=seed,
        taxon_prefix=taxon_prefix,
    )


def simulate_brownian_traits(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float | None = None,
    sigma_squared: float | None = None,
    seed: int = 1,
) -> ContinuousTraitSimulationReport:
    from .continuous import simulate_brownian_traits as simulate_brownian_traits_impl

    return simulate_brownian_traits_impl(
        tree_path,
        root_state=root_state,
        sigma=sigma,
        sigma_squared=sigma_squared,
        seed=seed,
    )


def simulate_ou_traits(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float = 1.0,
    alpha: float = 1.0,
    theta: float = 0.0,
    seed: int = 1,
) -> ContinuousTraitSimulationReport:
    from .continuous import simulate_ou_traits as simulate_ou_traits_impl

    return simulate_ou_traits_impl(
        tree_path,
        root_state=root_state,
        sigma=sigma,
        alpha=alpha,
        theta=theta,
        seed=seed,
    )


def simulate_early_burst_traits(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float = 1.0,
    rate_change: float = 1.0,
    seed: int = 1,
) -> ContinuousTraitSimulationReport:
    from .continuous import (
        simulate_early_burst_traits as simulate_early_burst_traits_impl,
    )

    return simulate_early_burst_traits_impl(
        tree_path,
        root_state=root_state,
        sigma=sigma,
        rate_change=rate_change,
        seed=seed,
    )


def simulate_speciational_traits(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float | None = None,
    sigma_squared: float | None = None,
    seed: int = 1,
) -> ContinuousTraitSimulationReport:
    from .continuous import (
        simulate_speciational_traits as simulate_speciational_traits_impl,
    )

    return simulate_speciational_traits_impl(
        tree_path,
        root_state=root_state,
        sigma=sigma,
        sigma_squared=sigma_squared,
        seed=seed,
    )


def simulate_brownian_trait_collection(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float | None = None,
    sigma_squared: float | None = None,
    replicates: int = 128,
    seed: int = 1,
) -> ContinuousTraitSimulationCollectionReport:
    from .continuous import (
        simulate_brownian_trait_collection as simulate_brownian_trait_collection_impl,
    )

    return simulate_brownian_trait_collection_impl(
        tree_path,
        root_state=root_state,
        sigma=sigma,
        sigma_squared=sigma_squared,
        replicates=replicates,
        seed=seed,
    )


def simulate_speciational_trait_collection(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float | None = None,
    sigma_squared: float | None = None,
    replicates: int = 128,
    seed: int = 1,
) -> ContinuousTraitSimulationCollectionReport:
    from .continuous import (
        simulate_speciational_trait_collection as simulate_speciational_trait_collection_impl,
    )

    return simulate_speciational_trait_collection_impl(
        tree_path,
        root_state=root_state,
        sigma=sigma,
        sigma_squared=sigma_squared,
        replicates=replicates,
        seed=seed,
    )


def simulate_correlated_brownian_traits(
    tree_path: Path,
    *,
    trait_names: list[str] | tuple[str, ...],
    evolutionary_covariance_matrix: list[list[float]]
    | tuple[tuple[float, ...], ...]
    | None = None,
    evolutionary_correlation_matrix: list[list[float]]
    | tuple[tuple[float, ...], ...]
    | None = None,
    trait_standard_deviations: list[float] | tuple[float, ...] | None = None,
    root_states: list[float] | tuple[float, ...] | None = None,
    seed: int = 1,
) -> CorrelatedContinuousTraitSimulationReport:
    from .correlated import (
        simulate_correlated_brownian_traits as simulate_correlated_brownian_traits_impl,
    )

    return simulate_correlated_brownian_traits_impl(
        tree_path,
        trait_names=trait_names,
        evolutionary_covariance_matrix=evolutionary_covariance_matrix,
        evolutionary_correlation_matrix=evolutionary_correlation_matrix,
        trait_standard_deviations=trait_standard_deviations,
        root_states=root_states,
        seed=seed,
    )


def simulate_correlated_brownian_trait_collection(
    tree_path: Path,
    *,
    trait_names: list[str] | tuple[str, ...],
    evolutionary_covariance_matrix: list[list[float]]
    | tuple[tuple[float, ...], ...]
    | None = None,
    evolutionary_correlation_matrix: list[list[float]]
    | tuple[tuple[float, ...], ...]
    | None = None,
    trait_standard_deviations: list[float] | tuple[float, ...] | None = None,
    root_states: list[float] | tuple[float, ...] | None = None,
    replicates: int = 128,
    seed: int = 1,
) -> CorrelatedContinuousTraitSimulationCollectionReport:
    from .correlated import (
        simulate_correlated_brownian_trait_collection as simulate_correlated_brownian_trait_collection_impl,
    )

    return simulate_correlated_brownian_trait_collection_impl(
        tree_path,
        trait_names=trait_names,
        evolutionary_covariance_matrix=evolutionary_covariance_matrix,
        evolutionary_correlation_matrix=evolutionary_correlation_matrix,
        trait_standard_deviations=trait_standard_deviations,
        root_states=root_states,
        replicates=replicates,
        seed=seed,
    )


def simulate_discrete_traits(
    tree_path: Path,
    *,
    states: list[str],
    transition_rate: float = 1.0,
    root_state: str | None = None,
    seed: int = 1,
) -> DiscreteTraitSimulationReport:
    from .discrete_traits import (
        simulate_discrete_traits as simulate_discrete_traits_impl,
    )

    return simulate_discrete_traits_impl(
        tree_path,
        states=states,
        transition_rate=transition_rate,
        root_state=root_state,
        seed=seed,
    )


def simulate_dna_alignment(
    tree_path: Path,
    *,
    sequence_length: int,
    substitution_rate: float = 1.0,
    seed: int = 1,
) -> AlignmentSimulationReport:
    from .alignment import simulate_dna_alignment as simulate_dna_alignment_impl

    return simulate_dna_alignment_impl(
        tree_path,
        sequence_length=sequence_length,
        substitution_rate=substitution_rate,
        seed=seed,
    )


def simulate_protein_alignment(
    tree_path: Path,
    *,
    sequence_length: int,
    substitution_rate: float = 1.0,
    seed: int = 1,
) -> AlignmentSimulationReport:
    from .alignment import (
        simulate_protein_alignment as simulate_protein_alignment_impl,
    )

    return simulate_protein_alignment_impl(
        tree_path,
        sequence_length=sequence_length,
        substitution_rate=substitution_rate,
        seed=seed,
    )


def write_tree_set(path: Path, trees: list[PhyloTree]) -> Path:
    from .trees import write_tree_set as write_tree_set_impl

    return write_tree_set_impl(path, trees)


def write_simulated_tree(path: Path, tree: PhyloTree) -> Path:
    from .trees import write_simulated_tree as write_simulated_tree_impl

    return write_simulated_tree_impl(path, tree)


def write_tree_simulation_record_table(
    path: Path, report: TreeSimulationReport
) -> Path:
    from .trees import (
        write_tree_simulation_record_table as write_tree_simulation_record_table_impl,
    )

    return write_tree_simulation_record_table_impl(path, report)


def write_tree_simulation_envelope_table(
    path: Path,
    report: TreeSimulationReport,
) -> Path:
    from .trees import (
        write_tree_simulation_envelope_table as write_tree_simulation_envelope_table_impl,
    )

    return write_tree_simulation_envelope_table_impl(path, report)


def write_continuous_trait_table(
    path: Path, report: ContinuousTraitSimulationReport
) -> Path:
    from .continuous import write_continuous_trait_table as write_continuous_trait_table_impl

    return write_continuous_trait_table_impl(path, report)


def write_continuous_trait_collection_table(
    path: Path,
    report: ContinuousTraitSimulationCollectionReport,
) -> Path:
    from .continuous import (
        write_continuous_trait_collection_table as write_continuous_trait_collection_table_impl,
    )

    return write_continuous_trait_collection_table_impl(path, report)


def write_continuous_trait_collection_summary_table(
    path: Path,
    report: ContinuousTraitSimulationCollectionReport,
) -> Path:
    from .continuous import (
        write_continuous_trait_collection_summary_table as write_continuous_trait_collection_summary_table_impl,
    )

    return write_continuous_trait_collection_summary_table_impl(path, report)


def write_correlated_continuous_trait_table(
    path: Path,
    report: CorrelatedContinuousTraitSimulationReport,
) -> Path:
    from .correlated import (
        write_correlated_continuous_trait_table as write_correlated_continuous_trait_table_impl,
    )

    return write_correlated_continuous_trait_table_impl(path, report)


def write_correlated_continuous_trait_collection_table(
    path: Path,
    report: CorrelatedContinuousTraitSimulationCollectionReport,
) -> Path:
    from .correlated import (
        write_correlated_continuous_trait_collection_table as write_correlated_continuous_trait_collection_table_impl,
    )

    return write_correlated_continuous_trait_collection_table_impl(path, report)


def write_correlated_continuous_trait_collection_summary_table(
    path: Path,
    report: CorrelatedContinuousTraitSimulationCollectionReport,
) -> Path:
    from .correlated import (
        write_correlated_continuous_trait_collection_summary_table as write_correlated_continuous_trait_collection_summary_table_impl,
    )

    return write_correlated_continuous_trait_collection_summary_table_impl(path, report)


def write_discrete_trait_table(
    path: Path, report: DiscreteTraitSimulationReport
) -> Path:
    from .discrete_traits import (
        write_discrete_trait_table as write_discrete_trait_table_impl,
    )

    return write_discrete_trait_table_impl(path, report)


def write_discrete_history_tip_truth_table(
    path: Path,
    report: DiscreteHistorySimulationCollectionReport,
) -> Path:
    from .discrete_histories import (
        write_discrete_history_tip_truth_table as write_discrete_history_tip_truth_table_impl,
    )

    return write_discrete_history_tip_truth_table_impl(path, report)


def write_discrete_history_node_truth_table(
    path: Path,
    report: DiscreteHistorySimulationCollectionReport,
) -> Path:
    from .discrete_histories import (
        write_discrete_history_node_truth_table as write_discrete_history_node_truth_table_impl,
    )

    return write_discrete_history_node_truth_table_impl(path, report)


def write_discrete_history_branch_truth_table(
    path: Path,
    report: DiscreteHistorySimulationCollectionReport,
) -> Path:
    from .discrete_histories import (
        write_discrete_history_branch_truth_table as write_discrete_history_branch_truth_table_impl,
    )

    return write_discrete_history_branch_truth_table_impl(path, report)


def write_discrete_history_event_table(
    path: Path,
    report: DiscreteHistorySimulationCollectionReport,
) -> Path:
    from .discrete_histories import (
        write_discrete_history_event_table as write_discrete_history_event_table_impl,
    )

    return write_discrete_history_event_table_impl(path, report)


def write_discrete_history_segment_table(
    path: Path,
    report: DiscreteHistorySimulationCollectionReport,
) -> Path:
    from .discrete_histories import (
        write_discrete_history_segment_table as write_discrete_history_segment_table_impl,
    )

    return write_discrete_history_segment_table_impl(path, report)


def write_discrete_history_summary_table(
    path: Path,
    report: DiscreteHistorySimulationCollectionReport,
) -> Path:
    from .discrete_histories import (
        write_discrete_history_summary_table as write_discrete_history_summary_table_impl,
    )

    return write_discrete_history_summary_table_impl(path, report)


def write_simulated_alignment(path: Path, report: AlignmentSimulationReport) -> Path:
    from .alignment import write_simulated_alignment as write_simulated_alignment_impl

    return write_simulated_alignment_impl(path, report)


def validate_geiger_sim_char_reference_examples():
    from .geiger_sim_char_reference import (
        validate_geiger_sim_char_reference_examples as validate_geiger_sim_char_reference_examples_impl,
    )

    return validate_geiger_sim_char_reference_examples_impl()
