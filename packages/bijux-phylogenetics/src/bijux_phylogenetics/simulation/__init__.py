"""Curated public API for simulation workflows and generated truth ledgers."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, Any

from .contracts import (
    AlignmentSimulationReport,
    CoalescentSkylineSummaryRow,
    CoalescentWaitingTimeSummaryRow,
    ContinuousTraitSimulationCollectionReport,
    ContinuousTraitSimulationReport,
    ContinuousTraitSimulationSummaryRow,
    CorrelatedContinuousTraitSimulationCollectionReport,
    CorrelatedContinuousTraitSimulationReport,
    DiscreteHistoryRateRow,
    DiscreteHistorySimulationCollectionReport,
    DiscreteHistorySummaryRow,
    DiscreteTraitSimulationReport,
    MultispeciesCoalescentBranchRow,
    MultispeciesCoalescentEventRow,
    MultispeciesCoalescentReport,
    MultispeciesCoalescentSampleRow,
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

_PUBLIC_NAME_TO_MODULE = {
    "simulate_birth_death_trees": ".trees",
    "simulate_coalescent_tree": ".trees",
    "simulate_coalescent_trees": ".trees",
    "simulate_multispecies_coalescent_gene_tree": ".trees",
    "simulate_random_tree": ".trees",
    "simulate_random_trees": ".trees",
    "write_coalescent_skyline_table": ".trees",
    "write_coalescent_waiting_time_table": ".trees",
    "write_multispecies_coalescent_branch_table": ".trees",
    "write_multispecies_coalescent_event_table": ".trees",
    "write_simulated_tree": ".trees",
    "write_tree_set": ".trees",
    "write_tree_simulation_envelope_table": ".trees",
    "write_tree_simulation_record_table": ".trees",
    "simulate_brownian_trait_collection": ".continuous",
    "simulate_brownian_traits": ".continuous",
    "simulate_early_burst_traits": ".continuous",
    "simulate_ou_traits": ".continuous",
    "simulate_speciational_trait_collection": ".continuous",
    "simulate_speciational_traits": ".continuous",
    "write_continuous_trait_collection_summary_table": ".continuous",
    "write_continuous_trait_collection_table": ".continuous",
    "write_continuous_trait_table": ".continuous",
    "simulate_correlated_brownian_trait_collection": ".continuous",
    "simulate_correlated_brownian_traits": ".continuous",
    "write_correlated_continuous_trait_collection_summary_table": ".continuous",
    "write_correlated_continuous_trait_collection_table": ".continuous",
    "write_correlated_continuous_trait_table": ".continuous",
    "simulate_discrete_traits": ".discrete",
    "write_discrete_trait_table": ".discrete",
    "simulate_discrete_histories": ".discrete",
    "write_discrete_history_branch_truth_table": ".discrete",
    "write_discrete_history_event_table": ".discrete",
    "write_discrete_history_node_truth_table": ".discrete",
    "write_discrete_history_segment_table": ".discrete",
    "write_discrete_history_summary_table": ".discrete",
    "write_discrete_history_tip_truth_table": ".discrete",
    "simulate_dna_alignment": ".alignment",
    "simulate_protein_alignment": ".alignment",
    "write_simulated_alignment": ".alignment",
    "validate_geiger_sim_char_reference_examples": ".reference",
}

__all__ = [
    "AlignmentSimulationReport",
    "CoalescentSkylineSummaryRow",
    "CoalescentWaitingTimeSummaryRow",
    "ContinuousTraitSimulationCollectionReport",
    "ContinuousTraitSimulationReport",
    "ContinuousTraitSimulationSummaryRow",
    "CorrelatedContinuousTraitSimulationCollectionReport",
    "CorrelatedContinuousTraitSimulationReport",
    "DiscreteHistoryRateRow",
    "DiscreteHistorySimulationCollectionReport",
    "DiscreteHistorySummaryRow",
    "DiscreteTraitSimulationReport",
    "MultispeciesCoalescentBranchRow",
    "MultispeciesCoalescentEventRow",
    "MultispeciesCoalescentReport",
    "MultispeciesCoalescentSampleRow",
    "SimulatedContinuousNode",
    "SimulatedContinuousTrait",
    "SimulatedCorrelatedContinuousTrait",
    "SimulatedDiscreteBranchHistory",
    "SimulatedDiscreteNode",
    "SimulatedDiscreteStateSegment",
    "SimulatedDiscreteTrait",
    "SimulatedDiscreteTransitionEvent",
    "SimulatedTreeRecord",
    "TreeSimulationEnvelopeMetric",
    "TreeSimulationReport",
    *_PUBLIC_NAME_TO_MODULE,
]

if TYPE_CHECKING:
    from .alignment import (
        simulate_dna_alignment as simulate_dna_alignment,
    )
    from .alignment import (
        simulate_protein_alignment as simulate_protein_alignment,
    )
    from .alignment import (
        write_simulated_alignment as write_simulated_alignment,
    )
    from .continuous import (
        simulate_brownian_trait_collection as simulate_brownian_trait_collection,
    )
    from .continuous import (
        simulate_brownian_traits as simulate_brownian_traits,
    )
    from .continuous import (
        simulate_correlated_brownian_trait_collection as simulate_correlated_brownian_trait_collection,
    )
    from .continuous import (
        simulate_correlated_brownian_traits as simulate_correlated_brownian_traits,
    )
    from .continuous import (
        simulate_early_burst_traits as simulate_early_burst_traits,
    )
    from .continuous import (
        simulate_ou_traits as simulate_ou_traits,
    )
    from .continuous import (
        simulate_speciational_trait_collection as simulate_speciational_trait_collection,
    )
    from .continuous import (
        simulate_speciational_traits as simulate_speciational_traits,
    )
    from .continuous import (
        write_continuous_trait_collection_summary_table as write_continuous_trait_collection_summary_table,
    )
    from .continuous import (
        write_continuous_trait_collection_table as write_continuous_trait_collection_table,
    )
    from .continuous import (
        write_continuous_trait_table as write_continuous_trait_table,
    )
    from .continuous import (
        write_correlated_continuous_trait_collection_summary_table as write_correlated_continuous_trait_collection_summary_table,
    )
    from .continuous import (
        write_correlated_continuous_trait_collection_table as write_correlated_continuous_trait_collection_table,
    )
    from .continuous import (
        write_correlated_continuous_trait_table as write_correlated_continuous_trait_table,
    )
    from .discrete import (
        simulate_discrete_histories as simulate_discrete_histories,
    )
    from .discrete import (
        simulate_discrete_traits as simulate_discrete_traits,
    )
    from .discrete import (
        write_discrete_history_branch_truth_table as write_discrete_history_branch_truth_table,
    )
    from .discrete import (
        write_discrete_history_event_table as write_discrete_history_event_table,
    )
    from .discrete import (
        write_discrete_history_node_truth_table as write_discrete_history_node_truth_table,
    )
    from .discrete import (
        write_discrete_history_segment_table as write_discrete_history_segment_table,
    )
    from .discrete import (
        write_discrete_history_summary_table as write_discrete_history_summary_table,
    )
    from .discrete import (
        write_discrete_history_tip_truth_table as write_discrete_history_tip_truth_table,
    )
    from .discrete import (
        write_discrete_trait_table as write_discrete_trait_table,
    )
    from .reference import (
        validate_geiger_sim_char_reference_examples as validate_geiger_sim_char_reference_examples,
    )
    from .trees import (
        simulate_birth_death_trees as simulate_birth_death_trees,
    )
    from .trees import (
        simulate_coalescent_tree as simulate_coalescent_tree,
    )
    from .trees import (
        simulate_coalescent_trees as simulate_coalescent_trees,
    )
    from .trees import (
        simulate_multispecies_coalescent_gene_tree as simulate_multispecies_coalescent_gene_tree,
    )
    from .trees import (
        simulate_random_tree as simulate_random_tree,
    )
    from .trees import (
        simulate_random_trees as simulate_random_trees,
    )
    from .trees import (
        write_coalescent_skyline_table as write_coalescent_skyline_table,
    )
    from .trees import (
        write_coalescent_waiting_time_table as write_coalescent_waiting_time_table,
    )
    from .trees import (
        write_multispecies_coalescent_branch_table as write_multispecies_coalescent_branch_table,
    )
    from .trees import (
        write_multispecies_coalescent_event_table as write_multispecies_coalescent_event_table,
    )
    from .trees import (
        write_simulated_tree as write_simulated_tree,
    )
    from .trees import (
        write_tree_set as write_tree_set,
    )
    from .trees import (
        write_tree_simulation_envelope_table as write_tree_simulation_envelope_table,
    )
    from .trees import (
        write_tree_simulation_record_table as write_tree_simulation_record_table,
    )


def __getattr__(name: str) -> Any:
    if name not in _PUBLIC_NAME_TO_MODULE:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(_PUBLIC_NAME_TO_MODULE[name], __name__)
    value = getattr(module, name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
