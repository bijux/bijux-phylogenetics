from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..tree_sets import (
    CladeFrequencyReport,
    ConsensusTreeReport,
    GeneTreeQuartetConcordanceRow,
    TreeSetProcessingSummary,
    TreeSetReport,
    TreeSetWorkflowBudgetReport,
)


@dataclass(frozen=True, slots=True)
class TreeDistanceDistributionRow:
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    pair_count: int
    frequency: float


@dataclass(slots=True)
class PosteriorTopologyDiversityReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    rooted_topology_count: int
    dominant_topology_frequency: float
    effective_topology_count: float
    pair_count: int
    mean_robinson_foulds_distance: float
    mean_normalized_robinson_foulds_distance: float
    maximum_robinson_foulds_distance: int
    maximum_normalized_robinson_foulds_distance: float
    unstable_clade_count: int
    rf_distribution: list[TreeDistanceDistributionRow]


@dataclass(frozen=True, slots=True)
class TreeTopologyCluster:
    rooted_topology_id: str
    tree_indices: list[int]
    tree_count: int
    frequency: float
    representative_index: int
    representative_newick: str


@dataclass(slots=True)
class TreeTopologyClusterReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    rooted_topology_count: int
    clusters: list[TreeTopologyCluster]


@dataclass(frozen=True, slots=True)
class TaxonPlacementSignature:
    signature: str
    tree_count: int
    frequency: float


@dataclass(frozen=True, slots=True)
class UnstableTaxon:
    taxon: str
    unique_placements: int
    dominant_frequency: float
    instability_score: float
    placements: list[TaxonPlacementSignature]


@dataclass(slots=True)
class UnstableTaxaReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    taxa: list[UnstableTaxon]


@dataclass(frozen=True, slots=True)
class RogueTaxonScoreRow:
    taxon: str
    rank: int
    mean_terminal_branch_length: float | None
    baseline_consensus_resolution: float
    pruned_consensus_resolution: float
    consensus_resolution_delta: float
    baseline_mean_support_percent: float
    pruned_mean_support_percent: float
    mean_support_percent_delta: float
    baseline_mean_normalized_robinson_foulds: float
    pruned_mean_normalized_robinson_foulds: float
    normalized_robinson_foulds_stability_delta: float
    baseline_rooted_topology_count: int
    pruned_rooted_topology_count: int
    rooted_topology_count_delta: int
    baseline_dominant_topology_frequency: float
    pruned_dominant_topology_frequency: float
    dominant_topology_frequency_delta: float
    pruned_consensus_newick: str


@dataclass(slots=True)
class RogueTaxonDetectionReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    consensus_threshold: float
    ranking_objective: str
    baseline_consensus_newick: str
    baseline_consensus_resolution: float
    baseline_mean_support_percent: float
    baseline_mean_normalized_robinson_foulds: float
    baseline_rooted_topology_count: int
    baseline_dominant_topology_frequency: float
    rows: list[RogueTaxonScoreRow]


@dataclass(frozen=True, slots=True)
class UnstableClade:
    clade: str
    tree_count: int
    frequency: float
    conflict_count: int
    instability_score: float
    support_classification: str
    conflicting_clades: list[str]


@dataclass(slots=True)
class UnstableCladeReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    clades: list[UnstableClade]


@dataclass(frozen=True, slots=True)
class BootstrapUnstableBranch:
    clade: str
    bootstrap_tree_count: int
    bootstrap_frequency: float
    bootstrap_support_percent: float
    conflict_count: int
    instability_score: float
    support_classification: str
    conflicting_clades: list[str]


@dataclass(slots=True)
class BootstrapTreeSetSummaryReport:
    path: Path
    consensus_threshold: float
    robust_support_threshold: float
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    summary: TreeSetReport
    clade_frequencies: CladeFrequencyReport
    consensus: ConsensusTreeReport
    diversity: PosteriorTopologyDiversityReport
    unstable_clades: UnstableCladeReport
    unstable_branch_count: int
    unstable_branches: list[BootstrapUnstableBranch]
    warnings: list[str]


@dataclass(slots=True)
class BootstrapTreeSetArtifactReport:
    input_path: Path
    out_dir: Path
    prefix: str
    summary_report: BootstrapTreeSetSummaryReport
    budget_report: TreeSetWorkflowBudgetReport
    output_paths: dict[str, Path]


@dataclass(frozen=True, slots=True)
class CladeFrequencyDelta:
    clade: str
    left_frequency: float
    right_frequency: float
    delta: float


@dataclass(slots=True)
class PosteriorTreeSetComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_tree_count: int
    right_tree_count: int
    left_rooted_topology_count: int
    right_rooted_topology_count: int
    shared_rooted_topology_count: int
    mean_between_set_robinson_foulds: float
    mean_between_set_normalized_robinson_foulds: float
    clade_frequency_deltas: list[CladeFrequencyDelta]


@dataclass(frozen=True, slots=True)
class PosteriorTopologicalDiversitySummary:
    tree_count: int
    rooted_topology_count: int
    dominant_topology_frequency: float
    effective_topology_count: float
    mean_within_set_robinson_foulds: float
    mean_within_set_normalized_robinson_foulds: float


@dataclass(slots=True)
class PosteriorTopologicalDiversityComparisonReport:
    left_path: Path
    right_path: Path
    left_summary: PosteriorTopologicalDiversitySummary
    right_summary: PosteriorTopologicalDiversitySummary
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class PosteriorTopologyMode:
    rooted_topology_id: str
    representative_index: int
    representative_newick: str
    tree_indices: list[int]
    tree_count: int
    frequency: float


@dataclass(slots=True)
class PosteriorTopologyMultimodalityReport:
    path: Path
    tree_count: int
    rooted_topology_count: int
    dominant_mode_frequency: float
    mode_count: int
    multimodal: bool
    modes: list[PosteriorTopologyMode]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class CladeCredibilityConflict:
    left_clade: str
    left_frequency: float
    right_clade: str
    right_frequency: float
    combined_frequency: float


@dataclass(slots=True)
class CladeCredibilityConflictReport:
    path: Path
    tree_count: int
    credibility_threshold: float
    high_credibility_clade_count: int
    conflict_count: int
    conflicts: list[CladeCredibilityConflict]


@dataclass(frozen=True, slots=True)
class GeneTreeConflictReferenceTree:
    selection_method: str
    rooted_topology_id: str
    frequency: float
    newick: str


@dataclass(slots=True)
class GeneTreeConflictQuartetSummary:
    branch_count: int
    total_quartet_count: int
    concordant_quartet_count: int
    discordant_first_quartet_count: int
    discordant_second_quartet_count: int
    uninformative_quartet_count: int
    informative_quartet_count: int
    rows: list[GeneTreeQuartetConcordanceRow]


@dataclass(slots=True)
class GeneTreeConflictSummaryReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    credibility_threshold: float
    rogue_consensus_threshold: float
    reference_tree: GeneTreeConflictReferenceTree
    clade_frequencies: CladeFrequencyReport
    quartet_concordance: GeneTreeConflictQuartetSummary
    rogue_taxa: RogueTaxonDetectionReport
    clade_conflicts: CladeCredibilityConflictReport


@dataclass(slots=True)
class GeneTreeConflictArtifactReport:
    input_path: Path
    out_dir: Path
    prefix: str
    summary_report: GeneTreeConflictSummaryReport
    output_paths: dict[str, Path]


@dataclass(frozen=True, slots=True)
class UncertaintyAwareCladeConclusion:
    clade: str
    frequency: float
    conclusion: str
    rationale: str


@dataclass(slots=True)
class UncertaintyAwareConclusionSummaryReport:
    path: Path
    tree_count: int
    robust_clade_count: int
    uncertain_clade_count: int
    conflicting_clade_count: int
    robust_clades: list[UncertaintyAwareCladeConclusion]
    uncertain_clades: list[UncertaintyAwareCladeConclusion]
    conflicting_clades: list[UncertaintyAwareCladeConclusion]


@dataclass(frozen=True, slots=True)
class BootstrapPosteriorCladeComparison:
    clade: str
    bootstrap_support: float | None
    posterior_frequency: float | None
    absolute_delta: float | None
    agreement: str


@dataclass(slots=True)
class BootstrapPosteriorSupportComparisonReport:
    bootstrap_tree_path: Path
    posterior_tree_set_path: Path
    posterior_tree_count: int
    shared_taxa: list[str]
    high_conflict_clade_count: int
    topology_mismatch_detected: bool
    topology_mismatch_clade_count: int
    rows: list[BootstrapPosteriorCladeComparison]


@dataclass(frozen=True, slots=True)
class TreeSetBenchmarkRow:
    tree_count: int
    taxon_count: int
    replicate: int
    elapsed_seconds: float
    peak_memory_bytes: int
    rooted_topology_count: int
    unstable_taxon_count: int
    unstable_clade_count: int
    robust_clade_count: int


@dataclass(slots=True)
class TreeSetScalingBenchmarkReport:
    tree_counts: list[int]
    taxon_counts: list[int]
    rows: list[TreeSetBenchmarkRow]


@dataclass(slots=True)
class TreeSetStorageRiskReport:
    path: Path
    file_size_bytes: int
    file_size_megabytes: float
    tree_count: int
    rooted_topology_count: int
    shared_taxon_count: int
    mean_bytes_per_tree: float
    risk_level: str
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class TreeSetThinningSensitivityRow:
    thinning_interval: int
    retained_tree_count: int
    retained_fraction: float
    rooted_topology_count: int
    shared_rooted_topology_count: int
    dominant_topology_frequency: float
    dominant_topology_delta: float
    robust_clade_count: int
    uncertain_clade_count: int
    conflicting_clade_count: int
    warnings: list[str]


@dataclass(slots=True)
class TreeSetThinningSensitivityReport:
    path: Path
    original_tree_count: int
    original_rooted_topology_count: int
    original_dominant_topology_frequency: float
    rows: list[TreeSetThinningSensitivityRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class ConsensusThresholdSensitivityRow:
    threshold: float
    informative_clade_count: int
    rooted_topology_id: str
    consensus_newick: str
    warnings: list[str]


@dataclass(slots=True)
class ConsensusThresholdSensitivityReport:
    path: Path
    tree_count: int
    rows: list[ConsensusThresholdSensitivityRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class TreeSetMaturityGateCheck:
    name: str
    satisfied: bool
    details: str


@dataclass(slots=True)
class TreeSetMaturityGateReport:
    path: Path
    decision: str
    checks: list[TreeSetMaturityGateCheck]
    warnings: list[str]
