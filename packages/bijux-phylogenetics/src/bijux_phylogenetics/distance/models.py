from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.compare.topology import (
    BranchLengthComparisonReport,
    TreeComparisonReport,
)

DistanceModel = str
GapHandlingMode = str
AmbiguityPolicy = str
MissingDistancePolicy = str


@dataclass(frozen=True, slots=True)
class GeneticDistanceModelParameters:
    """Alignment-wide parameters used by composition-aware DNA distance models."""

    informative_base_count: int
    base_frequency_a: float
    base_frequency_c: float
    base_frequency_g: float
    base_frequency_t: float
    purine_frequency: float
    pyrimidine_frequency: float
    f81_limit: float
    tn93_ag_coefficient: float | None
    tn93_ct_coefficient: float | None
    tn93_transversion_coefficient: float | None


@dataclass(frozen=True, slots=True)
class PairwiseGeneticDistance:
    """One pairwise genetic distance entry for an aligned dataset."""

    left_identifier: str
    right_identifier: str
    distance: float | None
    comparable_sites: int
    mismatch_sites: float
    transition_sites: float
    ag_transition_sites: float
    ct_transition_sites: float
    transversion_sites: float
    ambiguity_sites: int
    skipped_sites: int
    saturated: bool
    saturation_reason: str | None


@dataclass(slots=True)
class GeneticDistanceMatrix:
    """Deterministic pairwise genetic distance matrix for one alignment."""

    path: Path
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    inferred_alphabet: str
    alignment_length: int
    identifiers: list[str]
    model_parameters: GeneticDistanceModelParameters | None
    warnings: list[str]
    pairs: list[PairwiseGeneticDistance]


@dataclass(slots=True)
class DistanceTreeBuildReport:
    """Explicit report for a distance-based tree build."""

    alignment_path: Path
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    method: str
    method_policy: DistanceTreeMethodPolicy
    taxon_count: int
    pair_count: int
    assumptions: DistanceMethodAssumptionReport
    missing_distance_policy_report: MissingDistancePolicyReport


@dataclass(frozen=True, slots=True)
class MinimumEvolutionBranchFit:
    """One fitted branch length for a fixed-topology minimum-evolution score."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    fitted_branch_length: float


@dataclass(slots=True)
class MinimumEvolutionScoreReport:
    """Minimum-evolution score and fitted branch lengths for one fixed topology."""

    taxa: list[str]
    pair_count: int
    branch_count: int
    minimum_evolution_score: float
    total_fitted_branch_length: float
    negative_branch_count: int
    branch_fits: list[MinimumEvolutionBranchFit]


@dataclass(frozen=True, slots=True)
class FitchMargoliashBranchFit:
    """One fitted branch length from a Fitch-Margoliash weighted least-squares fit."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    fitted_branch_length: float


@dataclass(slots=True)
class FitchMargoliashFitReport:
    """Weighted least-squares branch fit and RSS summary for one fixed topology."""

    taxa: list[str]
    pair_count: int
    branch_count: int
    weighting_power: float
    residual_sum_squares: float
    weighted_residual_sum_squares: float
    matrix_rank: int
    negative_branch_count: int
    branch_fits: list[FitchMargoliashBranchFit]


@dataclass(frozen=True, slots=True)
class OrdinaryLeastSquaresBranchFit:
    """One fitted branch length from an ordinary least-squares fixed-topology fit."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    fitted_branch_length: float


@dataclass(slots=True)
class OrdinaryLeastSquaresFitReport:
    """Unweighted fixed-topology least-squares fit with residual matrix diagnostics."""

    taxa: list[str]
    pair_count: int
    branch_count: int
    residual_sum_squares: float
    matrix_rank: int
    condition_number: float
    negative_branch_count: int
    fitted_distance_matrix: list[list[float]]
    residual_matrix: list[list[float]]
    branch_fits: list[OrdinaryLeastSquaresBranchFit]


@dataclass(frozen=True, slots=True)
class NonnegativeLeastSquaresBranchFit:
    """One fitted branch length from a nonnegative least-squares fixed-topology fit."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    fitted_branch_length: float


@dataclass(frozen=True, slots=True)
class NonnegativeLeastSquaresActiveConstraint:
    """One branch pinned to zero by the active nonnegative constraint set."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]


@dataclass(slots=True)
class NonnegativeLeastSquaresFitReport:
    """Fixed-topology NNLS fit with explicit active zero constraints."""

    taxa: list[str]
    pair_count: int
    branch_count: int
    residual_sum_squares: float
    condition_number: float
    active_constraint_count: int
    fitted_distance_matrix: list[list[float]]
    residual_matrix: list[list[float]]
    branch_fits: list[NonnegativeLeastSquaresBranchFit]
    active_constraints: list[NonnegativeLeastSquaresActiveConstraint]


@dataclass(frozen=True, slots=True)
class PatristicResidualRow:
    """One observed-versus-tree distance residual for a unique taxon pair."""

    left_identifier: str
    right_identifier: str
    observed_distance: float
    fitted_distance: float
    residual: float
    absolute_residual: float
    rank: int


@dataclass(slots=True)
class PatristicResidualDiagnosticsReport:
    """Observed-versus-tree patristic residual diagnostics for one distance matrix."""

    matrix_path: Path | None
    tree_path: Path | None
    taxa: list[str]
    pair_count: int
    residual_sum_squares: float
    max_absolute_residual: float
    rows: list[PatristicResidualRow]


@dataclass(frozen=True, slots=True)
class DistanceTaxonInfluenceRow:
    """One leave-one-taxon-out distance-tree influence observation."""

    taxon: str
    retained_taxa: list[str]
    raw_missing_pair_count: int
    baseline_residual_sum_squares: float
    leave_one_out_residual_sum_squares: float
    residual_sum_squares_improvement: float
    baseline_rooted_robinson_foulds_distance: int
    leave_one_out_rooted_robinson_foulds_distance: int
    rooted_robinson_foulds_improvement: int
    baseline_rooted_normalized_robinson_foulds: float
    leave_one_out_rooted_normalized_robinson_foulds: float
    rooted_normalized_robinson_foulds_improvement: float
    topology_improved: bool
    residual_improved: bool
    influence_rank: int


@dataclass(slots=True)
class DistanceTaxonInfluenceReport:
    """Rank leave-one-out taxon effects on one distance-tree fit and reference match."""

    source_path: Path
    source_kind: str
    reference_tree_path: Path
    method: str
    missing_distance_policy: MissingDistancePolicy
    taxa: list[str]
    baseline_residual_sum_squares: float
    baseline_rooted_robinson_foulds_distance: int
    baseline_rooted_normalized_robinson_foulds: float
    rows: list[DistanceTaxonInfluenceRow]


@dataclass(frozen=True, slots=True)
class DistanceTaxonJackknifeRow:
    """One leave-one-taxon-out rebuilt-tree jackknife observation."""

    removed_taxon: str
    retained_taxa: list[str]
    pruned_baseline_tree_newick: str
    rebuilt_tree_newick: str
    pruned_baseline_residual_sum_squares: float
    rebuilt_residual_sum_squares: float
    residual_sum_squares_change: float
    rooted_robinson_foulds_distance: int
    rooted_normalized_robinson_foulds: float
    reference_only_clades: list[str]
    rebuilt_only_clades: list[str]
    affected_clades: list[str]
    topology_changed: bool


@dataclass(slots=True)
class DistanceTaxonJackknifeReport:
    """Leave-one-taxon-out topology jackknife over one imported distance matrix."""

    source_path: Path
    source_kind: str
    method: str
    missing_distance_policy: MissingDistancePolicy
    taxa: list[str]
    baseline_tree_newick: str
    baseline_residual_sum_squares: float
    rows: list[DistanceTaxonJackknifeRow]


@dataclass(frozen=True, slots=True)
class DistanceMethodComparisonRow:
    """One distance-tree method scored against the same imported distance matrix."""

    method: str
    tree_newick: str
    patristic_residual_sum_squares: float
    balanced_minimum_evolution_score: float
    ordinary_least_squares_residual_sum_squares: float
    ordinary_least_squares_negative_branch_count: int
    assumption_warnings: list[str]


@dataclass(frozen=True, slots=True)
class DistanceMethodRfRow:
    """One pairwise rooted RF comparison across two distance-tree methods."""

    left_method: str
    right_method: str
    rooted_robinson_foulds_distance: int
    rooted_normalized_robinson_foulds: float


@dataclass(frozen=True, slots=True)
class DistanceMethodWarningRow:
    """One explicit matrix-level or method-level assumption warning."""

    warning_rank: int
    scope: str
    method: str | None
    warning: str


@dataclass(slots=True)
class DistanceMethodComparisonReport:
    """Compare owned distance-tree methods and fixed-topology scores on one matrix."""

    source_path: Path
    source_kind: str
    missing_distance_policy: MissingDistancePolicy
    taxa: list[str]
    compared_methods: list[str]
    rows: list[DistanceMethodComparisonRow]
    rf_rows: list[DistanceMethodRfRow]
    warning_rows: list[DistanceMethodWarningRow]


@dataclass(frozen=True, slots=True)
class DistanceUltrametricityViolation:
    """One taxon triple whose pairwise distances violate the three-point condition."""

    left_identifier: str
    middle_identifier: str
    right_identifier: str
    left_middle_distance: float
    left_right_distance: float
    middle_right_distance: float
    second_largest_distance: float
    largest_distance: float
    violation: float


@dataclass(slots=True)
class DistanceUltrametricityDiagnosticsReport:
    """Three-point ultrametricity diagnostics for one pairwise distance matrix."""

    source_path: Path | None
    source_kind: str
    taxon_count: int
    defined_pair_count: int
    tested_triple_count: int
    skipped_triple_count: int
    tolerance: float
    ultrametric: bool
    max_violation: float
    violating_triples: list[DistanceUltrametricityViolation]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class DistanceFourPointViolation:
    """One quartet whose pairwise distances violate the four-point condition."""

    first_identifier: str
    second_identifier: str
    third_identifier: str
    fourth_identifier: str
    split_ab_cd_sum: float
    split_ac_bd_sum: float
    split_ad_bc_sum: float
    best_split: str
    violation_magnitude: float


@dataclass(slots=True)
class DistanceAdditivityDiagnosticsReport:
    """Four-point additivity diagnostics for one pairwise distance matrix."""

    source_path: Path | None
    source_kind: str
    taxon_count: int
    defined_pair_count: int
    tested_quartet_count: int
    skipped_quartet_count: int
    tolerance: float
    additive: bool
    max_violation: float
    violating_quartets: list[DistanceFourPointViolation]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class BalancedMinimumEvolutionNniTraceRow:
    """One deterministic event row in a rooted BME NNI search trace."""

    event_index: int
    event_kind: str
    iteration: int
    score_before: float | None
    score_after: float
    score_delta: float | None
    tree_before_newick: str | None
    tree_after_newick: str
    pivot_branch_id: str | None
    sibling_clade_id: str | None
    exchanged_clade_id: str | None
    stopping_reason: str | None


@dataclass(slots=True)
class BalancedMinimumEvolutionNniSearchReport:
    """Complete rooted NNI hill-climb report over one NJ or BIONJ starting tree."""

    algorithm: str
    matrix_path: Path | None
    start_method: str
    taxon_count: int
    pair_count: int
    start_tree_newick: str
    start_score: float
    final_tree_newick: str
    final_score: float
    accepted_move_count: int
    evaluated_neighbor_count: int
    stopping_reason: str
    trace_rows: list[BalancedMinimumEvolutionNniTraceRow]


@dataclass(slots=True)
class DistanceTreeTopologyComparison:
    """Topology comparison between NJ and UPGMA trees built from one alignment."""

    alignment_path: Path
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    shared_taxa: list[str]
    nj_informative_clades: int
    upgma_informative_clades: int
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    topology_equal: bool
    same_unrooted_topology: bool
    same_taxa_different_rooting: bool


@dataclass(frozen=True, slots=True)
class UPGMAUltrametricViolation:
    """One taxon triple whose distances are inconsistent with an ultrametric clock."""

    left_identifier: str
    middle_identifier: str
    right_identifier: str
    smallest_distance: float
    middle_distance: float
    largest_distance: float
    deviation: float


@dataclass(slots=True)
class DistanceMethodAssumptionReport:
    """Audit whether an alignment or matrix respects core distance-tree assumptions."""

    source_path: Path
    source_kind: str
    taxon_count: int
    pair_count: int
    nj_assumptions: list[str]
    upgma_assumptions: list[str]
    ultrametric_compatible: bool
    ultrametric_tolerance: float
    upgma_ultrametric_violations: list[UPGMAUltrametricViolation]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class ImportedDistanceEntry:
    """One directional entry from an imported long-form distance matrix table."""

    left_identifier: str
    right_identifier: str
    distance: float
    comparable_sites: int | None


@dataclass(frozen=True, slots=True)
class DistanceMatrixAsymmetry:
    """Two directional entries that disagree numerically."""

    left_identifier: str
    right_identifier: str
    left_to_right_distance: float
    right_to_left_distance: float


@dataclass(frozen=True, slots=True)
class NonMetricDistanceObservation:
    """One triangle-inequality violation within an imported distance matrix."""

    left_identifier: str
    middle_identifier: str
    right_identifier: str
    direct_distance: float
    indirect_distance: float


@dataclass(slots=True)
class ImportedDistanceMatrixReport:
    """Validation report for an imported long-form distance matrix table."""

    path: Path
    identifiers: list[str]
    pair_count: int
    complete: bool
    zero_diagonal: bool
    symmetric: bool
    nonnegative: bool
    missing_pairs: list[str]
    diagonal_problems: list[str]
    negative_distance_pairs: list[str]
    asymmetric_pairs: list[DistanceMatrixAsymmetry]
    nonmetric_observations: list[NonMetricDistanceObservation]
    warnings: list[str]


@dataclass(slots=True)
class ImportedDistanceTreeBuildReport:
    """Explicit report for building a tree from an imported distance matrix."""

    matrix_path: Path
    method: str
    method_policy: DistanceTreeMethodPolicy
    taxon_count: int
    pair_count: int
    assumptions: DistanceMethodAssumptionReport
    missing_distance_policy_report: MissingDistancePolicyReport


@dataclass(slots=True)
class ImportedDistanceMatrixQualityReport:
    """Diagnostics over an imported long-form distance matrix."""

    validation: ImportedDistanceMatrixReport
    saturated_pairs: list[SaturatedDistancePair]
    low_information_pairs: list[LowInformationPair]
    low_information_pair_cutoff: int | None
    saturation_audit_scale: str
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class DistanceReferenceObservation:
    """One reference example used to verify core distance calculations."""

    case: str
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    left_identifier: str
    right_identifier: str
    expected_distance: float | None
    observed_distance: float | None
    comparable_sites: int
    expected_ambiguity_sites: int | None
    observed_ambiguity_sites: int
    passed: bool


@dataclass(slots=True)
class DistanceReferenceValidationReport:
    """Validation of built-in reference distance examples."""

    observations: list[DistanceReferenceObservation]
    tree_observations: list[DistanceTreeReferenceObservation]
    all_passed: bool


@dataclass(frozen=True, slots=True)
class DistanceTreeReferenceObservation:
    """One reference example used to validate distance-tree clustering."""

    case: str
    method: str
    matrix_path: Path
    expected_clades: list[str]
    observed_clades: list[str]
    passed: bool


@dataclass(frozen=True, slots=True)
class SaturatedDistancePair:
    """One pair that reaches an undefined or unreliable correction regime."""

    left_identifier: str
    right_identifier: str
    distance: float | None
    comparable_sites: int
    reason: str


@dataclass(frozen=True, slots=True)
class MissingDistanceImputation:
    """One missing pair imputed under an explicit distance-policy rule."""

    left_identifier: str
    right_identifier: str
    imputed_distance: float
    policy: MissingDistancePolicy
    rationale: str


@dataclass(slots=True)
class MissingDistancePolicyReport:
    """Resolution report for one explicit missing-distance policy decision."""

    policy: MissingDistancePolicy
    taxon_count: int
    requested_pair_count: int
    missing_pairs: list[str]
    imputed_rows: list[MissingDistanceImputation]
    unresolved_pairs: list[str]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class DistanceSaturationWarning:
    """One pair-level saturation warning reviewed before tree inference."""

    left_identifier: str
    right_identifier: str
    distance: float | None
    comparable_sites: int
    warning_kind: str
    reason: str
    blocks_tree_inference: bool


@dataclass(slots=True)
class DistanceSaturationDiagnosticsReport:
    """Pair-level saturation diagnostics for one computed distance matrix."""

    alignment_path: Path
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    taxon_count: int
    pair_count: int
    blocking_warning_count: int
    warning_rows: list[DistanceSaturationWarning]
    warnings: list[str]
    blocks_tree_inference: bool


@dataclass(frozen=True, slots=True)
class DistanceOutlierPair:
    """One pair whose distance is unusually large relative to the dataset."""

    left_identifier: str
    right_identifier: str
    distance: float
    note: str


@dataclass(frozen=True, slots=True)
class LowInformationPair:
    """One pair with too few comparable sites for robust interpretation."""

    left_identifier: str
    right_identifier: str
    comparable_sites: int
    note: str


@dataclass(slots=True)
class DistanceMethodAssessment:
    """Decision about whether the computed matrix is suitable for distance methods."""

    decision: str
    reasons: list[str]


@dataclass(frozen=True, slots=True)
class DistanceTreeMethodPolicy:
    """Stable support policy for one distance-tree method surface."""

    method: str
    supported: bool
    reference_surface: str | None
    support_scope: str
    summary: str
    limitations: list[str]


@dataclass(slots=True)
class DistanceMatrixQualityReport:
    """Diagnostics over a computed distance matrix."""

    path: Path
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    inferred_alphabet: str
    taxon_count: int
    pair_count: int
    saturated_pairs: list[SaturatedDistancePair]
    high_distance_outliers: list[DistanceOutlierPair]
    low_information_pairs: list[LowInformationPair]
    assumptions: DistanceMethodAssumptionReport
    warnings: list[str]
    method_assessment: DistanceMethodAssessment


@dataclass(frozen=True, slots=True)
class DistanceBootstrapSupportRow:
    """One clade support row across bootstrap replicate trees."""

    clade: str
    tree_count: int
    frequency: float


@dataclass(frozen=True, slots=True)
class DistanceBootstrapReplicateRow:
    """One bootstrap replicate draw with its rebuilt distance-tree outcome."""

    replicate_index: int
    sampled_site_indices: list[int]
    tree_newick: str


@dataclass(slots=True)
class DistanceBootstrapReport:
    """Bootstrap summary for a distance-based tree-building workflow."""

    alignment_path: Path
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    method: str
    replicates: int
    seed: int
    tree_count: int
    consensus_newick: str
    replicate_rows: list[DistanceBootstrapReplicateRow]
    support: list[DistanceBootstrapSupportRow]


@dataclass(slots=True)
class DistanceBootstrapSupportSummary:
    """Reviewer-facing summary over bootstrap clade support frequencies."""

    alignment_path: Path
    method: str
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    replicates: int
    clade_count: int
    minimum_frequency: float | None
    maximum_frequency: float | None
    median_frequency: float | None
    weak_clade_count: int
    warnings: list[str]


@dataclass(slots=True)
class DistanceTreeReferenceComparisonReport:
    """Compare one built distance tree against a reviewer-supplied reference tree."""

    alignment_path: Path
    reference_tree_path: Path
    method: str
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    topology: TreeComparisonReport
    branch_lengths: BranchLengthComparisonReport
    warnings: list[str]


@dataclass(slots=True)
class DistanceModelComparisonRow:
    """One supported distance model summarized over the same alignment."""

    model: DistanceModel
    defined_pair_count: int
    saturated_pair_count: int
    low_information_pair_count: int
    mean_distance: float | None
    maximum_distance: float | None
    decision: str
    reasons: list[str]


@dataclass(slots=True)
class DistanceModelComparisonReport:
    """Comparison of all supported distance models for one alignment."""

    alignment_path: Path
    inferred_alphabet: str
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    rows: list[DistanceModelComparisonRow]
    warnings: list[str]


@dataclass(slots=True)
class DistanceGapPolicyDeltaRow:
    """One taxon pair whose distance changed across gap-handling policies."""

    left_identifier: str
    right_identifier: str
    pairwise_distance: float | None
    complete_distance: float | None
    pairwise_comparable_sites: int
    complete_comparable_sites: int
    distance_delta: float | None
    comparable_site_delta: int


@dataclass(slots=True)
class DistanceGapPolicySensitivityReport:
    """Summarize how pairwise versus complete deletion changes the same analysis."""

    alignment_path: Path
    model: DistanceModel
    ambiguity_policy: AmbiguityPolicy
    changed_pair_count: int
    pair_count: int
    rows: list[DistanceGapPolicyDeltaRow]
    warnings: list[str]


@dataclass(slots=True)
class DistanceMethodMaturityCheck:
    """One explicit maturity criterion for distance-analysis surfaces."""

    name: str
    satisfied: bool
    details: str


@dataclass(slots=True)
class DistanceMethodMaturityGateReport:
    """High-level maturity gate over one distance-analysis workflow."""

    alignment_path: Path
    method: str
    method_policy: DistanceTreeMethodPolicy
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    decision: str
    checks: list[DistanceMethodMaturityCheck]
    warnings: list[str]


@dataclass(slots=True)
class DistanceMethodReport:
    """Structured machine-readable report for one distance-analysis workflow."""

    alignment_path: Path
    method: str
    method_policy: DistanceTreeMethodPolicy
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    matrix: GeneticDistanceMatrix
    quality: DistanceMatrixQualityReport
    assumptions: DistanceMethodAssumptionReport
    reference_validation: DistanceReferenceValidationReport
    built_tree_newick: str
    alternative_tree_newick: str
    topology_comparison: DistanceTreeTopologyComparison
    bootstrap_summary: DistanceBootstrapSupportSummary
    model_comparison: DistanceModelComparisonReport
    gap_policy_sensitivity: DistanceGapPolicySensitivityReport
    maturity_gate: DistanceMethodMaturityGateReport


@dataclass(slots=True)
class DistanceReproducibilityBundleReport:
    """Reproducibility bundle written for one distance-analysis run."""

    out_dir: Path
    alignment_path: Path
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    method: str
    replicates: int
    files: list[Path]
