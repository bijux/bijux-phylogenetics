from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.benchmark import (
    BenchmarkObservation,
    LargeAlignmentScalingWorkflowBenchmark,
    LargeTreeScalingWorkflowBenchmark,
    LargeTreeSetScalingWorkflowBenchmark,
    WorkflowPracticalLimitEntry,
)


@dataclass(frozen=True, slots=True)
class CorpusDatasetCase:
    """One checked-in dataset case used in a benchmark corpus."""

    name: str
    tree_path: Path
    metadata_path: Path
    traits_path: Path
    alignment_path: Path | None = None
    tip_dates_path: Path | None = None
    calibration_path: Path | None = None
    required_allowed_analyses: tuple[str, ...] = ()
    forbidden_blockers: tuple[str, ...] = ()
    allowed_warning_prefixes: tuple[str, ...] = ()


@dataclass(slots=True)
class CorpusDatasetCaseResult:
    """Observed result for one evaluated dataset case."""

    name: str
    passed: bool
    readiness_decision: str
    analysis_taxa: list[str]
    allowed_analyses: list[str]
    blocked_analyses: list[str]
    blockers: list[str]
    warnings: list[str]
    notes: list[str]
    observed_code: str | None = None


@dataclass(slots=True)
class BenchmarkCorpusReport:
    """Reviewer-facing summary for one benchmark corpus."""

    goal_id: int
    corpus: str
    passed: bool
    case_count: int
    passed_case_count: int
    failed_case_count: int
    cases: list[CorpusDatasetCaseResult]
    limitations: list[str]


@dataclass(slots=True)
class RegressionDatasetCaseResult:
    """Observed-versus-expected summary for one regression dataset case."""

    name: str
    passed: bool
    expected: dict[str, object]
    observed: dict[str, object]
    notes: list[str]


@dataclass(slots=True)
class RegressionDatasetCorpusReport:
    """Stable biological summary snapshots tracked across releases."""

    goal_id: int
    corpus: str
    passed: bool
    case_count: int
    passed_case_count: int
    failed_case_count: int
    cases: list[RegressionDatasetCaseResult]
    limitations: list[str]


@dataclass(slots=True)
class MethodAccuracyRow:
    """One validation surface summarized for accuracy, error, and coverage."""

    surface: str
    accuracy: float
    passed_count: int
    failed_count: int
    coverage_count: int
    bias_notes: list[str]
    error_notes: list[str]


@dataclass(slots=True)
class MethodAccuracyDashboard:
    """Goal 246 dashboard across the main checked-in validation surfaces."""

    goal_id: int
    rows: list[MethodAccuracyRow]
    limitations: list[str]


@dataclass(slots=True)
class BenchmarkDashboardRow:
    """One workflow scaling curve in runtime or memory dashboards."""

    workflow: str
    scaling_axis: str
    observations: list[BenchmarkObservation]


@dataclass(slots=True)
class RuntimeBenchmarkDashboard:
    """Goal 247 runtime scaling summary across major benchmark axes."""

    goal_id: int
    rows: list[BenchmarkDashboardRow]
    limitations: list[str]


@dataclass(slots=True)
class MemoryBenchmarkDashboard:
    """Goal 248 memory scaling summary across major benchmark axes."""

    goal_id: int
    rows: list[BenchmarkDashboardRow]
    limitations: list[str]


@dataclass(slots=True)
class LargeTreeScalingBenchmarkDashboard:
    """Goal 221 scaling summary for large-tree reviewer workflows."""

    goal_id: int
    workflows: list[LargeTreeScalingWorkflowBenchmark]
    limitations: list[str]


@dataclass(slots=True)
class LargeAlignmentScalingBenchmarkDashboard:
    """Goal 222 scaling summary for large-alignment reviewer workflows."""

    goal_id: int
    workflows: list[LargeAlignmentScalingWorkflowBenchmark]
    limitations: list[str]


@dataclass(slots=True)
class LargeTreeSetScalingBenchmarkDashboard:
    """Goal 223 scaling summary for large-tree-set reviewer workflows."""

    goal_id: int
    workflows: list[LargeTreeSetScalingWorkflowBenchmark]
    limitations: list[str]


@dataclass(slots=True)
class WorkflowPracticalLimitDashboard:
    """Goal 224 tested practical limits for governed workflow surfaces."""

    goal_id: int
    entries: list[WorkflowPracticalLimitEntry]
    limitations: list[str]


@dataclass(slots=True)
class MethodLimitationEntry:
    """One method with its assumptions, invalid inputs, and current trust status."""

    method: str
    status: str
    validated_by: list[str]
    assumptions: list[str]
    invalid_inputs: list[str]
    limitations: list[str]


@dataclass(slots=True)
class MethodLimitationRegistry:
    """Goal 250 registry of method assumptions and current trust boundaries."""

    goal_id: int
    entries: list[MethodLimitationEntry]
    limitations: list[str]


@dataclass(slots=True)
class ScientificValidationClaim:
    """One reviewer-facing claim bucketed by validation confidence."""

    status: str
    claim: str
    evidence: list[str]


@dataclass(slots=True)
class ScientificValidationReport:
    """Goal 249 summary of validated, unvalidated, experimental, and unsafe claims."""

    goal_id: int
    claims: list[ScientificValidationClaim]
    limitations: list[str]


@dataclass(slots=True)
class SimulationReproducibilityCase:
    """One seeded simulation surface checked for exact repeatability."""

    surface: str
    passed: bool
    digest: str
    notes: list[str]


@dataclass(slots=True)
class SimulationReproducibilityReport:
    """Goal 251 seeded simulation reproducibility validation."""

    goal_id: int
    passed: bool
    cases: list[SimulationReproducibilityCase]
    limitations: list[str]


@dataclass(slots=True)
class InvariantMixtureModelCoverageCase:
    """One eligible finite-state model checked under an active invariant mixture."""

    model_name: str
    surface_model_name: str
    invariant_proportion: float
    site_count: int
    pattern_count: int
    fixed_rate_log_likelihood: float
    invariant_log_likelihood: float
    likelihood_delta: float
    invariant_component_likelihood: float
    variable_component_likelihood: float
    likelihood_changed: bool
    invariant_component_active: bool
    passed: bool
    notes: list[str]


@dataclass(slots=True)
class InvariantMixtureFamilyCoverageEntry:
    """One eligible finite-state family summarized across its governed model cases."""

    family_name: str
    supported_model_names: list[str]
    covered_model_names: list[str]
    model_count: int
    passed_model_count: int
    passed: bool
    model_cases: list[InvariantMixtureModelCoverageCase]
    notes: list[str]


@dataclass(slots=True)
class InvariantMixtureFamilyCoverageReport:
    """Reviewer-facing proof that invariant mass is active across eligible models."""

    goal_id: int
    passed: bool
    family_count: int
    model_count: int
    passed_family_count: int
    passed_model_count: int
    families: list[InvariantMixtureFamilyCoverageEntry]
    limitations: list[str]


@dataclass(slots=True)
class GammaInvariantMixtureModelCoverageCase:
    """One eligible finite-state model checked under one active combined gamma-plus-invariant mixture."""

    model_name: str
    surface_model_name: str
    alpha: float
    category_count: int
    invariant_proportion: float
    site_count: int
    pattern_count: int
    gamma_log_likelihood: float
    gamma_invariant_log_likelihood: float
    likelihood_delta: float
    invariant_component_likelihood: float
    variable_component_likelihood: float
    category_likelihood_count: int
    likelihood_changed: bool
    invariant_component_active: bool
    variable_component_active: bool
    boundary_warning_count: int
    hit_lower_invariant_proportion_boundary: bool
    hit_upper_invariant_proportion_boundary: bool
    boundary_warnings: list[str]
    passed: bool
    notes: list[str]


@dataclass(slots=True)
class GammaInvariantMixtureFamilyCoverageEntry:
    """One eligible finite-state family summarized across its governed gamma-plus-invariant cases."""

    family_name: str
    supported_model_names: list[str]
    covered_model_names: list[str]
    model_count: int
    passed_model_count: int
    passed: bool
    model_cases: list[GammaInvariantMixtureModelCoverageCase]
    notes: list[str]


@dataclass(slots=True)
class GammaInvariantMixtureFamilyCoverageReport:
    """Reviewer-facing proof that combined gamma and invariant mass stay active across eligible models."""

    goal_id: int
    passed: bool
    family_count: int
    model_count: int
    passed_family_count: int
    passed_model_count: int
    families: list[GammaInvariantMixtureFamilyCoverageEntry]
    limitations: list[str]


@dataclass(slots=True)
class SharedContinuousFixtureCatalogEntry:
    """One governed continuous-trait fixture used by live `geiger` parity."""

    fixture_id: str
    tree_fixture_id: str
    trait_table_fixture_id: str
    trait_name: str
    standard_error_trait_name: str | None
    taxon_column: str
    taxon_count: int | None
    trait_row_count: int
    supported_model_names: list[str]
    covered_model_names: list[str]
    live_case_count: int
    validation_expectation: str
    geiger_reference_expectation: str
    feature_tags: list[str]
    notes: str


@dataclass(slots=True)
class SharedContinuousFixtureCatalogReport:
    """Goal 252 governed inventory for live `geiger::fitContinuous` fixtures."""

    goal_id: int
    fixture_count: int
    live_case_count: int
    fixtures: list[SharedContinuousFixtureCatalogEntry]
    limitations: list[str]


@dataclass(slots=True)
class SharedDiscreteFixtureCatalogEntry:
    """One governed discrete-trait fixture used by live `geiger` parity."""

    fixture_id: str
    tree_fixture_id: str
    trait_table_fixture_id: str
    trait_name: str
    taxon_column: str
    trait_kind: str
    taxon_count: int | None
    trait_row_count: int
    supported_model_names: list[str]
    covered_model_names: list[str]
    live_case_count: int
    validation_expectation: str
    geiger_reference_expectation: str
    feature_tags: list[str]
    transition_matrix_metadata: dict[str, object] | None
    notes: str


@dataclass(slots=True)
class SharedDiscreteFixtureCatalogReport:
    """Goal 253 governed inventory for live `geiger::fitDiscrete` fixtures."""

    goal_id: int
    fixture_count: int
    live_case_count: int
    fixtures: list[SharedDiscreteFixtureCatalogEntry]
    limitations: list[str]


@dataclass(slots=True)
class SubstitutionParameterPriorCoverageEntry:
    """One governed substitution-prior target checked for score sensitivity."""

    target_name: str
    surface_name: str
    prior_family: str
    contribution_count: int
    covered_names: list[str]
    baseline_log_prior: float
    comparison_log_prior: float
    log_prior_delta: float
    changed: bool
    passed: bool
    notes: list[str]


@dataclass(slots=True)
class SubstitutionParameterPriorCoverageReport:
    """Proof that governed substitution-prior targets change prior scores."""

    goal_id: int
    passed: bool
    target_count: int
    passed_target_count: int
    entries: list[SubstitutionParameterPriorCoverageEntry]
    limitations: list[str]
