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
