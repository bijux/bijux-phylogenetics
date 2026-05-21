from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.validation.validation_corpus.benchmark_dashboards import (
    build_large_alignment_scaling_benchmark_dashboard,
    build_large_tree_scaling_benchmark_dashboard,
    build_large_tree_set_scaling_benchmark_dashboard,
    build_memory_benchmark_dashboard,
    build_method_accuracy_dashboard,
    build_runtime_benchmark_dashboard,
    build_workflow_practical_limit_dashboard,
)
from bijux_phylogenetics.validation.validation_corpus.dataset_corpora import (
    build_broken_benchmark_corpus,
    build_clean_benchmark_corpus,
    build_messy_benchmark_corpus,
)
from bijux_phylogenetics.validation.validation_corpus.presentation import (
    write_validation_corpus_json as _write_validation_corpus_json,
)
from bijux_phylogenetics.validation.validation_corpus.regression_corpus import (
    build_regression_dataset_corpus,
)
from bijux_phylogenetics.validation.validation_corpus.reproducibility import (
    validate_simulation_reproducibility as _validate_simulation_reproducibility,
)
from bijux_phylogenetics.validation.validation_corpus.scientific_review import (
    build_method_limitation_registry as _build_method_limitation_registry,
    build_scientific_validation_report as _build_scientific_validation_report,
)
from bijux_phylogenetics.validation.validation_corpus.contracts import (
    BenchmarkCorpusReport,
    BenchmarkDashboardRow,
    CorpusDatasetCase,
    CorpusDatasetCaseResult,
    LargeAlignmentScalingBenchmarkDashboard,
    LargeTreeScalingBenchmarkDashboard,
    LargeTreeSetScalingBenchmarkDashboard,
    MemoryBenchmarkDashboard,
    MethodAccuracyDashboard,
    MethodAccuracyRow,
    MethodLimitationEntry,
    MethodLimitationRegistry,
    RegressionDatasetCaseResult,
    RegressionDatasetCorpusReport,
    RuntimeBenchmarkDashboard,
    ScientificValidationClaim,
    ScientificValidationReport,
    SimulationReproducibilityCase,
    SimulationReproducibilityReport,
    WorkflowPracticalLimitDashboard,
)

_CORPUS_EXPORT_TYPES = (
    BenchmarkCorpusReport,
    BenchmarkDashboardRow,
    CorpusDatasetCase,
    CorpusDatasetCaseResult,
    LargeAlignmentScalingBenchmarkDashboard,
    LargeTreeScalingBenchmarkDashboard,
    LargeTreeSetScalingBenchmarkDashboard,
    MemoryBenchmarkDashboard,
    MethodAccuracyDashboard,
    MethodAccuracyRow,
    MethodLimitationEntry,
    MethodLimitationRegistry,
    RegressionDatasetCaseResult,
    RegressionDatasetCorpusReport,
    RuntimeBenchmarkDashboard,
    ScientificValidationClaim,
    ScientificValidationReport,
    SimulationReproducibilityCase,
    SimulationReproducibilityReport,
    WorkflowPracticalLimitDashboard,
)

_CORPUS_EXPORT_BUILDERS = (
    build_broken_benchmark_corpus,
    build_clean_benchmark_corpus,
    build_messy_benchmark_corpus,
    build_regression_dataset_corpus,
    build_method_accuracy_dashboard,
    build_runtime_benchmark_dashboard,
    build_memory_benchmark_dashboard,
    build_large_tree_scaling_benchmark_dashboard,
    build_large_alignment_scaling_benchmark_dashboard,
    build_large_tree_set_scaling_benchmark_dashboard,
    build_workflow_practical_limit_dashboard,
)

def build_method_limitation_registry() -> MethodLimitationRegistry:
    """Enumerate major method families with explicit assumptions and trust boundaries."""
    return _build_method_limitation_registry()


def build_scientific_validation_report(
    *, fixtures_root: Path | None = None
) -> ScientificValidationReport:
    """Separate validated, unvalidated, experimental, and unsafe claims for reviewers."""
    return _build_scientific_validation_report(fixtures_root=fixtures_root)


def write_validation_corpus_json(path: Path, report: object) -> Path:
    """Write a validation-corpus or dashboard report as deterministic JSON."""
    return _write_validation_corpus_json(path, report)


def validate_simulation_reproducibility(
    *, fixtures_root: Path | None = None
) -> SimulationReproducibilityReport:
    """Verify that repeated simulations with the same seed produce identical structured results."""
    return _validate_simulation_reproducibility(fixtures_root=fixtures_root)
