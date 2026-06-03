from __future__ import annotations

from .benchmark_dashboards import (
    build_large_alignment_scaling_benchmark_dashboard,
    build_large_tree_scaling_benchmark_dashboard,
    build_large_tree_set_scaling_benchmark_dashboard,
    build_memory_benchmark_dashboard,
    build_method_accuracy_dashboard,
    build_runtime_benchmark_dashboard,
    build_workflow_practical_limit_dashboard,
)
from .contracts import (
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
from .dataset_corpora import (
    build_broken_benchmark_corpus,
    build_clean_benchmark_corpus,
    build_messy_benchmark_corpus,
)
from .presentation import write_validation_corpus_json
from .regression_corpus import build_regression_dataset_corpus
from .reproducibility import validate_simulation_reproducibility
from .scientific_review import (
    build_method_limitation_registry,
    build_scientific_validation_report,
)

__all__ = [
    "build_broken_benchmark_corpus",
    "build_clean_benchmark_corpus",
    "build_large_alignment_scaling_benchmark_dashboard",
    "build_large_tree_scaling_benchmark_dashboard",
    "build_large_tree_set_scaling_benchmark_dashboard",
    "build_memory_benchmark_dashboard",
    "build_messy_benchmark_corpus",
    "build_method_accuracy_dashboard",
    "build_method_limitation_registry",
    "build_regression_dataset_corpus",
    "build_runtime_benchmark_dashboard",
    "build_scientific_validation_report",
    "build_workflow_practical_limit_dashboard",
    "BenchmarkCorpusReport",
    "BenchmarkDashboardRow",
    "CorpusDatasetCase",
    "CorpusDatasetCaseResult",
    "LargeAlignmentScalingBenchmarkDashboard",
    "LargeTreeScalingBenchmarkDashboard",
    "LargeTreeSetScalingBenchmarkDashboard",
    "MemoryBenchmarkDashboard",
    "MethodAccuracyDashboard",
    "MethodAccuracyRow",
    "MethodLimitationEntry",
    "MethodLimitationRegistry",
    "RegressionDatasetCaseResult",
    "RegressionDatasetCorpusReport",
    "RuntimeBenchmarkDashboard",
    "ScientificValidationClaim",
    "ScientificValidationReport",
    "SimulationReproducibilityCase",
    "SimulationReproducibilityReport",
    "WorkflowPracticalLimitDashboard",
    "validate_simulation_reproducibility",
    "write_validation_corpus_json",
]
