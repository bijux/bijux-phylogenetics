from __future__ import annotations

from bijux_phylogenetics.benchmark import (
    benchmark_alignment_site_scaling,
    benchmark_large_alignment_scaling,
    benchmark_large_tree_scaling,
    benchmark_large_tree_set_scaling,
    benchmark_tree_comparison,
    benchmark_tree_set_consensus,
    benchmark_tree_validation,
    benchmark_workflow_practical_limits,
)
from bijux_phylogenetics.validation.reference import (
    build_core_workflow_validation_report,
)

from .contracts import (
    BenchmarkDashboardRow,
    LargeAlignmentScalingBenchmarkDashboard,
    LargeTreeScalingBenchmarkDashboard,
    LargeTreeSetScalingBenchmarkDashboard,
    MemoryBenchmarkDashboard,
    MethodAccuracyDashboard,
    MethodAccuracyRow,
    RuntimeBenchmarkDashboard,
    WorkflowPracticalLimitDashboard,
)
from .dataset_corpora import (
    build_broken_benchmark_corpus,
    build_clean_benchmark_corpus,
    build_messy_benchmark_corpus,
)
from .regression_corpus import build_regression_dataset_corpus


def build_method_accuracy_dashboard(*, fixtures_root=None) -> MethodAccuracyDashboard:
    """Summarize validation accuracy, error counts, and coverage across benchmark surfaces."""
    core = build_core_workflow_validation_report(fixtures_root=fixtures_root)
    clean = build_clean_benchmark_corpus(fixtures_root=fixtures_root)
    broken = build_broken_benchmark_corpus(fixtures_root=fixtures_root)
    messy = build_messy_benchmark_corpus(fixtures_root=fixtures_root)
    regression = build_regression_dataset_corpus(fixtures_root=fixtures_root)

    def row(
        surface: str,
        passed_count: int,
        failed_count: int,
        coverage_count: int,
        bias_notes: list[str],
        error_notes: list[str],
    ) -> MethodAccuracyRow:
        total = max(coverage_count, 1)
        return MethodAccuracyRow(
            surface=surface,
            accuracy=round(passed_count / total, 15),
            passed_count=passed_count,
            failed_count=failed_count,
            coverage_count=coverage_count,
            bias_notes=bias_notes,
            error_notes=error_notes,
        )

    rows = [
        row(
            "level1-reference-validation",
            core.passed_fixture_count,
            core.failed_fixture_count,
            core.total_fixture_count,
            core.limitations,
            [case.fixture_name for case in core.failure_gallery if not case.passed],
        ),
        row(
            "clean-benchmark-corpus",
            clean.passed_case_count,
            clean.failed_case_count,
            clean.case_count,
            clean.limitations,
            [case.name for case in clean.cases if not case.passed],
        ),
        row(
            "broken-benchmark-corpus",
            broken.passed_case_count,
            broken.failed_case_count,
            broken.case_count,
            broken.limitations,
            [case.name for case in broken.cases if not case.passed],
        ),
        row(
            "messy-benchmark-corpus",
            messy.passed_case_count,
            messy.failed_case_count,
            messy.case_count,
            messy.limitations,
            [case.name for case in messy.cases if not case.passed],
        ),
        row(
            "regression-dataset-corpus",
            regression.passed_case_count,
            regression.failed_case_count,
            regression.case_count,
            regression.limitations,
            [case.name for case in regression.cases if not case.passed],
        ),
    ]
    return MethodAccuracyDashboard(
        goal_id=246,
        rows=rows,
        limitations=[
            "accuracy currently summarizes checked-in fixture and corpus pass rates; it does not yet replace external software comparison studies",
        ],
    )


def build_runtime_benchmark_dashboard(
    *, replicates: int = 1
) -> RuntimeBenchmarkDashboard:
    """Summarize runtime scaling across taxa, sites, tree counts, and posterior-like samples."""
    rows = [
        BenchmarkDashboardRow(
            workflow="tree-validation",
            scaling_axis="taxa",
            observations=benchmark_tree_validation(replicates=replicates).observations,
        ),
        BenchmarkDashboardRow(
            workflow="tree-comparison",
            scaling_axis="taxa",
            observations=benchmark_tree_comparison(replicates=replicates).observations,
        ),
        BenchmarkDashboardRow(
            workflow="alignment-diagnostics",
            scaling_axis="sites",
            observations=benchmark_alignment_site_scaling(
                replicates=replicates
            ).observations,
        ),
        BenchmarkDashboardRow(
            workflow="tree-set-consensus",
            scaling_axis="posterior_samples",
            observations=benchmark_tree_set_consensus(
                replicates=replicates
            ).observations,
        ),
    ]
    return RuntimeBenchmarkDashboard(
        goal_id=247,
        rows=rows,
        limitations=[
            "runtime summaries measure local benchmark fixtures and should be re-run on target hardware before operational promises are made",
        ],
    )


def build_memory_benchmark_dashboard(
    *, replicates: int = 1
) -> MemoryBenchmarkDashboard:
    """Summarize peak memory scaling across the main benchmark axes."""
    rows = [
        BenchmarkDashboardRow(
            workflow="tree-validation",
            scaling_axis="taxa",
            observations=benchmark_tree_validation(replicates=replicates).observations,
        ),
        BenchmarkDashboardRow(
            workflow="tree-comparison",
            scaling_axis="taxa",
            observations=benchmark_tree_comparison(replicates=replicates).observations,
        ),
        BenchmarkDashboardRow(
            workflow="alignment-diagnostics",
            scaling_axis="sites",
            observations=benchmark_alignment_site_scaling(
                replicates=replicates
            ).observations,
        ),
        BenchmarkDashboardRow(
            workflow="tree-set-consensus",
            scaling_axis="posterior_samples",
            observations=benchmark_tree_set_consensus(
                replicates=replicates
            ).observations,
        ),
    ]
    return MemoryBenchmarkDashboard(
        goal_id=248,
        rows=rows,
        limitations=[
            "memory summaries capture Python-side peak allocations during benchmark runs and do not model every external-engine workflow",
        ],
    )


def build_large_tree_scaling_benchmark_dashboard(
    *,
    replicates: int = 1,
    tip_counts: list[int] | None = None,
) -> LargeTreeScalingBenchmarkDashboard:
    """Summarize realistic large-tree scaling across validation and review workflows."""
    report = benchmark_large_tree_scaling(
        replicates=replicates,
        tip_counts=tip_counts,
    )
    return LargeTreeScalingBenchmarkDashboard(
        goal_id=221,
        workflows=report.workflows,
        limitations=report.limitations,
    )


def build_large_alignment_scaling_benchmark_dashboard(
    *,
    replicates: int = 1,
    size_classes: list[tuple[str, int, int]] | None = None,
) -> LargeAlignmentScalingBenchmarkDashboard:
    """Summarize realistic large-alignment scaling across review workflows."""
    report = benchmark_large_alignment_scaling(
        replicates=replicates,
        size_classes=size_classes,
    )
    return LargeAlignmentScalingBenchmarkDashboard(
        goal_id=222,
        workflows=report.workflows,
        limitations=report.limitations,
    )


def build_large_tree_set_scaling_benchmark_dashboard(
    *,
    replicates: int = 1,
    size_classes: list[tuple[str, int, int]] | None = None,
) -> LargeTreeSetScalingBenchmarkDashboard:
    """Summarize realistic large-tree-set scaling across posterior review workflows."""
    report = benchmark_large_tree_set_scaling(
        replicates=replicates,
        size_classes=size_classes,
    )
    return LargeTreeSetScalingBenchmarkDashboard(
        goal_id=223,
        workflows=report.workflows,
        limitations=report.limitations,
    )


def build_workflow_practical_limit_dashboard(
    *,
    replicates: int = 1,
    tree_tip_counts: list[int] | None = None,
    alignment_size_classes: list[tuple[str, int, int]] | None = None,
    tree_set_size_classes: list[tuple[str, int, int]] | None = None,
    stress_tiers: list[str] | None = None,
) -> WorkflowPracticalLimitDashboard:
    """Summarize the largest governed workflow classes currently exercised in benchmark and stress lanes."""
    report = benchmark_workflow_practical_limits(
        replicates=replicates,
        tree_tip_counts=tree_tip_counts,
        alignment_size_classes=alignment_size_classes,
        tree_set_size_classes=tree_set_size_classes,
        stress_tiers=stress_tiers,
    )
    return WorkflowPracticalLimitDashboard(
        goal_id=224,
        entries=report.entries,
        limitations=report.limitations,
    )
