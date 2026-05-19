from __future__ import annotations

from pathlib import Path
import tempfile

from .fixtures import (
    comparative_stress_payload as _comparative_stress_payload,
    large_alignment_stress_payload as _large_alignment_stress_payload,
    resolve_stress_tier_config as _resolve_stress_tier_config,
    supermatrix_stress_payload as _supermatrix_stress_payload,
    table_generation_stress_payload as _table_generation_stress_payload,
    tree_set_stress_payload as _tree_set_stress_payload,
)
from .measurement import (
    max_peak_memory_bytes as _max_peak_memory_bytes,
    max_runtime_seconds as _max_runtime_seconds,
    measure_stress_workload as _measure_stress_workload,
)
from .models import (
    AlignmentDiagnosticsBenchmarkReport as AlignmentDiagnosticsBenchmarkReport,
    AlignmentSiteBenchmarkReport as AlignmentSiteBenchmarkReport,
    BenchmarkObservation as BenchmarkObservation,
    LargeAlignmentScalingBenchmarkReport as LargeAlignmentScalingBenchmarkReport,
    LargeAlignmentScalingObservation as LargeAlignmentScalingObservation,
    LargeAlignmentScalingWorkflowBenchmark as LargeAlignmentScalingWorkflowBenchmark,
    LargeDatasetStressObservation as LargeDatasetStressObservation,
    LargeDatasetStressSuiteReport as LargeDatasetStressSuiteReport,
    LargeTreeScalingBenchmarkReport as LargeTreeScalingBenchmarkReport,
    LargeTreeScalingWorkflowBenchmark as LargeTreeScalingWorkflowBenchmark,
    LargeTreeSetScalingBenchmarkReport as LargeTreeSetScalingBenchmarkReport,
    LargeTreeSetScalingObservation as LargeTreeSetScalingObservation,
    LargeTreeSetScalingWorkflowBenchmark as LargeTreeSetScalingWorkflowBenchmark,
    TreeComparisonBenchmarkReport as TreeComparisonBenchmarkReport,
    TreeSetConsensusBenchmarkReport as TreeSetConsensusBenchmarkReport,
    TreeValidationBenchmarkReport as TreeValidationBenchmarkReport,
    WorkflowPracticalLimitEntry as WorkflowPracticalLimitEntry,
    WorkflowPracticalLimitReport as WorkflowPracticalLimitReport,
    _StressObservationPayload as _StressObservationPayload,
    _StressTierConfig as _StressTierConfig,
)
from .model_fitting import (
    LargeTreeModelFittingBenchmarkBundle as LargeTreeModelFittingBenchmarkBundle,
    LargeTreeModelFittingBenchmarkReport as LargeTreeModelFittingBenchmarkReport,
    LargeTreeModelFittingObservation as LargeTreeModelFittingObservation,
    LargeTreeModelFittingThreshold as LargeTreeModelFittingThreshold,
    benchmark_large_tree_model_fitting as benchmark_large_tree_model_fitting,
    write_large_tree_model_fitting_bundle as write_large_tree_model_fitting_bundle,
    write_large_tree_model_fitting_observation_table as write_large_tree_model_fitting_observation_table,
    write_large_tree_model_fitting_summary_table as write_large_tree_model_fitting_summary_table,
)
from .review_benchmarks import (
    benchmark_alignment_diagnostics as benchmark_alignment_diagnostics,
    benchmark_alignment_site_scaling as benchmark_alignment_site_scaling,
    benchmark_tree_comparison as benchmark_tree_comparison,
    benchmark_tree_set_consensus as benchmark_tree_set_consensus,
    benchmark_tree_validation as benchmark_tree_validation,
)
from .scaling_benchmarks import (
    benchmark_large_alignment_scaling as benchmark_large_alignment_scaling,
    benchmark_large_tree_scaling as benchmark_large_tree_scaling,
    benchmark_large_tree_set_scaling as benchmark_large_tree_set_scaling,
)
from .real_dataset_macroevolution import (
    RealDatasetMacroevolutionAlignmentReviewRow as RealDatasetMacroevolutionAlignmentReviewRow,
    RealDatasetMacroevolutionBenchmarkBundle as RealDatasetMacroevolutionBenchmarkBundle,
    RealDatasetMacroevolutionBenchmarkDemoResult as RealDatasetMacroevolutionBenchmarkDemoResult,
    RealDatasetMacroevolutionBenchmarkReport as RealDatasetMacroevolutionBenchmarkReport,
    RealDatasetMacroevolutionModelRow as RealDatasetMacroevolutionModelRow,
    RealDatasetMacroevolutionParityRow as RealDatasetMacroevolutionParityRow,
    RealDatasetMacroevolutionSummaryRow as RealDatasetMacroevolutionSummaryRow,
    benchmark_real_dataset_macroevolution as benchmark_real_dataset_macroevolution,
    run_real_dataset_macroevolution_benchmark_demo as run_real_dataset_macroevolution_benchmark_demo,
    write_geiger_real_dataset_reference_payload_table as write_geiger_real_dataset_reference_payload_table,
    write_real_dataset_macroevolution_alignment_review_table as write_real_dataset_macroevolution_alignment_review_table,
    write_real_dataset_macroevolution_bundle as write_real_dataset_macroevolution_bundle,
    write_real_dataset_macroevolution_model_table as write_real_dataset_macroevolution_model_table,
    write_real_dataset_macroevolution_parity_table as write_real_dataset_macroevolution_parity_table,
    write_real_dataset_macroevolution_summary_table as write_real_dataset_macroevolution_summary_table,
)


def benchmark_large_dataset_stress_suite(
    *,
    tier: str = "small",
) -> LargeDatasetStressSuiteReport:
    """Benchmark large owned workloads across one governed stress tier."""
    config = _resolve_stress_tier_config(tier)
    observations: list[LargeDatasetStressObservation] = []
    limitations = [
        "resource peaks are measured with python tracemalloc where possible and reuse stage-level engine memory observations when an owned workflow already records them",
        "timeout_seconds is a workload budget recorded for review; only engine-backed workflows enforce it internally during execution",
    ]
    with tempfile.TemporaryDirectory(prefix=f"bijux-stress-{config.tier}-") as tmpdir:
        root = Path(tmpdir)
        workloads = [
            lambda: _large_alignment_stress_payload(
                root=root / "alignment", config=config
            ),
            lambda: _supermatrix_stress_payload(
                root=root / "supermatrix", config=config
            ),
            lambda: _tree_set_stress_payload(root=root / "tree-set", config=config),
            lambda: _comparative_stress_payload(
                root=root / "comparative", config=config
            ),
            lambda: _table_generation_stress_payload(
                root=root / "tables", config=config
            ),
        ]
        for workload in workloads:
            payload, runtime_seconds, peak_memory_bytes, memory_observation_kind = (
                _measure_stress_workload(workload)
            )
            observations.append(
                LargeDatasetStressObservation(
                    workload=payload.workload,
                    tier=config.tier,
                    timeout_seconds=config.timeout_seconds,
                    input_size_bytes=payload.input_size_bytes,
                    sequence_count=payload.sequence_count,
                    alignment_length=payload.alignment_length,
                    tree_count=payload.tree_count,
                    taxon_count=payload.taxon_count,
                    locus_count=payload.locus_count,
                    runtime_seconds=round(runtime_seconds, 15),
                    peak_memory_bytes=peak_memory_bytes,
                    memory_observation_kind=memory_observation_kind,
                    output_row_count=payload.output_row_count,
                    notes=payload.notes,
                )
            )
    return LargeDatasetStressSuiteReport(
        tier=config.tier,
        observations=observations,
        limitations=limitations,
    )
def _large_tree_limit_entries(
    report: LargeTreeScalingBenchmarkReport,
) -> list[WorkflowPracticalLimitEntry]:
    entries: list[WorkflowPracticalLimitEntry] = []
    for workflow in report.workflows:
        entries.append(
            WorkflowPracticalLimitEntry(
                workflow=workflow.workflow,
                evidence_source="large-tree-scaling",
                tested_taxon_limit=max(report.tip_counts),
                tested_site_limit=None,
                tested_tree_limit=2 if workflow.workflow == "tree-comparison" else 1,
                tested_posterior_size=None,
                max_runtime_seconds=_max_runtime_seconds(workflow.observations),
                max_peak_memory_bytes=_max_peak_memory_bytes(workflow.observations),
                memory_observation_kind="python-tracemalloc",
                notes=workflow.notes,
            )
        )
    return entries


def _large_alignment_limit_entries(
    report: LargeAlignmentScalingBenchmarkReport,
) -> list[WorkflowPracticalLimitEntry]:
    entries: list[WorkflowPracticalLimitEntry] = []
    for workflow in report.workflows:
        entries.append(
            WorkflowPracticalLimitEntry(
                workflow=workflow.workflow,
                evidence_source="large-alignment-scaling",
                tested_taxon_limit=max(report.sequence_counts),
                tested_site_limit=max(report.alignment_lengths),
                tested_tree_limit=None,
                tested_posterior_size=None,
                max_runtime_seconds=_max_runtime_seconds(workflow.observations),
                max_peak_memory_bytes=_max_peak_memory_bytes(workflow.observations),
                memory_observation_kind="python-tracemalloc",
                notes=workflow.notes,
            )
        )
    return entries


def _large_tree_set_limit_entries(
    report: LargeTreeSetScalingBenchmarkReport,
) -> list[WorkflowPracticalLimitEntry]:
    entries: list[WorkflowPracticalLimitEntry] = []
    for workflow in report.workflows:
        entries.append(
            WorkflowPracticalLimitEntry(
                workflow=workflow.workflow,
                evidence_source="large-tree-set-scaling",
                tested_taxon_limit=max(report.tip_counts),
                tested_site_limit=None,
                tested_tree_limit=max(report.tree_counts),
                tested_posterior_size=max(report.tree_counts),
                max_runtime_seconds=_max_runtime_seconds(workflow.observations),
                max_peak_memory_bytes=_max_peak_memory_bytes(workflow.observations),
                memory_observation_kind="python-tracemalloc",
                notes=workflow.notes,
            )
        )
    return entries


def _stress_limit_entries(
    reports: list[LargeDatasetStressSuiteReport],
) -> list[WorkflowPracticalLimitEntry]:
    rows_by_workload: dict[str, list[LargeDatasetStressObservation]] = {}
    for report in reports:
        for row in report.observations:
            rows_by_workload.setdefault(row.workload, []).append(row)

    entries: list[WorkflowPracticalLimitEntry] = []
    for workload in sorted(rows_by_workload):
        rows = rows_by_workload[workload]
        tiers = sorted({row.tier for row in rows})
        max_runtime = round(max(row.runtime_seconds for row in rows), 15)
        max_peak_memory = max(row.peak_memory_bytes for row in rows)
        memory_kinds = sorted(
            {row.memory_observation_kind for row in rows if row.memory_observation_kind}
        )
        notes: list[str] = []
        for row in rows:
            for note in row.notes:
                if note not in notes:
                    notes.append(note)
        notes.append(
            "tested through governed stress tiers: " + ", ".join(tiers)
        )
        entries.append(
            WorkflowPracticalLimitEntry(
                workflow=workload,
                evidence_source="stress-suite",
                tested_taxon_limit=max(
                    row.taxon_count for row in rows if row.taxon_count is not None
                )
                if any(row.taxon_count is not None for row in rows)
                else None,
                tested_site_limit=max(
                    row.alignment_length
                    for row in rows
                    if row.alignment_length is not None
                )
                if any(row.alignment_length is not None for row in rows)
                else None,
                tested_tree_limit=max(
                    row.tree_count for row in rows if row.tree_count is not None
                )
                if any(row.tree_count is not None for row in rows)
                else None,
                tested_posterior_size=max(
                    row.tree_count
                    for row in rows
                    if row.workload == "posterior-tree-set-consensus"
                    and row.tree_count is not None
                )
                if any(
                    row.workload == "posterior-tree-set-consensus"
                    and row.tree_count is not None
                    for row in rows
                )
                else None,
                max_runtime_seconds=max_runtime,
                max_peak_memory_bytes=max_peak_memory,
                memory_observation_kind=(
                    None if not memory_kinds else ",".join(memory_kinds)
                ),
                notes=notes,
            )
        )
    return entries


def benchmark_workflow_practical_limits(
    *,
    replicates: int = 1,
    tree_tip_counts: list[int] | None = None,
    alignment_size_classes: list[tuple[str, int, int]] | None = None,
    tree_set_size_classes: list[tuple[str, int, int]] | None = None,
    stress_tiers: list[str] | None = None,
) -> WorkflowPracticalLimitReport:
    """Report the largest governed workflow classes currently exercised in benchmark and stress lanes."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    tiers = ["heavy"] if stress_tiers is None else list(stress_tiers)
    if not tiers:
        raise ValueError("stress_tiers must contain at least one governed tier")

    tree_report = benchmark_large_tree_scaling(
        replicates=replicates,
        tip_counts=tree_tip_counts,
    )
    alignment_report = benchmark_large_alignment_scaling(
        replicates=replicates,
        size_classes=alignment_size_classes,
    )
    tree_set_report = benchmark_large_tree_set_scaling(
        replicates=replicates,
        size_classes=tree_set_size_classes,
    )
    stress_reports = [
        benchmark_large_dataset_stress_suite(tier=tier) for tier in tiers
    ]

    entries = [
        *_large_tree_limit_entries(tree_report),
        *_large_alignment_limit_entries(alignment_report),
        *_large_tree_set_limit_entries(tree_set_report),
        *_stress_limit_entries(stress_reports),
    ]
    limitations: list[str] = []
    for report in (tree_report, alignment_report, tree_set_report):
        for item in report.limitations:
            if item not in limitations:
                limitations.append(item)
    for report in stress_reports:
        for item in report.limitations:
            if item not in limitations:
                limitations.append(item)
    limitations.append(
        "practical limits report tested maxima from governed benchmark and stress lanes; it does not claim untested workflows or hardware-specific guarantees"
    )
    return WorkflowPracticalLimitReport(
        replicates=replicates,
        stress_tiers=tiers,
        entries=entries,
        limitations=limitations,
    )
