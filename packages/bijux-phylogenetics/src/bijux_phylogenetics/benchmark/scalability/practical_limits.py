from __future__ import annotations

from .._measurement import max_peak_memory_bytes, max_runtime_seconds
from ..contracts import (
    LargeAlignmentScalingBenchmarkReport,
    LargeDatasetStressObservation,
    LargeDatasetStressSuiteReport,
    LargeTreeScalingBenchmarkReport,
    LargeTreeSetScalingBenchmarkReport,
    WorkflowPracticalLimitEntry,
    WorkflowPracticalLimitReport,
)
from .scaling import (
    benchmark_large_alignment_scaling,
    benchmark_large_tree_scaling,
    benchmark_large_tree_set_scaling,
)
from .stress import benchmark_large_dataset_stress_suite


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
                max_runtime_seconds=max_runtime_seconds(workflow.observations),
                max_peak_memory_bytes=max_peak_memory_bytes(workflow.observations),
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
                max_runtime_seconds=max_runtime_seconds(workflow.observations),
                max_peak_memory_bytes=max_peak_memory_bytes(workflow.observations),
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
                max_runtime_seconds=max_runtime_seconds(workflow.observations),
                max_peak_memory_bytes=max_peak_memory_bytes(workflow.observations),
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
        notes.append("tested through governed stress tiers: " + ", ".join(tiers))
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


def _collect_limitations(*reports: object) -> list[str]:
    limitations: list[str] = []
    for report in reports:
        for item in report.limitations:
            if item not in limitations:
                limitations.append(item)
    return limitations


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
    stress_reports = [benchmark_large_dataset_stress_suite(tier=tier) for tier in tiers]

    entries = [
        *_large_tree_limit_entries(tree_report),
        *_large_alignment_limit_entries(alignment_report),
        *_large_tree_set_limit_entries(tree_set_report),
        *_stress_limit_entries(stress_reports),
    ]
    limitations = _collect_limitations(
        tree_report,
        alignment_report,
        tree_set_report,
        *stress_reports,
    )
    limitations.append(
        "practical limits report tested maxima from governed benchmark and stress lanes; it does not claim untested workflows or hardware-specific guarantees"
    )
    return WorkflowPracticalLimitReport(
        replicates=replicates,
        stress_tiers=tiers,
        entries=entries,
        limitations=limitations,
    )
