from __future__ import annotations

from .contracts import (
    RabiesMethodSensitivitySlurmAssumptionRow,
    RabiesMethodSensitivitySlurmPlanningReport,
)
from .interfaces import WorkflowReportLike
from .sizing import build_job_plan_row


def build_rabies_method_sensitivity_slurm_planning_report(
    report: WorkflowReportLike,
) -> RabiesMethodSensitivitySlurmPlanningReport:
    """Estimate one schedulable Slurm job per declared workflow variant."""
    task_records = {record.variant_id: record for record in report.task_records}
    rows = tuple(
        build_job_plan_row(
            dataset_id=report.dataset.dataset_id,
            workflow_prefix=report.dataset.workflow_prefix,
            taxon_count=report.dataset.taxon_count,
            variant_run=variant_run,
            task_record=task_records[variant_run.config.variant_id],
            iqtree_threads=report.iqtree_threads,
            bootstrap_replicates=report.bootstrap_replicates,
        )
        for variant_run in report.variant_runs
    )
    assumptions = (
        RabiesMethodSensitivitySlurmAssumptionRow(
            assumption_id="observed-output-footprint",
            parameter="observed_output_bytes",
            value="sum of files currently written under each variant output root",
            rationale=(
                "Output-size and scratch estimates are anchored to the real bytes "
                "written by the current governed workflow run instead of a synthetic "
                "placeholder volume."
            ),
        ),
        RabiesMethodSensitivitySlurmAssumptionRow(
            assumption_id="alignment-mode-cpu-bump",
            parameter="estimated_cpus_per_task",
            value="ginsi variants reserve at least 2 CPUs; auto variants follow iqtree_threads",
            rationale=(
                "The ginsi alignment mode is materially heavier than the auto mode, "
                "so the planner allocates an extra CPU even when the governed IQ-TREE "
                "thread count stays at 1."
            ),
        ),
        RabiesMethodSensitivitySlurmAssumptionRow(
            assumption_id="bootstrap-driven-wallclock",
            parameter="estimated_wallclock_minutes",
            value="20-minute floor plus linear scaling with trimmed sites, taxa, and bootstrap count",
            rationale=(
                "The batch wallclock estimate should stay conservative when the "
                "bootstrap burden rises, even on this compact rabies panel."
            ),
        ),
        RabiesMethodSensitivitySlurmAssumptionRow(
            assumption_id="minimum-memory-floor",
            parameter="estimated_memory_mib",
            value="1024 MiB minimum, rounded to 256 MiB blocks",
            rationale=(
                "The dataset is small, but Slurm requests still need enough room for "
                "alignment, inference, logging, and report materialization overhead."
            ),
        ),
        RabiesMethodSensitivitySlurmAssumptionRow(
            assumption_id="scratch-buffer",
            parameter="estimated_scratch_mib",
            value="observed output footprint inflated 64x with a 256 MiB minimum",
            rationale=(
                "Scratch space should cover intermediate alignments, trees, manifests, "
                "and logs without pretending that the final output footprint is the "
                "whole temporary working set."
            ),
        ),
    )
    return RabiesMethodSensitivitySlurmPlanningReport(
        dataset_id=report.dataset.dataset_id,
        workflow_prefix=report.dataset.workflow_prefix,
        job_count=len(rows),
        bootstrap_replicates=report.bootstrap_replicates,
        iqtree_threads=report.iqtree_threads,
        total_estimated_core_hours=round(
            sum(row.estimated_core_hours for row in rows), 2
        ),
        maximum_estimated_memory_mib=max(row.estimated_memory_mib for row in rows),
        maximum_estimated_wallclock_minutes=max(
            row.estimated_wallclock_minutes for row in rows
        ),
        total_estimated_scratch_mib=sum(row.estimated_scratch_mib for row in rows),
        total_estimated_output_mib=sum(row.estimated_output_mib for row in rows),
        assumptions=assumptions,
        rows=rows,
    )
