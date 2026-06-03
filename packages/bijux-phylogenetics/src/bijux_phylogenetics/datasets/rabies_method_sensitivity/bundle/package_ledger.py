from __future__ import annotations

from pathlib import Path

from ..models import (
    RabiesMethodSensitivityCladeRow,
    RabiesMethodSensitivityConclusionRow,
    RabiesMethodSensitivityPanelWorkflowReport,
    RabiesMethodSensitivityPreprocessingComparisonRow,
)
from .shared import _format_float, _format_optional_float, _write_tsv


def _write_workflow_summary_table(
    path: Path, report: RabiesMethodSensitivityPanelWorkflowReport
) -> Path:
    serious_conflicts = [
        variant.inference_comparison.conclusion_summary.serious_conflict_count
        for variant in report.variant_runs
    ]
    rows = [
        [
            "dataset_id",
            "variant_count",
            "stable_clade_count",
            "changed_clade_count",
            "preprocessing_change_pair_count",
            "rooted_engine_change_variant_count",
            "serious_conflict_variant_count",
            "maximum_serious_conflict_count",
        ],
        [
            report.dataset.dataset_id,
            str(len(report.variant_runs)),
            str(len(report.stable_clade_rows)),
            str(len(report.changed_clade_rows)),
            str(
                sum(
                    1
                    for row in report.preprocessing_comparison_rows
                    if row.robinson_foulds_distance > 0
                    or row.same_taxa_different_rooting
                )
            ),
            str(
                sum(
                    1
                    for variant in report.variant_runs
                    if variant.rooted_engine_comparison.robinson_foulds_distance > 0
                    or variant.rooted_engine_comparison.same_taxa_different_rooting
                )
            ),
            str(sum(1 for value in serious_conflicts if value > 0)),
            str(max(serious_conflicts)),
        ],
    ]
    return _write_tsv(path, rows)


def _write_variant_summary_table(
    path: Path, report: RabiesMethodSensitivityPanelWorkflowReport
) -> Path:
    rows = [
        [
            "variant_id",
            "alignment_mode",
            "trimming_mode",
            "trim_gap_threshold",
            "alignment_length",
            "trimmed_alignment_length",
            "selected_model",
            "minimum_support",
            "maximum_support",
            "stable_clade_count",
            "engine_specific_clade_count",
            "serious_conflict_count",
            "rooted_engine_rf_distance",
            "rooted_engine_same_taxa_different_rooting",
        ]
    ]
    for variant in report.variant_runs:
        summary = variant.inference_comparison.conclusion_summary
        rows.append(
            [
                variant.config.variant_id,
                variant.config.alignment_mode,
                variant.config.trimming_mode,
                _format_float(variant.config.trim_gap_threshold),
                str(variant.alignment_length),
                str(variant.trimmed_alignment_length),
                variant.inference_comparison.selected_model,
                _format_optional_float(
                    variant.inference_comparison.iqtree_support_workflow.bootstrap_support_summary.minimum_support
                    if variant.inference_comparison.iqtree_support_workflow.bootstrap_support_summary
                    is not None
                    else None
                ),
                _format_optional_float(
                    variant.inference_comparison.iqtree_support_workflow.bootstrap_support_summary.maximum_support
                    if variant.inference_comparison.iqtree_support_workflow.bootstrap_support_summary
                    is not None
                    else None
                ),
                str(summary.stable_clade_count),
                str(summary.engine_specific_clade_count),
                str(summary.serious_conflict_count),
                str(variant.rooted_engine_comparison.robinson_foulds_distance),
                str(
                    variant.rooted_engine_comparison.same_taxa_different_rooting
                ).lower(),
            ]
        )
    return _write_tsv(path, rows)


def _write_preprocessing_comparison_table(
    path: Path,
    rows: tuple[RabiesMethodSensitivityPreprocessingComparisonRow, ...],
) -> Path:
    rendered = [
        [
            "left_variant_id",
            "right_variant_id",
            "comparison_axis",
            "robinson_foulds_distance",
            "normalized_robinson_foulds",
            "same_taxa_different_rooting",
        ]
    ]
    for row in rows:
        rendered.append(
            [
                row.left_variant_id,
                row.right_variant_id,
                row.comparison_axis,
                str(row.robinson_foulds_distance),
                _format_float(row.normalized_robinson_foulds),
                str(row.same_taxa_different_rooting).lower(),
            ]
        )
    return _write_tsv(path, rendered)


def _write_clade_table(
    path: Path, rows: tuple[RabiesMethodSensitivityCladeRow, ...]
) -> Path:
    rendered = [
        [
            "split_id",
            "conclusion_class",
            "evidence_class",
            "occurrence_count",
            "variant_count",
            "detail",
        ]
    ]
    for row in rows:
        rendered.append(
            [
                row.split_id,
                row.conclusion_class,
                row.evidence_class,
                str(row.occurrence_count),
                str(row.variant_count),
                row.detail,
            ]
        )
    return _write_tsv(path, rendered)


def _write_conclusion_summary_table(
    path: Path, rows: tuple[RabiesMethodSensitivityConclusionRow, ...]
) -> Path:
    rendered = [
        [
            "conclusion_id",
            "method_axis",
            "stability_status",
            "claim",
            "evidence",
            "caution",
        ]
    ]
    for row in rows:
        rendered.append(
            [
                row.conclusion_id,
                row.method_axis,
                row.stability_status,
                row.claim,
                row.evidence,
                row.caution,
            ]
        )
    return _write_tsv(path, rendered)


def _write_parallel_execution_summary_table(
    path: Path, report: RabiesMethodSensitivityPanelWorkflowReport
) -> Path:
    rows = [
        [
            "variant_id",
            "label",
            "execution_mode",
            "status",
            "log_path",
            "error_code",
        ]
    ]
    for task in report.task_records:
        rows.append(
            [
                task.variant_id,
                task.label,
                task.execution_mode,
                task.status,
                Path("parallel-logs", task.log_path.name).as_posix(),
                "" if task.error_code is None else task.error_code,
            ]
        )
    return _write_tsv(path, rows)
