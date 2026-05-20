from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .formatting import _format_number, _format_optional_bool
from .models import ContinuousModeRecoveryReport


def write_continuous_mode_recovery_summary_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write one integrated recovery summary row per governed simulation case."""
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "label",
            "tree_path",
            "generating_model",
            "expected_selected_model",
            "bijux_selected_model",
            "geiger_selected_model",
            "bijux_selection_matches_expectation",
            "geiger_selection_matches_expectation",
            "parameter_row_count",
            "parameter_comparison_row_count",
            "bijux_parameter_pass_count",
            "geiger_parameter_pass_count",
            "parameter_closer_to_truth_count_bijux",
            "parameter_closer_to_truth_count_geiger",
            "expected_warning_count",
            "expected_warning_kinds_present",
            "warning_count",
            "notes",
        ],
        rows=[
            {
                "case_id": case.scenario.case_id,
                "label": case.scenario.label,
                "tree_path": str(case.tree_path),
                "generating_model": case.scenario.generating_model,
                "expected_selected_model": case.scenario.expected_selected_model or "",
                "bijux_selected_model": case.selected_model or "",
                "geiger_selected_model": case.geiger_selected_model or "",
                "bijux_selection_matches_expectation": _format_optional_bool(
                    case.selection_matches_expectation
                ),
                "geiger_selection_matches_expectation": _format_optional_bool(
                    case.geiger_selection_matches_expectation
                ),
                "parameter_row_count": str(len(case.parameter_rows)),
                "parameter_comparison_row_count": str(
                    len(case.parameter_comparison_rows)
                ),
                "bijux_parameter_pass_count": str(
                    sum(
                        1
                        for row in case.parameter_rows
                        if row.recovery_engine == "bijux" and row.within_tolerance
                    )
                ),
                "geiger_parameter_pass_count": str(
                    sum(
                        1
                        for row in case.parameter_rows
                        if row.recovery_engine == "geiger" and row.within_tolerance
                    )
                ),
                "parameter_closer_to_truth_count_bijux": str(
                    sum(
                        1
                        for row in case.parameter_comparison_rows
                        if row.closer_engine == "bijux"
                    )
                ),
                "parameter_closer_to_truth_count_geiger": str(
                    sum(
                        1
                        for row in case.parameter_comparison_rows
                        if row.closer_engine == "geiger"
                    )
                ),
                "expected_warning_count": str(
                    len(case.scenario.expected_warning_kinds)
                ),
                "expected_warning_kinds_present": str(
                    case.expected_warning_kinds_present
                ).lower(),
                "warning_count": str(len(case.warning_rows)),
                "notes": case.scenario.notes,
            }
            for case in report.case_reports
        ],
    )


def write_continuous_mode_recovery_parameter_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write all truth-versus-fit parameter-recovery rows across the cases."""
    rows = [
        {
            "case_id": row.case_id,
            "generating_model": row.generating_model,
            "recovery_engine": row.recovery_engine,
            "fitted_model": row.fitted_model,
            "fit_status": row.fit_status,
            "parameter": row.parameter,
            "true_value": _format_number(row.true_value),
            "estimated_value": _format_number(row.estimated_value),
            "absolute_error": _format_number(row.absolute_error),
            "relative_error": _format_number(row.relative_error),
            "tolerance": _format_number(row.tolerance),
            "within_tolerance": str(row.within_tolerance).lower(),
            "interpretation": row.interpretation,
        }
        for case in report.case_reports
        for row in case.parameter_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "generating_model",
            "recovery_engine",
            "fitted_model",
            "fit_status",
            "parameter",
            "true_value",
            "estimated_value",
            "absolute_error",
            "relative_error",
            "tolerance",
            "within_tolerance",
            "interpretation",
        ],
        rows=rows,
    )


def write_continuous_mode_recovery_parameter_comparison_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write paired Bijux-versus-geiger parameter-recovery comparisons."""
    rows = [
        {
            "case_id": row.case_id,
            "generating_model": row.generating_model,
            "parameter": row.parameter,
            "true_value": _format_number(row.true_value),
            "bijux_estimated_value": _format_number(row.bijux_estimated_value),
            "geiger_estimated_value": _format_number(row.geiger_estimated_value),
            "bijux_absolute_error": _format_number(row.bijux_absolute_error),
            "geiger_absolute_error": _format_number(row.geiger_absolute_error),
            "bijux_within_tolerance": str(row.bijux_within_tolerance).lower(),
            "geiger_within_tolerance": str(row.geiger_within_tolerance).lower(),
            "closer_engine": row.closer_engine,
            "tolerance": _format_number(row.tolerance),
            "interpretation": row.interpretation,
        }
        for case in report.case_reports
        for row in case.parameter_comparison_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "generating_model",
            "parameter",
            "true_value",
            "bijux_estimated_value",
            "geiger_estimated_value",
            "bijux_absolute_error",
            "geiger_absolute_error",
            "bijux_within_tolerance",
            "geiger_within_tolerance",
            "closer_engine",
            "tolerance",
            "interpretation",
        ],
        rows=rows,
    )


def write_continuous_mode_recovery_model_choice_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write candidate-model review rows for both recovery engines."""
    rows = [
        {
            "case_id": row.case_id,
            "generating_model": row.generating_model,
            "recovery_engine": row.recovery_engine,
            "expected_selected_model": row.expected_selected_model or "",
            "model": row.model,
            "parameter_count": str(row.parameter_count),
            "log_likelihood": _format_number(row.log_likelihood),
            "aic": _format_number(row.aic),
            "aicc": _format_number(row.aicc),
            "selected": str(row.selected).lower(),
        }
        for case in report.case_reports
        for row in case.model_choice_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "generating_model",
            "recovery_engine",
            "expected_selected_model",
            "model",
            "parameter_count",
            "log_likelihood",
            "aic",
            "aicc",
            "selected",
        ],
        rows=rows,
    )


def write_continuous_mode_recovery_execution_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write fit and model-comparison execution status rows for each engine."""
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "recovery_engine",
            "operation",
            "fitted_model",
            "fit_status",
            "selected_model",
            "optimizer_name",
            "converged",
            "hit_lower_parameter_boundary",
            "hit_upper_parameter_boundary",
            "warning_count",
            "failure_reason",
        ],
        rows=[
            {
                "case_id": row.case_id,
                "recovery_engine": row.recovery_engine,
                "operation": row.operation,
                "fitted_model": row.fitted_model,
                "fit_status": row.fit_status,
                "selected_model": row.selected_model or "",
                "optimizer_name": row.optimizer_name or "",
                "converged": _format_optional_bool(row.converged),
                "hit_lower_parameter_boundary": _format_optional_bool(
                    row.hit_lower_parameter_boundary
                ),
                "hit_upper_parameter_boundary": _format_optional_bool(
                    row.hit_upper_parameter_boundary
                ),
                "warning_count": str(row.warning_count),
                "failure_reason": row.failure_reason or "",
            }
            for case in report.case_reports
            for row in case.execution_rows
        ],
    )


def write_continuous_mode_recovery_warning_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write the governed warning rows observed across recovery cases."""
    return write_taxon_rows(
        path,
        columns=["case_id", "recovery_engine", "fitted_model", "kind", "message"],
        rows=[
            {
                "case_id": row.case_id,
                "recovery_engine": row.recovery_engine,
                "fitted_model": row.fitted_model,
                "kind": row.kind,
                "message": row.message,
            }
            for case in report.case_reports
            for row in case.warning_rows
        ],
    )
