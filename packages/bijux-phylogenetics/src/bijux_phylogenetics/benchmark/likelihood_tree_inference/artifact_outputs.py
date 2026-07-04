from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .contracts import (
    MaximumLikelihoodWrapperBenchmarkReport,
    NativeMaximumLikelihoodAccuracyBenchmarkReport,
    NativeMaximumLikelihoodMemoryBenchmarkReport,
    NativeMaximumLikelihoodSpeedBenchmarkReport,
    NativeMaximumLikelihoodSuiteReport,
)


def write_native_maximum_likelihood_accuracy_summary_table(
    path: Path,
    report: NativeMaximumLikelihoodAccuracyBenchmarkReport,
) -> Path:
    """Write the stable summary counts for one native ML accuracy benchmark run."""
    return write_taxon_rows(
        path,
        columns=[
            "observation_count",
            "wrapper_observation_count",
            "model_names",
            "search_methods",
            "topology_recovery_case_count",
            "parameter_recovery_case_count",
            "parameter_not_applicable_case_count",
            "likelihood_quality_case_count",
            "wrapper_supported_case_count",
            "wrapper_clear_case_count",
            "all_truth_cases_clear",
            "all_supported_wrapper_cases_clear",
            "limitations",
        ],
        rows=[
            {
                "observation_count": len(report.observations),
                "wrapper_observation_count": len(report.wrapper_observations),
                "model_names": " | ".join(report.model_names),
                "search_methods": " | ".join(report.search_methods),
                "topology_recovery_case_count": report.topology_recovery_case_count,
                "parameter_recovery_case_count": report.parameter_recovery_case_count,
                "parameter_not_applicable_case_count": (
                    report.parameter_not_applicable_case_count
                ),
                "likelihood_quality_case_count": report.likelihood_quality_case_count,
                "wrapper_supported_case_count": report.wrapper_supported_case_count,
                "wrapper_clear_case_count": report.wrapper_clear_case_count,
                "all_truth_cases_clear": report.all_truth_cases_clear,
                "all_supported_wrapper_cases_clear": (
                    report.all_supported_wrapper_cases_clear
                ),
                "limitations": " | ".join(report.limitations),
            }
        ],
    )


def write_native_maximum_likelihood_accuracy_observation_table(
    path: Path,
    report: NativeMaximumLikelihoodAccuracyBenchmarkReport,
) -> Path:
    """Write one stable observation row per governed native ML accuracy lane."""
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "model_name",
            "search_method",
            "taxon_count",
            "site_count",
            "simulation_seed",
            "selected_model_name",
            "topology_recovered",
            "topology_recovery_status",
            "robinson_foulds_distance",
            "normalized_robinson_foulds",
            "parameter_recovery_status",
            "parameter_truth",
            "estimated_parameter_values",
            "parameter_absolute_deltas",
            "true_topology_log_likelihood",
            "inferred_log_likelihood",
            "likelihood_shortfall",
            "likelihood_quality_status",
        ],
        rows=[
            {
                "case_id": row.case_id,
                "model_name": row.model_name,
                "search_method": row.search_method,
                "taxon_count": row.taxon_count,
                "site_count": row.site_count,
                "simulation_seed": row.simulation_seed,
                "selected_model_name": row.selected_model_name,
                "topology_recovered": row.topology_recovered,
                "topology_recovery_status": row.topology_recovery_status,
                "robinson_foulds_distance": row.robinson_foulds_distance,
                "normalized_robinson_foulds": format(
                    row.normalized_robinson_foulds,
                    ".15g",
                ),
                "parameter_recovery_status": row.parameter_recovery_status,
                "parameter_truth": json.dumps(row.parameter_truth, sort_keys=True),
                "estimated_parameter_values": json.dumps(
                    row.estimated_parameter_values,
                    sort_keys=True,
                ),
                "parameter_absolute_deltas": json.dumps(
                    row.parameter_absolute_deltas,
                    sort_keys=True,
                ),
                "true_topology_log_likelihood": format(
                    row.true_topology_log_likelihood,
                    ".15g",
                ),
                "inferred_log_likelihood": format(
                    row.inferred_log_likelihood,
                    ".15g",
                ),
                "likelihood_shortfall": format(row.likelihood_shortfall, ".15g"),
                "likelihood_quality_status": row.likelihood_quality_status,
            }
            for row in report.observations
        ],
    )


def write_native_maximum_likelihood_accuracy_wrapper_table(
    path: Path,
    report: NativeMaximumLikelihoodAccuracyBenchmarkReport,
) -> Path:
    """Write one stable wrapper-evidence row per native ML accuracy benchmark run."""
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "wrapper_engine",
            "status",
            "supported",
            "blocking",
            "topology_equal",
            "native_log_likelihood",
            "wrapper_log_likelihood",
            "log_likelihood_delta",
            "support_split_count",
            "notes",
        ],
        rows=[
            {
                "case_id": row.case_id,
                "wrapper_engine": row.wrapper_engine,
                "status": row.status,
                "supported": row.supported,
                "blocking": row.blocking,
                "topology_equal": row.topology_equal,
                "native_log_likelihood": (
                    ""
                    if row.native_log_likelihood is None
                    else format(row.native_log_likelihood, ".15g")
                ),
                "wrapper_log_likelihood": (
                    ""
                    if row.wrapper_log_likelihood is None
                    else format(row.wrapper_log_likelihood, ".15g")
                ),
                "log_likelihood_delta": (
                    ""
                    if row.log_likelihood_delta is None
                    else format(row.log_likelihood_delta, ".15g")
                ),
                "support_split_count": len(row.support_comparisons),
                "notes": " | ".join(row.notes),
            }
            for row in report.wrapper_observations
        ],
    )


def write_native_maximum_likelihood_accuracy_run_json(
    path: Path,
    report: NativeMaximumLikelihoodAccuracyBenchmarkReport,
) -> Path:
    """Write the machine-readable native ML accuracy benchmark payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_native_maximum_likelihood_memory_summary_table(
    path: Path,
    report: NativeMaximumLikelihoodMemoryBenchmarkReport,
) -> Path:
    """Write the stable summary counts for one native ML memory benchmark run."""
    return write_taxon_rows(
        path,
        columns=[
            "replicates",
            "observation_count",
            "trend_row_count",
            "size_class_labels",
            "model_names",
            "search_methods",
            "stable_memory_case_count",
            "elevated_memory_case_count",
            "memory_blowup_case_count",
            "limitations",
        ],
        rows=[
            {
                "replicates": report.replicates,
                "observation_count": len(report.observations),
                "trend_row_count": len(report.trend_rows),
                "size_class_labels": " | ".join(report.size_class_labels),
                "model_names": " | ".join(report.model_names),
                "search_methods": " | ".join(report.search_methods),
                "stable_memory_case_count": report.stable_memory_case_count,
                "elevated_memory_case_count": report.elevated_memory_case_count,
                "memory_blowup_case_count": report.memory_blowup_case_count,
                "limitations": " | ".join(report.limitations),
            }
        ],
    )


def write_native_maximum_likelihood_memory_observation_table(
    path: Path,
    report: NativeMaximumLikelihoodMemoryBenchmarkReport,
) -> Path:
    """Write one stable observation row per native ML memory benchmark lane."""
    return write_taxon_rows(
        path,
        columns=[
            "label",
            "size_class",
            "model_name",
            "search_method",
            "taxon_count",
            "site_count",
            "runtime_seconds",
            "peak_memory_bytes",
            "memory_observation_kind",
            "selected_model_name",
            "best_run_source_label",
            "search_algorithm",
            "final_log_likelihood",
        ],
        rows=[
            {
                "label": row.label,
                "size_class": row.size_class,
                "model_name": row.model_name,
                "search_method": row.search_method,
                "taxon_count": row.taxon_count,
                "site_count": row.site_count,
                "runtime_seconds": format(row.runtime_seconds, ".15g"),
                "peak_memory_bytes": row.peak_memory_bytes,
                "memory_observation_kind": row.memory_observation_kind,
                "selected_model_name": row.selected_model_name,
                "best_run_source_label": row.best_run_source_label,
                "search_algorithm": row.search_algorithm,
                "final_log_likelihood": format(row.final_log_likelihood, ".15g"),
            }
            for row in report.observations
        ],
    )


def write_native_maximum_likelihood_memory_trend_table(
    path: Path,
    report: NativeMaximumLikelihoodMemoryBenchmarkReport,
) -> Path:
    """Write one stable peak-memory growth row per native ML model and search lane."""
    return write_taxon_rows(
        path,
        columns=[
            "model_name",
            "search_method",
            "baseline_size_class",
            "comparison_size_class",
            "baseline_taxon_count",
            "comparison_taxon_count",
            "baseline_site_count",
            "comparison_site_count",
            "baseline_peak_memory_bytes",
            "comparison_peak_memory_bytes",
            "peak_memory_delta_bytes",
            "peak_memory_ratio",
            "memory_growth_status",
        ],
        rows=[
            {
                "model_name": row.model_name,
                "search_method": row.search_method,
                "baseline_size_class": row.baseline_size_class,
                "comparison_size_class": row.comparison_size_class,
                "baseline_taxon_count": row.baseline_taxon_count,
                "comparison_taxon_count": row.comparison_taxon_count,
                "baseline_site_count": row.baseline_site_count,
                "comparison_site_count": row.comparison_site_count,
                "baseline_peak_memory_bytes": row.baseline_peak_memory_bytes,
                "comparison_peak_memory_bytes": row.comparison_peak_memory_bytes,
                "peak_memory_delta_bytes": row.peak_memory_delta_bytes,
                "peak_memory_ratio": format(row.peak_memory_ratio, ".15g"),
                "memory_growth_status": row.memory_growth_status,
            }
            for row in report.trend_rows
        ],
    )


def write_native_maximum_likelihood_memory_run_json(
    path: Path,
    report: NativeMaximumLikelihoodMemoryBenchmarkReport,
) -> Path:
    """Write the machine-readable native ML memory benchmark payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_native_maximum_likelihood_speed_summary_table(
    path: Path,
    report: NativeMaximumLikelihoodSpeedBenchmarkReport,
) -> Path:
    """Write the stable summary counts for one native ML speed benchmark run."""
    return write_taxon_rows(
        path,
        columns=[
            "replicates",
            "observation_count",
            "trend_row_count",
            "size_class_labels",
            "model_names",
            "search_methods",
            "limitations",
        ],
        rows=[
            {
                "replicates": report.replicates,
                "observation_count": len(report.observations),
                "trend_row_count": len(report.trend_rows),
                "size_class_labels": " | ".join(report.size_class_labels),
                "model_names": " | ".join(report.model_names),
                "search_methods": " | ".join(report.search_methods),
                "limitations": " | ".join(report.limitations),
            }
        ],
    )


def write_native_maximum_likelihood_speed_observation_table(
    path: Path,
    report: NativeMaximumLikelihoodSpeedBenchmarkReport,
) -> Path:
    """Write one stable observation row per native ML speed benchmark lane."""
    return write_taxon_rows(
        path,
        columns=[
            "label",
            "size_class",
            "model_name",
            "search_method",
            "taxon_count",
            "site_count",
            "runtime_seconds",
            "peak_memory_bytes",
            "selected_model_name",
            "best_run_source_label",
            "search_algorithm",
            "final_log_likelihood",
        ],
        rows=[
            {
                "label": row.label,
                "size_class": row.size_class,
                "model_name": row.model_name,
                "search_method": row.search_method,
                "taxon_count": row.taxon_count,
                "site_count": row.site_count,
                "runtime_seconds": format(row.runtime_seconds, ".15g"),
                "peak_memory_bytes": row.peak_memory_bytes,
                "selected_model_name": row.selected_model_name,
                "best_run_source_label": row.best_run_source_label,
                "search_algorithm": row.search_algorithm,
                "final_log_likelihood": format(row.final_log_likelihood, ".15g"),
            }
            for row in report.observations
        ],
    )


def write_native_maximum_likelihood_speed_trend_table(
    path: Path,
    report: NativeMaximumLikelihoodSpeedBenchmarkReport,
) -> Path:
    """Write one stable runtime trend row per native ML model and search lane."""
    return write_taxon_rows(
        path,
        columns=[
            "model_name",
            "search_method",
            "baseline_size_class",
            "comparison_size_class",
            "baseline_taxon_count",
            "comparison_taxon_count",
            "baseline_site_count",
            "comparison_site_count",
            "baseline_runtime_seconds",
            "comparison_runtime_seconds",
            "runtime_delta_seconds",
            "runtime_ratio",
        ],
        rows=[
            {
                "model_name": row.model_name,
                "search_method": row.search_method,
                "baseline_size_class": row.baseline_size_class,
                "comparison_size_class": row.comparison_size_class,
                "baseline_taxon_count": row.baseline_taxon_count,
                "comparison_taxon_count": row.comparison_taxon_count,
                "baseline_site_count": row.baseline_site_count,
                "comparison_site_count": row.comparison_site_count,
                "baseline_runtime_seconds": format(
                    row.baseline_runtime_seconds,
                    ".15g",
                ),
                "comparison_runtime_seconds": format(
                    row.comparison_runtime_seconds,
                    ".15g",
                ),
                "runtime_delta_seconds": format(row.runtime_delta_seconds, ".15g"),
                "runtime_ratio": format(row.runtime_ratio, ".15g"),
            }
            for row in report.trend_rows
        ],
    )


def write_native_maximum_likelihood_speed_run_json(
    path: Path,
    report: NativeMaximumLikelihoodSpeedBenchmarkReport,
) -> Path:
    """Write the machine-readable native ML speed benchmark payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_native_maximum_likelihood_suite_summary_table(
    path: Path,
    report: NativeMaximumLikelihoodSuiteReport,
) -> Path:
    """Write the stable summary counts for one unified native ML benchmark suite."""
    return write_taxon_rows(
        path,
        columns=[
            "observation_count",
            "summary_row_count",
            "native_advantage_case_count",
            "native_bug_case_count",
            "expected_assumption_difference_case_count",
            "unsupported_case_count",
            "all_supported_cases_clear",
            "limitations",
        ],
        rows=[
            {
                "observation_count": len(report.observations),
                "summary_row_count": len(report.summary_rows),
                "native_advantage_case_count": report.native_advantage_case_count,
                "native_bug_case_count": report.native_bug_case_count,
                "expected_assumption_difference_case_count": (
                    report.expected_assumption_difference_case_count
                ),
                "unsupported_case_count": report.unsupported_case_count,
                "all_supported_cases_clear": report.all_supported_cases_clear,
                "limitations": " | ".join(report.limitations),
            }
        ],
    )


def write_native_maximum_likelihood_suite_observation_table(
    path: Path,
    report: NativeMaximumLikelihoodSuiteReport,
) -> Path:
    """Write one stable observation row per unified native ML benchmark case."""
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "dimension",
            "status",
            "native_surface",
            "comparison_surface",
            "measurement_unit",
            "native_value",
            "comparison_value",
            "delta_value",
            "notes",
        ],
        rows=[
            {
                "case_id": row.case_id,
                "dimension": row.dimension,
                "status": row.status,
                "native_surface": row.native_surface,
                "comparison_surface": row.comparison_surface,
                "measurement_unit": row.measurement_unit or "",
                "native_value": _format_optional_float(row.native_value),
                "comparison_value": _format_optional_float(row.comparison_value),
                "delta_value": _format_optional_float(row.delta_value),
                "notes": " | ".join(row.notes),
            }
            for row in report.observations
        ],
    )


def write_native_maximum_likelihood_suite_run_json(
    path: Path,
    report: NativeMaximumLikelihoodSuiteReport,
) -> Path:
    """Write the machine-readable native ML benchmark suite payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_maximum_likelihood_wrapper_comparison_summary_table(
    path: Path,
    report: MaximumLikelihoodWrapperBenchmarkReport,
) -> Path:
    """Write the stable summary counts for one maximum-likelihood wrapper benchmark run."""
    return write_taxon_rows(
        path,
        columns=[
            "case_count",
            "supported_case_count",
            "tolerance_match_case_count",
            "expected_model_assumption_difference_case_count",
            "unsupported_case_count",
            "native_bug_case_count",
            "blocking_case_count",
            "limitations",
        ],
        rows=[
            {
                "case_count": report.case_count,
                "supported_case_count": report.supported_case_count,
                "tolerance_match_case_count": report.tolerance_match_case_count,
                "expected_model_assumption_difference_case_count": (
                    report.expected_model_assumption_difference_case_count
                ),
                "unsupported_case_count": report.unsupported_case_count,
                "native_bug_case_count": report.native_bug_case_count,
                "blocking_case_count": report.blocking_case_count,
                "limitations": " | ".join(report.limitations),
            }
        ],
    )


def write_maximum_likelihood_wrapper_comparison_observation_table(
    path: Path,
    report: MaximumLikelihoodWrapperBenchmarkReport,
) -> Path:
    """Write one stable observation row per maximum-likelihood wrapper benchmark case."""
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "wrapper_engine",
            "native_surface",
            "wrapper_surface",
            "comparison_policy",
            "status",
            "supported",
            "blocking",
            "native_runtime_seconds",
            "wrapper_runtime_seconds",
            "runtime_delta_seconds",
            "native_peak_memory_bytes",
            "wrapper_peak_memory_bytes",
            "peak_memory_delta_bytes",
            "native_log_likelihood",
            "wrapper_log_likelihood",
            "log_likelihood_delta",
            "topology_equal",
            "robinson_foulds_distance",
            "support_split_count",
            "max_support_percent_delta",
            "notes",
        ],
        rows=[
            {
                "case_id": row.case_id,
                "wrapper_engine": row.wrapper_engine,
                "native_surface": row.native_surface,
                "wrapper_surface": row.wrapper_surface,
                "comparison_policy": row.comparison_policy,
                "status": row.status,
                "supported": row.supported,
                "blocking": row.blocking,
                "native_runtime_seconds": _format_optional_float(
                    row.native_runtime_seconds
                ),
                "wrapper_runtime_seconds": _format_optional_float(
                    row.wrapper_runtime_seconds
                ),
                "runtime_delta_seconds": _format_optional_float(
                    row.runtime_delta_seconds
                ),
                "native_peak_memory_bytes": _format_optional_int(
                    row.native_peak_memory_bytes
                ),
                "wrapper_peak_memory_bytes": _format_optional_int(
                    row.wrapper_peak_memory_bytes
                ),
                "peak_memory_delta_bytes": _format_optional_int(
                    row.peak_memory_delta_bytes
                ),
                "native_log_likelihood": _format_optional_float(
                    row.native_log_likelihood
                ),
                "wrapper_log_likelihood": _format_optional_float(
                    row.wrapper_log_likelihood
                ),
                "log_likelihood_delta": _format_optional_float(
                    row.log_likelihood_delta
                ),
                "topology_equal": _format_optional_bool(row.topology_equal),
                "robinson_foulds_distance": _format_optional_int(
                    row.robinson_foulds_distance
                ),
                "support_split_count": len(row.support_comparisons),
                "max_support_percent_delta": _format_optional_float(
                    _max_support_percent_delta(row)
                ),
                "notes": " | ".join(row.notes),
            }
            for row in report.observations
        ],
    )


def write_maximum_likelihood_wrapper_comparison_run_json(
    path: Path,
    report: MaximumLikelihoodWrapperBenchmarkReport,
) -> Path:
    """Write the machine-readable maximum-likelihood wrapper benchmark payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "case_count": report.case_count,
                "supported_case_count": report.supported_case_count,
                "tolerance_match_case_count": report.tolerance_match_case_count,
                "expected_model_assumption_difference_case_count": (
                    report.expected_model_assumption_difference_case_count
                ),
                "unsupported_case_count": report.unsupported_case_count,
                "native_bug_case_count": report.native_bug_case_count,
                "blocking_case_count": report.blocking_case_count,
                "all_supported_cases_clear": report.all_supported_cases_clear,
                "limitations": report.limitations,
                "summary_rows": [asdict(row) for row in report.summary_rows],
                "observations": [asdict(row) for row in report.observations],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return path


def _format_optional_bool(value: bool | None) -> str:
    if value is None:
        return ""
    return "True" if value else "False"


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".15g")


def _format_optional_int(value: int | None) -> str:
    if value is None:
        return ""
    return str(value)


def _max_support_percent_delta(row) -> float | None:
    deltas = [
        abs(comparison.support_percent_delta)
        for comparison in row.support_comparisons
        if comparison.support_percent_delta is not None
    ]
    if not deltas:
        return None
    return max(deltas)
