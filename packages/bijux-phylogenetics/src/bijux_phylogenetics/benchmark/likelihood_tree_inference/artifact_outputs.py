from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .contracts import MaximumLikelihoodWrapperBenchmarkReport


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
