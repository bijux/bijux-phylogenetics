from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .contracts import CorrelatedTraitEvolutionReport
from .statistics import _format_optional


def write_correlated_trait_summary_table(
    path: Path,
    report: CorrelatedTraitEvolutionReport,
) -> Path:
    """Write one summary ledger for correlated trait evolution."""
    return write_taxon_rows(
        path,
        columns=[
            "analysis_kind",
            "left_trait",
            "right_trait",
            "taxon_column",
            "tree_taxon_count",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "observation_row_count",
            "association_measure_name",
            "association_measure_value",
            "evolutionary_covariance",
            "evolutionary_correlation",
            "lower_95_confidence_interval",
            "upper_95_confidence_interval",
            "independent_parameter_count",
            "independent_log_likelihood",
            "independent_aic",
            "correlated_parameter_count",
            "correlated_log_likelihood",
            "correlated_aic",
            "better_model",
            "likelihood_ratio_statistic",
            "likelihood_ratio_degrees_of_freedom",
            "likelihood_ratio_p_value",
            "likelihood_ratio_p_value_method",
            "left_root_estimate",
            "right_root_estimate",
            "left_state_order",
            "right_state_order",
            "joint_state_count",
            "warning_count",
        ],
        rows=[
            {
                "analysis_kind": report.analysis_kind,
                "left_trait": report.left_trait,
                "right_trait": report.right_trait,
                "taxon_column": report.taxon_column,
                "tree_taxon_count": report.tree_taxon_count,
                "analyzed_taxon_count": len(report.analyzed_taxa),
                "excluded_taxon_count": len(report.excluded_taxa),
                "observation_row_count": len(report.observation_rows),
                "association_measure_name": report.association_measure_name,
                "association_measure_value": format(
                    report.association_measure_value, ".15g"
                ),
                "evolutionary_covariance": _format_optional(
                    report.evolutionary_covariance
                ),
                "evolutionary_correlation": _format_optional(
                    report.evolutionary_correlation
                ),
                "lower_95_confidence_interval": _format_optional(
                    report.lower_95_confidence_interval
                ),
                "upper_95_confidence_interval": _format_optional(
                    report.upper_95_confidence_interval
                ),
                "independent_parameter_count": report.independent_parameter_count,
                "independent_log_likelihood": format(
                    report.independent_log_likelihood, ".15g"
                ),
                "independent_aic": format(report.independent_aic, ".15g"),
                "correlated_parameter_count": report.correlated_parameter_count,
                "correlated_log_likelihood": format(
                    report.correlated_log_likelihood, ".15g"
                ),
                "correlated_aic": format(report.correlated_aic, ".15g"),
                "better_model": report.better_model,
                "likelihood_ratio_statistic": format(
                    report.likelihood_ratio_statistic, ".15g"
                ),
                "likelihood_ratio_degrees_of_freedom": (
                    report.likelihood_ratio_degrees_of_freedom
                ),
                "likelihood_ratio_p_value": format(
                    report.likelihood_ratio_p_value, ".15g"
                ),
                "likelihood_ratio_p_value_method": (
                    report.likelihood_ratio_p_value_method
                ),
                "left_root_estimate": _format_optional(report.left_root_estimate),
                "right_root_estimate": _format_optional(report.right_root_estimate),
                "left_state_order": ",".join(report.left_state_order),
                "right_state_order": ",".join(report.right_state_order),
                "joint_state_count": len(report.joint_state_counts),
                "warning_count": len(report.warnings),
            }
        ],
    )


def write_correlated_trait_comparison_table(
    path: Path,
    report: CorrelatedTraitEvolutionReport,
) -> Path:
    """Write one independent-versus-correlated model comparison ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "model_kind",
            "model_description",
            "parameter_count",
            "log_likelihood",
            "aic",
            "delta_aic",
            "selected",
        ],
        rows=[
            {
                "model_kind": row.model_kind,
                "model_description": row.model_description,
                "parameter_count": row.parameter_count,
                "log_likelihood": format(row.log_likelihood, ".15g"),
                "aic": format(row.aic, ".15g"),
                "delta_aic": format(row.delta_aic, ".15g"),
                "selected": str(row.selected).lower(),
            }
            for row in report.comparison_rows
        ],
    )


def write_correlated_trait_observation_table(
    path: Path,
    report: CorrelatedTraitEvolutionReport,
) -> Path:
    """Write one detailed observation ledger for a coupling analysis."""
    return write_taxon_rows(
        path,
        columns=[
            "row_kind",
            "label",
            "taxon",
            "left_taxa",
            "right_taxa",
            "left_numeric_value",
            "right_numeric_value",
            "expected_variance",
            "left_state",
            "right_state",
            "joint_state",
        ],
        rows=[
            {
                "row_kind": row.row_kind,
                "label": row.label,
                "taxon": row.taxon or "",
                "left_taxa": ",".join(row.left_taxa),
                "right_taxa": ",".join(row.right_taxa),
                "left_numeric_value": _format_optional(row.left_numeric_value),
                "right_numeric_value": _format_optional(row.right_numeric_value),
                "expected_variance": _format_optional(row.expected_variance),
                "left_state": row.left_state or "",
                "right_state": row.right_state or "",
                "joint_state": row.joint_state or "",
            }
            for row in report.observation_rows
        ],
    )


def write_correlated_trait_exclusion_table(
    path: Path,
    report: CorrelatedTraitEvolutionReport,
) -> Path:
    """Write one excluded-taxon ledger for correlated trait evolution."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason", "missing_traits"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
                "missing_traits": ",".join(row.missing_traits),
            }
            for row in report.excluded_taxa
        ],
    )
