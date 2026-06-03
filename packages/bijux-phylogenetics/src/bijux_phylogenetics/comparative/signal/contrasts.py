from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative._math import (
    student_t_quantile,
    student_t_two_sided_p_value,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .core import (
    IndependentContrastReport,
    compute_phylogenetic_independent_contrasts,
)


@dataclass(slots=True)
class IndependentContrastRegressionRow:
    """One node-level row for PIC regression through the origin."""

    node: str
    predictor_contrast: float
    response_contrast: float
    fitted_response_contrast: float
    residual: float
    leverage_fraction: float


@dataclass(slots=True)
class IndependentContrastRegressionReport:
    """Regression-through-origin over matched phylogenetic independent contrasts."""

    tree_path: Path
    traits_path: Path
    response_trait: str
    predictor_trait: str
    contrast_count: int
    slope: float
    standard_error: float
    test_statistic: float
    p_value: float
    lower_95_confidence_interval: float
    upper_95_confidence_interval: float
    residual_sum_of_squares: float
    r_squared_through_origin: float
    rows: list[IndependentContrastRegressionRow]


def summarize_independent_contrast_regression(
    tree_path: Path,
    traits_path: Path,
    *,
    response_trait: str,
    predictor_trait: str,
    taxon_column: str | None = None,
) -> IndependentContrastRegressionReport:
    """Regress one PIC response trait on one PIC predictor trait through the origin."""
    response_report = compute_phylogenetic_independent_contrasts(
        tree_path,
        traits_path,
        trait=response_trait,
        taxon_column=taxon_column,
    )
    predictor_report = compute_phylogenetic_independent_contrasts(
        tree_path,
        traits_path,
        trait=predictor_trait,
        taxon_column=taxon_column,
    )
    response_by_node = {row.node: row.contrast for row in response_report.contrasts}
    predictor_by_node = {row.node: row.contrast for row in predictor_report.contrasts}
    if response_by_node.keys() != predictor_by_node.keys():
        raise ComparativeMethodError(
            "independent contrast regression requires matching internal nodes across both traits"
        )
    ordered_nodes = [row.node for row in response_report.contrasts]
    predictor_values = [predictor_by_node[node] for node in ordered_nodes]
    response_values = [response_by_node[node] for node in ordered_nodes]
    sum_xx = sum(value * value for value in predictor_values)
    if sum_xx == 0.0:
        raise ComparativeMethodError(
            "independent contrast regression requires non-zero predictor contrasts"
        )
    slope = (
        sum(
            predictor * response
            for predictor, response in zip(
                predictor_values, response_values, strict=True
            )
        )
        / sum_xx
    )
    fitted_values = [slope * value for value in predictor_values]
    residuals = [
        response - fitted
        for response, fitted in zip(response_values, fitted_values, strict=True)
    ]
    residual_sum_of_squares = sum(value * value for value in residuals)
    total_sum_of_squares = sum(value * value for value in response_values)
    r_squared = (
        1.0 - (residual_sum_of_squares / total_sum_of_squares)
        if total_sum_of_squares
        else 1.0
    )
    degrees_of_freedom = len(ordered_nodes) - 1
    if degrees_of_freedom <= 0:
        raise ComparativeMethodError(
            "independent contrast regression requires at least two contrasts"
        )
    residual_variance = residual_sum_of_squares / degrees_of_freedom
    standard_error = (residual_variance / sum_xx) ** 0.5
    test_statistic = slope / standard_error if standard_error else 0.0
    p_value = student_t_two_sided_p_value(test_statistic, degrees_of_freedom)
    critical_value = student_t_quantile(0.975, degrees_of_freedom)
    interval_radius = critical_value * standard_error
    leverage_denominator = sum_xx
    rows = [
        IndependentContrastRegressionRow(
            node=node,
            predictor_contrast=predictor,
            response_contrast=response,
            fitted_response_contrast=fitted,
            residual=residual,
            leverage_fraction=(predictor * predictor) / leverage_denominator,
        )
        for node, predictor, response, fitted, residual in zip(
            ordered_nodes,
            predictor_values,
            response_values,
            fitted_values,
            residuals,
            strict=True,
        )
    ]
    return IndependentContrastRegressionReport(
        tree_path=tree_path,
        traits_path=traits_path,
        response_trait=response_trait,
        predictor_trait=predictor_trait,
        contrast_count=len(ordered_nodes),
        slope=slope,
        standard_error=standard_error,
        test_statistic=test_statistic,
        p_value=p_value,
        lower_95_confidence_interval=slope - interval_radius,
        upper_95_confidence_interval=slope + interval_radius,
        residual_sum_of_squares=residual_sum_of_squares,
        r_squared_through_origin=r_squared,
        rows=rows,
    )


def write_independent_contrast_table(
    path: Path, report: IndependentContrastReport
) -> Path:
    """Write one flat phylogenetic independent-contrast ledger as TSV or CSV."""
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "node_id",
            "node",
            "left_taxa",
            "right_taxa",
            "contrast",
            "expected_variance",
            "ancestral_value",
            "root_estimate",
        ],
        rows=[
            {
                "trait": report.trait,
                "node_id": row.node_id,
                "node": row.node,
                "left_taxa": ",".join(row.left_taxa),
                "right_taxa": ",".join(row.right_taxa),
                "contrast": format(row.contrast, ".15g"),
                "expected_variance": format(row.expected_variance, ".15g"),
                "ancestral_value": format(row.ancestral_value, ".15g"),
                "root_estimate": format(report.root_estimate, ".15g"),
            }
            for row in report.contrasts
        ],
    )


def write_independent_contrast_regression_table(
    path: Path, report: IndependentContrastRegressionReport
) -> Path:
    """Write one regression-through-origin ledger over matched PIC values."""
    return write_taxon_rows(
        path,
        columns=[
            "response_trait",
            "predictor_trait",
            "node",
            "predictor_contrast",
            "response_contrast",
            "fitted_response_contrast",
            "residual",
            "leverage_fraction",
            "slope",
            "standard_error",
            "test_statistic",
            "p_value",
            "lower_95_confidence_interval",
            "upper_95_confidence_interval",
            "residual_sum_of_squares",
            "r_squared_through_origin",
        ],
        rows=[
            {
                "response_trait": report.response_trait,
                "predictor_trait": report.predictor_trait,
                "node": row.node,
                "predictor_contrast": format(row.predictor_contrast, ".15g"),
                "response_contrast": format(row.response_contrast, ".15g"),
                "fitted_response_contrast": format(
                    row.fitted_response_contrast, ".15g"
                ),
                "residual": format(row.residual, ".15g"),
                "leverage_fraction": format(row.leverage_fraction, ".15g"),
                "slope": format(report.slope, ".15g"),
                "standard_error": format(report.standard_error, ".15g"),
                "test_statistic": format(report.test_statistic, ".15g"),
                "p_value": format(report.p_value, ".15g"),
                "lower_95_confidence_interval": format(
                    report.lower_95_confidence_interval, ".15g"
                ),
                "upper_95_confidence_interval": format(
                    report.upper_95_confidence_interval, ".15g"
                ),
                "residual_sum_of_squares": format(
                    report.residual_sum_of_squares, ".15g"
                ),
                "r_squared_through_origin": format(
                    report.r_squared_through_origin, ".15g"
                ),
            }
            for row in report.rows
        ],
    )
