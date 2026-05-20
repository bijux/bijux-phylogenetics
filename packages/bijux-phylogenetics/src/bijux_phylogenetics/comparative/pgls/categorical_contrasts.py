from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative.pgls.design import inspect_pgls_inputs
from bijux_phylogenetics.comparative.pgls.fitting import run_pgls
from bijux_phylogenetics.comparative.pgls.models import (
    ComparativeFormulaSpecification,
    PGLSCoefficient,
)
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree


@dataclass(slots=True)
class PGLSCategoricalContrastRow:
    """One categorical level interpreted against the fitted comparative model."""

    predictor: str
    source_column: str
    encoding_scheme: str
    level: str
    baseline_level: str | None
    is_reference_level: bool
    coefficient_name: str | None
    coefficient_estimate: float | None
    standard_error: float | None
    test_statistic: float | None
    p_value: float | None
    lower_95_confidence_interval: float | None
    upper_95_confidence_interval: float | None
    observed_taxon_count: int
    missing_category_taxa: list[str]


@dataclass(slots=True)
class PGLSCategoricalContrastReport:
    """Reviewer-facing summary of categorical predictor encoding and coefficients."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    formula: ComparativeFormulaSpecification
    categorical_predictor_count: int
    rows: list[PGLSCategoricalContrastRow]


def summarize_pgls_categorical_contrasts(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> PGLSCategoricalContrastReport:
    """Summarize how categorical predictors are encoded and estimated in one PGLS fit."""
    input_report = inspect_pgls_inputs(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
    )
    model_report = run_pgls(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    coefficient_by_name = {
        coefficient.name: coefficient for coefficient in model_report.coefficients
    }
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    tree = load_tree(tree_path)
    rows: list[PGLSCategoricalContrastRow] = []
    for predictor_report in input_report.predictors:
        if predictor_report.kind == "numeric":
            continue
        source_column = predictor_report.source_column or predictor_report.name
        missing_category_taxa = sorted(
            [
                taxon
                for taxon in tree.tip_names
                if taxon not in rows_by_taxon
                or not rows_by_taxon[taxon].get(source_column, "").strip()
            ]
        )
        encoding_scheme = (
            "reference-level"
            if input_report.formula.include_intercept
            else "full-indicator"
        )
        if input_report.formula.include_intercept:
            rows.append(
                PGLSCategoricalContrastRow(
                    predictor=predictor_report.name,
                    source_column=source_column,
                    encoding_scheme=encoding_scheme,
                    level=predictor_report.reference_level or "",
                    baseline_level=predictor_report.reference_level,
                    is_reference_level=True,
                    coefficient_name=None,
                    coefficient_estimate=None,
                    standard_error=None,
                    test_statistic=None,
                    p_value=None,
                    lower_95_confidence_interval=None,
                    upper_95_confidence_interval=None,
                    observed_taxon_count=(predictor_report.level_counts or {}).get(
                        predictor_report.reference_level or "", 0
                    ),
                    missing_category_taxa=missing_category_taxa,
                )
            )
        for coefficient_name in predictor_report.encoded_columns or []:
            level = coefficient_name.removeprefix(
                f"{predictor_report.name}["
            ).removesuffix("]")
            coefficient = coefficient_by_name[coefficient_name]
            rows.append(
                _build_contrast_row(
                    predictor=predictor_report.name,
                    source_column=source_column,
                    encoding_scheme=encoding_scheme,
                    level=level,
                    baseline_level=predictor_report.reference_level,
                    coefficient_name=coefficient_name,
                    coefficient=coefficient,
                    observed_taxon_count=(predictor_report.level_counts or {}).get(
                        level, 0
                    ),
                    missing_category_taxa=missing_category_taxa,
                )
            )
    return PGLSCategoricalContrastReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        formula=input_report.formula,
        categorical_predictor_count=len(input_report.categorical_predictors),
        rows=rows,
    )


def write_pgls_categorical_contrast_table(
    path: Path, report: PGLSCategoricalContrastReport
) -> Path:
    """Write categorical contrast rows as CSV or TSV."""
    return write_taxon_rows(
        path,
        columns=[
            "predictor",
            "source_column",
            "encoding_scheme",
            "level",
            "baseline_level",
            "is_reference_level",
            "coefficient_name",
            "coefficient_estimate",
            "standard_error",
            "test_statistic",
            "p_value",
            "lower_95_confidence_interval",
            "upper_95_confidence_interval",
            "observed_taxon_count",
            "missing_category_taxa",
        ],
        rows=[
            {
                "predictor": row.predictor,
                "source_column": row.source_column,
                "encoding_scheme": row.encoding_scheme,
                "level": row.level,
                "baseline_level": row.baseline_level or "",
                "is_reference_level": "true" if row.is_reference_level else "false",
                "coefficient_name": row.coefficient_name or "",
                "coefficient_estimate": _format_optional_float(
                    row.coefficient_estimate
                ),
                "standard_error": _format_optional_float(row.standard_error),
                "test_statistic": _format_optional_float(row.test_statistic),
                "p_value": _format_optional_float(row.p_value),
                "lower_95_confidence_interval": _format_optional_float(
                    row.lower_95_confidence_interval
                ),
                "upper_95_confidence_interval": _format_optional_float(
                    row.upper_95_confidence_interval
                ),
                "observed_taxon_count": str(row.observed_taxon_count),
                "missing_category_taxa": ",".join(row.missing_category_taxa),
            }
            for row in report.rows
        ],
    )


def _build_contrast_row(
    *,
    predictor: str,
    source_column: str,
    encoding_scheme: str,
    level: str,
    baseline_level: str | None,
    coefficient_name: str,
    coefficient: PGLSCoefficient,
    observed_taxon_count: int,
    missing_category_taxa: list[str],
) -> PGLSCategoricalContrastRow:
    return PGLSCategoricalContrastRow(
        predictor=predictor,
        source_column=source_column,
        encoding_scheme=encoding_scheme,
        level=level,
        baseline_level=baseline_level,
        is_reference_level=False,
        coefficient_name=coefficient_name,
        coefficient_estimate=coefficient.estimate,
        standard_error=coefficient.standard_error,
        test_statistic=coefficient.test_statistic,
        p_value=coefficient.p_value,
        lower_95_confidence_interval=coefficient.lower_95_confidence_interval,
        upper_95_confidence_interval=coefficient.upper_95_confidence_interval,
        observed_taxon_count=observed_taxon_count,
        missing_category_taxa=missing_category_taxa,
    )


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".15g")
