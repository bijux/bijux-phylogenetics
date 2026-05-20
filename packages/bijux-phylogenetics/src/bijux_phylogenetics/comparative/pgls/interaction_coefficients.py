from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative.pgls.design import inspect_pgls_inputs
from bijux_phylogenetics.comparative.pgls.fitting import run_pgls
from bijux_phylogenetics.comparative.pgls.models import (
    PGLSInputReport,
    PGLSPredictorClassification,
    PGLSResult,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.runtime.errors import ComparativeMethodError


@dataclass(slots=True)
class PGLSInteractionCoefficientRow:
    """One reviewer-facing interaction coefficient from a fitted PGLS model."""

    interaction_term: str
    interaction_kind: str
    coefficient_name: str
    component_terms: list[str]
    component_source_columns: list[str]
    component_kinds: list[str]
    component_columns: list[str]
    component_levels: list[str | None]
    omitted_reference_levels: list[str]
    estimate: float
    standard_error: float
    test_statistic: float
    p_value: float
    lower_95_confidence_interval: float
    upper_95_confidence_interval: float


@dataclass(slots=True)
class PGLSInteractionCoefficientReport:
    """Explicit interaction-coefficient review surface for one fitted PGLS model."""

    tree_path: Path
    traits_path: Path
    response: str
    formula: str
    lambda_value: float
    interaction_term_count: int
    rows: list[PGLSInteractionCoefficientRow]


def summarize_pgls_interaction_coefficients(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> PGLSInteractionCoefficientReport:
    """Summarize fitted PGLS interaction coefficients in an interpretable ledger."""
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
    return _build_interaction_report(input_report, model_report)


def write_pgls_interaction_coefficient_table(
    path: Path, report: PGLSInteractionCoefficientReport
) -> Path:
    """Write one flat interaction-coefficient ledger as TSV or CSV."""
    return write_taxon_rows(
        path,
        columns=[
            "interaction_term",
            "interaction_kind",
            "coefficient_name",
            "component_terms",
            "component_source_columns",
            "component_kinds",
            "component_columns",
            "component_levels",
            "omitted_reference_levels",
            "estimate",
            "standard_error",
            "test_statistic",
            "p_value",
            "lower_95_confidence_interval",
            "upper_95_confidence_interval",
        ],
        rows=[
            {
                "interaction_term": row.interaction_term,
                "interaction_kind": row.interaction_kind,
                "coefficient_name": row.coefficient_name,
                "component_terms": ";".join(row.component_terms),
                "component_source_columns": ";".join(row.component_source_columns),
                "component_kinds": ";".join(row.component_kinds),
                "component_columns": ";".join(row.component_columns),
                "component_levels": ";".join(
                    level or "" for level in row.component_levels
                ),
                "omitted_reference_levels": ";".join(row.omitted_reference_levels),
                "estimate": format(row.estimate, ".15g"),
                "standard_error": format(row.standard_error, ".15g"),
                "test_statistic": format(row.test_statistic, ".15g"),
                "p_value": format(row.p_value, ".15g"),
                "lower_95_confidence_interval": format(
                    row.lower_95_confidence_interval, ".15g"
                ),
                "upper_95_confidence_interval": format(
                    row.upper_95_confidence_interval, ".15g"
                ),
            }
            for row in report.rows
        ],
    )


def _build_interaction_report(
    input_report: PGLSInputReport, model_report: PGLSResult
) -> PGLSInteractionCoefficientReport:
    predictor_by_name = {row.name: row for row in input_report.predictors}
    coefficient_by_name = {row.name: row for row in model_report.coefficients}
    rows: list[PGLSInteractionCoefficientRow] = []
    for interaction_audit in input_report.formula_audit.interaction_terms:
        component_reports = [
            predictor_by_name[term] for term in interaction_audit.component_terms
        ]
        interaction_kind = "-by-".join(
            _interaction_component_kind(report) for report in component_reports
        )
        omitted_reference_levels = [
            f"{report.name}={report.reference_level}"
            for report in component_reports
            if report.kind != "numeric" and report.reference_level is not None
        ]
        for encoded_column in interaction_audit.encoded_columns:
            coefficient = coefficient_by_name.get(encoded_column)
            if coefficient is None:
                raise ComparativeMethodError(
                    f"interaction coefficient '{encoded_column}' is missing from the fitted PGLS model"
                )
            component_columns = encoded_column.split(":")
            if len(component_columns) != len(component_reports):
                raise ComparativeMethodError(
                    f"interaction coefficient '{encoded_column}' does not match the expected encoded interaction structure"
                )
            rows.append(
                PGLSInteractionCoefficientRow(
                    interaction_term=interaction_audit.term,
                    interaction_kind=interaction_kind,
                    coefficient_name=coefficient.name,
                    component_terms=list(interaction_audit.component_terms),
                    component_source_columns=[
                        report.source_column or report.name
                        for report in component_reports
                    ],
                    component_kinds=[
                        _interaction_component_kind(report)
                        for report in component_reports
                    ],
                    component_columns=component_columns,
                    component_levels=[
                        _extract_component_level(report, encoded_component)
                        for report, encoded_component in zip(
                            component_reports, component_columns, strict=True
                        )
                    ],
                    omitted_reference_levels=omitted_reference_levels,
                    estimate=coefficient.estimate,
                    standard_error=coefficient.standard_error,
                    test_statistic=coefficient.test_statistic,
                    p_value=coefficient.p_value,
                    lower_95_confidence_interval=coefficient.lower_95_confidence_interval,
                    upper_95_confidence_interval=coefficient.upper_95_confidence_interval,
                )
            )
    return PGLSInteractionCoefficientReport(
        tree_path=model_report.tree_path,
        traits_path=model_report.traits_path,
        response=model_report.response,
        formula=model_report.formula.formula,
        lambda_value=model_report.lambda_value,
        interaction_term_count=len(input_report.formula_audit.interaction_terms),
        rows=rows,
    )


def _interaction_component_kind(report: PGLSPredictorClassification) -> str:
    return "continuous" if report.kind == "numeric" else "categorical"


def _extract_component_level(
    report: PGLSPredictorClassification, encoded_component: str
) -> str | None:
    if report.kind == "numeric":
        return None
    prefix = f"{report.name}["
    if not encoded_component.startswith(prefix) or not encoded_component.endswith("]"):
        raise ComparativeMethodError(
            f"encoded interaction component '{encoded_component}' does not match categorical predictor '{report.name}'"
        )
    return encoded_component[len(prefix) : -1]
