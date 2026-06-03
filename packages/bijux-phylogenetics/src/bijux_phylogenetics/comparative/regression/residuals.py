from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative.pgls import inspect_pgls_inputs, run_pgls
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

_OUTLIER_THRESHOLD = 2.0


@dataclass(slots=True)
class PhylogeneticResidualExclusion:
    """One taxon excluded before phylogenetic residual review."""

    taxon: str
    reason: str
    details: str


@dataclass(slots=True)
class PhylogeneticResidualCoefficientRow:
    """One fitted regression coefficient used to derive phylogenetic residuals."""

    name: str
    estimate: float
    standard_error: float
    p_value: float
    lower_95_confidence_interval: float
    upper_95_confidence_interval: float


@dataclass(slots=True)
class PhylogeneticResidualTaxonRow:
    """One taxon-level phylogenetic residual row."""

    taxon: str
    input_order: int
    tree_tip_label: str
    observed_value: float
    fitted_value: float
    residual: float
    standardized_residual: float
    abs_standardized_residual: float
    leverage: float
    outlier: bool


@dataclass(slots=True)
class PhylogeneticResidualReport:
    """Reviewer-facing phylogenetic residual summary for one response-predictor pair."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    response: str
    predictor: str
    method: str
    tree_taxon_count: int
    analyzed_taxa: list[str]
    analyzed_taxon_count: int
    excluded_taxa: list[PhylogeneticResidualExclusion]
    lambda_value: float
    lambda_estimation_mode: str
    log_likelihood: float
    aic: float
    coefficient_rows: list[PhylogeneticResidualCoefficientRow]
    taxon_rows: list[PhylogeneticResidualTaxonRow]
    outlier_threshold: float
    outlier_taxa: list[str]
    max_abs_standardized_residual: float | None
    warnings: list[str]


def summarize_phylogenetic_residuals(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str,
    predictor: str,
    taxon_column: str | None = None,
    method: str = "lambda",
) -> PhylogeneticResidualReport:
    """Summarize phylogenetically corrected fitted values and residuals for one predictor."""
    if method not in {"brownian", "lambda"}:
        raise ComparativeMethodError(
            "phylogenetic residual method must be 'brownian' or 'lambda'"
        )
    tree = load_tree(tree_path)
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    input_report = inspect_pgls_inputs(
        tree_path,
        traits_path,
        response=response,
        predictors=[predictor],
        taxon_column=taxon_column,
    )
    if not input_report.ready:
        raise ComparativeMethodError("; ".join(input_report.blockers))
    if (
        len(input_report.predictors) != 1
        or input_report.predictors[0].kind != "numeric"
    ):
        raise ComparativeMethodError(
            "phylogenetic residual review requires exactly one numeric predictor"
        )
    lambda_value: float | str = 1.0 if method == "brownian" else "estimate"
    model = run_pgls(
        tree_path,
        traits_path,
        response=response,
        predictors=[predictor],
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    fitted_by_taxon = {row.taxon: row for row in model.diagnostics.fitted_observed_rows}
    leverage_by_taxon = {row.taxon: row for row in model.diagnostics.leverage_rows}
    outlier_taxa = {row.taxon for row in model.diagnostics.outlier_taxa}
    input_order_taxa = [
        row[table.taxon_column]
        for row in table.rows
        if row[table.taxon_column] in fitted_by_taxon
    ]
    taxon_rows: list[PhylogeneticResidualTaxonRow] = []
    for input_order, taxon in enumerate(input_order_taxa, start=1):
        fitted_row = fitted_by_taxon[taxon]
        leverage_row = leverage_by_taxon[taxon]
        taxon_rows.append(
            PhylogeneticResidualTaxonRow(
                taxon=taxon,
                input_order=input_order,
                tree_tip_label=taxon,
                observed_value=fitted_row.observed,
                fitted_value=fitted_row.fitted,
                residual=fitted_row.residual,
                standardized_residual=leverage_row.standardized_residual,
                abs_standardized_residual=abs(leverage_row.standardized_residual),
                leverage=leverage_row.leverage,
                outlier=taxon in outlier_taxa,
            )
        )
    taxon_rows.sort(key=lambda row: row.input_order)
    exclusions = _build_exclusions(
        table=table,
        tree_taxa=set(tree.tip_names),
        analyzed_taxa=set(model.taxa),
        response=response,
        predictor=predictor,
    )
    warnings = list(input_report.warnings)
    if outlier_taxa:
        warnings.append(
            "one or more taxa have unusually large phylogenetically corrected residuals"
        )
    coefficient_rows = [
        PhylogeneticResidualCoefficientRow(
            name=coefficient.name,
            estimate=coefficient.estimate,
            standard_error=coefficient.standard_error,
            p_value=coefficient.p_value,
            lower_95_confidence_interval=coefficient.lower_95_confidence_interval,
            upper_95_confidence_interval=coefficient.upper_95_confidence_interval,
        )
        for coefficient in model.coefficients
    ]
    return PhylogeneticResidualReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        response=response,
        predictor=predictor,
        method=method,
        tree_taxon_count=tree.tip_count,
        analyzed_taxa=input_order_taxa,
        analyzed_taxon_count=len(input_order_taxa),
        excluded_taxa=exclusions,
        lambda_value=model.lambda_value,
        lambda_estimation_mode=model.lambda_fit.mode,
        log_likelihood=model.log_likelihood,
        aic=model.aic,
        coefficient_rows=coefficient_rows,
        taxon_rows=taxon_rows,
        outlier_threshold=_OUTLIER_THRESHOLD,
        outlier_taxa=[row.taxon for row in taxon_rows if row.outlier],
        max_abs_standardized_residual=(
            max((row.abs_standardized_residual for row in taxon_rows), default=None)
        ),
        warnings=warnings,
    )


def write_phylogenetic_residual_summary_table(
    path: Path,
    report: PhylogeneticResidualReport,
) -> Path:
    """Write one summary ledger for phylogenetic residual review."""
    return write_taxon_rows(
        path,
        columns=[
            "response",
            "predictor",
            "method",
            "taxon_column",
            "tree_taxon_count",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "lambda_value",
            "lambda_estimation_mode",
            "log_likelihood",
            "aic",
            "outlier_threshold",
            "outlier_count",
            "max_abs_standardized_residual",
        ],
        rows=[
            {
                "response": report.response,
                "predictor": report.predictor,
                "method": report.method,
                "taxon_column": report.taxon_column,
                "tree_taxon_count": str(report.tree_taxon_count),
                "analyzed_taxon_count": str(report.analyzed_taxon_count),
                "excluded_taxon_count": str(len(report.excluded_taxa)),
                "lambda_value": format(report.lambda_value, ".15g"),
                "lambda_estimation_mode": report.lambda_estimation_mode,
                "log_likelihood": format(report.log_likelihood, ".15g"),
                "aic": format(report.aic, ".15g"),
                "outlier_threshold": format(report.outlier_threshold, ".15g"),
                "outlier_count": str(len(report.outlier_taxa)),
                "max_abs_standardized_residual": _format_optional(
                    report.max_abs_standardized_residual
                ),
            }
        ],
    )


def write_phylogenetic_residual_taxon_table(
    path: Path,
    report: PhylogeneticResidualReport,
) -> Path:
    """Write one taxon-level residual ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "input_order",
            "tree_tip_label",
            "observed_value",
            "fitted_value",
            "residual",
            "standardized_residual",
            "abs_standardized_residual",
            "leverage",
            "outlier",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "input_order": str(row.input_order),
                "tree_tip_label": row.tree_tip_label,
                "observed_value": format(row.observed_value, ".15g"),
                "fitted_value": format(row.fitted_value, ".15g"),
                "residual": format(row.residual, ".15g"),
                "standardized_residual": format(row.standardized_residual, ".15g"),
                "abs_standardized_residual": format(
                    row.abs_standardized_residual, ".15g"
                ),
                "leverage": format(row.leverage, ".15g"),
                "outlier": "true" if row.outlier else "false",
            }
            for row in report.taxon_rows
        ],
    )


def write_phylogenetic_residual_coefficient_table(
    path: Path,
    report: PhylogeneticResidualReport,
) -> Path:
    """Write one coefficient ledger for phylogenetic residual review."""
    return write_taxon_rows(
        path,
        columns=[
            "name",
            "estimate",
            "standard_error",
            "p_value",
            "lower_95_confidence_interval",
            "upper_95_confidence_interval",
        ],
        rows=[
            {
                "name": row.name,
                "estimate": format(row.estimate, ".15g"),
                "standard_error": format(row.standard_error, ".15g"),
                "p_value": format(row.p_value, ".15g"),
                "lower_95_confidence_interval": format(
                    row.lower_95_confidence_interval, ".15g"
                ),
                "upper_95_confidence_interval": format(
                    row.upper_95_confidence_interval, ".15g"
                ),
            }
            for row in report.coefficient_rows
        ],
    )


def write_phylogenetic_residual_exclusion_table(
    path: Path,
    report: PhylogeneticResidualReport,
) -> Path:
    """Write one excluded-taxa ledger for phylogenetic residual review."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason", "details"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
                "details": row.details,
            }
            for row in report.excluded_taxa
        ],
    )


def _build_exclusions(
    *,
    table,
    tree_taxa: set[str],
    analyzed_taxa: set[str],
    response: str,
    predictor: str,
) -> list[PhylogeneticResidualExclusion]:
    excluded: list[PhylogeneticResidualExclusion] = []
    seen_taxa: set[str] = set()
    for row in table.rows:
        taxon = row[table.taxon_column]
        if taxon in seen_taxa:
            continue
        seen_taxa.add(taxon)
        if taxon not in tree_taxa:
            excluded.append(
                PhylogeneticResidualExclusion(
                    taxon=taxon,
                    reason="absent_from_tree",
                    details="taxon is present in the trait table but absent from the tree",
                )
            )
            continue
        if taxon in analyzed_taxa:
            continue
        missing_columns = [
            column for column in (response, predictor) if not row.get(column, "")
        ]
        reason = "missing_value" if missing_columns else "excluded_before_fit"
        details = (
            f"taxon is missing required value(s): {', '.join(missing_columns)}"
            if missing_columns
            else "taxon was excluded before residual fitting"
        )
        excluded.append(
            PhylogeneticResidualExclusion(
                taxon=taxon,
                reason=reason,
                details=details,
            )
        )
    for taxon in sorted(tree_taxa - seen_taxa):
        excluded.append(
            PhylogeneticResidualExclusion(
                taxon=taxon,
                reason="missing_from_trait_table",
                details="taxon is present in the tree but absent from the trait table",
            )
        )
    return excluded


def _format_optional(value: float | None) -> str:
    return "" if value is None else format(value, ".15g")
