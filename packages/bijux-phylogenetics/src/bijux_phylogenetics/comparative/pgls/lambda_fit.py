from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.pgls.fitting import run_pgls
from bijux_phylogenetics.comparative.pgls.models import (
    PGLSLambdaFitReport,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows


def summarize_pgls_lambda_fit(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> PGLSLambdaFitReport:
    """Summarize the Pagel lambda fit surface used by one PGLS model."""
    return run_pgls(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    ).lambda_fit


def write_pgls_lambda_profile_table(path: Path, report: PGLSLambdaFitReport) -> Path:
    """Write one lambda likelihood-profile table as TSV or CSV."""
    return write_taxon_rows(
        path,
        columns=[
            "mode",
            "lambda_value",
            "log_likelihood",
            "delta_log_likelihood",
            "within_95_confidence_interval",
            "profile_lower_95_confidence_interval",
            "profile_upper_95_confidence_interval",
            "null_log_likelihood",
            "brownian_log_likelihood",
        ],
        rows=[
            {
                "mode": report.mode,
                "lambda_value": format(row.lambda_value, ".15g"),
                "log_likelihood": format(row.log_likelihood, ".15g"),
                "delta_log_likelihood": format(row.delta_log_likelihood, ".15g"),
                "within_95_confidence_interval": str(
                    row.within_95_confidence_interval
                ).lower(),
                "profile_lower_95_confidence_interval": (
                    ""
                    if report.lower_95_confidence_interval is None
                    else format(report.lower_95_confidence_interval, ".15g")
                ),
                "profile_upper_95_confidence_interval": (
                    ""
                    if report.upper_95_confidence_interval is None
                    else format(report.upper_95_confidence_interval, ".15g")
                ),
                "null_log_likelihood": format(report.null_log_likelihood, ".15g"),
                "brownian_log_likelihood": format(
                    report.brownian_log_likelihood, ".15g"
                ),
            }
            for row in report.profile_rows
        ],
    )
