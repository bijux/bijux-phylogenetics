from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .fitting import run_pgls
from .models import ComparativeHypothesisTestRow, ComparativeMultipleTestingReport


def run_pgls_multiple_testing(
    tree_path: Path,
    traits_path: Path,
    *,
    responses: list[str],
    predictors: list[str],
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
    adjustment_method: str = "benjamini-hochberg",
) -> ComparativeMultipleTestingReport:
    """Run repeated PGLS fits and adjust coefficient p-values across the family of tests."""
    if not responses:
        raise ComparativeMethodError(
            "multiple-testing analysis requires at least one response trait"
        )
    if adjustment_method != "benjamini-hochberg":
        raise ComparativeMethodError(
            "supported multiple-testing adjustment is 'benjamini-hochberg'"
        )
    rows: list[ComparativeHypothesisTestRow] = []
    for response in responses:
        report = run_pgls(
            tree_path,
            traits_path,
            response=response,
            predictors=predictors,
            taxon_column=taxon_column,
            lambda_value=lambda_value,
        )
        for coefficient in report.coefficients:
            if coefficient.name == "intercept":
                continue
            rows.append(
                ComparativeHypothesisTestRow(
                    response=response,
                    term=coefficient.name,
                    estimate=coefficient.estimate,
                    p_value=coefficient.p_value,
                    adjusted_p_value=coefficient.p_value,
                    significant=False,
                )
            )
    adjusted = _benjamini_hochberg_adjustment([row.p_value for row in rows])
    for row, adjusted_p_value in zip(rows, adjusted, strict=True):
        row.adjusted_p_value = adjusted_p_value
        row.significant = adjusted_p_value <= 0.05
    raw_significant_count = sum(1 for row in rows if row.p_value <= 0.05)
    adjusted_significant_count = sum(1 for row in rows if row.significant)
    return ComparativeMultipleTestingReport(
        tree_path=tree_path,
        traits_path=traits_path,
        responses=list(responses),
        predictors=list(predictors),
        adjustment_method=adjustment_method,
        family_size=len(rows),
        raw_significant_count=raw_significant_count,
        adjusted_significant_count=adjusted_significant_count,
        rows=rows,
    )


def _benjamini_hochberg_adjustment(p_values: list[float]) -> list[float]:
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    adjusted = [1.0] * len(p_values)
    running_minimum = 1.0
    total = len(p_values)
    for rank, (index, p_value) in enumerate(reversed(indexed), start=1):
        denominator = total - rank + 1
        candidate = min(1.0, (p_value * total) / denominator)
        running_minimum = min(running_minimum, candidate)
        adjusted[index] = running_minimum
    return adjusted
