from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .core import (
    BlombergKReport,
    PagelLambdaReport,
    PhylogeneticSignalInputAudit,
    PhylogeneticSignalTestReport,
    compute_blombergs_k,
    compute_phylogenetic_signal_test,
    estimate_pagels_lambda,
)


@dataclass(slots=True)
class PhylogeneticSignalSummaryReport:
    """Reviewer-facing summary over comparative phylogenetic signal metrics."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    input_audit: PhylogeneticSignalInputAudit
    blombergs_k: BlombergKReport
    pagels_lambda: PagelLambdaReport
    signal_test: PhylogeneticSignalTestReport
    lambda_likelihood_ratio_statistic: float
    lambda_likelihood_ratio_p_value: float
    lambda_p_value_method: str


def summarize_phylogenetic_signal(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    permutations: int = 199,
    seed: int = 1,
) -> PhylogeneticSignalSummaryReport:
    """Summarize K, lambda, and permutation evidence for one numeric trait."""
    blomberg = compute_blombergs_k(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    lambda_report = estimate_pagels_lambda(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    signal_test = compute_phylogenetic_signal_test(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        permutations=permutations,
        seed=seed,
    )
    return PhylogeneticSignalSummaryReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=blomberg.taxon_count,
        input_audit=blomberg.input_audit,
        blombergs_k=blomberg,
        pagels_lambda=lambda_report,
        signal_test=signal_test,
        lambda_likelihood_ratio_statistic=lambda_report.likelihood_ratio_statistic,
        lambda_likelihood_ratio_p_value=lambda_report.likelihood_ratio_p_value,
        lambda_p_value_method=lambda_report.p_value_method,
    )


def write_phylogenetic_signal_summary_table(
    path: Path, report: PhylogeneticSignalSummaryReport
) -> Path:
    """Write one flat summary table for a phylogenetic signal workflow."""
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_count",
            "blombergs_k",
            "blombergs_generalized_mean",
            "blombergs_observed_mean_square",
            "blombergs_phylogenetic_mean_square",
            "blombergs_expected_mean_square_ratio",
            "signal_permutation_p_value",
            "pagels_lambda",
            "lambda_log_likelihood",
            "lambda_null_log_likelihood",
            "lambda_brownian_log_likelihood",
            "lambda_likelihood_ratio_statistic",
            "lambda_likelihood_ratio_p_value",
            "lambda_p_value_method",
            "permutations",
            "permuted_k_at_or_above_observed",
            "signal_null_k_minimum",
            "signal_null_k_mean",
            "signal_null_k_maximum",
        ],
        rows=[
            {
                "trait": report.trait,
                "taxon_count": report.taxon_count,
                "blombergs_k": format(report.blombergs_k.k, ".15g"),
                "blombergs_generalized_mean": format(
                    report.blombergs_k.generalized_mean, ".15g"
                ),
                "blombergs_observed_mean_square": format(
                    report.blombergs_k.observed_mean_square, ".15g"
                ),
                "blombergs_phylogenetic_mean_square": format(
                    report.blombergs_k.phylogenetic_mean_square, ".15g"
                ),
                "blombergs_expected_mean_square_ratio": format(
                    report.blombergs_k.expected_mean_square_ratio, ".15g"
                ),
                "signal_permutation_p_value": format(
                    report.signal_test.p_value, ".15g"
                ),
                "pagels_lambda": format(report.pagels_lambda.lambda_value, ".15g"),
                "lambda_log_likelihood": format(
                    report.pagels_lambda.log_likelihood, ".15g"
                ),
                "lambda_null_log_likelihood": format(
                    report.pagels_lambda.null_log_likelihood, ".15g"
                ),
                "lambda_brownian_log_likelihood": format(
                    report.pagels_lambda.brownian_log_likelihood, ".15g"
                ),
                "lambda_likelihood_ratio_statistic": format(
                    report.lambda_likelihood_ratio_statistic, ".15g"
                ),
                "lambda_likelihood_ratio_p_value": format(
                    report.lambda_likelihood_ratio_p_value, ".15g"
                ),
                "lambda_p_value_method": report.lambda_p_value_method,
                "permutations": report.signal_test.permutations,
                "permuted_k_at_or_above_observed": (
                    report.signal_test.permuted_k_at_or_above_observed
                ),
                "signal_null_k_minimum": format(
                    report.signal_test.null_distribution_minimum, ".15g"
                ),
                "signal_null_k_mean": format(
                    report.signal_test.null_distribution_mean, ".15g"
                ),
                "signal_null_k_maximum": format(
                    report.signal_test.null_distribution_maximum, ".15g"
                ),
            }
        ],
    )


def write_phylogenetic_signal_permutation_table(
    path: Path, report: PhylogeneticSignalSummaryReport
) -> Path:
    """Write one permutation ledger for the Blomberg-K signal test."""
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "observed_k",
            "estimated_lambda",
            "permutations",
            "signal_permutation_p_value",
            "permutation_index",
            "permuted_k",
            "at_or_above_observed",
        ],
        rows=[
            {
                "trait": report.trait,
                "observed_k": format(report.signal_test.observed_k, ".15g"),
                "estimated_lambda": format(report.signal_test.estimated_lambda, ".15g"),
                "permutations": report.signal_test.permutations,
                "signal_permutation_p_value": format(
                    report.signal_test.p_value, ".15g"
                ),
                "permutation_index": row.permutation_index,
                "permuted_k": format(row.permuted_k, ".15g"),
                "at_or_above_observed": str(row.at_or_above_observed).lower(),
            }
            for row in report.signal_test.permutation_rows
        ],
    )
