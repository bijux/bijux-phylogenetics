from __future__ import annotations

import csv
from functools import cache
import json
import math
from pathlib import Path

from bijux_phylogenetics.ancestral import (
    reconstruct_continuous_evolutionary_mode_states,
)
from bijux_phylogenetics.comparative import (
    compare_continuous_evolutionary_modes,
    fit_continuous_evolutionary_mode,
    rescale_tree_early_burst,
    rescale_tree_ornstein_uhlenbeck,
)
from bijux_phylogenetics.comparative.pgls import run_pgls
from bijux_phylogenetics.comparative.signal import estimate_pagels_lambda

from .definitions import STUDY_ID, STUDY_ONE_REFERENCE_ROOT


def study_root(repo_root: Path) -> Path:
    return Path(repo_root) / "evidence-book" / "studies" / STUDY_ID


def source_reference_paths(repo_root: Path) -> tuple[Path, Path]:
    root = Path(repo_root) / STUDY_ONE_REFERENCE_ROOT
    return root / "reference_trimmed_primatetree.nwk", root / "reference_primate.csv"


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


@cache
def load_r_reference_results(repo_root: Path) -> dict[str, object]:
    return read_json(study_root(repo_root) / "reference" / "reference_results.json")


def rounded(value: float) -> float:
    return float(format(float(value), ".15g"))


def rounded_display(value: float, digits: int = 4) -> float:
    return round(float(value), digits)


def diagnostic_summary_from_series(
    fitted_values: list[float],
    residuals: list[float],
) -> dict[str, float | int]:
    residual_mean = sum(residuals) / len(residuals)
    residual_variance = sum((value - residual_mean) ** 2 for value in residuals) / max(
        1,
        len(residuals) - 1,
    )
    residual_sd = math.sqrt(max(residual_variance, 0.0))
    if residual_sd == 0.0:
        standardized = [0.0 for _ in residuals]
    else:
        standardized = [(value - residual_mean) / residual_sd for value in residuals]
    mean_fitted = sum(fitted_values) / len(fitted_values)
    abs_residuals = [abs(value) for value in residuals]
    mean_abs_residual = sum(abs_residuals) / len(abs_residuals)
    fitted_ss = sum((value - mean_fitted) ** 2 for value in fitted_values)
    abs_residual_ss = sum((value - mean_abs_residual) ** 2 for value in abs_residuals)
    if fitted_ss == 0.0 or abs_residual_ss == 0.0:
        abs_residual_fitted_correlation = 0.0
    else:
        covariance = sum(
            (left - mean_fitted) * (right - mean_abs_residual)
            for left, right in zip(fitted_values, abs_residuals, strict=True)
        ) / len(fitted_values)
        abs_residual_fitted_correlation = covariance / math.sqrt(
            fitted_ss / len(fitted_values) * abs_residual_ss / len(fitted_values)
        )

    ordered_residuals = sorted(residuals)
    quantiles = []
    for index in range(len(ordered_residuals)):
        probability = (index + 0.5) / len(ordered_residuals)
        quantiles.append(inverse_normal_cdf(probability))
    qq_correlation = pearson_correlation(quantiles, ordered_residuals)
    return {
        "residual_mean": rounded(residual_mean),
        "residual_variance": rounded(residual_variance),
        "residual_sd": rounded(residual_sd),
        "max_abs_z_residual": rounded(max(abs(value) for value in standardized)),
        "abs_residual_fitted_correlation": rounded(abs_residual_fitted_correlation),
        "qq_correlation": rounded(qq_correlation),
        "outlier_count_abs_z_ge_2": sum(
            1 for value in standardized if abs(value) >= 2.0
        ),
    }


def inverse_normal_cdf(probability: float) -> float:
    return math.sqrt(2.0) * inverse_error(2.0 * probability - 1.0)


def inverse_error(value: float) -> float:
    # Winitzki approximation is sufficient for reviewer-facing QQ summaries.
    if value <= -1.0:
        return float("-inf")
    if value >= 1.0:
        return float("inf")
    a = 0.147
    signed = 1.0 if value >= 0.0 else -1.0
    ln = math.log(1.0 - value * value)
    term = (2.0 / (math.pi * a)) + (ln / 2.0)
    return signed * math.sqrt(math.sqrt(term * term - (ln / a)) - term)


def pearson_correlation(left: list[float], right: list[float]) -> float:
    mean_left = sum(left) / len(left)
    mean_right = sum(right) / len(right)
    left_ss = sum((value - mean_left) ** 2 for value in left)
    right_ss = sum((value - mean_right) ** 2 for value in right)
    if left_ss == 0.0 or right_ss == 0.0:
        return 0.0
    covariance = sum(
        (left_value - mean_left) * (right_value - mean_right)
        for left_value, right_value in zip(left, right, strict=True)
    )
    return covariance / math.sqrt(left_ss * right_ss)


def r_squared(observed: list[float], fitted: list[float]) -> float:
    mean_observed = sum(observed) / len(observed)
    total = sum((value - mean_observed) ** 2 for value in observed)
    residual = sum(
        (value - fit) ** 2 for value, fit in zip(observed, fitted, strict=True)
    )
    if total == 0.0:
        return 1.0
    return 1.0 - (residual / total)


def ordered_trait_values(
    traits_path: Path,
    taxa: list[str],
    *,
    trait: str,
    taxon_column: str,
) -> list[float]:
    with traits_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    values_by_taxon = {
        row[taxon_column]: float(row[trait])
        for row in rows
        if row.get(taxon_column) and row.get(trait)
    }
    return [rounded_display(values_by_taxon[taxon]) for taxon in taxa]


def tree_rescaling_payload(report: object) -> dict[str, object]:
    branch_rows = []
    for row in report.branch_rows:
        descendant_taxa: object = list(row.descendant_taxa)
        if len(row.descendant_taxa) == 1:
            descendant_taxa = row.descendant_taxa[0]
        branch_rows.append(
            {
                "node": row.node,
                "descendant_taxa": descendant_taxa,
                "branch_length": rounded_display(row.transformed_branch_length),
                "parent_depth": rounded_display(row.parent_depth),
                "child_depth": rounded_display(row.child_depth),
            }
        )
    return {
        "branch_count": len(branch_rows),
        "total_branch_length": rounded_display(report.transformed_total_branch_length),
        "branch_rows": branch_rows,
    }


def continuous_mode_fit_payload(
    report: object,
    *,
    parameter_key: str | None,
    parameter_count: int,
    tip_values: list[float] | None = None,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "root_state": rounded_display(report.root_state),
        "rate": rounded_display(report.rate),
        "log_likelihood": rounded_display(report.log_likelihood),
        "aic": rounded_display(report.aic),
        "parameter_count": parameter_count,
    }
    if tip_values is not None:
        payload["tip_values"] = tip_values
    if parameter_key is not None and report.parameter_value is not None:
        payload[parameter_key] = rounded_display(report.parameter_value)
    return payload


def likelihood_ratio_payload(report: object) -> dict[str, object]:
    return {
        "statistic": rounded_display(report.statistic),
        "p_value": rounded(report.p_value),
    }


def ancestral_reconstruction_payload(
    report: object,
    *,
    parameter_key: str | None = None,
) -> dict[str, object]:
    internal_rows = [
        estimate for estimate in report.reconstruction.estimates if not estimate.is_tip
    ]
    rows = [
        {
            "node_index": index,
            "node": estimate.node,
            "estimate": rounded_display(estimate.estimate),
        }
        for index, estimate in enumerate(internal_rows, start=1)
    ]
    payload: dict[str, object] = {
        "node_count": len(rows),
        "first_five_estimates": [row["estimate"] for row in rows[:5]],
        "recent_five_estimates": [row["estimate"] for row in rows[-5:]],
        "rows": rows,
    }
    if parameter_key is not None and report.parameter_value is not None:
        payload[parameter_key] = rounded_display(report.parameter_value)
    return payload


@cache
def load_python_results(repo_root: Path) -> dict[str, object]:
    baseline = baseline_gls_report(repo_root)
    fixed_reference = fixed_reference_lambda_pgls_report(repo_root)
    estimated = estimated_lambda_pgls_report(repo_root)
    signal = signal_report(repo_root)
    brownian_fit, ou_fit, early_burst_fit = continuous_mode_fit_reports(repo_root)
    mode_comparison = mode_comparison_report(repo_root)
    transformed_tree_reports = transformed_tree_reports_for_repo(repo_root)
    brownian_ancestral, early_burst_ancestral = ancestral_reconstruction_reports(
        repo_root
    )
    tip_values = ordered_trait_values(
        source_reference_paths(repo_root)[1],
        brownian_fit.taxa,
        trait="longevity",
        taxon_column="species",
    )
    return {
        "source_contract": source_contract_payload(baseline),
        "baseline_gls": baseline_gls_payload(baseline),
        "fixed_reference_lambda_pgls": fixed_reference_lambda_pgls_payload(
            fixed_reference
        ),
        "estimated_lambda_pgls": estimated_lambda_pgls_payload(estimated),
        "signal_test": signal_test_payload(signal),
        "tree_rescaling": tree_rescaling_payloads(transformed_tree_reports),
        "continuous_mode_fits": continuous_mode_fit_payloads(
            brownian_fit,
            ou_fit,
            early_burst_fit,
            tip_values=tip_values,
        ),
        "likelihood_ratio_tests": likelihood_ratio_test_payloads(mode_comparison),
        "ancestral_reconstruction": ancestral_reconstruction_payloads(
            brownian_ancestral,
            early_burst_ancestral,
        ),
        "coverage_boundaries": coverage_boundary_payload(),
    }


@cache
def baseline_gls_report(repo_root: Path):
    tree_path, traits_path = source_reference_paths(repo_root)
    return run_pgls(
        tree_path,
        traits_path,
        response="longevity",
        predictors=["social_group_size"],
        taxon_column="species",
        lambda_value=0.0,
    )


@cache
def fixed_reference_lambda_pgls_report(repo_root: Path):
    tree_path, traits_path = source_reference_paths(repo_root)
    reference_results = load_r_reference_results(repo_root)
    lambda_value = float(
        reference_results["fixed_reference_lambda_pgls"]["lambda_value"]
    )
    return run_pgls(
        tree_path,
        traits_path,
        response="longevity",
        predictors=["social_group_size"],
        taxon_column="species",
        lambda_value=lambda_value,
    )


@cache
def estimated_lambda_pgls_report(repo_root: Path):
    tree_path, traits_path = source_reference_paths(repo_root)
    return run_pgls(
        tree_path,
        traits_path,
        response="longevity",
        predictors=["social_group_size"],
        taxon_column="species",
        lambda_value="estimate",
    )


@cache
def signal_report(repo_root: Path):
    tree_path, traits_path = source_reference_paths(repo_root)
    return estimate_pagels_lambda(
        tree_path,
        traits_path,
        trait="longevity",
        taxon_column="species",
    )


@cache
def transformed_tree_reports_for_repo(repo_root: Path) -> dict[str, object]:
    tree_path, _ = source_reference_paths(repo_root)
    return {
        "ou_alpha_1": rescale_tree_ornstein_uhlenbeck(tree_path, alpha=1.0),
        "ou_alpha_10": rescale_tree_ornstein_uhlenbeck(tree_path, alpha=10.0),
        "early_burst_2": rescale_tree_early_burst(tree_path, rate_change=2.0),
        "late_burst_minus_2": rescale_tree_early_burst(tree_path, rate_change=-2.0),
    }


@cache
def continuous_mode_fit_reports(repo_root: Path):
    tree_path, traits_path = source_reference_paths(repo_root)
    brownian_fit = fit_continuous_evolutionary_mode(
        tree_path,
        traits_path,
        trait="longevity",
        taxon_column="species",
        mode="brownian",
    )
    ou_fit = fit_continuous_evolutionary_mode(
        tree_path,
        traits_path,
        trait="longevity",
        taxon_column="species",
        mode="ornstein-uhlenbeck",
        ou_bounds=(1e-6, 10.0),
    )
    early_burst_fit = fit_continuous_evolutionary_mode(
        tree_path,
        traits_path,
        trait="longevity",
        taxon_column="species",
        mode="early-burst",
        early_burst_bounds=(1e-6, 50.0),
    )
    return brownian_fit, ou_fit, early_burst_fit


@cache
def mode_comparison_report(repo_root: Path):
    tree_path, traits_path = source_reference_paths(repo_root)
    return compare_continuous_evolutionary_modes(
        tree_path,
        traits_path,
        trait="longevity",
        taxon_column="species",
        ou_bounds=(1e-6, 10.0),
        early_burst_bounds=(1e-6, 50.0),
    )


@cache
def ancestral_reconstruction_reports(repo_root: Path):
    tree_path, traits_path = source_reference_paths(repo_root)
    brownian_ancestral = reconstruct_continuous_evolutionary_mode_states(
        tree_path,
        traits_path,
        trait="longevity",
        taxon_column="species",
        mode="brownian",
    )
    early_burst_ancestral = reconstruct_continuous_evolutionary_mode_states(
        tree_path,
        traits_path,
        trait="longevity",
        taxon_column="species",
        mode="early-burst",
        rate_change=-2.0,
    )
    return brownian_ancestral, early_burst_ancestral


def source_contract_payload(baseline) -> dict[str, object]:
    return {
        "row_count": baseline.taxon_count,
        "tip_count": baseline.taxon_count,
        "predictor": "social_group_size",
        "response": "longevity",
    }


def _gls_like_payload(report) -> dict[str, object]:
    return {
        "coefficients": {
            row.name: rounded(row.estimate) for row in report.coefficients
        },
        "standard_errors": {
            row.name: rounded(row.standard_error) for row in report.coefficients
        },
        "p_values": {row.name: rounded(row.p_value) for row in report.coefficients},
        "log_likelihood": rounded(report.log_likelihood),
        "aic": rounded(report.aic),
        "r_squared": rounded(report.r_squared),
        "diagnostics": diagnostic_summary_from_series(
            report.fitted_values,
            report.residuals,
        ),
    }


def baseline_gls_payload(baseline) -> dict[str, object]:
    return _gls_like_payload(baseline)


def fixed_reference_lambda_pgls_payload(fixed_report) -> dict[str, object]:
    payload = _gls_like_payload(fixed_report)
    payload["lambda_value"] = rounded(fixed_report.lambda_value)
    return payload


def estimated_lambda_pgls_payload(estimated) -> dict[str, object]:
    payload = _gls_like_payload(estimated)
    payload["lambda_value"] = rounded(estimated.lambda_value)
    return payload


def signal_test_payload(signal) -> dict[str, object]:
    return {
        "estimated_lambda": rounded(signal.lambda_value),
        "estimated_log_likelihood": rounded(signal.log_likelihood),
        "null_log_likelihood": rounded(signal.null_log_likelihood),
        "likelihood_ratio": rounded(
            -2.0 * (signal.null_log_likelihood - signal.log_likelihood)
        ),
        "p_value": rounded(
            math.erfc(
                math.sqrt(
                    max(
                        0.0,
                        -2.0 * (signal.null_log_likelihood - signal.log_likelihood),
                    )
                    / 2.0
                )
            )
        ),
    }


def tree_rescaling_payloads(
    transformed_tree_reports: dict[str, object],
) -> dict[str, object]:
    return {
        key: tree_rescaling_payload(report)
        for key, report in transformed_tree_reports.items()
    }


def continuous_mode_fit_payloads(
    brownian_fit,
    ou_fit,
    early_burst_fit,
    *,
    tip_values: dict[str, float],
) -> dict[str, object]:
    return {
        "brownian": continuous_mode_fit_payload(
            brownian_fit,
            parameter_key=None,
            parameter_count=2,
            tip_values=tip_values,
        ),
        "ornstein_uhlenbeck": continuous_mode_fit_payload(
            ou_fit,
            parameter_key="alpha",
            parameter_count=3,
        ),
        "early_burst": continuous_mode_fit_payload(
            early_burst_fit,
            parameter_key="rate_change",
            parameter_count=3,
        ),
    }


def likelihood_ratio_test_payloads(mode_comparison) -> dict[str, object]:
    return {
        report.comparison_id.replace("-", "_"): likelihood_ratio_payload(report)
        for report in mode_comparison.likelihood_ratio_tests
    }


def ancestral_reconstruction_payloads(
    brownian_ancestral,
    early_burst_ancestral,
) -> dict[str, object]:
    return {
        "brownian": ancestral_reconstruction_payload(brownian_ancestral),
        "early_burst": ancestral_reconstruction_payload(
            early_burst_ancestral,
            parameter_key="rate_change",
        ),
    }


def coverage_boundary_payload() -> dict[str, object]:
    return {
        "uncovered_fragments": ["mode-linked-intercept-models"],
        "notes": [
            "The lecture corBlomberg likelihood sweep remains outside the current canonical runtime parity surface.",
            "The governed evidence closes transformed-tree, fitContinuous, likelihood-ratio, and ancestral-state parity without overstating the remaining intercept-mode boundary.",
        ],
    }
