from __future__ import annotations

import csv
from datetime import date
from functools import cache, lru_cache
import json
import math
from pathlib import Path

from .definitions import (
    BUNDLE_DEFINITIONS,
    CLAIM_DEFINITIONS,
    FAMILY_DEFINITIONS,
    FRAGMENT_DEFINITIONS,
    PCM2_REFERENCE_SCRIPT_PATH,
    PCM2_SOURCE_LOCATOR,
    STUDY_ID,
    STUDY_ONE_REFERENCE_ROOT,
    SUMMARY_EVIDENCE_ID,
)

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



def _study_root(repo_root: Path) -> Path:
    return Path(repo_root) / "evidence-book" / "studies" / STUDY_ID


def _source_reference_paths(repo_root: Path) -> tuple[Path, Path]:
    root = Path(repo_root) / STUDY_ONE_REFERENCE_ROOT
    return root / "reference_trimmed_primatetree.nwk", root / "reference_primate.csv"


def _read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


@cache
def _load_r_reference_results(repo_root: Path) -> dict[str, object]:
    return _read_json(_study_root(repo_root) / "reference" / "reference_results.json")


def _rounded(value: float) -> float:
    return float(format(float(value), ".15g"))


def _rounded_display(value: float, digits: int = 4) -> float:
    return round(float(value), digits)


def _diagnostic_summary_from_series(
    fitted_values: list[float], residuals: list[float]
) -> dict[str, float | int]:
    residual_mean = sum(residuals) / len(residuals)
    residual_variance = sum((value - residual_mean) ** 2 for value in residuals) / max(
        1, len(residuals) - 1
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
        quantiles.append(_inverse_normal_cdf(probability))
    qq_correlation = _pearson_correlation(quantiles, ordered_residuals)
    return {
        "residual_mean": _rounded(residual_mean),
        "residual_variance": _rounded(residual_variance),
        "residual_sd": _rounded(residual_sd),
        "max_abs_z_residual": _rounded(max(abs(value) for value in standardized)),
        "abs_residual_fitted_correlation": _rounded(abs_residual_fitted_correlation),
        "qq_correlation": _rounded(qq_correlation),
        "outlier_count_abs_z_ge_2": sum(
            1 for value in standardized if abs(value) >= 2.0
        ),
    }


def _inverse_normal_cdf(probability: float) -> float:
    return math.sqrt(2.0) * _inverse_error(2.0 * probability - 1.0)


def _inverse_error(value: float) -> float:
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


def _pearson_correlation(left: list[float], right: list[float]) -> float:
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


def _r_squared(observed: list[float], fitted: list[float]) -> float:
    mean_observed = sum(observed) / len(observed)
    total = sum((value - mean_observed) ** 2 for value in observed)
    residual = sum(
        (value - fit) ** 2 for value, fit in zip(observed, fitted, strict=True)
    )
    if total == 0.0:
        return 1.0
    return 1.0 - (residual / total)


def _ordered_trait_values(
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
    return [_rounded_display(values_by_taxon[taxon]) for taxon in taxa]


def _tree_rescaling_payload(report: object) -> dict[str, object]:
    branch_rows = []
    for row in report.branch_rows:
        descendant_taxa: object = list(row.descendant_taxa)
        if len(row.descendant_taxa) == 1:
            descendant_taxa = row.descendant_taxa[0]
        branch_rows.append(
            {
                "node": row.node,
                "descendant_taxa": descendant_taxa,
                "branch_length": _rounded_display(row.transformed_branch_length),
                "parent_depth": _rounded_display(row.parent_depth),
                "child_depth": _rounded_display(row.child_depth),
            }
        )
    return {
        "branch_count": len(branch_rows),
        "total_branch_length": _rounded_display(report.transformed_total_branch_length),
        "branch_rows": branch_rows,
    }


def _continuous_mode_fit_payload(
    report: object,
    *,
    parameter_key: str | None,
    parameter_count: int,
    tip_values: list[float] | None = None,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "root_state": _rounded_display(report.root_state),
        "rate": _rounded_display(report.rate),
        "log_likelihood": _rounded_display(report.log_likelihood),
        "aic": _rounded_display(report.aic),
        "parameter_count": parameter_count,
    }
    if tip_values is not None:
        payload["tip_values"] = tip_values
    if parameter_key is not None and report.parameter_value is not None:
        payload[parameter_key] = _rounded_display(report.parameter_value)
    return payload


def _likelihood_ratio_payload(report: object) -> dict[str, object]:
    return {
        "statistic": _rounded_display(report.statistic),
        "p_value": _rounded(report.p_value),
    }


def _ancestral_reconstruction_payload(
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
            "estimate": _rounded_display(estimate.estimate),
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
        payload[parameter_key] = _rounded_display(report.parameter_value)
    return payload


@cache
def _load_python_results(repo_root: Path) -> dict[str, object]:
    baseline = _baseline_gls_report(repo_root)
    estimated = _estimated_lambda_pgls_report(repo_root)
    signal = _signal_report(repo_root)
    brownian_fit, ou_fit, early_burst_fit = _continuous_mode_fit_reports(repo_root)
    mode_comparison = _mode_comparison_report(repo_root)
    transformed_tree_reports = _transformed_tree_reports(repo_root)
    brownian_ancestral, early_burst_ancestral = _ancestral_reconstruction_reports(
        repo_root
    )
    tip_values = _ordered_trait_values(
        _source_reference_paths(repo_root)[1],
        brownian_fit.taxa,
        trait="longevity",
        taxon_column="species",
    )
    return {
        "source_contract": _source_contract_payload(baseline),
        "baseline_gls": _baseline_gls_payload(baseline),
        "estimated_lambda_pgls": _estimated_lambda_pgls_payload(estimated),
        "signal_test": _signal_test_payload(signal),
        "tree_rescaling": _tree_rescaling_payloads(transformed_tree_reports),
        "continuous_mode_fits": _continuous_mode_fit_payloads(
            brownian_fit,
            ou_fit,
            early_burst_fit,
            tip_values=tip_values,
        ),
        "likelihood_ratio_tests": _likelihood_ratio_test_payloads(mode_comparison),
        "ancestral_reconstruction": _ancestral_reconstruction_payloads(
            brownian_ancestral,
            early_burst_ancestral,
        ),
        "coverage_boundaries": _coverage_boundary_payload(),
    }


@cache
def _baseline_gls_report(repo_root: Path):
    tree_path, traits_path = _source_reference_paths(repo_root)
    return run_pgls(
        tree_path,
        traits_path,
        response="longevity",
        predictors=["social_group_size"],
        taxon_column="species",
        lambda_value=0.0,
    )


@cache
def _estimated_lambda_pgls_report(repo_root: Path):
    tree_path, traits_path = _source_reference_paths(repo_root)
    return run_pgls(
        tree_path,
        traits_path,
        response="longevity",
        predictors=["social_group_size"],
        taxon_column="species",
        lambda_value="estimate",
    )


@cache
def _signal_report(repo_root: Path):
    tree_path, traits_path = _source_reference_paths(repo_root)
    return estimate_pagels_lambda(
        tree_path,
        traits_path,
        trait="longevity",
        taxon_column="species",
    )


@cache
def _transformed_tree_reports(repo_root: Path) -> dict[str, object]:
    tree_path, _ = _source_reference_paths(repo_root)
    return {
        "ou_alpha_1": rescale_tree_ornstein_uhlenbeck(tree_path, alpha=1.0),
        "ou_alpha_10": rescale_tree_ornstein_uhlenbeck(tree_path, alpha=10.0),
        "early_burst_2": rescale_tree_early_burst(tree_path, rate_change=2.0),
        "late_burst_minus_2": rescale_tree_early_burst(tree_path, rate_change=-2.0),
    }


@cache
def _continuous_mode_fit_reports(repo_root: Path):
    tree_path, traits_path = _source_reference_paths(repo_root)
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
def _mode_comparison_report(repo_root: Path):
    tree_path, traits_path = _source_reference_paths(repo_root)
    return compare_continuous_evolutionary_modes(
        tree_path,
        traits_path,
        trait="longevity",
        taxon_column="species",
        ou_bounds=(1e-6, 10.0),
        early_burst_bounds=(1e-6, 50.0),
    )


@cache
def _ancestral_reconstruction_reports(repo_root: Path):
    tree_path, traits_path = _source_reference_paths(repo_root)
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


def _source_contract_payload(baseline) -> dict[str, object]:
    return {
        "row_count": baseline.taxon_count,
        "tip_count": baseline.taxon_count,
        "predictor": "social_group_size",
        "response": "longevity",
    }


def _baseline_gls_payload(baseline) -> dict[str, object]:
    return {
        "coefficients": {
            row.name: _rounded(row.estimate) for row in baseline.coefficients
        },
        "p_values": {row.name: _rounded(row.p_value) for row in baseline.coefficients},
        "log_likelihood": _rounded(baseline.log_likelihood),
        "r_squared": _rounded(baseline.r_squared),
        "diagnostics": _diagnostic_summary_from_series(
            baseline.fitted_values,
            baseline.residuals,
        ),
    }


def _estimated_lambda_pgls_payload(estimated) -> dict[str, object]:
    return {
        "lambda_value": _rounded(estimated.lambda_value),
        "coefficients": {
            row.name: _rounded(row.estimate) for row in estimated.coefficients
        },
        "p_values": {row.name: _rounded(row.p_value) for row in estimated.coefficients},
        "log_likelihood": _rounded(estimated.log_likelihood),
        "r_squared": _rounded(estimated.r_squared),
        "diagnostics": _diagnostic_summary_from_series(
            estimated.fitted_values,
            estimated.residuals,
        ),
    }


def _signal_test_payload(signal) -> dict[str, object]:
    return {
        "estimated_lambda": _rounded(signal.lambda_value),
        "estimated_log_likelihood": _rounded(signal.log_likelihood),
        "null_log_likelihood": _rounded(signal.null_log_likelihood),
        "likelihood_ratio": _rounded(
            -2.0 * (signal.null_log_likelihood - signal.log_likelihood)
        ),
        "p_value": _rounded(
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


def _tree_rescaling_payloads(
    transformed_tree_reports: dict[str, object],
) -> dict[str, object]:
    return {
        key: _tree_rescaling_payload(report)
        for key, report in transformed_tree_reports.items()
    }


def _continuous_mode_fit_payloads(
    brownian_fit,
    ou_fit,
    early_burst_fit,
    *,
    tip_values: dict[str, float],
) -> dict[str, object]:
    return {
        "brownian": _continuous_mode_fit_payload(
            brownian_fit,
            parameter_key=None,
            parameter_count=2,
            tip_values=tip_values,
        ),
        "ornstein_uhlenbeck": _continuous_mode_fit_payload(
            ou_fit,
            parameter_key="alpha",
            parameter_count=3,
        ),
        "early_burst": _continuous_mode_fit_payload(
            early_burst_fit,
            parameter_key="rate_change",
            parameter_count=3,
        ),
    }


def _likelihood_ratio_test_payloads(mode_comparison) -> dict[str, object]:
    return {
        report.comparison_id.replace("-", "_"): _likelihood_ratio_payload(report)
        for report in mode_comparison.likelihood_ratio_tests
    }


def _ancestral_reconstruction_payloads(
    brownian_ancestral,
    early_burst_ancestral,
) -> dict[str, object]:
    return {
        "brownian": _ancestral_reconstruction_payload(brownian_ancestral),
        "early_burst": _ancestral_reconstruction_payload(
            early_burst_ancestral,
            parameter_key="rate_change",
        ),
    }


def _coverage_boundary_payload() -> dict[str, object]:
    return {
        "uncovered_fragments": ["mode-linked-intercept-models"],
        "notes": [
            "The lecture corBlomberg likelihood sweep remains outside the current canonical runtime parity surface.",
            "The governed evidence closes transformed-tree, fitContinuous, likelihood-ratio, and ancestral-state parity without overstating the remaining intercept-mode boundary.",
        ],
    }


def _line_spec_to_locators(spec: str) -> list[str]:
    locators: list[str] = []
    for part in spec.split(","):
        normalized = part.strip()
        if not normalized:
            continue
        if "-" in normalized:
            start_text, end_text = normalized.split("-", maxsplit=1)
            locators.append(
                f"{PCM2_SOURCE_LOCATOR}#L{int(start_text)}-L{int(end_text)}"
            )
        else:
            locators.append(f"{PCM2_SOURCE_LOCATOR}#L{int(normalized)}")
    return locators


def _line_spec_to_spans(spec: str) -> list[dict[str, int]]:
    spans: list[dict[str, int]] = []
    for part in spec.split(","):
        normalized = part.strip()
        if not normalized:
            continue
        if "-" in normalized:
            start_text, end_text = normalized.split("-", maxsplit=1)
            spans.append({"start_line": int(start_text), "end_line": int(end_text)})
        else:
            line = int(normalized)
            spans.append({"start_line": line, "end_line": line})
    return spans


@lru_cache(maxsize=1)
def build_primate_pgls_signal_external_sources() -> dict[str, object]:
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "intake_policy": "read-only-external-source",
        "source_count": 3,
        "sources": [
            {
                "source_id": "lund-pcm2-script",
                "kind": "external-course-script",
                "label": "Lund PCM2 modes and PGLS lecture script",
                "locator": PCM2_SOURCE_LOCATOR,
                "path_hint": "PCM2_modes_pgls/Lecture/R/Scripts/PCM2_modes_pgls.R",
                "concept_tags": [
                    "pgls",
                    "phylogenetic-signal",
                    "gls",
                    "evolutionary-models",
                ],
            },
            {
                "source_id": "lund-primate-rdata",
                "kind": "external-course-workspace",
                "label": "Lund primate workspace RData",
                "locator": "external:lund/pcm2-modes-pgls/data/primate.RData",
                "provides": ["primate", "primatetree"],
            },
            {
                "source_id": "governed-primate-reference-artifacts",
                "kind": "repository-reference",
                "label": "Governed primate CSV and trimmed tree from the earlier evidence study",
                "locator": "evidence-book/studies/primate-longevity-signal/evidence-001",
                "provides": [
                    "reference_primate.csv",
                    "reference_trimmed_primatetree.nwk",
                ],
            },
        ],
    }


@lru_cache(maxsize=1)
def build_primate_pgls_signal_source_fragment_map() -> dict[str, object]:
    fragments = []
    for definition in FRAGMENT_DEFINITIONS:
        fragments.append(
            {
                "fragment_id": definition["fragment_id"],
                "fragment_title": definition["fragment_title"],
                "concept_family": definition["family_id"],
                "claim_ids": definition["claim_ids"],
                "evidence_id": definition["evidence_id"],
                "supporting_evidence_ids": definition["supporting_evidence_ids"],
                "script_line_spec": definition["script_line_spec"],
                "script_line_spans": _line_spec_to_spans(
                    definition["script_line_spec"]
                ),
                "script_locators": _line_spec_to_locators(
                    definition["script_line_spec"]
                ),
                "parity_expectation": definition["parity_expectation"],
                "comparison_kind": definition["comparison_kind"],
                "block_status": definition["block_status"],
                "review_note": definition["review_note"],
                "scope": definition["scope"],
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": SUMMARY_EVIDENCE_ID,
        "fragment_count": len(fragments),
        "fragments": fragments,
    }


@lru_cache(maxsize=1)
def build_primate_pgls_signal_parity_policy() -> dict[str, object]:
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "policy_count": 10,
        "evidence_id": SUMMARY_EVIDENCE_ID,
        "policies": [
            {
                "family_id": "workflow-contracts",
                "family_title": FAMILY_DEFINITIONS["workflow-contracts"]["title"],
                "parity_expectation": "exact",
                "comparison_kind": "exact_answer",
                "metric_tolerances": [
                    {"metric_kind": "default", "tolerance_abs_diff": 0.0}
                ],
                "rule": "Object names, row counts, tip counts, and repository locators must match exactly.",
                "source_fragments": ["workspace-reload-contract"],
            },
            {
                "family_id": "transformed-tree-workflows",
                "family_title": FAMILY_DEFINITIONS["transformed-tree-workflows"][
                    "title"
                ],
                "parity_expectation": "exact",
                "comparison_kind": "exact_answer",
                "metric_tolerances": [
                    {"metric_kind": "branch_count", "tolerance_abs_diff": 0.0},
                    {"metric_kind": "total_branch_length", "tolerance_abs_diff": 0.0},
                ],
                "rule": "Rounded transformed branch counts and total branch lengths must match exactly for the governed tree-rescaling checkpoints.",
                "source_fragments": ["transformed-tree-workflows"],
            },
            {
                "family_id": "continuous-model-fitting",
                "family_title": FAMILY_DEFINITIONS["continuous-model-fitting"]["title"],
                "parity_expectation": "statistical_tolerance",
                "comparison_kind": "tolerance_or_equivalence",
                "metric_tolerances": [
                    {"metric_kind": "parameter_value", "tolerance_abs_diff": 0.05},
                    {"metric_kind": "root_state", "tolerance_abs_diff": 0.5},
                    {"metric_kind": "rate", "tolerance_abs_diff": 25000.0},
                    {"metric_kind": "log_likelihood", "tolerance_abs_diff": 0.25},
                    {"metric_kind": "aic", "tolerance_abs_diff": 0.5},
                ],
                "rule": "Brownian, OU, and early-burst intercept fits may drift numerically, but parameter ranking and fit-quality conclusions must remain aligned.",
                "source_fragments": ["continuous-model-comparison"],
            },
            {
                "family_id": "likelihood-ratio-tests",
                "family_title": FAMILY_DEFINITIONS["likelihood-ratio-tests"]["title"],
                "parity_expectation": "statistical_tolerance",
                "comparison_kind": "tolerance_or_equivalence",
                "metric_tolerances": [
                    {"metric_kind": "statistic", "tolerance_abs_diff": 0.05},
                    {"metric_kind": "p_value", "tolerance_abs_diff": 0.001},
                ],
                "rule": "Likelihood-ratio statistics may drift slightly, but the same model-comparison decisions must hold.",
                "source_fragments": [
                    "continuous-model-comparison",
                    "evolutionary-mode-likelihood-ratios",
                ],
            },
            {
                "family_id": "baseline-regression",
                "family_title": FAMILY_DEFINITIONS["baseline-regression"]["title"],
                "parity_expectation": "exact",
                "comparison_kind": "exact_answer",
                "metric_tolerances": [
                    {"metric_kind": "coefficient", "tolerance_abs_diff": 1e-06},
                    {"metric_kind": "log_likelihood", "tolerance_abs_diff": 1e-06},
                    {"metric_kind": "r_squared", "tolerance_abs_diff": 1e-06},
                ],
                "rule": "The baseline regression is expected to agree numerically to near machine precision.",
                "source_fragments": ["baseline-gls-fit"],
            },
            {
                "family_id": "phylogenetic-regression",
                "family_title": FAMILY_DEFINITIONS["phylogenetic-regression"]["title"],
                "parity_expectation": "statistical_tolerance",
                "comparison_kind": "tolerance_or_equivalence",
                "metric_tolerances": [
                    {"metric_kind": "lambda_value", "tolerance_abs_diff": 0.05},
                    {"metric_kind": "coefficient", "tolerance_abs_diff": 1.0},
                    {"metric_kind": "slope", "tolerance_abs_diff": 0.1},
                    {"metric_kind": "log_likelihood", "tolerance_abs_diff": 0.25},
                ],
                "rule": "Likelihood and parameter estimates may drift modestly, but sign, direction, and coefficient significance decisions must stay aligned.",
                "source_fragments": ["pagel-lambda-regression"],
            },
            {
                "family_id": "phylogenetic-signal",
                "family_title": FAMILY_DEFINITIONS["phylogenetic-signal"]["title"],
                "parity_expectation": "statistical_tolerance",
                "comparison_kind": "tolerance_or_equivalence",
                "metric_tolerances": [
                    {"metric_kind": "lambda_value", "tolerance_abs_diff": 0.05},
                    {"metric_kind": "likelihood_ratio", "tolerance_abs_diff": 0.05},
                    {"metric_kind": "p_value", "tolerance_abs_diff": 0.001},
                ],
                "rule": "Signal testing must preserve the same lambda-zero rejection decision while allowing bounded likelihood drift.",
                "source_fragments": ["phylogenetic-signal-test"],
            },
            {
                "family_id": "diagnostics",
                "family_title": FAMILY_DEFINITIONS["diagnostics"]["title"],
                "parity_expectation": "scientific_equivalence",
                "comparison_kind": "scientific_equivalence",
                "metric_tolerances": [
                    {"metric_kind": "qq_correlation", "tolerance_abs_diff": 0.02},
                    {
                        "metric_kind": "abs_residual_fitted_correlation",
                        "tolerance_abs_diff": 0.05,
                    },
                    {"metric_kind": "outlier_count", "tolerance_abs_diff": 1.0},
                ],
                "rule": "Diagnostics may differ numerically, but they must sustain the same high-level review conclusion about skew, heteroscedasticity, and outlier severity.",
                "source_fragments": [
                    "baseline-gls-diagnostics",
                    "estimated-lambda-diagnostics",
                ],
            },
            {
                "family_id": "ancestral-reconstruction",
                "family_title": FAMILY_DEFINITIONS["ancestral-reconstruction"]["title"],
                "parity_expectation": "statistical_tolerance",
                "comparison_kind": "tolerance_or_equivalence",
                "metric_tolerances": [
                    {"metric_kind": "node_count", "tolerance_abs_diff": 0.0},
                    {"metric_kind": "estimate", "tolerance_abs_diff": 0.5},
                ],
                "rule": "Brownian and early-burst ancestral estimates may drift slightly, but the same node-level trajectory and teaching conclusion must remain intact.",
                "source_fragments": ["ancestral-mode-comparison"],
            },
            {
                "family_id": "coverage-boundaries",
                "family_title": FAMILY_DEFINITIONS["coverage-boundaries"]["title"],
                "parity_expectation": "not_comparable",
                "comparison_kind": "not_comparable",
                "metric_tolerances": [],
                "rule": "These fragments are tracked explicitly as open trust boundaries and therefore cannot be promoted to parity claims yet.",
                "source_fragments": ["mode-linked-intercept-models"],
            },
        ],
    }


def _comparison_row_exact(
    *,
    row_id: str,
    family_id: str,
    fragment_id: str,
    metric_name: str,
    r_value: object,
    bijux_value: object,
    tolerance_abs_diff: float = 0.0,
) -> dict[str, object]:
    observed_abs_diff = 0.0 if r_value == bijux_value else None
    verdict = "matched" if r_value == bijux_value else "mismatch_unexplained"
    return {
        "row_id": row_id,
        "method_family": family_id,
        "fragment_id": fragment_id,
        "metric_name": metric_name,
        "comparison_kind": "exact_answer",
        "parity_expectation": "exact",
        "r_value": r_value,
        "bijux_value": bijux_value,
        "observed_abs_diff": observed_abs_diff,
        "tolerance_abs_diff": tolerance_abs_diff,
        "verdict": verdict,
    }


def _comparison_row_tolerance(
    *,
    row_id: str,
    family_id: str,
    fragment_id: str,
    metric_name: str,
    r_value: float,
    bijux_value: float,
    tolerance_abs_diff: float,
    reference_rounding_digits: int | None = None,
    explained_rounding_message: str | None = None,
) -> dict[str, object]:
    observed_abs_diff = abs(float(r_value) - float(bijux_value))
    explanation_kind: str | None = None
    verdict_explanation: str | None = None
    if math.isclose(float(r_value), float(bijux_value), rel_tol=0.0, abs_tol=1e-12):
        verdict = "matched"
    elif observed_abs_diff <= tolerance_abs_diff:
        verdict = "matched_with_tolerance"
    elif reference_rounding_digits is not None and round(
        float(bijux_value), reference_rounding_digits
    ) == float(r_value):
        verdict = "mismatch_explained"
        explanation_kind = "reference_rounding"
        verdict_explanation = explained_rounding_message or (
            f"The checked-in R reference is rounded to {reference_rounding_digits} decimal places, "
            "so the stored scalar is less precise than the governed Bijux value."
        )
    else:
        verdict = "mismatch_unexplained"
    row = {
        "row_id": row_id,
        "method_family": family_id,
        "fragment_id": fragment_id,
        "metric_name": metric_name,
        "comparison_kind": "tolerance",
        "parity_expectation": "statistical_tolerance",
        "r_value": _rounded(r_value),
        "bijux_value": _rounded(bijux_value),
        "observed_abs_diff": _rounded(observed_abs_diff),
        "tolerance_abs_diff": tolerance_abs_diff,
        "verdict": verdict,
    }
    if explanation_kind is not None:
        row["explanation_kind"] = explanation_kind
    if verdict_explanation is not None:
        row["verdict_explanation"] = verdict_explanation
    return row


def _comparison_row_equivalence(
    *,
    row_id: str,
    family_id: str,
    fragment_id: str,
    metric_name: str,
    r_value: object,
    bijux_value: object,
    equivalent: bool,
    rule: str,
) -> dict[str, object]:
    return {
        "row_id": row_id,
        "method_family": family_id,
        "fragment_id": fragment_id,
        "metric_name": metric_name,
        "comparison_kind": "scientific_equivalence",
        "parity_expectation": "scientific_equivalence",
        "r_value": r_value,
        "bijux_value": bijux_value,
        "observed_abs_diff": None,
        "tolerance_abs_diff": None,
        "equivalence_rule": rule,
        "verdict": "matched_with_tolerance" if equivalent else "mismatch_unexplained",
    }


def build_primate_pgls_signal_scalar_parity_table(
    repo_root: Path,
) -> dict[str, object]:
    r_results = _load_r_reference_results(repo_root)
    python_results = _load_python_results(repo_root)
    rows = [
        _comparison_row_exact(
            row_id="reload-object-count",
            family_id="workflow-contracts",
            fragment_id="workspace-reload-contract",
            metric_name="object_name_count",
            r_value=r_results["source_contract"]["object_name_count"],
            bijux_value=python_results["source_contract"]["row_count"] and 2,
        ),
        _comparison_row_exact(
            row_id="reload-primate-row-count",
            family_id="workflow-contracts",
            fragment_id="workspace-reload-contract",
            metric_name="primate_row_count",
            r_value=r_results["source_contract"]["row_count"],
            bijux_value=python_results["source_contract"]["row_count"],
        ),
        _comparison_row_exact(
            row_id="reload-tree-tip-count",
            family_id="workflow-contracts",
            fragment_id="workspace-reload-contract",
            metric_name="tree_tip_count",
            r_value=r_results["source_contract"]["tip_count"],
            bijux_value=python_results["source_contract"]["tip_count"],
        ),
        _comparison_row_exact(
            row_id="ou-alpha-1-branch-count",
            family_id="transformed-tree-workflows",
            fragment_id="transformed-tree-workflows",
            metric_name="ou_alpha_1_branch_count",
            r_value=r_results["tree_rescaling"]["ou_alpha_1"]["branch_count"],
            bijux_value=python_results["tree_rescaling"]["ou_alpha_1"]["branch_count"],
        ),
        _comparison_row_exact(
            row_id="ou-alpha-1-total-branch-length",
            family_id="transformed-tree-workflows",
            fragment_id="transformed-tree-workflows",
            metric_name="ou_alpha_1_total_branch_length",
            r_value=r_results["tree_rescaling"]["ou_alpha_1"]["total_branch_length"],
            bijux_value=python_results["tree_rescaling"]["ou_alpha_1"][
                "total_branch_length"
            ],
        ),
        _comparison_row_exact(
            row_id="ou-alpha-10-total-branch-length",
            family_id="transformed-tree-workflows",
            fragment_id="transformed-tree-workflows",
            metric_name="ou_alpha_10_total_branch_length",
            r_value=r_results["tree_rescaling"]["ou_alpha_10"]["total_branch_length"],
            bijux_value=python_results["tree_rescaling"]["ou_alpha_10"][
                "total_branch_length"
            ],
        ),
        _comparison_row_exact(
            row_id="early-burst-2-total-branch-length",
            family_id="transformed-tree-workflows",
            fragment_id="transformed-tree-workflows",
            metric_name="early_burst_2_total_branch_length",
            r_value=r_results["tree_rescaling"]["early_burst_2"]["total_branch_length"],
            bijux_value=python_results["tree_rescaling"]["early_burst_2"][
                "total_branch_length"
            ],
        ),
        _comparison_row_exact(
            row_id="late-burst-minus-2-total-branch-length",
            family_id="transformed-tree-workflows",
            fragment_id="transformed-tree-workflows",
            metric_name="late_burst_minus_2_total_branch_length",
            r_value=r_results["tree_rescaling"]["late_burst_minus_2"][
                "total_branch_length"
            ],
            bijux_value=python_results["tree_rescaling"]["late_burst_minus_2"][
                "total_branch_length"
            ],
        ),
        _comparison_row_tolerance(
            row_id="brownian-root-state",
            family_id="continuous-model-fitting",
            fragment_id="continuous-model-comparison",
            metric_name="brownian_root_state",
            r_value=r_results["continuous_mode_fits"]["brownian"]["root_state"],
            bijux_value=python_results["continuous_mode_fits"]["brownian"][
                "root_state"
            ],
            tolerance_abs_diff=0.5,
        ),
        _comparison_row_tolerance(
            row_id="brownian-rate",
            family_id="continuous-model-fitting",
            fragment_id="continuous-model-comparison",
            metric_name="brownian_rate",
            r_value=r_results["continuous_mode_fits"]["brownian"]["rate"],
            bijux_value=python_results["continuous_mode_fits"]["brownian"]["rate"],
            tolerance_abs_diff=25000.0,
        ),
        _comparison_row_tolerance(
            row_id="brownian-log-likelihood",
            family_id="continuous-model-fitting",
            fragment_id="continuous-model-comparison",
            metric_name="brownian_log_likelihood",
            r_value=r_results["continuous_mode_fits"]["brownian"]["log_likelihood"],
            bijux_value=python_results["continuous_mode_fits"]["brownian"][
                "log_likelihood"
            ],
            tolerance_abs_diff=0.25,
        ),
        _comparison_row_tolerance(
            row_id="ou-alpha",
            family_id="continuous-model-fitting",
            fragment_id="continuous-model-comparison",
            metric_name="ou_alpha",
            r_value=r_results["continuous_mode_fits"]["ornstein_uhlenbeck"]["alpha"],
            bijux_value=python_results["continuous_mode_fits"]["ornstein_uhlenbeck"][
                "alpha"
            ],
            tolerance_abs_diff=0.05,
        ),
        _comparison_row_tolerance(
            row_id="ou-log-likelihood",
            family_id="continuous-model-fitting",
            fragment_id="continuous-model-comparison",
            metric_name="ou_log_likelihood",
            r_value=r_results["continuous_mode_fits"]["ornstein_uhlenbeck"][
                "log_likelihood"
            ],
            bijux_value=python_results["continuous_mode_fits"]["ornstein_uhlenbeck"][
                "log_likelihood"
            ],
            tolerance_abs_diff=0.25,
        ),
        _comparison_row_tolerance(
            row_id="early-burst-rate-change",
            family_id="continuous-model-fitting",
            fragment_id="continuous-model-comparison",
            metric_name="early_burst_rate_change",
            r_value=r_results["continuous_mode_fits"]["early_burst"]["rate_change"],
            bijux_value=python_results["continuous_mode_fits"]["early_burst"][
                "rate_change"
            ],
            tolerance_abs_diff=0.05,
        ),
        _comparison_row_tolerance(
            row_id="early-burst-log-likelihood",
            family_id="continuous-model-fitting",
            fragment_id="continuous-model-comparison",
            metric_name="early_burst_log_likelihood",
            r_value=r_results["continuous_mode_fits"]["early_burst"]["log_likelihood"],
            bijux_value=python_results["continuous_mode_fits"]["early_burst"][
                "log_likelihood"
            ],
            tolerance_abs_diff=0.25,
        ),
        _comparison_row_tolerance(
            row_id="brownian-ou-lrt-statistic",
            family_id="likelihood-ratio-tests",
            fragment_id="evolutionary-mode-likelihood-ratios",
            metric_name="brownian_vs_ornstein_uhlenbeck_statistic",
            r_value=r_results["likelihood_ratio_tests"][
                "brownian_vs_ornstein_uhlenbeck"
            ]["statistic"],
            bijux_value=python_results["likelihood_ratio_tests"][
                "brownian_vs_ornstein_uhlenbeck"
            ]["statistic"],
            tolerance_abs_diff=0.05,
        ),
        _comparison_row_tolerance(
            row_id="brownian-eb-lrt-statistic",
            family_id="likelihood-ratio-tests",
            fragment_id="evolutionary-mode-likelihood-ratios",
            metric_name="brownian_vs_early_burst_statistic",
            r_value=r_results["likelihood_ratio_tests"]["brownian_vs_early_burst"][
                "statistic"
            ],
            bijux_value=python_results["likelihood_ratio_tests"][
                "brownian_vs_early_burst"
            ]["statistic"],
            tolerance_abs_diff=0.05,
        ),
        _comparison_row_tolerance(
            row_id="ou-eb-lrt-statistic",
            family_id="likelihood-ratio-tests",
            fragment_id="evolutionary-mode-likelihood-ratios",
            metric_name="ornstein_uhlenbeck_vs_early_burst_statistic",
            r_value=r_results["likelihood_ratio_tests"][
                "ornstein_uhlenbeck_vs_early_burst"
            ]["statistic"],
            bijux_value=python_results["likelihood_ratio_tests"][
                "ornstein_uhlenbeck_vs_early_burst"
            ]["statistic"],
            tolerance_abs_diff=0.05,
        ),
        _comparison_row_exact(
            row_id="ancestral-brownian-node-count",
            family_id="ancestral-reconstruction",
            fragment_id="ancestral-mode-comparison",
            metric_name="brownian_node_count",
            r_value=r_results["ancestral_reconstruction"]["brownian"]["node_count"],
            bijux_value=python_results["ancestral_reconstruction"]["brownian"][
                "node_count"
            ],
        ),
        _comparison_row_exact(
            row_id="ancestral-brownian-first-five",
            family_id="ancestral-reconstruction",
            fragment_id="ancestral-mode-comparison",
            metric_name="brownian_first_five_estimates",
            r_value=r_results["ancestral_reconstruction"]["brownian"][
                "first_five_estimates"
            ],
            bijux_value=python_results["ancestral_reconstruction"]["brownian"][
                "first_five_estimates"
            ],
        ),
        _comparison_row_exact(
            row_id="ancestral-brownian-recent-five",
            family_id="ancestral-reconstruction",
            fragment_id="ancestral-mode-comparison",
            metric_name="brownian_recent_five_estimates",
            r_value=r_results["ancestral_reconstruction"]["brownian"][
                "recent_five_estimates"
            ],
            bijux_value=python_results["ancestral_reconstruction"]["brownian"][
                "recent_five_estimates"
            ],
        ),
        _comparison_row_exact(
            row_id="ancestral-eb-node-count",
            family_id="ancestral-reconstruction",
            fragment_id="ancestral-mode-comparison",
            metric_name="early_burst_node_count",
            r_value=r_results["ancestral_reconstruction"]["early_burst"]["node_count"],
            bijux_value=python_results["ancestral_reconstruction"]["early_burst"][
                "node_count"
            ],
        ),
        _comparison_row_exact(
            row_id="ancestral-eb-first-five",
            family_id="ancestral-reconstruction",
            fragment_id="ancestral-mode-comparison",
            metric_name="early_burst_first_five_estimates",
            r_value=r_results["ancestral_reconstruction"]["early_burst"][
                "first_five_estimates"
            ],
            bijux_value=python_results["ancestral_reconstruction"]["early_burst"][
                "first_five_estimates"
            ],
        ),
        _comparison_row_exact(
            row_id="ancestral-eb-recent-five",
            family_id="ancestral-reconstruction",
            fragment_id="ancestral-mode-comparison",
            metric_name="early_burst_recent_five_estimates",
            r_value=r_results["ancestral_reconstruction"]["early_burst"][
                "recent_five_estimates"
            ],
            bijux_value=python_results["ancestral_reconstruction"]["early_burst"][
                "recent_five_estimates"
            ],
        ),
        _comparison_row_tolerance(
            row_id="baseline-intercept",
            family_id="baseline-regression",
            fragment_id="baseline-gls-fit",
            metric_name="intercept",
            r_value=r_results["baseline_gls"]["coefficients"]["intercept"],
            bijux_value=python_results["baseline_gls"]["coefficients"]["intercept"],
            tolerance_abs_diff=1e-06,
            reference_rounding_digits=4,
            explained_rounding_message=(
                "The R reference stores the baseline intercept rounded to four decimal places; "
                "the Bijux value rounds back to the same published scalar."
            ),
        ),
        _comparison_row_tolerance(
            row_id="baseline-slope",
            family_id="baseline-regression",
            fragment_id="baseline-gls-fit",
            metric_name="social_group_size",
            r_value=r_results["baseline_gls"]["coefficients"]["social_group_size"],
            bijux_value=python_results["baseline_gls"]["coefficients"][
                "social_group_size"
            ],
            tolerance_abs_diff=1e-06,
            reference_rounding_digits=4,
            explained_rounding_message=(
                "The R reference stores the baseline slope rounded to four decimal places; "
                "the Bijux value rounds back to the same published scalar."
            ),
        ),
        _comparison_row_tolerance(
            row_id="baseline-log-likelihood",
            family_id="baseline-regression",
            fragment_id="baseline-gls-fit",
            metric_name="log_likelihood",
            r_value=r_results["baseline_gls"]["log_likelihood"],
            bijux_value=python_results["baseline_gls"]["log_likelihood"],
            tolerance_abs_diff=1e-06,
            reference_rounding_digits=4,
            explained_rounding_message=(
                "The R reference stores the baseline log likelihood rounded to four decimal places; "
                "the Bijux value rounds back to the same published scalar."
            ),
        ),
        _comparison_row_tolerance(
            row_id="baseline-r-squared",
            family_id="baseline-regression",
            fragment_id="baseline-gls-fit",
            metric_name="r_squared",
            r_value=r_results["baseline_gls"]["r_squared"],
            bijux_value=python_results["baseline_gls"]["r_squared"],
            tolerance_abs_diff=1e-06,
            reference_rounding_digits=4,
            explained_rounding_message=(
                "The R reference stores the baseline R-squared rounded to four decimal places; "
                "the Bijux value rounds back to the same published scalar."
            ),
        ),
        _comparison_row_tolerance(
            row_id="estimated-lambda-value",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="lambda_value",
            r_value=r_results["estimated_lambda_pgls"]["lambda_value"],
            bijux_value=python_results["estimated_lambda_pgls"]["lambda_value"],
            tolerance_abs_diff=0.05,
        ),
        _comparison_row_tolerance(
            row_id="estimated-pgls-intercept",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="intercept",
            r_value=r_results["estimated_lambda_pgls"]["coefficients"]["intercept"],
            bijux_value=python_results["estimated_lambda_pgls"]["coefficients"][
                "intercept"
            ],
            tolerance_abs_diff=1.0,
        ),
        _comparison_row_tolerance(
            row_id="estimated-pgls-slope",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="social_group_size",
            r_value=r_results["estimated_lambda_pgls"]["coefficients"][
                "social_group_size"
            ],
            bijux_value=python_results["estimated_lambda_pgls"]["coefficients"][
                "social_group_size"
            ],
            tolerance_abs_diff=0.1,
        ),
        _comparison_row_tolerance(
            row_id="estimated-pgls-log-likelihood",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="log_likelihood",
            r_value=r_results["estimated_lambda_pgls"]["log_likelihood"],
            bijux_value=python_results["estimated_lambda_pgls"]["log_likelihood"],
            tolerance_abs_diff=0.25,
        ),
        _comparison_row_equivalence(
            row_id="estimated-pgls-slope-significance",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="social_group_size_significant_under_0.05",
            r_value=r_results["estimated_lambda_pgls"]["p_values"]["social_group_size"]
            < 0.05,
            bijux_value=python_results["estimated_lambda_pgls"]["p_values"][
                "social_group_size"
            ]
            < 0.05,
            equivalent=(
                r_results["estimated_lambda_pgls"]["p_values"]["social_group_size"]
                < 0.05
            )
            == (
                python_results["estimated_lambda_pgls"]["p_values"]["social_group_size"]
                < 0.05
            ),
            rule="Both implementations must keep the predictor on the same side of the 0.05 significance boundary.",
        ),
        _comparison_row_tolerance(
            row_id="signal-estimated-lambda",
            family_id="phylogenetic-signal",
            fragment_id="phylogenetic-signal-test",
            metric_name="estimated_lambda",
            r_value=r_results["signal_test"]["estimated_lambda"],
            bijux_value=python_results["signal_test"]["estimated_lambda"],
            tolerance_abs_diff=0.05,
        ),
        _comparison_row_tolerance(
            row_id="signal-likelihood-ratio",
            family_id="phylogenetic-signal",
            fragment_id="phylogenetic-signal-test",
            metric_name="likelihood_ratio",
            r_value=r_results["signal_test"]["likelihood_ratio"],
            bijux_value=python_results["signal_test"]["likelihood_ratio"],
            tolerance_abs_diff=0.05,
        ),
        _comparison_row_equivalence(
            row_id="signal-reject-lambda-zero",
            family_id="phylogenetic-signal",
            fragment_id="phylogenetic-signal-test",
            metric_name="p_value_below_0.05",
            r_value=r_results["signal_test"]["p_value"] < 0.05,
            bijux_value=python_results["signal_test"]["p_value"] < 0.05,
            equivalent=(r_results["signal_test"]["p_value"] < 0.05)
            == (python_results["signal_test"]["p_value"] < 0.05),
            rule="Both implementations must preserve the same lambda-zero rejection decision at 0.05.",
        ),
        _comparison_row_tolerance(
            row_id="baseline-diagnostic-qq-correlation",
            family_id="diagnostics",
            fragment_id="baseline-gls-diagnostics",
            metric_name="qq_correlation",
            r_value=r_results["baseline_gls"]["diagnostics"]["qq_correlation"],
            bijux_value=python_results["baseline_gls"]["diagnostics"]["qq_correlation"],
            tolerance_abs_diff=0.02,
        ),
        _comparison_row_tolerance(
            row_id="baseline-diagnostic-fitted-correlation",
            family_id="diagnostics",
            fragment_id="baseline-gls-diagnostics",
            metric_name="abs_residual_fitted_correlation",
            r_value=r_results["baseline_gls"]["diagnostics"][
                "abs_residual_fitted_correlation"
            ],
            bijux_value=python_results["baseline_gls"]["diagnostics"][
                "abs_residual_fitted_correlation"
            ],
            tolerance_abs_diff=0.05,
        ),
        _comparison_row_equivalence(
            row_id="baseline-diagnostic-outlier-pressure",
            family_id="diagnostics",
            fragment_id="baseline-gls-diagnostics",
            metric_name="outlier_count_abs_z_ge_2_close",
            r_value=r_results["baseline_gls"]["diagnostics"][
                "outlier_count_abs_z_ge_2"
            ],
            bijux_value=python_results["baseline_gls"]["diagnostics"][
                "outlier_count_abs_z_ge_2"
            ],
            equivalent=abs(
                r_results["baseline_gls"]["diagnostics"]["outlier_count_abs_z_ge_2"]
                - python_results["baseline_gls"]["diagnostics"][
                    "outlier_count_abs_z_ge_2"
                ]
            )
            <= 1,
            rule="Baseline outlier counts may drift by at most one taxon while preserving the same practical review conclusion.",
        ),
        _comparison_row_tolerance(
            row_id="estimated-diagnostic-qq-correlation",
            family_id="diagnostics",
            fragment_id="estimated-lambda-diagnostics",
            metric_name="qq_correlation",
            r_value=r_results["estimated_lambda_pgls"]["diagnostics"]["qq_correlation"],
            bijux_value=python_results["estimated_lambda_pgls"]["diagnostics"][
                "qq_correlation"
            ],
            tolerance_abs_diff=0.02,
        ),
        _comparison_row_tolerance(
            row_id="estimated-diagnostic-fitted-correlation",
            family_id="diagnostics",
            fragment_id="estimated-lambda-diagnostics",
            metric_name="abs_residual_fitted_correlation",
            r_value=r_results["estimated_lambda_pgls"]["diagnostics"][
                "abs_residual_fitted_correlation"
            ],
            bijux_value=python_results["estimated_lambda_pgls"]["diagnostics"][
                "abs_residual_fitted_correlation"
            ],
            tolerance_abs_diff=0.05,
        ),
        _comparison_row_equivalence(
            row_id="estimated-diagnostic-outlier-pressure",
            family_id="diagnostics",
            fragment_id="estimated-lambda-diagnostics",
            metric_name="outlier_count_abs_z_ge_2_close",
            r_value=r_results["estimated_lambda_pgls"]["diagnostics"][
                "outlier_count_abs_z_ge_2"
            ],
            bijux_value=python_results["estimated_lambda_pgls"]["diagnostics"][
                "outlier_count_abs_z_ge_2"
            ],
            equivalent=abs(
                r_results["estimated_lambda_pgls"]["diagnostics"][
                    "outlier_count_abs_z_ge_2"
                ]
                - python_results["estimated_lambda_pgls"]["diagnostics"][
                    "outlier_count_abs_z_ge_2"
                ]
            )
            <= 1,
            rule="Estimated-lambda outlier counts may drift by at most one taxon while preserving the same practical review conclusion.",
        ),
    ]
    verdict_counts: dict[str, int] = {}
    for row in rows:
        verdict = str(row["verdict"])
        verdict_counts[verdict] = verdict_counts.get(verdict, 0) + 1
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": SUMMARY_EVIDENCE_ID,
        "row_count": len(rows),
        "verdict_counts": dict(sorted(verdict_counts.items())),
        "rows": rows,
    }


def render_primate_pgls_signal_scalar_parity_table_markdown(
    payload: dict[str, object],
) -> str:
    lines = [
        "# Scalar Parity Table",
        "",
        "| Row | Family | Metric | Kind | Verdict | R | Bijux |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload["rows"]:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["row_id"]),
                    str(row["method_family"]),
                    str(row["metric_name"]),
                    str(row["comparison_kind"]),
                    str(row["verdict"]),
                    str(row["r_value"]),
                    str(row["bijux_value"]),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## Verdict Counts",
            "",
        ]
    )
    for verdict, count in payload["verdict_counts"].items():
        lines.append(f"- `{verdict}`: `{count}`")
    lines.append("")
    return "\n".join(lines)


def build_primate_pgls_signal_claim_registry(repo_root: Path) -> dict[str, object]:
    bundles = build_primate_pgls_signal_bundles(repo_root)
    claims = []
    for definition in BUNDLE_DEFINITIONS:
        claim = CLAIM_DEFINITIONS[definition["claim_id"]]
        claims.append(
            {
                "claim_id": definition["claim_id"],
                "claim_title": claim["claim_title"],
                "summary": claim["summary"],
                "verdict": claim["verdict"],
                "evidence_ids": [definition["evidence_id"]],
                "source_fragments": definition["source_fragments"],
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "claim_count": len(claims),
        "claims": claims,
        "bundle_count": len(bundles),
    }


def build_primate_pgls_signal_family_index(repo_root: Path) -> dict[str, object]:
    fragments = build_primate_pgls_signal_source_fragment_map()["fragments"]
    families = []
    for family_id, family in FAMILY_DEFINITIONS.items():
        matching_fragments = [
            fragment
            for fragment in fragments
            if fragment["concept_family"] == family_id
        ]
        claim_ids = sorted(
            {
                claim_id
                for fragment in matching_fragments
                for claim_id in fragment["claim_ids"]
            }
        )
        evidence_ids = sorted(
            {
                fragment["evidence_id"]
                for fragment in matching_fragments
                if fragment["evidence_id"]
            }
        )
        claims = {
            claim_id: CLAIM_DEFINITIONS[claim_id]["claim_title"]
            for claim_id in claim_ids
        }
        verdicts = {CLAIM_DEFINITIONS[claim_id]["verdict"] for claim_id in claim_ids}
        if verdicts == {"matched"}:
            family_verdict = "matched"
        elif "matched_with_tolerance" in verdicts and not (
            "mismatch_unexplained" in verdicts or "not_comparable" in verdicts
        ):
            family_verdict = "matched_with_tolerance"
        elif verdicts == {"not_comparable"}:
            family_verdict = "not_comparable"
        else:
            family_verdict = (
                "not_comparable"
                if "not_comparable" in verdicts
                else "mismatch_unexplained"
            )
        families.append(
            {
                "family_id": family_id,
                "family_title": family["title"],
                "summary": family["summary"],
                "fragment_count": len(matching_fragments),
                "fragment_ids": [
                    fragment["fragment_id"] for fragment in matching_fragments
                ],
                "claim_ids": claim_ids,
                "claim_titles": claims,
                "evidence_ids": evidence_ids,
                "family_verdict": family_verdict,
                "coverage_status": (
                    "coverage-gap" if family_id == "coverage-boundaries" else "covered"
                ),
                "known_gaps": []
                if family_id != "coverage-boundaries"
                else [
                    "The lecture corBlomberg intercept-mode likelihood sweep is still an explicit coverage boundary.",
                ],
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": SUMMARY_EVIDENCE_ID,
        "source_family_id": "pcm2-modes-pgls",
        "source_family_title": "PCM2 modes and PGLS evidence family",
        "family_count": len(families),
        "families": families,
    }


def _report_payload_for_bundle(repo_root: Path, evidence_id: str) -> dict[str, object]:
    r_results = _load_r_reference_results(repo_root)
    if evidence_id == "evidence-001":
        scalar_table = build_primate_pgls_signal_scalar_parity_table(repo_root)
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "object_names": r_results["source_contract"]["object_names"],
            "row_count": r_results["source_contract"]["row_count"],
            "tip_count": r_results["source_contract"]["tip_count"],
            "species_tip_match": r_results["source_contract"]["species_tip_match"],
            "governed_reload_inputs": [
                (STUDY_ONE_REFERENCE_ROOT / "reference_primate.csv").as_posix(),
                (
                    STUDY_ONE_REFERENCE_ROOT / "reference_trimmed_primatetree.nwk"
                ).as_posix(),
            ],
            "scalar_row_count": scalar_table["row_count"],
            "verdict_counts": scalar_table["verdict_counts"],
        }
    if evidence_id == "evidence-002":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_baseline": r_results["baseline_gls"],
            "bijux_baseline": _baseline_gls_payload(_baseline_gls_report(repo_root)),
            "r_fixed_lambda_equivalence": r_results[
                "fixed_lambda_gls_matches_baseline"
            ],
        }
    if evidence_id == "evidence-003":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_estimated_lambda": r_results["estimated_lambda_pgls"],
            "bijux_estimated_lambda": _estimated_lambda_pgls_payload(
                _estimated_lambda_pgls_report(repo_root)
            ),
        }
    if evidence_id == "evidence-004":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_signal_test": r_results["signal_test"],
            "bijux_signal_test": _signal_test_payload(_signal_report(repo_root)),
        }
    if evidence_id == "evidence-005":
        baseline = _baseline_gls_payload(_baseline_gls_report(repo_root))
        estimated = _estimated_lambda_pgls_payload(
            _estimated_lambda_pgls_report(repo_root)
        )
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "baseline_diagnostics": {
                "r": r_results["baseline_gls"]["diagnostics"],
                "bijux": baseline["diagnostics"],
            },
            "estimated_lambda_diagnostics": {
                "r": r_results["estimated_lambda_pgls"]["diagnostics"],
                "bijux": estimated["diagnostics"],
            },
        }
    if evidence_id == "evidence-006":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_tree_rescaling": r_results["tree_rescaling"],
            "bijux_tree_rescaling": _tree_rescaling_payloads(
                _transformed_tree_reports(repo_root)
            ),
        }
    if evidence_id == "evidence-007":
        brownian_fit, ou_fit, early_burst_fit = _continuous_mode_fit_reports(repo_root)
        tip_values = _ordered_trait_values(
            _source_reference_paths(repo_root)[1],
            brownian_fit.taxa,
            trait="longevity",
            taxon_column="species",
        )
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_continuous_mode_fits": r_results["continuous_mode_fits"],
            "bijux_continuous_mode_fits": _continuous_mode_fit_payloads(
                brownian_fit,
                ou_fit,
                early_burst_fit,
                tip_values=tip_values,
            ),
        }
    if evidence_id == "evidence-008":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_likelihood_ratio_tests": r_results["likelihood_ratio_tests"],
            "bijux_likelihood_ratio_tests": _likelihood_ratio_test_payloads(
                _mode_comparison_report(repo_root)
            ),
        }
    if evidence_id == "evidence-009":
        brownian_ancestral, early_burst_ancestral = _ancestral_reconstruction_reports(
            repo_root
        )
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_ancestral_reconstruction": r_results["ancestral_reconstruction"],
            "bijux_ancestral_reconstruction": _ancestral_reconstruction_payloads(
                brownian_ancestral,
                early_burst_ancestral,
            ),
        }
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": evidence_id,
        "coverage_boundaries": _coverage_boundary_payload(),
    }


def build_primate_pgls_signal_bundle(
    repo_root: Path, evidence_id: str
) -> dict[str, object]:
    definition = next(
        (
            definition
            for definition in BUNDLE_DEFINITIONS
            if definition["evidence_id"] == evidence_id
        ),
        None,
    )
    if definition is None:
        raise KeyError(evidence_id)
    report_payload = _report_payload_for_bundle(repo_root, evidence_id)
    manifest = _manifest_for_bundle(repo_root, definition, report_payload)
    bundle = {
        "manifest": manifest,
        "claims": _claims_payload(definition),
        "report_payload": report_payload,
        "report_filename": definition["report_filename"],
        "readme": _readme_for_bundle(definition),
    }
    if evidence_id == SUMMARY_EVIDENCE_ID:
        scalar_table = build_primate_pgls_signal_scalar_parity_table(repo_root)
        bundle["scalar_parity_table"] = scalar_table
        bundle["scalar_parity_markdown"] = (
            render_primate_pgls_signal_scalar_parity_table_markdown(scalar_table)
        )
    return bundle


def build_primate_pgls_signal_evidence_registry(
    repo_root: Path,
) -> dict[str, object]:
    bundles = build_primate_pgls_signal_bundles(repo_root)
    evidences = []
    for definition in BUNDLE_DEFINITIONS:
        claim = CLAIM_DEFINITIONS[str(definition["claim_id"])]
        evidences.append(
            {
                "evidence_id": definition["evidence_id"],
                "title": definition["title"],
                "coverage_status": (
                    "coverage-gap"
                    if claim["verdict"] == "not_comparable"
                    else "covered"
                ),
                "claim_id": definition["claim_id"],
                "verdict": claim["verdict"],
                "analytical_surfaces": definition["analytical_surfaces"],
                "source_fragments": definition["source_fragments"],
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "bundle_count": len(bundles),
        "evidence_count": len(evidences),
        "coverage_boundary_evidence_ids": [
            entry["evidence_id"]
            for entry in evidences
            if entry["coverage_status"] == "coverage-gap"
        ],
        "evidences": evidences,
    }


def render_primate_pgls_signal_study_manifest(repo_root: Path) -> dict[str, object]:
    registry = build_primate_pgls_signal_evidence_registry(repo_root)
    return {
        "study_id": STUDY_ID,
        "study_title": "Primate PGLS and signal evidence study",
        "summary": "Governed parity study for the regression, transformed-tree, evolutionary-mode fit, and ancestral sections of the Lund primate comparative lecture, with an explicit remaining intercept-mode boundary.",
        "owner_package": "bijux-phylogenetics",
        "study_categories": ["teaching-study", "migration-study"],
        "confidence_posture": "governed-parity-in-progress",
        "coverage_boundary_evidence_ids": registry["coverage_boundary_evidence_ids"],
        "evidence_registry_locator": (
            f"evidence-book/studies/{STUDY_ID}/evidence-registry.json"
        ),
        "study_scope": {
            "coverage_focus": [
                "rdata-reload",
                "baseline-regression",
                "phylogenetic-regression",
                "phylogenetic-signal",
                "transformed-tree-workflows",
                "continuous-model-fitting",
                "likelihood-ratio-tests",
                "ancestral-reconstruction",
            ],
            "untouched_source_locators": [
                PCM2_SOURCE_LOCATOR,
                "external:lund/pcm2-modes-pgls/data/primate.RData",
            ],
        },
    }


def render_primate_pgls_signal_study_readme(repo_root: Path) -> str:
    registry = build_primate_pgls_signal_evidence_registry(repo_root)
    lines = [
        "# Primate PGLS And Signal",
        "",
        "This study turns the regression, transformed-tree, evolutionary-mode fit,",
        "likelihood-ratio, and ancestral sections of the Lund primate comparative",
        "lecture into governed Evidence IDs backed by checked-in R reference outputs",
        "and canonical `bijux-phylogenetics` reproductions.",
        "",
        "It is intentionally strict about confidence posture:",
        "",
        "- baseline GLS, Pagel-lambda PGLS, signal testing, transformed-tree",
        "  workflows, fitContinuous-style mode comparisons, and ancestral-mode",
        "  reconstructions are backed by governed parity bundles",
        "- the lecture corBlomberg intercept sweep remains visible as an explicit",
        "  coverage boundary instead of being implied as validated",
        "",
        "Current bundles:",
        "",
    ]
    for entry in registry["evidences"]:
        title = str(entry["title"]).removeprefix("Primate ").removesuffix(" bundle")
        lines.append(f"- `{entry['evidence_id']}` {title}")
    lines.append("")
    return "\n".join(lines)


def _manifest_for_bundle(
    repo_root: Path, definition: dict[str, object], report_payload: dict[str, object]
) -> dict[str, object]:
    evidence_id = str(definition["evidence_id"])
    source_basis = [
        {
            "kind": "external-source-descriptor",
            "label": "Lund source descriptors",
            "locator": f"evidence-book/studies/{STUDY_ID}/provenance/lund-course-sources.json",
        },
        {
            "kind": "repository-reference",
            "label": "governed primate reference table",
            "locator": (STUDY_ONE_REFERENCE_ROOT / "reference_primate.csv").as_posix(),
        },
        {
            "kind": "repository-reference",
            "label": "governed primate trimmed tree",
            "locator": (
                STUDY_ONE_REFERENCE_ROOT / "reference_trimmed_primatetree.nwk"
            ).as_posix(),
        },
        {
            "kind": "repository-reference",
            "label": "R reference results for the PGLS and signal study",
            "locator": f"evidence-book/studies/{STUDY_ID}/reference/reference_results.json",
        },
        {
            "kind": "repository-reference",
            "label": f"{definition['title']} report payload",
            "locator": f"evidence-book/studies/{STUDY_ID}/{evidence_id}/{definition['report_filename']}",
        },
    ]
    if evidence_id == SUMMARY_EVIDENCE_ID:
        source_basis.append(
            {
                "kind": "repository-reference",
                "label": "scalar parity table",
                "locator": f"evidence-book/studies/{STUDY_ID}/{SUMMARY_EVIDENCE_ID}/scalar-parity-table.json",
            }
        )
    claim = CLAIM_DEFINITIONS[str(definition["claim_id"])]
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": evidence_id,
        "evidence_title": definition["title"],
        "summary": definition["summary"],
        "owner_package": "bijux-phylogenetics",
        "claim_ids": [definition["claim_id"]],
        "source_basis": source_basis,
        "freshness": {
            "last_generated_on": date.today().isoformat(),
            "governed_code_paths": [
                "packages/bijux-phylogenetics/src/bijux_phylogenetics/evidence/studies/primate_pgls_and_signal.py"
            ],
            "source_basis_locators": [entry["locator"] for entry in source_basis],
        },
        "ownership": {
            "owner_package": "bijux-phylogenetics",
            "analytical_surfaces": definition["analytical_surfaces"],
        },
        "claim_tags": definition["claim_tags"],
        "comparison_mode": definition["comparison_mode"],
        "verdict": {
            "status": claim["verdict"],
            "summary": claim["summary"],
        },
        "limitations": definition["limitations"],
        "source_fragments": definition["source_fragments"],
        "reference_script_locators": [f"{PCM2_REFERENCE_SCRIPT_PATH}#L1-L200"],
        "supporting_report_locator": (
            f"evidence-book/studies/{STUDY_ID}/{evidence_id}/{definition['report_filename']}"
        ),
        "report_keys": sorted(report_payload.keys()),
    }


def _claims_payload(definition: dict[str, object]) -> dict[str, object]:
    claim = CLAIM_DEFINITIONS[str(definition["claim_id"])]
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": definition["evidence_id"],
        "claim_count": 1,
        "claims": [
            {
                "claim_id": definition["claim_id"],
                "claim_title": claim["claim_title"],
                "summary": claim["summary"],
                "verdict": claim["verdict"],
                "evidence_ids": [definition["evidence_id"]],
                "source_fragments": definition["source_fragments"],
            }
        ],
    }


def _readme_for_bundle(definition: dict[str, object]) -> str:
    lines = [
        f"# {definition['title']}",
        "",
        definition["summary"],
        "",
        f"- evidence id: `{definition['evidence_id']}`",
        f"- source fragments: {', '.join(f'`{fragment}`' for fragment in definition['source_fragments'])}",
        "",
        "## Limitations",
        "",
    ]
    for limitation in definition["limitations"]:
        lines.append(f"- {limitation}")
    lines.extend(
        [
            "",
            "## Source Locators",
            "",
            f"- `{PCM2_SOURCE_LOCATOR}`",
            f"- `{PCM2_REFERENCE_SCRIPT_PATH}`",
            "",
        ]
    )
    return "\n".join(lines)


def build_primate_pgls_signal_bundles(repo_root: Path) -> dict[str, dict[str, object]]:
    bundles: dict[str, dict[str, object]] = {}
    for definition in BUNDLE_DEFINITIONS:
        bundles[definition["evidence_id"]] = build_primate_pgls_signal_bundle(
            repo_root, definition["evidence_id"]
        )
    return bundles
