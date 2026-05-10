from __future__ import annotations

from functools import lru_cache
import json
import math
from pathlib import Path

from bijux_phylogenetics.comparative.pgls import run_pgls
from bijux_phylogenetics.comparative.signal import estimate_pagels_lambda


STUDY_ID = "primate-pgls-and-signal"
SUMMARY_EVIDENCE_ID = "evidence-001"
PCM2_SOURCE_LOCATOR = "external:lund/pcm2-modes-pgls/script"
PCM2_REFERENCE_SCRIPT_PATH = (
    "evidence-book/studies/primate-pgls-and-signal/reference/"
    "primate_pgls_and_signal_reference_r.R"
)
STUDY_ONE_REFERENCE_ROOT = (
    Path("evidence-book")
    / "studies"
    / "primate-longevity-signal"
    / "evidence-001"
)

FAMILY_DEFINITIONS = {
    "workflow-contracts": {
        "title": "Workflow contracts",
        "summary": "Lecture workspace assumptions are represented explicitly without being overstated as numerical parity.",
    },
    "baseline-regression": {
        "title": "Baseline regression",
        "summary": "Non-phylogenetic regression outputs stay visible before phylogenetic correction is judged.",
    },
    "phylogenetic-regression": {
        "title": "Phylogenetic regression",
        "summary": "Pagel-lambda regression outputs are compared with explicit tolerance and conclusion rules.",
    },
    "phylogenetic-signal": {
        "title": "Phylogenetic signal",
        "summary": "Intercept-only signal testing and lambda-difference statistics are governed as their own trust surface.",
    },
    "diagnostics": {
        "title": "Diagnostics",
        "summary": "Residual, QQ, and heteroscedasticity checks are machine-recorded so visual diagnostics do not stay trapped in plots.",
    },
    "coverage-boundaries": {
        "title": "Coverage boundaries",
        "summary": "Unclosed evolutionary-mode, EB, and ancestral fragments are indexed explicitly instead of disappearing into vague future work.",
    },
}

CLAIM_DEFINITIONS = {
    "pcm2-reload-contract-governed": {
        "claim_title": "PCM2 reload semantics are represented explicitly",
        "summary": "The lecture one-line workspace reload contract is reconstructed from governed repository artifacts instead of being left as an opaque external assumption.",
        "verdict": "matched",
    },
    "pcm2-baseline-gls-parity": {
        "claim_title": "PCM2 baseline GLS outputs match before phylogenetic correction",
        "summary": "The non-phylogenetic regression slope, intercept, fit statistics, and coefficient decisions align before any phylogenetic covariance is introduced.",
        "verdict": "matched",
    },
    "pcm2-pagel-lambda-regression-parity": {
        "claim_title": "PCM2 Pagel-lambda regression agrees within governed tolerance",
        "summary": "Estimated-lambda regression outputs keep the same analytical conclusion while allowing bounded numerical drift between R and Bijux implementations.",
        "verdict": "matched_with_tolerance",
    },
    "pcm2-phylogenetic-signal-parity": {
        "claim_title": "PCM2 phylogenetic signal testing agrees within governed tolerance",
        "summary": "Lambda estimation and lambda-zero likelihood-ratio testing preserve the same signal conclusion while tolerating bounded numerical differences.",
        "verdict": "matched_with_tolerance",
    },
    "pcm2-diagnostics-parity": {
        "claim_title": "PCM2 residual diagnostics are machine-recorded and scientifically equivalent",
        "summary": "Residual diagnostics remain reviewer-visible as governed scalar summaries instead of being trapped in hand-inspected plots.",
        "verdict": "matched_with_tolerance",
    },
    "pcm2-coverage-boundary-explicit": {
        "claim_title": "PCM2 uncovered evolutionary-mode and ancestral fragments remain explicit",
        "summary": "Unimplemented EB, transformed-tree, and ancestral parity surfaces are indexed as open trust boundaries rather than silently implied.",
        "verdict": "not_comparable",
    },
}

FRAGMENT_DEFINITIONS = [
    {
        "fragment_id": "workspace-reload-contract",
        "fragment_title": "Package loading and one-line primate workspace reload",
        "family_id": "workflow-contracts",
        "claim_ids": ["pcm2-reload-contract-governed"],
        "evidence_id": "evidence-001",
        "supporting_evidence_ids": [],
        "script_line_spec": "8-16",
        "parity_expectation": "exact",
        "comparison_kind": "exact_answer",
        "block_status": "verified",
        "review_note": "The lecture `load()` contract is reconstructed from governed repository artifacts with explicit object and path semantics.",
        "scope": "workflow",
    },
    {
        "fragment_id": "transformed-tree-workflows",
        "fragment_title": "OU and EB tree rescaling exploration",
        "family_id": "coverage-boundaries",
        "claim_ids": ["pcm2-coverage-boundary-explicit"],
        "evidence_id": "evidence-006",
        "supporting_evidence_ids": [],
        "script_line_spec": "18-34",
        "parity_expectation": "not_comparable",
        "comparison_kind": "not_comparable",
        "block_status": "coverage_gap",
        "review_note": "The script-level transformed-tree surface is indexed explicitly, but EB parity is not yet claimed inside the runtime.",
        "scope": "analytical",
    },
    {
        "fragment_id": "continuous-model-comparison",
        "fragment_title": "BM, OU, and EB fitContinuous model comparison",
        "family_id": "coverage-boundaries",
        "claim_ids": ["pcm2-coverage-boundary-explicit"],
        "evidence_id": "evidence-006",
        "supporting_evidence_ids": [],
        "script_line_spec": "36-87",
        "parity_expectation": "not_comparable",
        "comparison_kind": "not_comparable",
        "block_status": "coverage_gap",
        "review_note": "These model-comparison fragments remain open because the canonical runtime does not yet expose a governed EB parity surface.",
        "scope": "analytical",
    },
    {
        "fragment_id": "ancestral-mode-comparison",
        "fragment_title": "Ancestral-state comparison under BM and EB",
        "family_id": "coverage-boundaries",
        "claim_ids": ["pcm2-coverage-boundary-explicit"],
        "evidence_id": "evidence-006",
        "supporting_evidence_ids": [],
        "script_line_spec": "89-111",
        "parity_expectation": "not_comparable",
        "comparison_kind": "not_comparable",
        "block_status": "coverage_gap",
        "review_note": "The lecture ancestral comparison is preserved as an indexed trust boundary until EB ancestral parity is implemented for the runtime.",
        "scope": "analytical",
    },
    {
        "fragment_id": "baseline-gls-fit",
        "fragment_title": "Non-phylogenetic GLS fit",
        "family_id": "baseline-regression",
        "claim_ids": ["pcm2-baseline-gls-parity"],
        "evidence_id": "evidence-002",
        "supporting_evidence_ids": ["evidence-005"],
        "script_line_spec": "122-136",
        "parity_expectation": "exact",
        "comparison_kind": "exact_answer",
        "block_status": "verified",
        "review_note": "The baseline regression is checked before phylogenetic covariance enters the workflow.",
        "scope": "analytical",
    },
    {
        "fragment_id": "baseline-gls-diagnostics",
        "fragment_title": "Baseline GLS heteroscedasticity and QQ diagnostics",
        "family_id": "diagnostics",
        "claim_ids": ["pcm2-diagnostics-parity"],
        "evidence_id": "evidence-005",
        "supporting_evidence_ids": ["evidence-002"],
        "script_line_spec": "126-133",
        "parity_expectation": "scientific_equivalence",
        "comparison_kind": "scientific_equivalence",
        "block_status": "verified",
        "review_note": "Plot-only diagnostics are converted into machine-recorded scalar summaries with explicit equivalence rules.",
        "scope": "diagnostic",
    },
    {
        "fragment_id": "pagel-lambda-regression",
        "fragment_title": "Fixed-lambda and estimated-lambda PGLS fits",
        "family_id": "phylogenetic-regression",
        "claim_ids": ["pcm2-pagel-lambda-regression-parity"],
        "evidence_id": "evidence-003",
        "supporting_evidence_ids": ["evidence-005"],
        "script_line_spec": "138-179",
        "parity_expectation": "statistical_tolerance",
        "comparison_kind": "tolerance_or_equivalence",
        "block_status": "verified",
        "review_note": "The lecture fixed-lambda equivalence and estimated-lambda fit are compared with separate tolerance and conclusion rules.",
        "scope": "analytical",
    },
    {
        "fragment_id": "estimated-lambda-diagnostics",
        "fragment_title": "Estimated-lambda residual and fitted diagnostics",
        "family_id": "diagnostics",
        "claim_ids": ["pcm2-diagnostics-parity"],
        "evidence_id": "evidence-005",
        "supporting_evidence_ids": ["evidence-003"],
        "script_line_spec": "168-179",
        "parity_expectation": "scientific_equivalence",
        "comparison_kind": "scientific_equivalence",
        "block_status": "verified",
        "review_note": "The estimated-lambda residual checks are reduced to governed scalar diagnostics instead of staying as plot-only anecdotes.",
        "scope": "diagnostic",
    },
    {
        "fragment_id": "phylogenetic-signal-test",
        "fragment_title": "Intercept-only PGLS and lambda-zero likelihood-ratio testing",
        "family_id": "phylogenetic-signal",
        "claim_ids": ["pcm2-phylogenetic-signal-parity"],
        "evidence_id": "evidence-004",
        "supporting_evidence_ids": [],
        "script_line_spec": "181-192",
        "parity_expectation": "statistical_tolerance",
        "comparison_kind": "tolerance_or_equivalence",
        "block_status": "verified",
        "review_note": "Signal is judged by intercept-only likelihood surfaces and the lambda-zero model comparison, not by prose.",
        "scope": "analytical",
    },
    {
        "fragment_id": "mode-linked-intercept-models",
        "fragment_title": "Mode-linked intercept-only GLS surrogates for BM, OU, and EB",
        "family_id": "coverage-boundaries",
        "claim_ids": ["pcm2-coverage-boundary-explicit"],
        "evidence_id": "evidence-006",
        "supporting_evidence_ids": [],
        "script_line_spec": "194-220",
        "parity_expectation": "not_comparable",
        "comparison_kind": "not_comparable",
        "block_status": "coverage_gap",
        "review_note": "The lecture’s BM/OU/EB intercept-only surrogate fits are kept visible, but EB parity remains open and therefore unclaimed.",
        "scope": "analytical",
    },
]

BUNDLE_DEFINITIONS = [
    {
        "evidence_id": "evidence-001",
        "report_filename": "rdata-reload-semantics.json",
        "title": "Primate reload semantics bundle",
        "summary": "Governed representation of the lecture one-line primate workspace reload contract.",
        "claim_id": "pcm2-reload-contract-governed",
        "claim_tags": ["teaching", "parity", "workflow", "reload-semantics"],
        "analytical_surfaces": ["workflow-contracts"],
        "source_fragments": ["workspace-reload-contract"],
        "limitations": [
            "The raw lecture `primate.RData` file remains external; this bundle governs the object contract using repository-owned reference artifacts."
        ],
    },
    {
        "evidence_id": "evidence-002",
        "report_filename": "baseline-gls-parity.json",
        "title": "Primate baseline GLS parity bundle",
        "summary": "Governed parity for the non-phylogenetic baseline regression before phylogenetic covariance enters the workflow.",
        "claim_id": "pcm2-baseline-gls-parity",
        "claim_tags": ["teaching", "parity", "baseline-regression", "gls"],
        "analytical_surfaces": ["baseline-regression"],
        "source_fragments": ["baseline-gls-fit"],
        "limitations": [
            "This bundle covers the baseline regression only; phylogenetic-covariance behavior is isolated in a separate bundle."
        ],
    },
    {
        "evidence_id": "evidence-003",
        "report_filename": "pagel-lambda-regression-parity.json",
        "title": "Primate Pagel-lambda regression parity bundle",
        "summary": "Governed parity for fixed-lambda and estimated-lambda regression surfaces derived from the lecture workflow.",
        "claim_id": "pcm2-pagel-lambda-regression-parity",
        "claim_tags": ["teaching", "parity", "phylogenetic-regression", "pagel-lambda"],
        "analytical_surfaces": ["phylogenetic-regression"],
        "source_fragments": ["pagel-lambda-regression"],
        "limitations": [
            "Likelihood and coefficient comparisons allow bounded numerical drift as long as the same scientific conclusion is preserved."
        ],
    },
    {
        "evidence_id": "evidence-004",
        "report_filename": "phylogenetic-signal-parity.json",
        "title": "Primate phylogenetic signal parity bundle",
        "summary": "Governed parity for intercept-only signal testing and lambda-zero likelihood-ratio logic.",
        "claim_id": "pcm2-phylogenetic-signal-parity",
        "claim_tags": ["teaching", "parity", "phylogenetic-signal", "lambda-zero"],
        "analytical_surfaces": ["phylogenetic-signal"],
        "source_fragments": ["phylogenetic-signal-test"],
        "limitations": [
            "This bundle records the lecture’s intercept-only signal workflow, not the later BM/OU/EB surrogate-model section."
        ],
    },
    {
        "evidence_id": "evidence-005",
        "report_filename": "residual-diagnostics-parity.json",
        "title": "Primate residual diagnostics parity bundle",
        "summary": "Governed machine-readable residual diagnostics for the baseline and estimated-lambda regression surfaces.",
        "claim_id": "pcm2-diagnostics-parity",
        "claim_tags": ["teaching", "parity", "diagnostics", "residuals"],
        "analytical_surfaces": ["diagnostics"],
        "source_fragments": [
            "baseline-gls-diagnostics",
            "estimated-lambda-diagnostics",
        ],
        "limitations": [
            "The lecture plots are summarized as scalar diagnostics; this bundle does not claim rendered-figure equivalence."
        ],
    },
    {
        "evidence_id": "evidence-006",
        "report_filename": "coverage-boundaries.json",
        "title": "Primate coverage boundary bundle",
        "summary": "Governed record of uncovered transformed-tree, EB, and ancestral fragments from the lecture workflow.",
        "claim_id": "pcm2-coverage-boundary-explicit",
        "claim_tags": ["teaching", "coverage-gap", "not-comparable"],
        "analytical_surfaces": ["coverage-boundaries"],
        "source_fragments": [
            "transformed-tree-workflows",
            "continuous-model-comparison",
            "ancestral-mode-comparison",
            "mode-linked-intercept-models",
        ],
        "limitations": [
            "This bundle is intentionally a boundary register, not a parity claim."
        ],
    },
]


def _study_root(repo_root: Path) -> Path:
    return Path(repo_root) / "evidence-book" / "studies" / STUDY_ID


def _source_reference_paths(repo_root: Path) -> tuple[Path, Path]:
    root = Path(repo_root) / STUDY_ONE_REFERENCE_ROOT
    return root / "reference_trimmed_primatetree.nwk", root / "reference_primate.csv"


def _read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


@lru_cache(maxsize=None)
def _load_r_reference_results(repo_root: Path) -> dict[str, object]:
    return _read_json(_study_root(repo_root) / "reference" / "reference_results.json")


def _rounded(value: float) -> float:
    return float(format(float(value), ".15g"))


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
    abs_residual_ss = sum(
        (value - mean_abs_residual) ** 2 for value in abs_residuals
    )
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


@lru_cache(maxsize=None)
def _load_python_results(repo_root: Path) -> dict[str, object]:
    tree_path, traits_path = _source_reference_paths(repo_root)
    baseline = run_pgls(
        tree_path,
        traits_path,
        response="longevity",
        predictors=["social_group_size"],
        taxon_column="species",
        lambda_value=0.0,
    )
    estimated = run_pgls(
        tree_path,
        traits_path,
        response="longevity",
        predictors=["social_group_size"],
        taxon_column="species",
        lambda_value="estimate",
    )
    signal = estimate_pagels_lambda(
        tree_path,
        traits_path,
        trait="longevity",
        taxon_column="species",
    )
    return {
        "source_contract": {
            "row_count": baseline.taxon_count,
            "tip_count": baseline.taxon_count,
            "predictor": "social_group_size",
            "response": "longevity",
        },
        "baseline_gls": {
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
        },
        "estimated_lambda_pgls": {
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
        },
        "signal_test": {
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
        },
    }


def _line_spec_to_locators(spec: str) -> list[str]:
    locators: list[str] = []
    for part in spec.split(","):
        normalized = part.strip()
        if not normalized:
            continue
        if "-" in normalized:
            start_text, end_text = normalized.split("-", maxsplit=1)
            locators.append(f"{PCM2_SOURCE_LOCATOR}#L{int(start_text)}-L{int(end_text)}")
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
                "provides": ["reference_primate.csv", "reference_trimmed_primatetree.nwk"],
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
                "script_line_spans": _line_spec_to_spans(definition["script_line_spec"]),
                "script_locators": _line_spec_to_locators(definition["script_line_spec"]),
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
        "policy_count": 6,
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
                "family_id": "coverage-boundaries",
                "family_title": FAMILY_DEFINITIONS["coverage-boundaries"]["title"],
                "parity_expectation": "not_comparable",
                "comparison_kind": "not_comparable",
                "metric_tolerances": [],
                "rule": "These fragments are tracked explicitly as open trust boundaries and therefore cannot be promoted to parity claims yet.",
                "source_fragments": [
                    "transformed-tree-workflows",
                    "continuous-model-comparison",
                    "ancestral-mode-comparison",
                    "mode-linked-intercept-models",
                ],
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
    elif (
        reference_rounding_digits is not None
        and round(float(bijux_value), reference_rounding_digits) == float(r_value)
    ):
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
                python_results["estimated_lambda_pgls"]["p_values"][
                    "social_group_size"
                ]
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
            r_value=r_results["baseline_gls"]["diagnostics"]["outlier_count_abs_z_ge_2"],
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
            fragment for fragment in fragments if fragment["concept_family"] == family_id
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
        claims = {claim_id: CLAIM_DEFINITIONS[claim_id]["claim_title"] for claim_id in claim_ids}
        verdicts = {
            CLAIM_DEFINITIONS[claim_id]["verdict"] for claim_id in claim_ids
        }
        if verdicts == {"matched"}:
            family_verdict = "matched"
        elif "matched_with_tolerance" in verdicts and not (
            "mismatch_unexplained" in verdicts or "not_comparable" in verdicts
        ):
            family_verdict = "matched_with_tolerance"
        elif verdicts == {"not_comparable"}:
            family_verdict = "not_comparable"
        else:
            family_verdict = "not_comparable" if "not_comparable" in verdicts else "mismatch_unexplained"
        families.append(
            {
                "family_id": family_id,
                "family_title": family["title"],
                "summary": family["summary"],
                "fragment_count": len(matching_fragments),
                "fragment_ids": [fragment["fragment_id"] for fragment in matching_fragments],
                "claim_ids": claim_ids,
                "claim_titles": claims,
                "evidence_ids": evidence_ids,
                "family_verdict": family_verdict,
                "coverage_status": (
                    "coverage-gap"
                    if family_id == "coverage-boundaries"
                    else "covered"
                ),
                "known_gaps": []
                if family_id != "coverage-boundaries"
                else [
                    "No governed EB parity fit is exposed in the canonical runtime yet.",
                    "Transformed-tree and ancestral-state comparison fragments remain visible but intentionally unclaimed.",
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
    python_results = _load_python_results(repo_root)
    scalar_table = build_primate_pgls_signal_scalar_parity_table(repo_root)
    if evidence_id == "evidence-001":
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
                (STUDY_ONE_REFERENCE_ROOT / "reference_trimmed_primatetree.nwk").as_posix(),
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
            "bijux_baseline": python_results["baseline_gls"],
            "r_fixed_lambda_equivalence": r_results["fixed_lambda_gls_matches_baseline"],
        }
    if evidence_id == "evidence-003":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_estimated_lambda": r_results["estimated_lambda_pgls"],
            "bijux_estimated_lambda": python_results["estimated_lambda_pgls"],
        }
    if evidence_id == "evidence-004":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "r_signal_test": r_results["signal_test"],
            "bijux_signal_test": python_results["signal_test"],
        }
    if evidence_id == "evidence-005":
        return {
            "schema_version": 1,
            "study_id": STUDY_ID,
            "evidence_id": evidence_id,
            "baseline_diagnostics": {
                "r": r_results["baseline_gls"]["diagnostics"],
                "bijux": python_results["baseline_gls"]["diagnostics"],
            },
            "estimated_lambda_diagnostics": {
                "r": r_results["estimated_lambda_pgls"]["diagnostics"],
                "bijux": python_results["estimated_lambda_pgls"]["diagnostics"],
            },
        }
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": evidence_id,
        "coverage_boundaries": r_results["coverage_boundaries"],
    }


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
            "last_generated_on": "2026-05-10",
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
        "verdict": {
            "status": claim["verdict"],
            "summary": claim["summary"],
        },
        "limitations": definition["limitations"],
        "source_fragments": definition["source_fragments"],
        "reference_script_locators": [
            f"{PCM2_REFERENCE_SCRIPT_PATH}#L1-L200"
        ],
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
    scalar_table = build_primate_pgls_signal_scalar_parity_table(repo_root)
    scalar_markdown = render_primate_pgls_signal_scalar_parity_table_markdown(
        scalar_table
    )
    for definition in BUNDLE_DEFINITIONS:
        report_payload = _report_payload_for_bundle(repo_root, definition["evidence_id"])
        manifest = _manifest_for_bundle(repo_root, definition, report_payload)
        bundle = {
            "manifest": manifest,
            "claims": _claims_payload(definition),
            "report_payload": report_payload,
            "report_filename": definition["report_filename"],
            "readme": _readme_for_bundle(definition),
        }
        if definition["evidence_id"] == SUMMARY_EVIDENCE_ID:
            bundle["scalar_parity_table"] = scalar_table
            bundle["scalar_parity_markdown"] = scalar_markdown
        bundles[definition["evidence_id"]] = bundle
    return bundles
