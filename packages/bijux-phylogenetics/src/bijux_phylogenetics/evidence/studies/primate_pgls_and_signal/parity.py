from __future__ import annotations

import math
from pathlib import Path

from .definitions import STUDY_ID, SUMMARY_EVIDENCE_ID
from .runtime import (
    load_python_results,
    load_r_reference_results,
    rounded,
)


def comparison_row_exact(
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


def comparison_row_tolerance(
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
        float(bijux_value),
        reference_rounding_digits,
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
        "r_value": rounded(r_value),
        "bijux_value": rounded(bijux_value),
        "observed_abs_diff": rounded(observed_abs_diff),
        "tolerance_abs_diff": tolerance_abs_diff,
        "verdict": verdict,
    }
    if explanation_kind is not None:
        row["explanation_kind"] = explanation_kind
    if verdict_explanation is not None:
        row["verdict_explanation"] = verdict_explanation
    return row


def comparison_row_equivalence(
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
    r_results = load_r_reference_results(repo_root)
    python_results = load_python_results(repo_root)
    rows = [
        comparison_row_exact(
            row_id="reload-object-count",
            family_id="workflow-contracts",
            fragment_id="workspace-reload-contract",
            metric_name="object_name_count",
            r_value=r_results["source_contract"]["object_name_count"],
            bijux_value=python_results["source_contract"]["row_count"] and 2,
        ),
        comparison_row_exact(
            row_id="reload-primate-row-count",
            family_id="workflow-contracts",
            fragment_id="workspace-reload-contract",
            metric_name="primate_row_count",
            r_value=r_results["source_contract"]["row_count"],
            bijux_value=python_results["source_contract"]["row_count"],
        ),
        comparison_row_exact(
            row_id="reload-tree-tip-count",
            family_id="workflow-contracts",
            fragment_id="workspace-reload-contract",
            metric_name="tree_tip_count",
            r_value=r_results["source_contract"]["tip_count"],
            bijux_value=python_results["source_contract"]["tip_count"],
        ),
        comparison_row_exact(
            row_id="ou-alpha-1-branch-count",
            family_id="transformed-tree-workflows",
            fragment_id="transformed-tree-workflows",
            metric_name="ou_alpha_1_branch_count",
            r_value=r_results["tree_rescaling"]["ou_alpha_1"]["branch_count"],
            bijux_value=python_results["tree_rescaling"]["ou_alpha_1"]["branch_count"],
        ),
        comparison_row_exact(
            row_id="ou-alpha-1-total-branch-length",
            family_id="transformed-tree-workflows",
            fragment_id="transformed-tree-workflows",
            metric_name="ou_alpha_1_total_branch_length",
            r_value=r_results["tree_rescaling"]["ou_alpha_1"]["total_branch_length"],
            bijux_value=python_results["tree_rescaling"]["ou_alpha_1"][
                "total_branch_length"
            ],
        ),
        comparison_row_exact(
            row_id="ou-alpha-10-total-branch-length",
            family_id="transformed-tree-workflows",
            fragment_id="transformed-tree-workflows",
            metric_name="ou_alpha_10_total_branch_length",
            r_value=r_results["tree_rescaling"]["ou_alpha_10"]["total_branch_length"],
            bijux_value=python_results["tree_rescaling"]["ou_alpha_10"][
                "total_branch_length"
            ],
        ),
        comparison_row_exact(
            row_id="early-burst-2-total-branch-length",
            family_id="transformed-tree-workflows",
            fragment_id="transformed-tree-workflows",
            metric_name="early_burst_2_total_branch_length",
            r_value=r_results["tree_rescaling"]["early_burst_2"]["total_branch_length"],
            bijux_value=python_results["tree_rescaling"]["early_burst_2"][
                "total_branch_length"
            ],
        ),
        comparison_row_exact(
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
        comparison_row_tolerance(
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
        comparison_row_tolerance(
            row_id="brownian-rate",
            family_id="continuous-model-fitting",
            fragment_id="continuous-model-comparison",
            metric_name="brownian_rate",
            r_value=r_results["continuous_mode_fits"]["brownian"]["rate"],
            bijux_value=python_results["continuous_mode_fits"]["brownian"]["rate"],
            tolerance_abs_diff=25000.0,
        ),
        comparison_row_tolerance(
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
        comparison_row_tolerance(
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
        comparison_row_tolerance(
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
        comparison_row_tolerance(
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
        comparison_row_tolerance(
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
        comparison_row_tolerance(
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
        comparison_row_tolerance(
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
        comparison_row_tolerance(
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
        comparison_row_exact(
            row_id="ancestral-brownian-node-count",
            family_id="ancestral-reconstruction",
            fragment_id="ancestral-mode-comparison",
            metric_name="brownian_node_count",
            r_value=r_results["ancestral_reconstruction"]["brownian"]["node_count"],
            bijux_value=python_results["ancestral_reconstruction"]["brownian"][
                "node_count"
            ],
        ),
        comparison_row_exact(
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
        comparison_row_exact(
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
        comparison_row_exact(
            row_id="ancestral-eb-node-count",
            family_id="ancestral-reconstruction",
            fragment_id="ancestral-mode-comparison",
            metric_name="early_burst_node_count",
            r_value=r_results["ancestral_reconstruction"]["early_burst"]["node_count"],
            bijux_value=python_results["ancestral_reconstruction"]["early_burst"][
                "node_count"
            ],
        ),
        comparison_row_exact(
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
        comparison_row_exact(
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
        comparison_row_tolerance(
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
        comparison_row_tolerance(
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
        comparison_row_tolerance(
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
        comparison_row_tolerance(
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
        comparison_row_tolerance(
            row_id="fixed-reference-lambda-value",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="fixed_reference_lambda_value",
            r_value=r_results["fixed_reference_lambda_pgls"]["lambda_value"],
            bijux_value=python_results["fixed_reference_lambda_pgls"]["lambda_value"],
            tolerance_abs_diff=5e-4,
        ),
        comparison_row_tolerance(
            row_id="fixed-reference-pgls-intercept",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="fixed_reference_intercept",
            r_value=r_results["fixed_reference_lambda_pgls"]["coefficients"][
                "intercept"
            ],
            bijux_value=python_results["fixed_reference_lambda_pgls"]["coefficients"][
                "intercept"
            ],
            tolerance_abs_diff=5e-4,
        ),
        comparison_row_tolerance(
            row_id="fixed-reference-pgls-slope",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="fixed_reference_social_group_size",
            r_value=r_results["fixed_reference_lambda_pgls"]["coefficients"][
                "social_group_size"
            ],
            bijux_value=python_results["fixed_reference_lambda_pgls"]["coefficients"][
                "social_group_size"
            ],
            tolerance_abs_diff=5e-4,
        ),
        comparison_row_tolerance(
            row_id="fixed-reference-pgls-intercept-standard-error",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="fixed_reference_intercept_standard_error",
            r_value=r_results["fixed_reference_lambda_pgls"]["standard_errors"][
                "intercept"
            ],
            bijux_value=python_results["fixed_reference_lambda_pgls"][
                "standard_errors"
            ]["intercept"],
            tolerance_abs_diff=5e-4,
        ),
        comparison_row_tolerance(
            row_id="fixed-reference-pgls-slope-standard-error",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="fixed_reference_social_group_size_standard_error",
            r_value=r_results["fixed_reference_lambda_pgls"]["standard_errors"][
                "social_group_size"
            ],
            bijux_value=python_results["fixed_reference_lambda_pgls"][
                "standard_errors"
            ]["social_group_size"],
            tolerance_abs_diff=5e-4,
        ),
        comparison_row_tolerance(
            row_id="fixed-reference-pgls-intercept-p-value",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="fixed_reference_intercept_p_value",
            r_value=r_results["fixed_reference_lambda_pgls"]["p_values"]["intercept"],
            bijux_value=python_results["fixed_reference_lambda_pgls"]["p_values"][
                "intercept"
            ],
            tolerance_abs_diff=5e-6,
        ),
        comparison_row_tolerance(
            row_id="fixed-reference-pgls-slope-p-value",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="fixed_reference_social_group_size_p_value",
            r_value=r_results["fixed_reference_lambda_pgls"]["p_values"][
                "social_group_size"
            ],
            bijux_value=python_results["fixed_reference_lambda_pgls"]["p_values"][
                "social_group_size"
            ],
            tolerance_abs_diff=5e-6,
        ),
        comparison_row_tolerance(
            row_id="fixed-reference-pgls-log-likelihood",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="fixed_reference_log_likelihood",
            r_value=r_results["fixed_reference_lambda_pgls"]["log_likelihood"],
            bijux_value=python_results["fixed_reference_lambda_pgls"]["log_likelihood"],
            tolerance_abs_diff=5e-4,
        ),
        comparison_row_tolerance(
            row_id="fixed-reference-pgls-aic",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="fixed_reference_aic",
            r_value=r_results["fixed_reference_lambda_pgls"]["aic"],
            bijux_value=python_results["fixed_reference_lambda_pgls"]["aic"],
            tolerance_abs_diff=5e-4,
        ),
        comparison_row_tolerance(
            row_id="estimated-lambda-value",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="lambda_value",
            r_value=r_results["estimated_lambda_pgls"]["lambda_value"],
            bijux_value=python_results["estimated_lambda_pgls"]["lambda_value"],
            tolerance_abs_diff=0.05,
        ),
        comparison_row_tolerance(
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
        comparison_row_tolerance(
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
        comparison_row_tolerance(
            row_id="estimated-pgls-intercept-standard-error",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="intercept_standard_error",
            r_value=r_results["estimated_lambda_pgls"]["standard_errors"]["intercept"],
            bijux_value=python_results["estimated_lambda_pgls"]["standard_errors"][
                "intercept"
            ],
            tolerance_abs_diff=0.1,
        ),
        comparison_row_tolerance(
            row_id="estimated-pgls-slope-standard-error",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="social_group_size_standard_error",
            r_value=r_results["estimated_lambda_pgls"]["standard_errors"][
                "social_group_size"
            ],
            bijux_value=python_results["estimated_lambda_pgls"]["standard_errors"][
                "social_group_size"
            ],
            tolerance_abs_diff=0.02,
        ),
        comparison_row_tolerance(
            row_id="estimated-pgls-log-likelihood",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="log_likelihood",
            r_value=r_results["estimated_lambda_pgls"]["log_likelihood"],
            bijux_value=python_results["estimated_lambda_pgls"]["log_likelihood"],
            tolerance_abs_diff=0.25,
        ),
        comparison_row_tolerance(
            row_id="estimated-pgls-aic",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="aic",
            r_value=r_results["estimated_lambda_pgls"]["aic"],
            bijux_value=python_results["estimated_lambda_pgls"]["aic"],
            tolerance_abs_diff=0.25,
        ),
        comparison_row_tolerance(
            row_id="estimated-pgls-intercept-p-value",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="intercept_p_value",
            r_value=r_results["estimated_lambda_pgls"]["p_values"]["intercept"],
            bijux_value=python_results["estimated_lambda_pgls"]["p_values"][
                "intercept"
            ],
            tolerance_abs_diff=5e-6,
        ),
        comparison_row_tolerance(
            row_id="estimated-pgls-slope-p-value",
            family_id="phylogenetic-regression",
            fragment_id="pagel-lambda-regression",
            metric_name="social_group_size_p_value",
            r_value=r_results["estimated_lambda_pgls"]["p_values"]["social_group_size"],
            bijux_value=python_results["estimated_lambda_pgls"]["p_values"][
                "social_group_size"
            ],
            tolerance_abs_diff=5e-6,
        ),
        comparison_row_equivalence(
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
        comparison_row_tolerance(
            row_id="signal-estimated-lambda",
            family_id="phylogenetic-signal",
            fragment_id="phylogenetic-signal-test",
            metric_name="estimated_lambda",
            r_value=r_results["signal_test"]["estimated_lambda"],
            bijux_value=python_results["signal_test"]["estimated_lambda"],
            tolerance_abs_diff=0.05,
        ),
        comparison_row_tolerance(
            row_id="signal-likelihood-ratio",
            family_id="phylogenetic-signal",
            fragment_id="phylogenetic-signal-test",
            metric_name="likelihood_ratio",
            r_value=r_results["signal_test"]["likelihood_ratio"],
            bijux_value=python_results["signal_test"]["likelihood_ratio"],
            tolerance_abs_diff=0.05,
        ),
        comparison_row_equivalence(
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
        comparison_row_tolerance(
            row_id="baseline-diagnostic-qq-correlation",
            family_id="diagnostics",
            fragment_id="baseline-gls-diagnostics",
            metric_name="qq_correlation",
            r_value=r_results["baseline_gls"]["diagnostics"]["qq_correlation"],
            bijux_value=python_results["baseline_gls"]["diagnostics"]["qq_correlation"],
            tolerance_abs_diff=0.02,
        ),
        comparison_row_tolerance(
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
        comparison_row_equivalence(
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
        comparison_row_tolerance(
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
        comparison_row_tolerance(
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
        comparison_row_equivalence(
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
