from __future__ import annotations

import json

from bijux_phylogenetics.evidence.provenance.method_tiers import MethodTierAssessment


def method_tier_summary_metrics(
    method_tier: MethodTierAssessment,
) -> list[tuple[str, str]]:
    return [
        ("Method Tier", method_tier.tier),
        ("Inference Mode", method_tier.inference_mode),
        (
            "Approximation",
            "none" if method_tier.approximation is None else method_tier.approximation,
        ),
    ]


def method_tier_section(method_tier: MethodTierAssessment) -> tuple[str, str]:
    lines = [
        f"tier: {method_tier.tier}",
        f"inference_mode: {method_tier.inference_mode}",
        f"summary: {method_tier.summary}",
    ]
    if method_tier.validation_basis:
        lines.append("validation_basis: " + "; ".join(method_tier.validation_basis))
    if method_tier.approximation is not None:
        lines.append(f"approximation: {method_tier.approximation}")
    if method_tier.warning is not None:
        lines.append(f"warning: {method_tier.warning}")
    return ("method-tier", "\n".join(lines))


def posterior_report_limitations(convergence_warnings: list[str]) -> list[str]:
    limitations = [
        "posterior clade support summarizes sampled trees under the fitted Bayesian model and should not be treated as direct proof of clade truth",
        "consensus and clade-frequency summaries can hide minority topologies, so interpretation should remain tied to convergence and tree-set dispersion checks",
        *convergence_warnings,
    ]
    return _deduplicate_limitations(limitations)


def calibration_audit_limitations(
    *,
    invalid_calibration_count: int,
    impossible_constraint_count: int,
    invalid_tip_count: int,
) -> list[str]:
    limitations = [
        "calibration and tip-date audits validate compatibility of the supplied constraints with the current tree and do not by themselves justify any dated-tree model choice",
        "a passing audit does not guarantee that downstream divergence-time estimates are robust to alternative calibrations, clock models, or taxon sampling decisions",
    ]
    if invalid_calibration_count:
        limitations.append(
            f"{invalid_calibration_count} calibration rows remain invalid and must be corrected before dated-tree interpretation is trusted"
        )
    if impossible_constraint_count:
        limitations.append(
            f"{impossible_constraint_count} impossible calibration constraints still conflict with the current topology"
        )
    if invalid_tip_count:
        limitations.append(
            f"{invalid_tip_count} tip-date rows remain invalid or mismatched to the current tree"
        )
    return _deduplicate_limitations(limitations)


def run_comparison_limitations(warnings: list[str]) -> list[str]:
    limitations = [
        "agreement between independent Bayesian runs only supports stability under the supplied model and priors and does not validate overall model adequacy",
        "parameter or topology differences across runs should block strong posterior interpretation until chain mixing, burn-in choice, and run configuration are reconciled",
        *warnings,
    ]
    return _deduplicate_limitations(limitations)


def ml_vs_bayesian_limitations(warnings: list[str]) -> list[str]:
    limitations = [
        "agreement or disagreement between maximum-likelihood and Bayesian summaries does not identify which inference framework is correct without external model checking",
        "topology and branch-length differences between ML and Bayesian trees should not be overinterpreted as biological rate shifts or timing evidence without checking model and taxon assumptions",
        *warnings,
    ]
    return _deduplicate_limitations(limitations)


def _deduplicate_limitations(limitations: list[str]) -> list[str]:
    normalized = []
    for item in limitations:
        text = item if isinstance(item, str) else json.dumps(item, sort_keys=True)
        text = text.strip()
        if text:
            normalized.append(text)
    return sorted(dict.fromkeys(normalized))
