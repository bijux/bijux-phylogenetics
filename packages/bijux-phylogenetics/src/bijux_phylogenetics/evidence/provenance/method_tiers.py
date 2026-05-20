from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

MethodTier = Literal["supported", "experimental", "advisory", "parser-only"]
InferenceMode = Literal["inference", "review-only", "parser-only"]


@dataclass(frozen=True, slots=True)
class MethodTierAssessment:
    """One durable trust classification for a user-facing analysis surface."""

    tier: MethodTier
    surface: str
    summary: str
    validation_basis: tuple[str, ...]
    inference_mode: InferenceMode
    approximation: str | None = None
    excluded_reference_surfaces: tuple[str, ...] = ()
    warning: str | None = None

    def __post_init__(self) -> None:
        if self.tier == "supported" and not any(
            basis.startswith(("reference-parity:", "real-engine-validation:"))
            for basis in self.validation_basis
        ):
            raise ValueError(
                "supported method tiers require reference parity or real-engine validation evidence"
            )
        if self.tier == "experimental" and self.warning is None:
            raise ValueError("experimental method tiers require a clear warning")
        if self.tier == "parser-only" and self.inference_mode != "parser-only":
            raise ValueError(
                "parser-only method tiers must declare parser-only inference mode"
            )


def fasta_to_tree_method_tier() -> MethodTierAssessment:
    return MethodTierAssessment(
        tier="supported",
        surface="fasta-to-tree",
        summary=(
            "This workflow runs governed external alignment and maximum-likelihood inference with real-engine validation coverage."
        ),
        validation_basis=(
            "real-engine-validation:fasta-to-tree",
            "artifact-schema:fasta-to-tree",
        ),
        inference_mode="inference",
    )


def comparative_report_method_tier() -> MethodTierAssessment:
    return MethodTierAssessment(
        tier="supported",
        surface="comparative-report-package",
        summary=(
            "This comparative package summarizes governed PGLS, signal, and trait-model surfaces backed by checked reference parity."
        ),
        validation_basis=(
            "reference-parity:pgls",
            "reference-parity:phylogenetic-signal",
        ),
        inference_mode="inference",
    )


def phylogenetic_logistic_method_tier(
    approximation_method: str,
) -> MethodTierAssessment:
    return MethodTierAssessment(
        tier="experimental",
        surface="phylogenetic-logistic",
        summary=(
            "This workflow fits an approximate phylogenetic logistic surface, should be treated as exploratory rather than publication-grade inference, and does not claim ape::compar.gee parity."
        ),
        validation_basis=(),
        inference_mode="inference",
        approximation=approximation_method,
        excluded_reference_surfaces=("ape::compar.gee",),
        warning=(
            "experimental method tier: phylogenetic logistic results use an approximate working-correlation fit and need cautious review before biological interpretation."
        ),
    )


def tree_report_method_tier() -> MethodTierAssessment:
    return MethodTierAssessment(
        tier="advisory",
        surface="tree-report-package",
        summary=(
            "This package audits and renders an existing tree for review; it does not perform phylogenetic inference."
        ),
        validation_basis=(
            "validation-surface:tree-audit",
            "support-reference:rendered-support-audit",
        ),
        inference_mode="review-only",
    )


def bayesian_report_method_tier(report_kind: str) -> MethodTierAssessment:
    return MethodTierAssessment(
        tier="parser-only",
        surface=report_kind,
        summary=(
            "This surface summarizes posterior artifacts produced by an external Bayesian engine and does not claim that Bijux performed the inference."
        ),
        validation_basis=(
            "parser-contract:bayesian-artifact-validation",
            "governed-reference:posterior-artifact-corpus",
        ),
        inference_mode="parser-only",
        warning=(
            "parser-only method tier: this report parses or audits external Bayesian outputs and does not itself infer posterior trees."
        ),
    )


def release_method_tier_inventory() -> list[MethodTierAssessment]:
    """Return the governed user-facing method inventory for release reporting."""
    return [
        fasta_to_tree_method_tier(),
        comparative_report_method_tier(),
        phylogenetic_logistic_method_tier("phylogenetic-working-correlation-gee"),
        tree_report_method_tier(),
        bayesian_report_method_tier("bayesian-posterior-report"),
        bayesian_report_method_tier("bayesian-diagnostics-report"),
        bayesian_report_method_tier("beast-calibration-report"),
        bayesian_report_method_tier("time-tree-readiness-report"),
        bayesian_report_method_tier("bayesian-run-comparison-report"),
    ]


def method_tier_metrics(assessment: MethodTierAssessment) -> dict[str, object]:
    return {
        "method_tier": assessment.tier,
        "method_surface": assessment.surface,
        "method_inference_mode": assessment.inference_mode,
        "method_validation_basis": list(assessment.validation_basis),
        "method_approximation": assessment.approximation,
        "method_excluded_reference_surfaces": list(
            assessment.excluded_reference_surfaces
        ),
    }


def method_tier_warnings(assessment: MethodTierAssessment) -> list[str]:
    warnings: list[str] = []
    if assessment.warning is not None:
        warnings.append(assessment.warning)
    for reference_surface in assessment.excluded_reference_surfaces:
        warnings.append(
            f"explicit non-claim: {assessment.surface} does not currently claim {reference_surface} parity and should not be interpreted as a drop-in {reference_surface} implementation."
        )
    return warnings
