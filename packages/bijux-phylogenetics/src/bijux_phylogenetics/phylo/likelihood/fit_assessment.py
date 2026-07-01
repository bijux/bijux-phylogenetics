from __future__ import annotations

from collections.abc import Sequence

from bijux_phylogenetics.phylo.likelihood.models import (
    LikelihoodFitAssessment,
    LikelihoodOptimizationBoundaryWarning,
)

_FLAT_WARNING_PREFIX = "flat-"
_WEAK_PATTERN_MULTIPLIER = 2


def build_likelihood_fit_assessment(
    *,
    site_count: int,
    pattern_count: int,
    optimized_parameter_count: int,
    converged: bool,
    boundary_warnings: Sequence[LikelihoodOptimizationBoundaryWarning],
) -> LikelihoodFitAssessment:
    """Classify one likelihood fit by optimizer behavior and information content."""
    warning_kinds = [warning.warning_kind for warning in boundary_warnings]
    if not converged:
        return LikelihoodFitAssessment(
            classification="failed",
            summary="optimizer did not converge to a trustworthy likelihood fit",
            optimized_parameter_count=optimized_parameter_count,
            boundary_warning_count=len(boundary_warnings),
            warning_kinds=warning_kinds,
            stable_conclusion_supported=False,
        )
    if _is_overparameterized(
        site_count=site_count,
        optimized_parameter_count=optimized_parameter_count,
    ):
        return LikelihoodFitAssessment(
            classification="overparameterized",
            summary=(
                "site count is too small for the fitted free-parameter count to support "
                "finite-sample interpretation"
            ),
            optimized_parameter_count=optimized_parameter_count,
            boundary_warning_count=len(boundary_warnings),
            warning_kinds=warning_kinds,
            stable_conclusion_supported=False,
        )
    if any(kind.startswith(_FLAT_WARNING_PREFIX) for kind in warning_kinds):
        return LikelihoodFitAssessment(
            classification="flat",
            summary="fit converged onto a flat likelihood boundary",
            optimized_parameter_count=optimized_parameter_count,
            boundary_warning_count=len(boundary_warnings),
            warning_kinds=warning_kinds,
            stable_conclusion_supported=False,
        )
    if boundary_warnings:
        return LikelihoodFitAssessment(
            classification="boundary",
            summary="fit converged only after hitting one or more optimizer boundaries",
            optimized_parameter_count=optimized_parameter_count,
            boundary_warning_count=len(boundary_warnings),
            warning_kinds=warning_kinds,
            stable_conclusion_supported=False,
        )
    if _has_weak_pattern_support(
        pattern_count=pattern_count,
        optimized_parameter_count=optimized_parameter_count,
    ):
        return LikelihoodFitAssessment(
            classification="weak",
            summary=(
                "pattern support is limited relative to the fitted free-parameter count"
            ),
            optimized_parameter_count=optimized_parameter_count,
            boundary_warning_count=0,
            warning_kinds=[],
            stable_conclusion_supported=False,
        )
    return LikelihoodFitAssessment(
        classification="stable",
        summary="fit converged away from boundaries with adequate pattern support",
        optimized_parameter_count=optimized_parameter_count,
        boundary_warning_count=0,
        warning_kinds=[],
        stable_conclusion_supported=True,
    )


def _is_overparameterized(
    *,
    site_count: int,
    optimized_parameter_count: int,
) -> bool:
    return optimized_parameter_count > 0 and site_count <= (
        optimized_parameter_count + 1
    )


def _has_weak_pattern_support(
    *,
    pattern_count: int,
    optimized_parameter_count: int,
) -> bool:
    return optimized_parameter_count > 0 and pattern_count <= (
        (_WEAK_PATTERN_MULTIPLIER * optimized_parameter_count) + 1
    )


__all__ = ["build_likelihood_fit_assessment"]
