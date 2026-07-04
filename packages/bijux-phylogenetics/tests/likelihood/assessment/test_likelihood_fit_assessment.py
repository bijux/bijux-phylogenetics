from __future__ import annotations

from bijux_phylogenetics.phylo.likelihood import (
    build_likelihood_fit_assessment,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    LikelihoodOptimizationBoundaryWarning,
)


def test_likelihood_fit_assessment_classifies_stable_fit() -> None:
    assessment = build_likelihood_fit_assessment(
        site_count=200,
        pattern_count=40,
        optimized_parameter_count=4,
        converged=True,
        boundary_warnings=[],
    )

    assert assessment.classification == "stable"
    assert assessment.stable_conclusion_supported is True


def test_likelihood_fit_assessment_classifies_weak_fit() -> None:
    assessment = build_likelihood_fit_assessment(
        site_count=200,
        pattern_count=5,
        optimized_parameter_count=2,
        converged=True,
        boundary_warnings=[],
    )

    assert assessment.classification == "weak"
    assert assessment.stable_conclusion_supported is False


def test_likelihood_fit_assessment_classifies_boundary_fit() -> None:
    assessment = build_likelihood_fit_assessment(
        site_count=200,
        pattern_count=40,
        optimized_parameter_count=2,
        converged=True,
        boundary_warnings=[_warning("branch-boundary")],
    )

    assert assessment.classification == "boundary"
    assert assessment.warning_kinds == ["branch-boundary"]


def test_likelihood_fit_assessment_classifies_overparameterized_fit() -> None:
    assessment = build_likelihood_fit_assessment(
        site_count=5,
        pattern_count=5,
        optimized_parameter_count=4,
        converged=True,
        boundary_warnings=[],
    )

    assert assessment.classification == "overparameterized"
    assert assessment.stable_conclusion_supported is False


def test_likelihood_fit_assessment_classifies_flat_fit() -> None:
    assessment = build_likelihood_fit_assessment(
        site_count=200,
        pattern_count=40,
        optimized_parameter_count=2,
        converged=True,
        boundary_warnings=[_warning("flat-alpha")],
    )

    assert assessment.classification == "flat"
    assert assessment.warning_kinds == ["flat-alpha"]


def test_likelihood_fit_assessment_classifies_failed_fit() -> None:
    assessment = build_likelihood_fit_assessment(
        site_count=200,
        pattern_count=40,
        optimized_parameter_count=2,
        converged=False,
        boundary_warnings=[],
    )

    assert assessment.classification == "failed"
    assert assessment.stable_conclusion_supported is False


def _warning(kind: str) -> LikelihoodOptimizationBoundaryWarning:
    return LikelihoodOptimizationBoundaryWarning(
        warning_kind=kind,
        affected_parameter="alpha",
        boundary_side="upper",
        observed_value=20.0,
        lower_bound=0.05,
        upper_bound=20.0,
        affected_branch_id=None,
        affected_branch_clade_id=None,
        message=f"{kind} warning",
    )
