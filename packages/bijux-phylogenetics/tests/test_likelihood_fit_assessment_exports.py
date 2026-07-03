from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    LikelihoodFitAssessment,
    build_likelihood_fit_assessment,
)
from bijux_phylogenetics.phylo.likelihood.fit_assessment import (
    build_likelihood_fit_assessment as build_likelihood_fit_assessment_impl,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    LikelihoodFitAssessment as LikelihoodFitAssessmentImpl,
)


def test_public_likelihood_exports_fit_assessment_surface() -> None:
    assert likelihood_api.LikelihoodFitAssessment is LikelihoodFitAssessment
    assert LikelihoodFitAssessment is LikelihoodFitAssessmentImpl
    assert (
        likelihood_api.build_likelihood_fit_assessment
        is build_likelihood_fit_assessment
    )
    assert build_likelihood_fit_assessment is build_likelihood_fit_assessment_impl
