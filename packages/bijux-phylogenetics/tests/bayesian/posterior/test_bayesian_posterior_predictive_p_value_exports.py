from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    PosteriorPredictivePValueReport,
    PosteriorPredictivePValueRow,
    summarize_posterior_predictive_p_values,
)
from bijux_phylogenetics.bayesian.posterior_predictive_p_values import (
    PosteriorPredictivePValueReport as PosteriorPredictivePValueReportImpl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_p_values import (
    PosteriorPredictivePValueRow as PosteriorPredictivePValueRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_p_values import (
    summarize_posterior_predictive_p_values as summarize_posterior_predictive_p_values_impl,
)


def test_bayesian_exports_posterior_predictive_p_value_surface() -> None:
    assert PosteriorPredictivePValueReport is PosteriorPredictivePValueReportImpl
    assert PosteriorPredictivePValueRow is PosteriorPredictivePValueRowImpl
    assert (
        summarize_posterior_predictive_p_values
        is summarize_posterior_predictive_p_values_impl
    )
