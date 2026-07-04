from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    PosteriorModelAveragedEstimateRow,
    PosteriorModelAveragingReport,
    PosteriorModelEstimateRow,
    PosteriorModelSupportRow,
    summarize_metropolis_hastings_model_averaged_estimates,
    summarize_posterior_model_averaged_estimates,
)
from bijux_phylogenetics.bayesian.posterior_model_averaging import (
    PosteriorModelAveragedEstimateRow as PosteriorModelAveragedEstimateRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_model_averaging import (
    PosteriorModelAveragingReport as PosteriorModelAveragingReportImpl,
)
from bijux_phylogenetics.bayesian.posterior_model_averaging import (
    PosteriorModelEstimateRow as PosteriorModelEstimateRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_model_averaging import (
    PosteriorModelSupportRow as PosteriorModelSupportRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_model_averaging import (
    summarize_metropolis_hastings_model_averaged_estimates as summarize_metropolis_hastings_model_averaged_estimates_impl,
)
from bijux_phylogenetics.bayesian.posterior_model_averaging import (
    summarize_posterior_model_averaged_estimates as summarize_posterior_model_averaged_estimates_impl,
)


def test_bayesian_exports_posterior_model_averaging_surface() -> None:
    assert PosteriorModelAveragedEstimateRow is PosteriorModelAveragedEstimateRowImpl
    assert PosteriorModelAveragingReport is PosteriorModelAveragingReportImpl
    assert PosteriorModelEstimateRow is PosteriorModelEstimateRowImpl
    assert PosteriorModelSupportRow is PosteriorModelSupportRowImpl
    assert (
        summarize_metropolis_hastings_model_averaged_estimates
        is summarize_metropolis_hastings_model_averaged_estimates_impl
    )
    assert (
        summarize_posterior_model_averaged_estimates
        is summarize_posterior_model_averaged_estimates_impl
    )
