from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    PosteriorContinuousTraitNodeSummaryRow,
    PosteriorContinuousTraitReport,
    PosteriorDiscreteTraitNodeSummaryRow,
    PosteriorDiscreteTraitReport,
    PosteriorDiscreteTraitStateProbabilityRow,
    summarize_brownian_continuous_trait_posterior_ancestral_states,
    summarize_continuous_trait_posterior_ancestral_states,
    summarize_discrete_trait_mk_posterior_ancestral_states,
    summarize_ornstein_uhlenbeck_continuous_trait_posterior_ancestral_states,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_traits import (
    PosteriorContinuousTraitNodeSummaryRow as PosteriorContinuousTraitNodeSummaryRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_traits import (
    PosteriorContinuousTraitReport as PosteriorContinuousTraitReportImpl,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_traits import (
    PosteriorDiscreteTraitNodeSummaryRow as PosteriorDiscreteTraitNodeSummaryRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_traits import (
    PosteriorDiscreteTraitReport as PosteriorDiscreteTraitReportImpl,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_traits import (
    PosteriorDiscreteTraitStateProbabilityRow as PosteriorDiscreteTraitStateProbabilityRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_traits import (
    summarize_brownian_continuous_trait_posterior_ancestral_states as summarize_brownian_continuous_trait_posterior_ancestral_states_impl,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_traits import (
    summarize_continuous_trait_posterior_ancestral_states as summarize_continuous_trait_posterior_ancestral_states_impl,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_traits import (
    summarize_discrete_trait_mk_posterior_ancestral_states as summarize_discrete_trait_mk_posterior_ancestral_states_impl,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_traits import (
    summarize_ornstein_uhlenbeck_continuous_trait_posterior_ancestral_states as summarize_ornstein_uhlenbeck_continuous_trait_posterior_ancestral_states_impl,
)


def test_bayesian_exports_posterior_ancestral_trait_surface() -> None:
    assert (
        PosteriorContinuousTraitNodeSummaryRow
        is PosteriorContinuousTraitNodeSummaryRowImpl
    )
    assert PosteriorContinuousTraitReport is PosteriorContinuousTraitReportImpl
    assert (
        PosteriorDiscreteTraitNodeSummaryRow is PosteriorDiscreteTraitNodeSummaryRowImpl
    )
    assert PosteriorDiscreteTraitReport is PosteriorDiscreteTraitReportImpl
    assert (
        PosteriorDiscreteTraitStateProbabilityRow
        is PosteriorDiscreteTraitStateProbabilityRowImpl
    )
    assert (
        summarize_brownian_continuous_trait_posterior_ancestral_states
        is summarize_brownian_continuous_trait_posterior_ancestral_states_impl
    )
    assert (
        summarize_continuous_trait_posterior_ancestral_states
        is summarize_continuous_trait_posterior_ancestral_states_impl
    )
    assert (
        summarize_discrete_trait_mk_posterior_ancestral_states
        is summarize_discrete_trait_mk_posterior_ancestral_states_impl
    )
    assert (
        summarize_ornstein_uhlenbeck_continuous_trait_posterior_ancestral_states
        is summarize_ornstein_uhlenbeck_continuous_trait_posterior_ancestral_states_impl
    )
