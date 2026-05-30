from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    CONTINUOUS_TRAIT_LOCATION_PRIOR_FAMILIES,
    ContinuousTraitLocationPriorModel,
    build_fixed_continuous_trait_location_prior,
    build_normal_continuous_trait_location_prior,
    evaluate_continuous_trait_location_log_prior,
)
from bijux_phylogenetics.bayesian.continuous_trait_location_priors import (
    CONTINUOUS_TRAIT_LOCATION_PRIOR_FAMILIES as CONTINUOUS_TRAIT_LOCATION_PRIOR_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.continuous_trait_location_priors import (
    ContinuousTraitLocationPriorModel as ContinuousTraitLocationPriorModelImpl,
)
from bijux_phylogenetics.bayesian.continuous_trait_location_priors import (
    build_fixed_continuous_trait_location_prior as build_fixed_continuous_trait_location_prior_impl,
)
from bijux_phylogenetics.bayesian.continuous_trait_location_priors import (
    build_normal_continuous_trait_location_prior as build_normal_continuous_trait_location_prior_impl,
)
from bijux_phylogenetics.bayesian.continuous_trait_location_priors import (
    evaluate_continuous_trait_location_log_prior as evaluate_continuous_trait_location_log_prior_impl,
)


def test_bayesian_exports_continuous_trait_location_prior_surface() -> None:
    assert (
        CONTINUOUS_TRAIT_LOCATION_PRIOR_FAMILIES
        == CONTINUOUS_TRAIT_LOCATION_PRIOR_FAMILIES_IMPL
    )
    assert ContinuousTraitLocationPriorModel is ContinuousTraitLocationPriorModelImpl
    assert (
        build_fixed_continuous_trait_location_prior
        is build_fixed_continuous_trait_location_prior_impl
    )
    assert (
        build_normal_continuous_trait_location_prior
        is build_normal_continuous_trait_location_prior_impl
    )
    assert (
        evaluate_continuous_trait_location_log_prior
        is evaluate_continuous_trait_location_log_prior_impl
    )
