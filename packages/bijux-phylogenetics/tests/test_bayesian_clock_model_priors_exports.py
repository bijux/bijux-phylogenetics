from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    CLOCK_MODEL_SCALAR_PRIOR_FAMILIES,
    ClockModelScalarPriorModel,
    build_exponential_clock_model_scalar_prior,
    build_fixed_clock_model_scalar_prior,
    build_gamma_clock_model_scalar_prior,
    build_lognormal_clock_model_scalar_prior,
    evaluate_clock_model_scalar_log_prior,
)
from bijux_phylogenetics.bayesian.clock_model_priors import (
    CLOCK_MODEL_SCALAR_PRIOR_FAMILIES as CLOCK_MODEL_SCALAR_PRIOR_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.clock_model_priors import (
    ClockModelScalarPriorModel as ClockModelScalarPriorModelImpl,
)
from bijux_phylogenetics.bayesian.clock_model_priors import (
    build_exponential_clock_model_scalar_prior as build_exponential_clock_model_scalar_prior_impl,
)
from bijux_phylogenetics.bayesian.clock_model_priors import (
    build_fixed_clock_model_scalar_prior as build_fixed_clock_model_scalar_prior_impl,
)
from bijux_phylogenetics.bayesian.clock_model_priors import (
    build_gamma_clock_model_scalar_prior as build_gamma_clock_model_scalar_prior_impl,
)
from bijux_phylogenetics.bayesian.clock_model_priors import (
    build_lognormal_clock_model_scalar_prior as build_lognormal_clock_model_scalar_prior_impl,
)
from bijux_phylogenetics.bayesian.clock_model_priors import (
    evaluate_clock_model_scalar_log_prior as evaluate_clock_model_scalar_log_prior_impl,
)


def test_bayesian_exports_clock_model_scalar_prior_surface() -> None:
    assert CLOCK_MODEL_SCALAR_PRIOR_FAMILIES == CLOCK_MODEL_SCALAR_PRIOR_FAMILIES_IMPL
    assert ClockModelScalarPriorModel is ClockModelScalarPriorModelImpl
    assert (
        build_exponential_clock_model_scalar_prior
        is build_exponential_clock_model_scalar_prior_impl
    )
    assert (
        build_fixed_clock_model_scalar_prior
        is build_fixed_clock_model_scalar_prior_impl
    )
    assert (
        build_gamma_clock_model_scalar_prior
        is build_gamma_clock_model_scalar_prior_impl
    )
    assert (
        build_lognormal_clock_model_scalar_prior
        is build_lognormal_clock_model_scalar_prior_impl
    )
    assert (
        evaluate_clock_model_scalar_log_prior
        is evaluate_clock_model_scalar_log_prior_impl
    )
