from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    CONTINUOUS_TRAIT_PRIOR_MODES,
    CONTINUOUS_TRAIT_PRIOR_TARGETS,
    CONTINUOUS_TRAIT_PROBABILITY_PRIOR_FAMILIES,
    CONTINUOUS_TRAIT_SCALAR_PRIOR_FAMILIES,
    ContinuousTraitModelPriorBundle,
    ContinuousTraitModelPriorEvaluationReport,
    ContinuousTraitModelPriorRow,
    ContinuousTraitProbabilityPriorModel,
    ContinuousTraitScalarPriorModel,
    build_beta_continuous_trait_probability_prior,
    build_continuous_trait_model_prior_bundle,
    build_exponential_continuous_trait_scalar_prior,
    build_fixed_continuous_trait_probability_prior,
    build_fixed_continuous_trait_scalar_prior,
    build_gamma_continuous_trait_scalar_prior,
    build_lognormal_continuous_trait_scalar_prior,
    evaluate_continuous_trait_model_log_prior,
    evaluate_continuous_trait_probability_log_prior,
    evaluate_continuous_trait_scalar_log_prior,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    CONTINUOUS_TRAIT_PRIOR_MODES as CONTINUOUS_TRAIT_PRIOR_MODES_IMPL,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    CONTINUOUS_TRAIT_PRIOR_TARGETS as CONTINUOUS_TRAIT_PRIOR_TARGETS_IMPL,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    CONTINUOUS_TRAIT_PROBABILITY_PRIOR_FAMILIES as CONTINUOUS_TRAIT_PROBABILITY_PRIOR_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    CONTINUOUS_TRAIT_SCALAR_PRIOR_FAMILIES as CONTINUOUS_TRAIT_SCALAR_PRIOR_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    ContinuousTraitModelPriorBundle as ContinuousTraitModelPriorBundleImpl,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    ContinuousTraitModelPriorEvaluationReport as ContinuousTraitModelPriorEvaluationReportImpl,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    ContinuousTraitModelPriorRow as ContinuousTraitModelPriorRowImpl,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    ContinuousTraitProbabilityPriorModel as ContinuousTraitProbabilityPriorModelImpl,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    ContinuousTraitScalarPriorModel as ContinuousTraitScalarPriorModelImpl,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    build_beta_continuous_trait_probability_prior as build_beta_continuous_trait_probability_prior_impl,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    build_continuous_trait_model_prior_bundle as build_continuous_trait_model_prior_bundle_impl,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    build_exponential_continuous_trait_scalar_prior as build_exponential_continuous_trait_scalar_prior_impl,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    build_fixed_continuous_trait_probability_prior as build_fixed_continuous_trait_probability_prior_impl,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    build_fixed_continuous_trait_scalar_prior as build_fixed_continuous_trait_scalar_prior_impl,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    build_gamma_continuous_trait_scalar_prior as build_gamma_continuous_trait_scalar_prior_impl,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    build_lognormal_continuous_trait_scalar_prior as build_lognormal_continuous_trait_scalar_prior_impl,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    evaluate_continuous_trait_model_log_prior as evaluate_continuous_trait_model_log_prior_impl,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    evaluate_continuous_trait_probability_log_prior as evaluate_continuous_trait_probability_log_prior_impl,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    evaluate_continuous_trait_scalar_log_prior as evaluate_continuous_trait_scalar_log_prior_impl,
)


def test_bayesian_exports_continuous_trait_model_prior_surface() -> None:
    assert CONTINUOUS_TRAIT_SCALAR_PRIOR_FAMILIES == (
        CONTINUOUS_TRAIT_SCALAR_PRIOR_FAMILIES_IMPL
    )
    assert CONTINUOUS_TRAIT_PROBABILITY_PRIOR_FAMILIES == (
        CONTINUOUS_TRAIT_PROBABILITY_PRIOR_FAMILIES_IMPL
    )
    assert CONTINUOUS_TRAIT_PRIOR_MODES == CONTINUOUS_TRAIT_PRIOR_MODES_IMPL
    assert CONTINUOUS_TRAIT_PRIOR_TARGETS == CONTINUOUS_TRAIT_PRIOR_TARGETS_IMPL
    assert ContinuousTraitScalarPriorModel is ContinuousTraitScalarPriorModelImpl
    assert (
        ContinuousTraitProbabilityPriorModel is ContinuousTraitProbabilityPriorModelImpl
    )
    assert ContinuousTraitModelPriorBundle is ContinuousTraitModelPriorBundleImpl
    assert (
        ContinuousTraitModelPriorEvaluationReport
        is ContinuousTraitModelPriorEvaluationReportImpl
    )
    assert ContinuousTraitModelPriorRow is ContinuousTraitModelPriorRowImpl
    assert (
        build_exponential_continuous_trait_scalar_prior
        is build_exponential_continuous_trait_scalar_prior_impl
    )
    assert (
        build_gamma_continuous_trait_scalar_prior
        is build_gamma_continuous_trait_scalar_prior_impl
    )
    assert (
        build_lognormal_continuous_trait_scalar_prior
        is build_lognormal_continuous_trait_scalar_prior_impl
    )
    assert (
        build_fixed_continuous_trait_scalar_prior
        is build_fixed_continuous_trait_scalar_prior_impl
    )
    assert (
        build_beta_continuous_trait_probability_prior
        is build_beta_continuous_trait_probability_prior_impl
    )
    assert (
        build_fixed_continuous_trait_probability_prior
        is build_fixed_continuous_trait_probability_prior_impl
    )
    assert (
        build_continuous_trait_model_prior_bundle
        is build_continuous_trait_model_prior_bundle_impl
    )
    assert (
        evaluate_continuous_trait_scalar_log_prior
        is evaluate_continuous_trait_scalar_log_prior_impl
    )
    assert (
        evaluate_continuous_trait_probability_log_prior
        is evaluate_continuous_trait_probability_log_prior_impl
    )
    assert (
        evaluate_continuous_trait_model_log_prior
        is evaluate_continuous_trait_model_log_prior_impl
    )
