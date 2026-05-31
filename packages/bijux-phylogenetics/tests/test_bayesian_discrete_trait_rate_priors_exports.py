from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    DISCRETE_TRAIT_RATE_PRIOR_FAMILIES,
    DISCRETE_TRAIT_RATE_PRIOR_MODELS,
    DiscreteTraitRatePriorEvaluationReport,
    DiscreteTraitRatePriorModel,
    DiscreteTraitRatePriorRow,
    build_exponential_discrete_trait_rate_prior,
    build_gamma_discrete_trait_rate_prior,
    build_lognormal_discrete_trait_rate_prior,
    evaluate_discrete_trait_rate_log_prior,
    evaluate_discrete_trait_rate_value_log_prior,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    DISCRETE_TRAIT_RATE_PRIOR_FAMILIES as DISCRETE_TRAIT_RATE_PRIOR_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    DISCRETE_TRAIT_RATE_PRIOR_MODELS as DISCRETE_TRAIT_RATE_PRIOR_MODELS_IMPL,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    DiscreteTraitRatePriorEvaluationReport as DiscreteTraitRatePriorEvaluationReportImpl,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    DiscreteTraitRatePriorModel as DiscreteTraitRatePriorModelImpl,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    DiscreteTraitRatePriorRow as DiscreteTraitRatePriorRowImpl,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    build_exponential_discrete_trait_rate_prior as build_exponential_discrete_trait_rate_prior_impl,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    build_gamma_discrete_trait_rate_prior as build_gamma_discrete_trait_rate_prior_impl,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    build_lognormal_discrete_trait_rate_prior as build_lognormal_discrete_trait_rate_prior_impl,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    evaluate_discrete_trait_rate_log_prior as evaluate_discrete_trait_rate_log_prior_impl,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    evaluate_discrete_trait_rate_value_log_prior as evaluate_discrete_trait_rate_value_log_prior_impl,
)


def test_bayesian_exports_discrete_trait_rate_prior_surface() -> None:
    assert DISCRETE_TRAIT_RATE_PRIOR_FAMILIES == DISCRETE_TRAIT_RATE_PRIOR_FAMILIES_IMPL
    assert DISCRETE_TRAIT_RATE_PRIOR_MODELS == DISCRETE_TRAIT_RATE_PRIOR_MODELS_IMPL
    assert DiscreteTraitRatePriorModel is DiscreteTraitRatePriorModelImpl
    assert (
        DiscreteTraitRatePriorEvaluationReport
        is DiscreteTraitRatePriorEvaluationReportImpl
    )
    assert DiscreteTraitRatePriorRow is DiscreteTraitRatePriorRowImpl
    assert (
        build_exponential_discrete_trait_rate_prior
        is build_exponential_discrete_trait_rate_prior_impl
    )
    assert (
        build_gamma_discrete_trait_rate_prior
        is build_gamma_discrete_trait_rate_prior_impl
    )
    assert (
        build_lognormal_discrete_trait_rate_prior
        is build_lognormal_discrete_trait_rate_prior_impl
    )
    assert (
        evaluate_discrete_trait_rate_log_prior
        is evaluate_discrete_trait_rate_log_prior_impl
    )
    assert (
        evaluate_discrete_trait_rate_value_log_prior
        is evaluate_discrete_trait_rate_value_log_prior_impl
    )
