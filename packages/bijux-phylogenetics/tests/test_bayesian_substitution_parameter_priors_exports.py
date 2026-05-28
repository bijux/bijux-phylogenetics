from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    POSITIVE_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES,
    PROBABILITY_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES,
    SIMPLEX_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES,
    SUBSTITUTION_PARAMETER_PRIOR_TARGETS,
    PositiveSubstitutionParameterPriorModel,
    ProbabilitySubstitutionParameterPriorModel,
    SimplexSubstitutionParameterPriorModel,
    SubstitutionParameterPriorBundle,
    SubstitutionParameterPriorEvaluationReport,
    SubstitutionParameterPriorRow,
    build_beta_probability_substitution_parameter_prior,
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_fixed_positive_substitution_parameter_prior,
    build_fixed_probability_substitution_parameter_prior,
    build_fixed_simplex_substitution_parameter_prior,
    build_gamma_positive_substitution_parameter_prior,
    build_lognormal_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
    evaluate_substitution_parameter_log_prior,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    POSITIVE_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES as POSITIVE_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    PROBABILITY_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES as PROBABILITY_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    SIMPLEX_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES as SIMPLEX_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    SUBSTITUTION_PARAMETER_PRIOR_TARGETS as SUBSTITUTION_PARAMETER_PRIOR_TARGETS_IMPL,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    PositiveSubstitutionParameterPriorModel as PositiveSubstitutionParameterPriorModelImpl,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    ProbabilitySubstitutionParameterPriorModel as ProbabilitySubstitutionParameterPriorModelImpl,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    SimplexSubstitutionParameterPriorModel as SimplexSubstitutionParameterPriorModelImpl,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    SubstitutionParameterPriorBundle as SubstitutionParameterPriorBundleImpl,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    SubstitutionParameterPriorEvaluationReport as SubstitutionParameterPriorEvaluationReportImpl,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    SubstitutionParameterPriorRow as SubstitutionParameterPriorRowImpl,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_beta_probability_substitution_parameter_prior as build_beta_probability_substitution_parameter_prior_impl,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_dirichlet_simplex_substitution_parameter_prior as build_dirichlet_simplex_substitution_parameter_prior_impl,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_exponential_positive_substitution_parameter_prior as build_exponential_positive_substitution_parameter_prior_impl,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_fixed_positive_substitution_parameter_prior as build_fixed_positive_substitution_parameter_prior_impl,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_fixed_probability_substitution_parameter_prior as build_fixed_probability_substitution_parameter_prior_impl,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_fixed_simplex_substitution_parameter_prior as build_fixed_simplex_substitution_parameter_prior_impl,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_gamma_positive_substitution_parameter_prior as build_gamma_positive_substitution_parameter_prior_impl,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_lognormal_positive_substitution_parameter_prior as build_lognormal_positive_substitution_parameter_prior_impl,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_substitution_parameter_prior_bundle as build_substitution_parameter_prior_bundle_impl,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    evaluate_substitution_parameter_log_prior as evaluate_substitution_parameter_log_prior_impl,
)


def test_bayesian_exports_substitution_parameter_prior_surface() -> None:
    assert (
        POSITIVE_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES
        == POSITIVE_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES_IMPL
    )
    assert (
        PROBABILITY_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES
        == PROBABILITY_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES_IMPL
    )
    assert (
        SIMPLEX_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES
        == SIMPLEX_SUBSTITUTION_PARAMETER_PRIOR_FAMILIES_IMPL
    )
    assert (
        SUBSTITUTION_PARAMETER_PRIOR_TARGETS
        == SUBSTITUTION_PARAMETER_PRIOR_TARGETS_IMPL
    )
    assert (
        PositiveSubstitutionParameterPriorModel
        is PositiveSubstitutionParameterPriorModelImpl
    )
    assert (
        ProbabilitySubstitutionParameterPriorModel
        is ProbabilitySubstitutionParameterPriorModelImpl
    )
    assert (
        SimplexSubstitutionParameterPriorModel
        is SimplexSubstitutionParameterPriorModelImpl
    )
    assert SubstitutionParameterPriorBundle is SubstitutionParameterPriorBundleImpl
    assert (
        SubstitutionParameterPriorEvaluationReport
        is SubstitutionParameterPriorEvaluationReportImpl
    )
    assert SubstitutionParameterPriorRow is SubstitutionParameterPriorRowImpl
    assert (
        build_beta_probability_substitution_parameter_prior
        is build_beta_probability_substitution_parameter_prior_impl
    )
    assert (
        build_dirichlet_simplex_substitution_parameter_prior
        is build_dirichlet_simplex_substitution_parameter_prior_impl
    )
    assert (
        build_exponential_positive_substitution_parameter_prior
        is build_exponential_positive_substitution_parameter_prior_impl
    )
    assert (
        build_fixed_positive_substitution_parameter_prior
        is build_fixed_positive_substitution_parameter_prior_impl
    )
    assert (
        build_fixed_probability_substitution_parameter_prior
        is build_fixed_probability_substitution_parameter_prior_impl
    )
    assert (
        build_fixed_simplex_substitution_parameter_prior
        is build_fixed_simplex_substitution_parameter_prior_impl
    )
    assert (
        build_gamma_positive_substitution_parameter_prior
        is build_gamma_positive_substitution_parameter_prior_impl
    )
    assert (
        build_lognormal_positive_substitution_parameter_prior
        is build_lognormal_positive_substitution_parameter_prior_impl
    )
    assert (
        build_substitution_parameter_prior_bundle
        is build_substitution_parameter_prior_bundle_impl
    )
    assert (
        evaluate_substitution_parameter_log_prior
        is evaluate_substitution_parameter_log_prior_impl
    )
