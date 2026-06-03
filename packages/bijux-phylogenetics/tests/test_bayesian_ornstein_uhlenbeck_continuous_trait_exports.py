from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    ORNSTEIN_UHLENBECK_CONTINUOUS_TRAIT_MODELS,
    OrnsteinUhlenbeckContinuousTraitIdentifiabilityWarning,
    OrnsteinUhlenbeckContinuousTraitModelDefinition,
    OrnsteinUhlenbeckContinuousTraitParameterSummary,
    OrnsteinUhlenbeckContinuousTraitPosteriorRow,
    OrnsteinUhlenbeckContinuousTraitProposalSchedule,
    OrnsteinUhlenbeckContinuousTraitRunReport,
    build_ornstein_uhlenbeck_continuous_trait_model_definition,
    build_ornstein_uhlenbeck_continuous_trait_proposal_schedule,
    run_ornstein_uhlenbeck_continuous_trait_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.ornstein_uhlenbeck_continuous_trait import (
    ORNSTEIN_UHLENBECK_CONTINUOUS_TRAIT_MODELS as ORNSTEIN_UHLENBECK_CONTINUOUS_TRAIT_MODELS_IMPL,
)
from bijux_phylogenetics.bayesian.ornstein_uhlenbeck_continuous_trait import (
    OrnsteinUhlenbeckContinuousTraitIdentifiabilityWarning as OrnsteinUhlenbeckContinuousTraitIdentifiabilityWarningImpl,
)
from bijux_phylogenetics.bayesian.ornstein_uhlenbeck_continuous_trait import (
    OrnsteinUhlenbeckContinuousTraitModelDefinition as OrnsteinUhlenbeckContinuousTraitModelDefinitionImpl,
)
from bijux_phylogenetics.bayesian.ornstein_uhlenbeck_continuous_trait import (
    OrnsteinUhlenbeckContinuousTraitParameterSummary as OrnsteinUhlenbeckContinuousTraitParameterSummaryImpl,
)
from bijux_phylogenetics.bayesian.ornstein_uhlenbeck_continuous_trait import (
    OrnsteinUhlenbeckContinuousTraitPosteriorRow as OrnsteinUhlenbeckContinuousTraitPosteriorRowImpl,
)
from bijux_phylogenetics.bayesian.ornstein_uhlenbeck_continuous_trait import (
    OrnsteinUhlenbeckContinuousTraitProposalSchedule as OrnsteinUhlenbeckContinuousTraitProposalScheduleImpl,
)
from bijux_phylogenetics.bayesian.ornstein_uhlenbeck_continuous_trait import (
    OrnsteinUhlenbeckContinuousTraitRunReport as OrnsteinUhlenbeckContinuousTraitRunReportImpl,
)
from bijux_phylogenetics.bayesian.ornstein_uhlenbeck_continuous_trait import (
    build_ornstein_uhlenbeck_continuous_trait_model_definition as build_ornstein_uhlenbeck_continuous_trait_model_definition_impl,
)
from bijux_phylogenetics.bayesian.ornstein_uhlenbeck_continuous_trait import (
    build_ornstein_uhlenbeck_continuous_trait_proposal_schedule as build_ornstein_uhlenbeck_continuous_trait_proposal_schedule_impl,
)
from bijux_phylogenetics.bayesian.ornstein_uhlenbeck_continuous_trait import (
    run_ornstein_uhlenbeck_continuous_trait_metropolis_hastings as run_ornstein_uhlenbeck_continuous_trait_metropolis_hastings_impl,
)


def test_bayesian_exports_ou_continuous_trait_surface() -> None:
    assert (
        ORNSTEIN_UHLENBECK_CONTINUOUS_TRAIT_MODELS
        == ORNSTEIN_UHLENBECK_CONTINUOUS_TRAIT_MODELS_IMPL
    )
    assert OrnsteinUhlenbeckContinuousTraitIdentifiabilityWarning is (
        OrnsteinUhlenbeckContinuousTraitIdentifiabilityWarningImpl
    )
    assert OrnsteinUhlenbeckContinuousTraitModelDefinition is (
        OrnsteinUhlenbeckContinuousTraitModelDefinitionImpl
    )
    assert OrnsteinUhlenbeckContinuousTraitParameterSummary is (
        OrnsteinUhlenbeckContinuousTraitParameterSummaryImpl
    )
    assert OrnsteinUhlenbeckContinuousTraitPosteriorRow is (
        OrnsteinUhlenbeckContinuousTraitPosteriorRowImpl
    )
    assert OrnsteinUhlenbeckContinuousTraitProposalSchedule is (
        OrnsteinUhlenbeckContinuousTraitProposalScheduleImpl
    )
    assert OrnsteinUhlenbeckContinuousTraitRunReport is (
        OrnsteinUhlenbeckContinuousTraitRunReportImpl
    )
    assert (
        build_ornstein_uhlenbeck_continuous_trait_model_definition
        is build_ornstein_uhlenbeck_continuous_trait_model_definition_impl
    )
    assert (
        build_ornstein_uhlenbeck_continuous_trait_proposal_schedule
        is build_ornstein_uhlenbeck_continuous_trait_proposal_schedule_impl
    )
    assert (
        run_ornstein_uhlenbeck_continuous_trait_metropolis_hastings
        is run_ornstein_uhlenbeck_continuous_trait_metropolis_hastings_impl
    )
