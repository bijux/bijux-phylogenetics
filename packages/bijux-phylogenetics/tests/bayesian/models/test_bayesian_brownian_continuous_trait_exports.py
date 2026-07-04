from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    BROWNIAN_CONTINUOUS_TRAIT_MODELS,
    BrownianContinuousTraitModelDefinition,
    BrownianContinuousTraitParameterSummary,
    BrownianContinuousTraitPosteriorRow,
    BrownianContinuousTraitProposalSchedule,
    BrownianContinuousTraitRunReport,
    build_brownian_continuous_trait_model_definition,
    build_brownian_continuous_trait_proposal_schedule,
    run_brownian_continuous_trait_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.brownian_continuous_trait import (
    BROWNIAN_CONTINUOUS_TRAIT_MODELS as BROWNIAN_CONTINUOUS_TRAIT_MODELS_IMPL,
)
from bijux_phylogenetics.bayesian.brownian_continuous_trait import (
    BrownianContinuousTraitModelDefinition as BrownianContinuousTraitModelDefinitionImpl,
)
from bijux_phylogenetics.bayesian.brownian_continuous_trait import (
    BrownianContinuousTraitParameterSummary as BrownianContinuousTraitParameterSummaryImpl,
)
from bijux_phylogenetics.bayesian.brownian_continuous_trait import (
    BrownianContinuousTraitPosteriorRow as BrownianContinuousTraitPosteriorRowImpl,
)
from bijux_phylogenetics.bayesian.brownian_continuous_trait import (
    BrownianContinuousTraitProposalSchedule as BrownianContinuousTraitProposalScheduleImpl,
)
from bijux_phylogenetics.bayesian.brownian_continuous_trait import (
    BrownianContinuousTraitRunReport as BrownianContinuousTraitRunReportImpl,
)
from bijux_phylogenetics.bayesian.brownian_continuous_trait import (
    build_brownian_continuous_trait_model_definition as build_brownian_continuous_trait_model_definition_impl,
)
from bijux_phylogenetics.bayesian.brownian_continuous_trait import (
    build_brownian_continuous_trait_proposal_schedule as build_brownian_continuous_trait_proposal_schedule_impl,
)
from bijux_phylogenetics.bayesian.brownian_continuous_trait import (
    run_brownian_continuous_trait_metropolis_hastings as run_brownian_continuous_trait_metropolis_hastings_impl,
)


def test_bayesian_exports_brownian_continuous_trait_surface() -> None:
    assert BROWNIAN_CONTINUOUS_TRAIT_MODELS == BROWNIAN_CONTINUOUS_TRAIT_MODELS_IMPL
    assert BrownianContinuousTraitModelDefinition is (
        BrownianContinuousTraitModelDefinitionImpl
    )
    assert BrownianContinuousTraitParameterSummary is (
        BrownianContinuousTraitParameterSummaryImpl
    )
    assert BrownianContinuousTraitPosteriorRow is (
        BrownianContinuousTraitPosteriorRowImpl
    )
    assert BrownianContinuousTraitProposalSchedule is (
        BrownianContinuousTraitProposalScheduleImpl
    )
    assert BrownianContinuousTraitRunReport is BrownianContinuousTraitRunReportImpl
    assert (
        build_brownian_continuous_trait_model_definition
        is build_brownian_continuous_trait_model_definition_impl
    )
    assert (
        build_brownian_continuous_trait_proposal_schedule
        is build_brownian_continuous_trait_proposal_schedule_impl
    )
    assert (
        run_brownian_continuous_trait_metropolis_hastings
        is run_brownian_continuous_trait_metropolis_hastings_impl
    )
