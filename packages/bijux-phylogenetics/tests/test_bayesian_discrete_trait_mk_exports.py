from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    DISCRETE_TRAIT_MK_MODELS,
    DISCRETE_TRAIT_MK_ROOT_PRIOR_MODES,
    DiscreteTraitMkModelDefinition,
    DiscreteTraitMkNodeStateSummary,
    DiscreteTraitMkPosteriorRow,
    DiscreteTraitMkProposalSchedule,
    DiscreteTraitMkRunReport,
    build_discrete_trait_mk_model_definition,
    build_discrete_trait_mk_proposal_schedule,
    run_discrete_trait_mk_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.discrete_trait_mk import (
    DISCRETE_TRAIT_MK_MODELS as DISCRETE_TRAIT_MK_MODELS_IMPL,
)
from bijux_phylogenetics.bayesian.discrete_trait_mk import (
    DISCRETE_TRAIT_MK_ROOT_PRIOR_MODES as DISCRETE_TRAIT_MK_ROOT_PRIOR_MODES_IMPL,
)
from bijux_phylogenetics.bayesian.discrete_trait_mk import (
    DiscreteTraitMkModelDefinition as DiscreteTraitMkModelDefinitionImpl,
)
from bijux_phylogenetics.bayesian.discrete_trait_mk import (
    DiscreteTraitMkNodeStateSummary as DiscreteTraitMkNodeStateSummaryImpl,
)
from bijux_phylogenetics.bayesian.discrete_trait_mk import (
    DiscreteTraitMkPosteriorRow as DiscreteTraitMkPosteriorRowImpl,
)
from bijux_phylogenetics.bayesian.discrete_trait_mk import (
    DiscreteTraitMkProposalSchedule as DiscreteTraitMkProposalScheduleImpl,
)
from bijux_phylogenetics.bayesian.discrete_trait_mk import (
    DiscreteTraitMkRunReport as DiscreteTraitMkRunReportImpl,
)
from bijux_phylogenetics.bayesian.discrete_trait_mk import (
    build_discrete_trait_mk_model_definition as build_discrete_trait_mk_model_definition_impl,
)
from bijux_phylogenetics.bayesian.discrete_trait_mk import (
    build_discrete_trait_mk_proposal_schedule as build_discrete_trait_mk_proposal_schedule_impl,
)
from bijux_phylogenetics.bayesian.discrete_trait_mk import (
    run_discrete_trait_mk_metropolis_hastings as run_discrete_trait_mk_metropolis_hastings_impl,
)


def test_bayesian_exports_discrete_trait_mk_surface() -> None:
    assert DISCRETE_TRAIT_MK_MODELS == DISCRETE_TRAIT_MK_MODELS_IMPL
    assert DISCRETE_TRAIT_MK_ROOT_PRIOR_MODES == DISCRETE_TRAIT_MK_ROOT_PRIOR_MODES_IMPL
    assert DiscreteTraitMkModelDefinition is DiscreteTraitMkModelDefinitionImpl
    assert DiscreteTraitMkNodeStateSummary is DiscreteTraitMkNodeStateSummaryImpl
    assert DiscreteTraitMkPosteriorRow is DiscreteTraitMkPosteriorRowImpl
    assert DiscreteTraitMkProposalSchedule is DiscreteTraitMkProposalScheduleImpl
    assert DiscreteTraitMkRunReport is DiscreteTraitMkRunReportImpl
    assert (
        build_discrete_trait_mk_model_definition
        is build_discrete_trait_mk_model_definition_impl
    )
    assert (
        build_discrete_trait_mk_proposal_schedule
        is build_discrete_trait_mk_proposal_schedule_impl
    )
    assert (
        run_discrete_trait_mk_metropolis_hastings
        is run_discrete_trait_mk_metropolis_hastings_impl
    )
