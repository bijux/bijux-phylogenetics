from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    MetropolisHastingsProposal,
    MetropolisHastingsRunReport,
    MetropolisHastingsStepRow,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsProposal as MetropolisHastingsProposalImpl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsRunReport as MetropolisHastingsRunReportImpl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsStepRow as MetropolisHastingsStepRowImpl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    run_metropolis_hastings_sampler as run_metropolis_hastings_sampler_impl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    score_bayesian_phylogenetic_state as score_bayesian_phylogenetic_state_impl,
)


def test_bayesian_exports_metropolis_hastings_surface() -> None:
    assert MetropolisHastingsProposal is MetropolisHastingsProposalImpl
    assert MetropolisHastingsStepRow is MetropolisHastingsStepRowImpl
    assert MetropolisHastingsRunReport is MetropolisHastingsRunReportImpl
    assert run_metropolis_hastings_sampler is run_metropolis_hastings_sampler_impl
    assert (
        score_bayesian_phylogenetic_state
        is score_bayesian_phylogenetic_state_impl
    )
