from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    CheckpointedMetropolisHastingsRunReport,
    MetropolisHastingsCheckpoint,
    MetropolisHastingsProposal,
    MetropolisHastingsRandomState,
    MetropolisHastingsRunReport,
    MetropolisHastingsStepRow,
    build_metropolis_hastings_checkpoint,
    build_metropolis_hastings_proposal,
    build_metropolis_hastings_random_state,
    deserialize_metropolis_hastings_checkpoint,
    deserialize_metropolis_hastings_checkpoint_json,
    resume_metropolis_hastings_sampler,
    run_checkpointed_metropolis_hastings_sampler,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
    serialize_metropolis_hastings_checkpoint,
    serialize_metropolis_hastings_checkpoint_json,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    CheckpointedMetropolisHastingsRunReport as CheckpointedMetropolisHastingsRunReportImpl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsCheckpoint as MetropolisHastingsCheckpointImpl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsProposal as MetropolisHastingsProposalImpl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsRandomState as MetropolisHastingsRandomStateImpl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsRunReport as MetropolisHastingsRunReportImpl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsStepRow as MetropolisHastingsStepRowImpl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    build_metropolis_hastings_checkpoint as build_metropolis_hastings_checkpoint_impl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    build_metropolis_hastings_proposal as build_metropolis_hastings_proposal_impl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    build_metropolis_hastings_random_state as build_metropolis_hastings_random_state_impl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    deserialize_metropolis_hastings_checkpoint as deserialize_metropolis_hastings_checkpoint_impl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    deserialize_metropolis_hastings_checkpoint_json as deserialize_metropolis_hastings_checkpoint_json_impl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    resume_metropolis_hastings_sampler as resume_metropolis_hastings_sampler_impl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    run_checkpointed_metropolis_hastings_sampler as run_checkpointed_metropolis_hastings_sampler_impl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    run_metropolis_hastings_sampler as run_metropolis_hastings_sampler_impl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    score_bayesian_phylogenetic_state as score_bayesian_phylogenetic_state_impl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    serialize_metropolis_hastings_checkpoint as serialize_metropolis_hastings_checkpoint_impl,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    serialize_metropolis_hastings_checkpoint_json as serialize_metropolis_hastings_checkpoint_json_impl,
)


def test_bayesian_exports_metropolis_hastings_surface() -> None:
    assert CheckpointedMetropolisHastingsRunReport is (
        CheckpointedMetropolisHastingsRunReportImpl
    )
    assert MetropolisHastingsCheckpoint is MetropolisHastingsCheckpointImpl
    assert MetropolisHastingsProposal is MetropolisHastingsProposalImpl
    assert MetropolisHastingsRandomState is MetropolisHastingsRandomStateImpl
    assert MetropolisHastingsStepRow is MetropolisHastingsStepRowImpl
    assert MetropolisHastingsRunReport is MetropolisHastingsRunReportImpl
    assert (
        build_metropolis_hastings_checkpoint
        is build_metropolis_hastings_checkpoint_impl
    )
    assert build_metropolis_hastings_proposal is build_metropolis_hastings_proposal_impl
    assert (
        build_metropolis_hastings_random_state
        is build_metropolis_hastings_random_state_impl
    )
    assert (
        deserialize_metropolis_hastings_checkpoint
        is deserialize_metropolis_hastings_checkpoint_impl
    )
    assert (
        deserialize_metropolis_hastings_checkpoint_json
        is deserialize_metropolis_hastings_checkpoint_json_impl
    )
    assert resume_metropolis_hastings_sampler is resume_metropolis_hastings_sampler_impl
    assert (
        run_checkpointed_metropolis_hastings_sampler
        is run_checkpointed_metropolis_hastings_sampler_impl
    )
    assert run_metropolis_hastings_sampler is run_metropolis_hastings_sampler_impl
    assert score_bayesian_phylogenetic_state is score_bayesian_phylogenetic_state_impl
    assert (
        serialize_metropolis_hastings_checkpoint
        is serialize_metropolis_hastings_checkpoint_impl
    )
    assert (
        serialize_metropolis_hastings_checkpoint_json
        is serialize_metropolis_hastings_checkpoint_json_impl
    )
