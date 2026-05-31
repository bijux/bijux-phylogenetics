from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood.search_replay import (
    replay_nucleotide_likelihood_nni_search_trace,
    replay_nucleotide_likelihood_spr_search_trace,
    replay_nucleotide_likelihood_tbr_search_trace,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodSearchTraceReplayReport,
    NucleotideLikelihoodSearchTraceReplayStep,
)


def test_public_likelihood_exports_search_trace_replay_surface() -> None:
    assert (
        likelihood_api.NucleotideLikelihoodSearchTraceReplayReport
        is NucleotideLikelihoodSearchTraceReplayReport
    )
    assert (
        likelihood_api.NucleotideLikelihoodSearchTraceReplayStep
        is NucleotideLikelihoodSearchTraceReplayStep
    )
    assert (
        likelihood_api.replay_nucleotide_likelihood_nni_search_trace
        is replay_nucleotide_likelihood_nni_search_trace
    )
    assert (
        likelihood_api.replay_nucleotide_likelihood_spr_search_trace
        is replay_nucleotide_likelihood_spr_search_trace
    )
    assert (
        likelihood_api.replay_nucleotide_likelihood_tbr_search_trace
        is replay_nucleotide_likelihood_tbr_search_trace
    )
