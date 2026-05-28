from __future__ import annotations

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    ParsimonyConsensusSummary,
    ParsimonyEqualBestConsensusReport,
    ParsimonyEqualBestTree,
)


def test_package_parsimony_gateway_exports_equal_best_consensus_contracts() -> None:
    assert parsimony_api.ParsimonyEqualBestTree is ParsimonyEqualBestTree
    assert (
        parsimony_api.ParsimonyEqualBestConsensusReport
        is ParsimonyEqualBestConsensusReport
    )
    assert parsimony_api.ParsimonyConsensusSummary is ParsimonyConsensusSummary
