from __future__ import annotations

from bijux_phylogenetics.phylo.likelihood import (
    SolvedCtmcStationaryDistribution,
    solve_ctmc_stationary_distribution,
    verify_ctmc_stationary_distribution,
)
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    SolvedCtmcStationaryDistribution as SolvedCtmcStationaryDistributionImpl,
)
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    solve_ctmc_stationary_distribution as solve_ctmc_stationary_distribution_impl,
)
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    verify_ctmc_stationary_distribution as verify_ctmc_stationary_distribution_impl,
)


def test_phylo_likelihood_exports_ctmc_stationary_distribution_surface() -> None:
    assert SolvedCtmcStationaryDistribution is SolvedCtmcStationaryDistributionImpl
    assert solve_ctmc_stationary_distribution is solve_ctmc_stationary_distribution_impl
    assert (
        verify_ctmc_stationary_distribution is verify_ctmc_stationary_distribution_impl
    )
