from __future__ import annotations

from bijux_phylogenetics.phylo.likelihood import (
    ValidatedCtmcRateMatrix,
    validate_ctmc_rate_matrix,
)
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    ValidatedCtmcRateMatrix as ValidatedCtmcRateMatrixImpl,
)
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    validate_ctmc_rate_matrix as validate_ctmc_rate_matrix_impl,
)


def test_phylo_likelihood_exports_ctmc_validation_surface() -> None:
    assert ValidatedCtmcRateMatrix is ValidatedCtmcRateMatrixImpl
    assert validate_ctmc_rate_matrix is validate_ctmc_rate_matrix_impl
