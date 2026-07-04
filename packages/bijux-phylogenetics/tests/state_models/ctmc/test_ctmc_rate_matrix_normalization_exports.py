from __future__ import annotations

from bijux_phylogenetics.phylo.likelihood import (
    compute_ctmc_expected_substitution_rate,
    normalize_ctmc_rate_matrix_by_expected_substitution_rate,
)
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    compute_ctmc_expected_substitution_rate as compute_ctmc_expected_substitution_rate_impl,
)
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    normalize_ctmc_rate_matrix_by_expected_substitution_rate as normalize_ctmc_rate_matrix_by_expected_substitution_rate_impl,
)


def test_phylo_likelihood_exports_ctmc_rate_matrix_normalization_surface() -> None:
    assert (
        compute_ctmc_expected_substitution_rate
        is compute_ctmc_expected_substitution_rate_impl
    )
    assert (
        normalize_ctmc_rate_matrix_by_expected_substitution_rate
        is normalize_ctmc_rate_matrix_by_expected_substitution_rate_impl
    )
