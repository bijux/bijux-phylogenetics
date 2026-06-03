from __future__ import annotations

from bijux_phylogenetics.phylo.likelihood import (
    CodonCtmcTreeLikelihoodReport,
    CodonStateSpace,
    build_equal_rate_codon_ctmc_rate_matrix,
    evaluate_codon_ctmc_tree_likelihood,
    evaluate_codon_ctmc_tree_likelihood_from_alignment,
    resolve_codon_state_space,
    validate_codon_frequency_vector,
)
from bijux_phylogenetics.phylo.likelihood.codon import (
    evaluate_codon_ctmc_tree_likelihood as evaluate_codon_ctmc_tree_likelihood_impl,
)
from bijux_phylogenetics.phylo.likelihood.codon import (
    evaluate_codon_ctmc_tree_likelihood_from_alignment as evaluate_codon_ctmc_tree_likelihood_from_alignment_impl,
)
from bijux_phylogenetics.phylo.likelihood.codon_states import (
    CodonStateSpace as CodonStateSpaceImpl,
)
from bijux_phylogenetics.phylo.likelihood.codon_states import (
    build_equal_rate_codon_ctmc_rate_matrix as build_equal_rate_codon_ctmc_rate_matrix_impl,
)
from bijux_phylogenetics.phylo.likelihood.codon_states import (
    resolve_codon_state_space as resolve_codon_state_space_impl,
)
from bijux_phylogenetics.phylo.likelihood.codon_states import (
    validate_codon_frequency_vector as validate_codon_frequency_vector_impl,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    CodonCtmcTreeLikelihoodReport as CodonCtmcTreeLikelihoodReportImpl,
)


def test_phylo_likelihood_exports_codon_ctmc_surface() -> None:
    assert CodonStateSpace is CodonStateSpaceImpl
    assert CodonCtmcTreeLikelihoodReport is CodonCtmcTreeLikelihoodReportImpl
    assert resolve_codon_state_space is resolve_codon_state_space_impl
    assert validate_codon_frequency_vector is validate_codon_frequency_vector_impl
    assert (
        build_equal_rate_codon_ctmc_rate_matrix
        is build_equal_rate_codon_ctmc_rate_matrix_impl
    )
    assert (
        evaluate_codon_ctmc_tree_likelihood is evaluate_codon_ctmc_tree_likelihood_impl
    )
    assert (
        evaluate_codon_ctmc_tree_likelihood_from_alignment
        is evaluate_codon_ctmc_tree_likelihood_from_alignment_impl
    )
