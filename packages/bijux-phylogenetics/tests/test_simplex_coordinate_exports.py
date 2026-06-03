from __future__ import annotations

from bijux_phylogenetics.phylo.likelihood import (
    DNA_EXCHANGEABILITY_LABELS,
    SimplexCoordinateParameterization,
    evaluate_gtr_tree_likelihood_from_unconstrained_exchangeabilities,
    evaluate_gtr_tree_likelihood_from_unconstrained_exchangeabilities_from_alignment,
    parameterize_dna_base_frequency_simplex,
    parameterize_dna_exchangeability_simplex,
    parameterize_named_simplex,
    resolve_anchor_normalized_dna_exchangeabilities_from_unconstrained,
    resolve_dna_base_frequencies_from_unconstrained,
    resolve_dna_base_frequency_simplex_from_unconstrained,
    resolve_dna_exchangeability_simplex_from_unconstrained,
    resolve_named_simplex_from_unconstrained,
)
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    DNA_EXCHANGEABILITY_LABELS as dna_exchangeability_labels_impl,
)
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    parameterize_dna_base_frequency_simplex as parameterize_dna_base_frequency_simplex_impl,
)
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    parameterize_dna_exchangeability_simplex as parameterize_dna_exchangeability_simplex_impl,
)
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    resolve_anchor_normalized_dna_exchangeabilities_from_unconstrained as resolve_anchor_normalized_dna_exchangeabilities_from_unconstrained_impl,
)
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    resolve_dna_base_frequencies_from_unconstrained as resolve_dna_base_frequencies_from_unconstrained_impl,
)
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    resolve_dna_base_frequency_simplex_from_unconstrained as resolve_dna_base_frequency_simplex_from_unconstrained_impl,
)
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    resolve_dna_exchangeability_simplex_from_unconstrained as resolve_dna_exchangeability_simplex_from_unconstrained_impl,
)
from bijux_phylogenetics.phylo.likelihood.gtr import (
    evaluate_gtr_tree_likelihood_from_unconstrained_exchangeabilities as evaluate_gtr_tree_likelihood_from_unconstrained_exchangeabilities_impl,
)
from bijux_phylogenetics.phylo.likelihood.gtr import (
    evaluate_gtr_tree_likelihood_from_unconstrained_exchangeabilities_from_alignment as evaluate_gtr_tree_likelihood_from_unconstrained_exchangeabilities_from_alignment_impl,
)
from bijux_phylogenetics.phylo.likelihood.simplex_coordinates import (
    SimplexCoordinateParameterization as simplex_coordinate_parameterization_impl,
)
from bijux_phylogenetics.phylo.likelihood.simplex_coordinates import (
    parameterize_named_simplex as parameterize_named_simplex_impl,
)
from bijux_phylogenetics.phylo.likelihood.simplex_coordinates import (
    resolve_named_simplex_from_unconstrained as resolve_named_simplex_from_unconstrained_impl,
)


def test_phylo_likelihood_exports_simplex_coordinate_surface() -> None:
    assert DNA_EXCHANGEABILITY_LABELS is dna_exchangeability_labels_impl
    assert SimplexCoordinateParameterization is simplex_coordinate_parameterization_impl
    assert parameterize_named_simplex is parameterize_named_simplex_impl
    assert (
        parameterize_dna_base_frequency_simplex
        is parameterize_dna_base_frequency_simplex_impl
    )
    assert (
        parameterize_dna_exchangeability_simplex
        is parameterize_dna_exchangeability_simplex_impl
    )
    assert (
        resolve_named_simplex_from_unconstrained
        is resolve_named_simplex_from_unconstrained_impl
    )
    assert (
        resolve_dna_base_frequency_simplex_from_unconstrained
        is resolve_dna_base_frequency_simplex_from_unconstrained_impl
    )
    assert (
        resolve_dna_base_frequencies_from_unconstrained
        is resolve_dna_base_frequencies_from_unconstrained_impl
    )
    assert (
        resolve_dna_exchangeability_simplex_from_unconstrained
        is resolve_dna_exchangeability_simplex_from_unconstrained_impl
    )
    assert (
        resolve_anchor_normalized_dna_exchangeabilities_from_unconstrained
        is resolve_anchor_normalized_dna_exchangeabilities_from_unconstrained_impl
    )
    assert (
        evaluate_gtr_tree_likelihood_from_unconstrained_exchangeabilities
        is evaluate_gtr_tree_likelihood_from_unconstrained_exchangeabilities_impl
    )
    assert (
        evaluate_gtr_tree_likelihood_from_unconstrained_exchangeabilities_from_alignment
        is evaluate_gtr_tree_likelihood_from_unconstrained_exchangeabilities_from_alignment_impl
    )
