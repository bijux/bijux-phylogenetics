from __future__ import annotations

from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood.codon_finite_state import (
    build_codon_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.discrete_finite_state import (
    build_discrete_mk_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    clear_shared_finite_state_transition_cache,
    compile_finite_state_likelihood_context,
    shared_finite_state_transition_cache_entry_count,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_finite_state import (
    build_selected_nucleotide_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.protein_finite_state import (
    build_empirical_protein_finite_state_surface,
    build_protein_poisson_finite_state_surface,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_selected_nucleotide_surface_reuses_shared_transition_cache_for_identical_model_state() -> (
    None
):
    clear_shared_finite_state_transition_cache()
    tree = load_tree(fixture("trees", "hky85_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta")
    )
    surface = build_selected_nucleotide_finite_state_surface(
        records,
        model_name="hky85",
        owner_name="family transition cache state",
        kappa=4.0,
        base_frequencies={"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
    )

    first_context = compile_finite_state_likelihood_context(tree, surface)
    cache_entry_count_after_first = shared_finite_state_transition_cache_entry_count()
    second_context = compile_finite_state_likelihood_context(tree, surface)

    assert first_context.transition_model_state is not None
    assert second_context.transition_model_state == first_context.transition_model_state
    assert cache_entry_count_after_first > 0
    assert (
        shared_finite_state_transition_cache_entry_count()
        == cache_entry_count_after_first
    )


def test_empirical_protein_surface_distinguishes_rate_matrix_state_in_shared_cache() -> (
    None
):
    clear_shared_finite_state_transition_cache()
    tree = load_tree(fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta")
    )

    first_context = compile_finite_state_likelihood_context(
        tree,
        build_empirical_protein_finite_state_surface(
            records,
            rate_matrix=_empirical_rate_matrix(
                boosted_pairs={("A", "C"): 0.45, ("D", "E"): 0.55},
            ),
            root_prior=_biased_protein_root_prior(),
            matrix_label="polar-a",
        ),
    )
    cache_entry_count_after_first = shared_finite_state_transition_cache_entry_count()
    second_context = compile_finite_state_likelihood_context(
        tree,
        build_empirical_protein_finite_state_surface(
            records,
            rate_matrix=_empirical_rate_matrix(
                boosted_pairs={("A", "C"): 0.20, ("D", "E"): 0.80},
            ),
            root_prior=_biased_protein_root_prior(),
            matrix_label="polar-b",
        ),
    )

    assert first_context.transition_model_state is not None
    assert second_context.transition_model_state is not None
    assert first_context.transition_model_state != second_context.transition_model_state
    assert shared_finite_state_transition_cache_entry_count() > cache_entry_count_after_first


def test_codon_and_protein_surfaces_keep_shared_transition_model_state_separate() -> (
    None
):
    clear_shared_finite_state_transition_cache()
    protein_tree = load_tree(
        fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk")
    )
    protein_records = load_fasta_alignment(
        fixture("alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta")
    )
    codon_tree = load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk"))
    codon_records = load_fasta_alignment(
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta")
    )

    protein_context = compile_finite_state_likelihood_context(
        protein_tree,
        build_protein_poisson_finite_state_surface(protein_records),
    )
    codon_context = compile_finite_state_likelihood_context(
        codon_tree,
        build_codon_finite_state_surface(codon_records),
    )

    assert protein_context.transition_model_state is not None
    assert codon_context.transition_model_state is not None
    assert protein_context.transition_model_state != codon_context.transition_model_state
    assert protein_context.transition_cache_key_by_node_id
    assert codon_context.transition_cache_key_by_node_id


def test_discrete_surface_root_prior_changes_shared_transition_model_state() -> None:
    clear_shared_finite_state_transition_cache()
    tree = load_tree(fixture("trees", "felsenstein_two_tip_tree.nwk"))
    rate_matrix = numpy.array(
        [
            [-1.0, 0.6, 0.4],
            [0.3, -0.8, 0.5],
            [0.2, 0.7, -0.9],
        ],
        dtype=float,
    )

    first_context = compile_finite_state_likelihood_context(
        tree,
        build_discrete_mk_finite_state_surface(
            state_order=["north", "south", "west"],
            rate_matrix=rate_matrix,
            root_prior=numpy.array([0.6, 0.2, 0.2], dtype=float),
            model_name="unordered discrete Mk",
        ),
    )
    second_context = compile_finite_state_likelihood_context(
        tree,
        build_discrete_mk_finite_state_surface(
            state_order=["north", "south", "west"],
            rate_matrix=rate_matrix,
            root_prior=numpy.array([0.2, 0.6, 0.2], dtype=float),
            model_name="unordered discrete Mk",
        ),
    )

    assert first_context.transition_model_state is not None
    assert second_context.transition_model_state is not None
    assert first_context.transition_model_state != second_context.transition_model_state


def _biased_protein_root_prior() -> numpy.ndarray:
    prior = numpy.full(20, 0.02, dtype=float)
    prior[0] = 0.18
    prior[1] = 0.10
    prior[2] = 0.14
    prior[3] = 0.12
    prior[4] = 0.06
    return prior / float(prior.sum())


def _empirical_rate_matrix(
    *,
    boosted_pairs: dict[tuple[str, str], float],
) -> numpy.ndarray:
    state_order = (
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    )
    state_index = {state: index for index, state in enumerate(state_order)}
    rate_matrix = numpy.full((len(state_order), len(state_order)), 0.02, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
    for (left_state, right_state), rate in boosted_pairs.items():
        left_index = state_index[left_state]
        right_index = state_index[right_state]
        rate_matrix[left_index, right_index] = rate
        rate_matrix[right_index, left_index] = rate
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    return rate_matrix
