from __future__ import annotations

from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import (
    build_codon_finite_state_surface,
    build_empirical_protein_finite_state_surface,
    build_protein_poisson_finite_state_surface,
    compile_finite_state_likelihood_context,
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
)
from tests.support.three_tip_pruning_reference import (
    build_three_tip_pruning_reference,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    return FIXTURES / "trees" / name


@pytest.mark.parametrize(
    (
        "family_name",
        "surface_builder",
        "states_by_taxon",
        "expected_root_prior_factory",
    ),
    [
        (
            "protein-poisson",
            lambda records: build_protein_poisson_finite_state_surface(records),
            {"A": "A", "B": "C", "C": "D"},
            lambda: numpy.full(20, 1.0 / 20.0, dtype=float),
        ),
        (
            "protein-empirical",
            lambda records: build_empirical_protein_finite_state_surface(
                records,
                rate_matrix=_compact_polar_rate_matrix(),
                root_prior=_biased_root_prior(),
                matrix_label="compact-polar",
            ),
            {"A": "A", "B": "C", "C": "D"},
            lambda: _biased_root_prior(),
        ),
        (
            "codon",
            lambda records: build_codon_finite_state_surface(records),
            {"A": "AAA", "B": "AAG", "C": "GCT"},
            lambda: None,
        ),
    ],
)
def test_protein_and_codon_surfaces_match_three_tip_reference(
    family_name: str,
    surface_builder,
    states_by_taxon: dict[str, str],
    expected_root_prior_factory,
) -> None:
    surface = surface_builder(_records_from_states(states_by_taxon))
    expected_root_prior = expected_root_prior_factory()

    if expected_root_prior is not None:
        assert numpy.allclose(
            surface.root_prior,
            expected_root_prior,
            rtol=0.0,
            atol=1e-12,
        ), family_name
    _assert_surface_matches_three_tip_reference(
        surface=surface,
        states_by_taxon=states_by_taxon,
    )


def _assert_surface_matches_three_tip_reference(
    *,
    surface,
    states_by_taxon: dict[str, str],
) -> None:
    tree = load_tree(fixture("felsenstein_three_tip_tree.nwk"))
    context = compile_finite_state_likelihood_context(tree, surface)
    leaf_vector_by_taxon = {
        taxon: surface.leaf_likelihood_vector_for_state(taxon, state)
        for taxon, state in states_by_taxon.items()
    }
    reference = build_three_tip_pruning_reference(
        tree=tree,
        root_prior=surface.root_prior,
        leaf_vector_by_taxon=leaf_vector_by_taxon,
        transition_matrix_by_node_id=context.transition_matrix_by_node_id,
    )
    pruning_pass = postorder_conditional_likelihoods(
        tree,
        state_count=surface.state_count,
        leaf_likelihood=lambda node: leaf_vector_by_taxon[node.name or ""],
        transition_matrix_for_child=lambda child: context.transition_matrix_by_node_id[
            child.node_id
        ],
    )
    log_likelihood = log_likelihood_from_root_prior(
        tree,
        pruning_pass,
        root_prior=surface.root_prior,
    )

    assert numpy.allclose(
        pruning_pass.conditional_for_node(reference.internal_child),
        reference.internal_normalized,
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.allclose(
        pruning_pass.conditional_for_node(tree.root),
        reference.root_normalized,
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.isclose(
        pruning_pass.subtree_log_scaling_for_node(reference.internal_child),
        reference.internal_log_scale,
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.isclose(
        pruning_pass.subtree_log_scaling_for_node(tree.root),
        reference.root_log_scale,
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.isclose(log_likelihood, reference.log_likelihood, rtol=0.0, atol=1e-12)


def _records_from_states(states_by_taxon: dict[str, str]) -> list[AlignmentRecord]:
    return [
        AlignmentRecord(identifier=taxon, sequence=states_by_taxon[taxon])
        for taxon in ["A", "B", "C"]
    ]


def _compact_polar_rate_matrix() -> numpy.ndarray:
    rate_matrix = numpy.full((20, 20), 0.02, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
    state_index = _protein_state_index()
    for (left_state, right_state), rate in {
        ("A", "C"): 0.45,
        ("C", "D"): 0.35,
        ("D", "E"): 0.55,
        ("A", "E"): 0.20,
    }.items():
        left_index = state_index[left_state]
        right_index = state_index[right_state]
        rate_matrix[left_index, right_index] = rate
        rate_matrix[right_index, left_index] = rate
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    return rate_matrix


def _biased_root_prior() -> numpy.ndarray:
    prior = numpy.full(20, 0.02, dtype=float)
    state_index = _protein_state_index()
    prior[state_index["A"]] = 0.18
    prior[state_index["C"]] = 0.10
    prior[state_index["D"]] = 0.14
    prior[state_index["E"]] = 0.12
    prior[state_index["F"]] = 0.06
    return prior / float(prior.sum())


def _protein_state_index() -> dict[str, int]:
    return {
        state: index
        for index, state in enumerate(
            (
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
        )
    }
