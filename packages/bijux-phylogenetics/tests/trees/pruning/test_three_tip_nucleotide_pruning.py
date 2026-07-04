from __future__ import annotations

from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import (
    compile_finite_state_likelihood_context,
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_finite_state import (
    build_selected_nucleotide_finite_state_surface,
)
from tests.support.three_tip_pruning_reference import (
    build_three_tip_pruning_reference,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    return FIXTURES / "trees" / name


@pytest.mark.parametrize(
    ("model_name", "surface_kwargs", "expected_root_prior"),
    [
        ("jc69", {}, numpy.full(4, 0.25, dtype=float)),
        ("k80", {"kappa": 4.0}, numpy.full(4, 0.25, dtype=float)),
        (
            "f81",
            {"base_frequencies": {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}},
            numpy.array([0.4, 0.1, 0.2, 0.3], dtype=float),
        ),
        (
            "hky85",
            {
                "kappa": 4.0,
                "base_frequencies": {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
            },
            numpy.array([0.4, 0.1, 0.2, 0.3], dtype=float),
        ),
        (
            "gtr",
            {
                "base_frequencies": {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
                "exchangeabilities": {
                    "AC": 1.0,
                    "AG": 4.5,
                    "AT": 0.8,
                    "CG": 1.6,
                    "CT": 2.4,
                    "GT": 3.1,
                },
            },
            numpy.array([0.4, 0.1, 0.2, 0.3], dtype=float),
        ),
    ],
)
def test_selected_nucleotide_surface_matches_three_tip_reference(
    model_name: str,
    surface_kwargs: dict[str, object],
    expected_root_prior: numpy.ndarray,
) -> None:
    states_by_taxon = {"A": "A", "B": "C", "C": "G"}
    surface = build_selected_nucleotide_finite_state_surface(
        _records_from_states(states_by_taxon),
        model_name=model_name,
        owner_name="three-tip pruning fixture",
        **surface_kwargs,
    )

    assert numpy.allclose(surface.root_prior, expected_root_prior, rtol=0.0, atol=1e-12)
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
