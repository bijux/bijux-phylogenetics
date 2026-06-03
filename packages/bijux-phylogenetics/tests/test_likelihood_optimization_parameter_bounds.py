from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
import bijux_phylogenetics.phylo.likelihood.empirical as empirical_likelihood
import bijux_phylogenetics.phylo.likelihood.gtr as gtr_likelihood
import bijux_phylogenetics.phylo.likelihood.jc69 as jc69_likelihood
import bijux_phylogenetics.phylo.likelihood.k80 as k80_likelihood
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_k80_kappa_optimization_rejects_initial_value_outside_bounds_before_evaluation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tree = load_tree(fixture("trees", "k80_kappa_optimization_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "k80_kappa_optimization_alignment_2_taxa.fasta")
    )

    def fail_if_called(*args: object, **kwargs: object) -> object:
        raise AssertionError("k80 likelihood evaluator should not run")

    monkeypatch.setattr(
        k80_likelihood,
        "_evaluate_k80_tree_likelihood_from_patterns",
        fail_if_called,
    )

    with pytest.raises(
        ValueError,
        match=r"K80 kappa optimization requires 'kappa' to lie within \[0.05, 20.0\]",
    ):
        k80_likelihood.optimize_k80_kappa(
            tree,
            records,
            initial_kappa=21.0,
            lower_kappa_bound=0.05,
            upper_kappa_bound=20.0,
        )


def test_gtr_exchangeability_optimization_rejects_initial_rate_outside_bounds_before_evaluation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tree = load_tree(
        fixture("trees", "gtr_exchangeability_optimization_tree_2_taxa.nwk")
    )
    records = load_fasta_alignment(
        fixture("alignments", "gtr_exchangeability_optimization_alignment_2_taxa.fasta")
    )

    def fail_if_called(*args: object, **kwargs: object) -> object:
        raise AssertionError("gtr likelihood evaluator should not run")

    monkeypatch.setattr(
        gtr_likelihood,
        "_evaluate_gtr_tree_likelihood_from_patterns",
        fail_if_called,
    )

    with pytest.raises(
        ValueError,
        match=r"GTR exchangeability optimization requires 'AG' to lie within \[0.05, 20.0\]",
    ):
        gtr_likelihood.optimize_gtr_exchangeabilities(
            tree,
            records,
            initial_exchangeabilities={
                "AC": 1.0,
                "AG": 21.0,
                "AT": 1.0,
                "CG": 1.0,
                "CT": 1.0,
                "GT": 1.0,
            },
            lower_exchangeability_bound=0.05,
            upper_exchangeability_bound=20.0,
        )


def test_jc69_branch_optimization_rejects_starting_branch_outside_bounds_before_evaluation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tree = load_tree(fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta")
    )
    for _parent, child in tree.iter_edges():
        child.branch_length = 5.5

    def fail_if_called(*args: object, **kwargs: object) -> object:
        raise AssertionError("jc69 likelihood evaluator should not run")

    monkeypatch.setattr(
        jc69_likelihood,
        "_evaluate_jc69_tree_likelihood_from_patterns",
        fail_if_called,
    )

    with pytest.raises(
        InvalidBranchLengthError,
        match=r"JC69 branch-length optimization requires 'branch length' to lie within \[1e-06, 5.0\]",
    ):
        jc69_likelihood.optimize_jc69_branch_lengths(
            tree,
            records,
            lower_branch_length_bound=1e-6,
            upper_branch_length_bound=5.0,
        )


def test_empirical_branch_optimization_rejects_invalid_alpha_before_objective(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_if_called(*args: object, **kwargs: object) -> object:
        raise AssertionError("empirical branch objective should not run")

    monkeypatch.setattr(
        empirical_likelihood,
        "_evaluate_empirical_protein_branch_optimization_objective",
        fail_if_called,
    )

    with pytest.raises(
        ValueError,
        match="discrete-gamma alpha must be a finite positive value",
    ):
        empirical_likelihood.optimize_empirical_protein_branch_lengths_from_alignment(
            fixture(
                "trees", "empirical_protein_branch_optimization_start_tree_2_taxa.nwk"
            ),
            fixture(
                "alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"
            ),
            rate_matrix=_compact_polar_rate_matrix(),
            likelihood_model="discrete-gamma",
            alpha=0.0,
            root_prior=_biased_root_prior(),
        )


def _compact_polar_rate_matrix():
    return empirical_rate_matrix_fixture()


def _biased_root_prior():
    return empirical_root_prior_fixture()


def empirical_rate_matrix_fixture():
    import numpy

    state_order = _protein_state_order()
    state_index = {state: index for index, state in enumerate(state_order)}
    rate_matrix = numpy.full((len(state_order), len(state_order)), 0.02, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
    boosted_pairs = {
        ("A", "C"): 0.45,
        ("C", "D"): 0.35,
        ("D", "E"): 0.55,
        ("A", "E"): 0.20,
    }
    for (left_state, right_state), rate in boosted_pairs.items():
        left_index = state_index[left_state]
        right_index = state_index[right_state]
        rate_matrix[left_index, right_index] = rate
        rate_matrix[right_index, left_index] = rate
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(rate_matrix[row_index, :].sum())
    return rate_matrix


def empirical_root_prior_fixture():
    import numpy

    prior = numpy.full(20, 0.02, dtype=float)
    state_index = {state: index for index, state in enumerate(_protein_state_order())}
    prior[state_index["A"]] = 0.18
    prior[state_index["C"]] = 0.10
    prior[state_index["D"]] = 0.14
    prior[state_index["E"]] = 0.12
    prior[state_index["F"]] = 0.06
    return prior / float(prior.sum())


def _protein_state_order() -> tuple[str, ...]:
    return (
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
