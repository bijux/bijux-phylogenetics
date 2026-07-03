from __future__ import annotations

import math
from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import codon as codon_module
from bijux_phylogenetics.phylo.likelihood import empirical as empirical_module
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_codon_ctmc_tree_likelihood,
    evaluate_empirical_protein_tree_likelihood,
    evaluate_f81_tree_likelihood,
    evaluate_gtr_tree_likelihood,
    evaluate_hky85_tree_likelihood,
    evaluate_jc69_tree_likelihood,
    evaluate_k80_tree_likelihood,
    evaluate_protein_poisson_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood import f81 as f81_module
from bijux_phylogenetics.phylo.likelihood import gtr as gtr_module
from bijux_phylogenetics.phylo.likelihood import hky85 as hky85_module
from bijux_phylogenetics.phylo.likelihood import jc69 as jc69_module
from bijux_phylogenetics.phylo.likelihood import k80 as k80_module
from bijux_phylogenetics.phylo.likelihood import poisson as poisson_module

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


@pytest.mark.parametrize(
    ("module", "evaluator", "tree_name", "alignment_name", "kwargs"),
    [
        (
            jc69_module,
            evaluate_jc69_tree_likelihood,
            "jc69_likelihood_tree_2_taxa.nwk",
            "jc69_likelihood_alignment_2_taxa.fasta",
            {},
        ),
        (
            k80_module,
            evaluate_k80_tree_likelihood,
            "k80_likelihood_tree_2_taxa.nwk",
            "k80_likelihood_alignment_2_taxa.fasta",
            {"kappa": 4.0},
        ),
        (
            f81_module,
            evaluate_f81_tree_likelihood,
            "f81_likelihood_tree_2_taxa.nwk",
            "f81_likelihood_alignment_2_taxa.fasta",
            {
                "base_frequencies": {"A": 0.55, "C": 0.15, "G": 0.15, "T": 0.15},
            },
        ),
        (
            hky85_module,
            evaluate_hky85_tree_likelihood,
            "hky85_likelihood_tree_2_taxa.nwk",
            "hky85_likelihood_alignment_2_taxa.fasta",
            {
                "kappa": 4.0,
                "base_frequencies": {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
            },
        ),
        (
            gtr_module,
            evaluate_gtr_tree_likelihood,
            "gtr_likelihood_tree_2_taxa.nwk",
            "gtr_likelihood_alignment_2_taxa.fasta",
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
        ),
        (
            codon_module,
            evaluate_codon_ctmc_tree_likelihood,
            "codon_likelihood_tree_2_taxa.nwk",
            "codon_likelihood_alignment_2_taxa.fasta",
            {},
        ),
    ],
)
def test_public_sequence_likelihood_reports_route_through_shared_tree_helper(
    monkeypatch: pytest.MonkeyPatch,
    module,
    evaluator,
    tree_name: str,
    alignment_name: str,
    kwargs: dict[str, object],
) -> None:
    tree = load_tree(fixture("trees", tree_name))
    records = load_fasta_alignment(fixture("alignments", alignment_name))
    routed_calls = 0
    original = module.evaluate_finite_state_tree_log_likelihood_from_patterns

    def wrapped(*args, **inner_kwargs):
        nonlocal routed_calls
        routed_calls += 1
        return original(*args, **inner_kwargs)

    monkeypatch.setattr(
        module,
        "evaluate_finite_state_tree_log_likelihood_from_patterns",
        wrapped,
    )

    report = evaluator(tree, records, **kwargs)

    assert routed_calls == 1
    assert math.isfinite(report.log_likelihood)


@pytest.mark.parametrize(
    ("module", "evaluator", "tree_name", "alignment_name", "kwargs_factory"),
    [
        (
            poisson_module,
            evaluate_protein_poisson_tree_likelihood,
            "protein_poisson_likelihood_tree_2_taxa.nwk",
            "protein_poisson_likelihood_alignment_2_taxa.fasta",
            dict,
        ),
        (
            empirical_module,
            evaluate_empirical_protein_tree_likelihood,
            "empirical_protein_likelihood_tree_2_taxa.nwk",
            "empirical_protein_likelihood_alignment_2_taxa.fasta",
            lambda: {
                "rate_matrix": _compact_polar_rate_matrix(),
                "root_prior": _biased_root_prior(),
                "matrix_label": "compact-polar",
            },
        ),
    ],
)
def test_public_protein_reports_route_through_shared_site_row_helper(
    monkeypatch: pytest.MonkeyPatch,
    module,
    evaluator,
    tree_name: str,
    alignment_name: str,
    kwargs_factory,
) -> None:
    tree = load_tree(fixture("trees", tree_name))
    records = load_fasta_alignment(fixture("alignments", alignment_name))
    kwargs = kwargs_factory()
    routed_calls = 0
    original = module.evaluate_finite_state_site_log_likelihood_rows_from_patterns

    def wrapped(*args, **inner_kwargs):
        nonlocal routed_calls
        routed_calls += 1
        return original(*args, **inner_kwargs)

    monkeypatch.setattr(
        module,
        "evaluate_finite_state_site_log_likelihood_rows_from_patterns",
        wrapped,
    )

    report = evaluator(tree, records, **kwargs)

    assert routed_calls == 1
    assert report.site_log_likelihoods
    assert math.isfinite(report.log_likelihood)


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
