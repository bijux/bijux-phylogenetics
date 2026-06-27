from __future__ import annotations

import math
from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    compress_alignment_site_patterns_from_records,
    evaluate_empirical_protein_tree_likelihood,
    evaluate_protein_poisson_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    compile_finite_state_likelihood_context,
    evaluate_finite_state_tree_log_likelihood_from_patterns,
)
from bijux_phylogenetics.phylo.likelihood import protein_finite_state as protein_surface_module
from bijux_phylogenetics.phylo.likelihood.protein_finite_state import (
    build_empirical_protein_finite_state_surface,
    build_protein_poisson_finite_state_surface,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_protein_poisson_surface_matches_existing_tree_report() -> None:
    tree = load_tree(fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta")
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    surface = build_protein_poisson_finite_state_surface(records)
    context = compile_finite_state_likelihood_context(tree, surface)

    protocol_log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )
    family_report = evaluate_protein_poisson_tree_likelihood(tree, records)

    assert math.isclose(
        protocol_log_likelihood,
        family_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_protein_poisson_surface_respects_ambiguity_vector_policy() -> None:
    tree = load_tree(fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture(
            "alignments",
            "protein_poisson_likelihood_alignment_ambiguity_2_taxa.fasta",
        )
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    surface = build_protein_poisson_finite_state_surface(
        records,
        ambiguity_policy="ambiguity-vector",
    )
    context = compile_finite_state_likelihood_context(tree, surface)

    protocol_log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )
    family_report = evaluate_protein_poisson_tree_likelihood(
        tree,
        records,
        ambiguity_policy="ambiguity-vector",
    )

    assert math.isclose(
        protocol_log_likelihood,
        family_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_protein_poisson_surface_uses_native_transition_surface(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tree = load_tree(fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta")
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    surface = build_protein_poisson_finite_state_surface(records)

    def _boom(*_args: object, **_kwargs: object) -> numpy.ndarray:
        raise AssertionError("protein Poisson surface should not call the generic CTMC engine")

    monkeypatch.setattr(protein_surface_module, "transition_probability_matrix", _boom)
    context = compile_finite_state_likelihood_context(tree, surface)

    assert math.isfinite(
        evaluate_finite_state_tree_log_likelihood_from_patterns(
            context,
            compressed_patterns,
        )
    )


def test_empirical_protein_surface_matches_existing_tree_report() -> None:
    tree = load_tree(fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta")
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    rate_matrix = _compact_polar_rate_matrix()
    root_prior = _biased_root_prior()
    surface = build_empirical_protein_finite_state_surface(
        records,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        matrix_label="compact-polar",
    )
    context = compile_finite_state_likelihood_context(tree, surface)

    protocol_log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )
    family_report = evaluate_empirical_protein_tree_likelihood(
        tree,
        records,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        matrix_label="compact-polar",
    )

    assert math.isclose(
        protocol_log_likelihood,
        family_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_empirical_protein_surface_respects_ambiguity_vector_policy() -> None:
    tree = load_tree(fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture(
            "alignments",
            "protein_poisson_likelihood_alignment_ambiguity_2_taxa.fasta",
        )
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    surface = build_empirical_protein_finite_state_surface(
        records,
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        ambiguity_policy="ambiguity-vector",
    )
    context = compile_finite_state_likelihood_context(tree, surface)

    protocol_log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )
    family_report = evaluate_empirical_protein_tree_likelihood(
        tree,
        records,
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        ambiguity_policy="ambiguity-vector",
    )

    assert math.isclose(
        protocol_log_likelihood,
        family_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def _compact_polar_rate_matrix() -> numpy.ndarray:
    return _build_empirical_rate_matrix(
        boosted_pairs={
            ("A", "C"): 0.45,
            ("C", "D"): 0.35,
            ("D", "E"): 0.55,
            ("A", "E"): 0.20,
        }
    )


def _biased_root_prior() -> numpy.ndarray:
    prior = numpy.full(20, 0.02, dtype=float)
    state_index = _protein_state_index()
    prior[state_index["A"]] = 0.18
    prior[state_index["C"]] = 0.10
    prior[state_index["D"]] = 0.14
    prior[state_index["E"]] = 0.12
    prior[state_index["F"]] = 0.06
    return prior / float(prior.sum())


def _build_empirical_rate_matrix(
    *,
    boosted_pairs: dict[tuple[str, str], float],
) -> numpy.ndarray:
    state_order = _protein_state_order()
    state_index = _protein_state_index()
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


def _protein_state_index() -> dict[str, int]:
    return {state: index for index, state in enumerate(_protein_state_order())}
