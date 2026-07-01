from __future__ import annotations

from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
import bijux_phylogenetics.phylo.likelihood.empirical as empirical_module
import bijux_phylogenetics.phylo.likelihood.hky85 as hky85_module
import bijux_phylogenetics.phylo.likelihood.jc69 as jc69_module
import bijux_phylogenetics.phylo.likelihood.k80 as k80_module
import bijux_phylogenetics.phylo.likelihood.substitution_model_selection as substitution_model_selection_module
from bijux_phylogenetics.phylo.likelihood.dna import (
    normalize_unambiguous_dna_records,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    compress_alignment_site_patterns_from_records,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_jc69_branch_length_optimization_routes_through_shared_likelihood_search(
    monkeypatch,
) -> None:
    observed_calls: list[tuple[tuple[float, float], str, str]] = []
    original_search = jc69_module.run_bounded_likelihood_search

    def observe_search(**kwargs):
        observed_calls.append(
            (
                (kwargs["lower_bound"], kwargs["upper_bound"]),
                kwargs["parameter_name"],
                kwargs["owner_name"],
            )
        )
        return original_search(**kwargs)

    monkeypatch.setattr(jc69_module, "run_bounded_likelihood_search", observe_search)

    report = jc69_module.optimize_jc69_branch_lengths_from_alignment(
        fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"),
        fixture("alignments", "jc69_branch_optimization_alignment_2_taxa.fasta"),
    )

    assert report.converged is True
    assert observed_calls
    assert all(call[0] == (1e-06, 5.0) for call in observed_calls)
    assert all(call[1] == "branch_length" for call in observed_calls)
    assert all(call[2] == "JC69 branch-length fit" for call in observed_calls)


def test_k80_kappa_optimization_routes_through_shared_likelihood_search(
    monkeypatch,
) -> None:
    observed_calls: list[tuple[tuple[float, float], str, str]] = []
    original_search = k80_module.run_bounded_likelihood_search

    def observe_search(**kwargs):
        observed_calls.append(
            (
                (kwargs["lower_bound"], kwargs["upper_bound"]),
                kwargs["parameter_name"],
                kwargs["owner_name"],
            )
        )
        return original_search(**kwargs)

    monkeypatch.setattr(k80_module, "run_bounded_likelihood_search", observe_search)

    report = k80_module.optimize_k80_kappa_from_alignment(
        fixture("trees", "k80_kappa_optimization_tree_2_taxa.nwk"),
        fixture("alignments", "k80_kappa_optimization_alignment_2_taxa.fasta"),
    )

    assert report.converged is True
    assert observed_calls == [((0.05, 20.0), "kappa", "K80 kappa fit")]


def test_hky85_kappa_optimization_routes_through_shared_likelihood_search(
    monkeypatch,
) -> None:
    observed_calls: list[tuple[tuple[float, float], str, str]] = []
    original_search = hky85_module.run_bounded_likelihood_search

    def observe_search(**kwargs):
        observed_calls.append(
            (
                (kwargs["lower_bound"], kwargs["upper_bound"]),
                kwargs["parameter_name"],
                kwargs["owner_name"],
            )
        )
        return original_search(**kwargs)

    monkeypatch.setattr(hky85_module, "run_bounded_likelihood_search", observe_search)

    report = hky85_module.optimize_hky85_kappa_from_alignment(
        fixture("trees", "hky85_kappa_optimization_tree_2_taxa.nwk"),
        fixture("alignments", "hky85_kappa_optimization_alignment_2_taxa.fasta"),
    )

    assert report.converged is True
    assert observed_calls == [((0.05, 20.0), "kappa", "HKY85 kappa fit")]


def test_empirical_invariant_mixture_routes_through_shared_likelihood_search(
    monkeypatch,
) -> None:
    observed_calls: list[tuple[tuple[float, float], str, str]] = []
    original_search = empirical_module.run_bounded_likelihood_search

    def observe_search(**kwargs):
        observed_calls.append(
            (
                (kwargs["lower_bound"], kwargs["upper_bound"]),
                kwargs["parameter_name"],
                kwargs["owner_name"],
            )
        )
        return original_search(**kwargs)

    monkeypatch.setattr(
        empirical_module,
        "run_bounded_likelihood_search",
        observe_search,
    )

    report = (
        empirical_module.optimize_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment(
            fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
            fixture(
                "alignments",
                "empirical_protein_invariant_mixture_alignment_2_taxa.fasta",
            ),
            rate_matrix=_compact_polar_rate_matrix(),
            root_prior=_biased_root_prior(),
            matrix_label="compact-polar",
            initial_invariant_proportion=0.1,
        )
    )

    assert report.converged is True
    assert observed_calls == [
        (
            (0.0, 0.95),
            "invariant_proportion",
            "empirical protein invariant-proportion fit",
        )
    ]


def test_empirical_gamma_invariant_mixture_routes_through_shared_likelihood_search(
    monkeypatch,
) -> None:
    observed_calls: list[tuple[tuple[float, float], str, str]] = []
    original_search = empirical_module.run_bounded_likelihood_search

    def observe_search(**kwargs):
        observed_calls.append(
            (
                (kwargs["lower_bound"], kwargs["upper_bound"]),
                kwargs["parameter_name"],
                kwargs["owner_name"],
            )
        )
        return original_search(**kwargs)

    monkeypatch.setattr(
        empirical_module,
        "run_bounded_likelihood_search",
        observe_search,
    )

    report = (
        empirical_module.optimize_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
            fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
            fixture(
                "alignments",
                "empirical_protein_invariant_mixture_alignment_2_taxa.fasta",
            ),
            rate_matrix=_compact_polar_rate_matrix(),
            root_prior=_biased_root_prior(),
            alpha=0.8,
            category_count=4,
            matrix_label="compact-polar",
            initial_invariant_proportion=0.1,
        )
    )

    assert report.converged is True
    assert observed_calls == [
        (
            (0.0, 0.95),
            "invariant_proportion",
            "empirical protein gamma-invariant fit",
        )
    ]


def test_gamma_candidate_routes_alpha_through_shared_coordinate_search(
    monkeypatch,
) -> None:
    observed_parameter_sets: list[set[str]] = []
    observed_owner_names: list[str] = []
    original_search = (
        substitution_model_selection_module.run_bounded_coordinate_likelihood_search
    )

    def observe_search(**kwargs):
        observed_parameter_sets.append(set(kwargs["initial_values"]))
        observed_owner_names.append(kwargs["owner_name"])
        return original_search(**kwargs)

    monkeypatch.setattr(
        substitution_model_selection_module,
        "run_bounded_coordinate_likelihood_search",
        observe_search,
    )

    tree = load_tree(fixture("trees", "k80_kappa_optimization_tree_2_taxa.nwk"))
    records = normalize_unambiguous_dna_records(
        load_fasta_alignment(
            fixture("alignments", "k80_kappa_optimization_alignment_2_taxa.fasta")
        ),
        model_name="bounded optimizer routing test",
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    empirical_base_frequencies = (
        substitution_model_selection_module.estimate_empirical_dna_base_frequencies(
            records
        )
    )

    row = substitution_model_selection_module._fit_substitution_model_candidate(
        tree,
        records,
        compressed_patterns,
        empirical_base_frequencies=empirical_base_frequencies,
        candidate_model_name="JC69+G",
        category_count=4,
        max_coordinate_passes=8,
        improvement_tolerance=1e-6,
    )

    assert row.fit_succeeded is True
    assert observed_parameter_sets == [{"alpha"}]
    assert observed_owner_names == ["substitution model rate-heterogeneity fit"]


def _compact_polar_rate_matrix() -> numpy.ndarray:
    state_order = tuple("ACDEFGHIKLMNPQRSTVWY")
    state_index = {state: index for index, state in enumerate(state_order)}
    rate_matrix = numpy.full((len(state_order), len(state_order)), 0.02, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
    for left_state, right_state, rate in (
        ("A", "C", 0.45),
        ("C", "D", 0.35),
        ("D", "E", 0.55),
        ("A", "E", 0.20),
    ):
        left_index = state_index[left_state]
        right_index = state_index[right_state]
        rate_matrix[left_index, right_index] = rate
        rate_matrix[right_index, left_index] = rate
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    return rate_matrix


def _biased_root_prior() -> numpy.ndarray:
    prior = numpy.full(20, 0.02, dtype=float)
    state_index = {state: index for index, state in enumerate("ACDEFGHIKLMNPQRSTVWY")}
    prior[state_index["A"]] = 0.18
    prior[state_index["C"]] = 0.10
    prior[state_index["D"]] = 0.14
    prior[state_index["E"]] = 0.12
    prior[state_index["F"]] = 0.06
    return prior / float(prior.sum())
