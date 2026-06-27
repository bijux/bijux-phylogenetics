from __future__ import annotations

import math
from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.ancestral.discrete.likelihood.rate_matrix import (
    rate_matrix_from_log_parameters,
)
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
    EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
    build_equal_rate_codon_ctmc_rate_matrix,
    compress_alignment_site_patterns_from_records,
    compute_ctmc_expected_substitution_rate,
    evaluate_codon_ctmc_tree_likelihood,
    evaluate_empirical_protein_tree_likelihood,
    evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
    evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_alignment,
    evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment,
    evaluate_finite_state_site_log_likelihood,
    evaluate_finite_state_tree_log_likelihood,
    evaluate_protein_poisson_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.codon import (
    compress_codon_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.codon_finite_state import (
    build_codon_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.discrete_finite_state import (
    build_discrete_mk_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import FiniteStateLikelihoodSurface
from bijux_phylogenetics.phylo.likelihood.poisson import (
    protein_poisson_rate_matrix,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.protein import (
    PROTEIN_STATE_ORDER,
    UNIFORM_PROTEIN_ROOT_PRIOR,
)
from bijux_phylogenetics.phylo.likelihood.protein_finite_state import (
    build_empirical_protein_finite_state_surface,
    build_protein_poisson_finite_state_surface,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_protein_poisson_declares_expected_rate_one_and_preserves_scale_compensated_likelihood() -> (
    None
):
    tree = load_tree(fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta")
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    surface = build_protein_poisson_finite_state_surface(records)
    report = evaluate_protein_poisson_tree_likelihood(tree, records)
    rate_matrix = protein_poisson_rate_matrix()

    assert surface.normalization_policy == EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY
    assert report.normalization_policy == EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY
    assert math.isclose(
        compute_ctmc_expected_substitution_rate(
            rate_matrix,
            UNIFORM_PROTEIN_ROOT_PRIOR,
            state_labels=PROTEIN_STATE_ORDER,
        ),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.log_likelihood,
        evaluate_finite_state_tree_log_likelihood(
            tree,
            compressed_patterns,
            surface=_scale_compensated_surface(
                surface,
                rate_matrix=rate_matrix,
                scale_factor=9.0,
            ),
        ),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_empirical_protein_declares_caller_scale_and_preserves_scale_compensated_likelihood() -> (
    None
):
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
    report = evaluate_empirical_protein_tree_likelihood(
        tree,
        records,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        matrix_label="compact-polar",
    )

    assert surface.normalization_policy == CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY
    assert report.normalization_policy == CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY
    assert math.isclose(
        report.log_likelihood,
        evaluate_finite_state_tree_log_likelihood(
            tree,
            compressed_patterns,
            surface=_scale_compensated_surface(
                surface,
                rate_matrix=rate_matrix,
                scale_factor=5.0,
            ),
        ),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_empirical_protein_rate_heterogeneity_reports_preserve_caller_scale_policy() -> (
    None
):
    tree_path = fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments", "empirical_protein_invariant_mixture_alignment_2_taxa.fasta"
    )
    kwargs = {
        "rate_matrix": _compact_polar_rate_matrix(),
        "root_prior": _biased_root_prior(),
        "matrix_label": "compact-polar",
    }

    gamma_report = evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_alignment(
        tree_path,
        alignment_path,
        alpha=0.8,
        category_count=4,
        **kwargs,
    )
    invariant_report = evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment(
        tree_path,
        alignment_path,
        invariant_proportion=0.35,
        **kwargs,
    )
    combined_report = evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
        tree_path,
        alignment_path,
        alpha=0.8,
        category_count=4,
        invariant_proportion=0.35,
        **kwargs,
    )

    assert gamma_report.normalization_policy == CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY
    assert (
        invariant_report.normalization_policy
        == CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY
    )
    assert (
        combined_report.normalization_policy
        == CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY
    )


def test_codon_ctmc_declares_expected_rate_one_and_preserves_scale_compensated_likelihood() -> (
    None
):
    tree = load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta")
    )
    compressed_patterns = compress_codon_site_patterns_from_records(records)
    surface = build_codon_finite_state_surface(records)
    report = evaluate_codon_ctmc_tree_likelihood(tree, records)
    state_space, rate_matrix, frequencies, _ = build_equal_rate_codon_ctmc_rate_matrix()

    assert surface.normalization_policy == EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY
    assert report.normalization_policy == EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY
    assert math.isclose(
        compute_ctmc_expected_substitution_rate(
            rate_matrix,
            frequencies,
            state_labels=state_space.state_order,
        ),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.log_likelihood,
        evaluate_finite_state_tree_log_likelihood(
            tree,
            compressed_patterns,
            surface=_scale_compensated_surface(
                surface,
                rate_matrix=rate_matrix,
                scale_factor=11.0,
            ),
        ),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


@pytest.mark.parametrize(
    ("model_name", "state_order", "states_by_taxon", "rate_matrix"),
    [
        (
            "unordered discrete Mk",
            ["north", "south", "west"],
            {"A": "north", "B": "south"},
            rate_matrix_from_log_parameters(
                numpy.log(numpy.array([0.8], dtype=float)),
                state_order=["north", "south", "west"],
                model="equal-rates",
                state_ordering="unordered",
                allowed_transition_pairs={
                    (0, 1),
                    (0, 2),
                    (1, 0),
                    (1, 2),
                    (2, 0),
                    (2, 1),
                },
            ),
        ),
        (
            "ordered discrete Mk",
            ["north", "south", "island"],
            {"A": "north", "B": "island"},
            rate_matrix_from_log_parameters(
                numpy.log(numpy.array([0.6], dtype=float)),
                state_order=["north", "south", "island"],
                model="equal-rates",
                state_ordering="ordered",
                allowed_transition_pairs={(0, 1), (1, 0), (1, 2), (2, 1)},
            ),
        ),
        (
            "irreversible discrete Mk",
            ["present", "absent"],
            {"A": "present", "B": "absent"},
            rate_matrix_from_log_parameters(
                numpy.log(numpy.array([0.4], dtype=float)),
                state_order=["present", "absent"],
                model="all-rates-different",
                state_ordering="unordered",
                allowed_transition_pairs={(0, 1)},
            ),
        ),
    ],
)
def test_discrete_mk_surfaces_declare_caller_scale_and_preserve_scale_compensated_likelihood(
    model_name: str,
    state_order: list[str],
    states_by_taxon: dict[str, str],
    rate_matrix: numpy.ndarray,
) -> None:
    tree = load_tree(fixture("trees", "felsenstein_two_tip_tree.nwk"))
    root_prior = numpy.full(len(state_order), 1.0 / len(state_order), dtype=float)
    surface = build_discrete_mk_finite_state_surface(
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        model_name=model_name,
    )
    scaled_surface = _scale_compensated_surface(
        surface,
        rate_matrix=rate_matrix,
        scale_factor=13.0,
    )
    states = tuple(states_by_taxon[taxon] for taxon in ["A", "B"])
    context = _compile_context(tree, surface)

    assert surface.normalization_policy == CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY
    assert math.isclose(
        evaluate_finite_state_site_log_likelihood(
            context,
            states,
            taxon_order=["A", "B"],
        ),
        evaluate_finite_state_site_log_likelihood(_compile_context(tree, scaled_surface), states, taxon_order=["A", "B"]),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def _compile_context(tree, surface: FiniteStateLikelihoodSurface):
    from bijux_phylogenetics.phylo.likelihood.finite_state import (
        compile_finite_state_likelihood_context,
    )

    return compile_finite_state_likelihood_context(tree, surface)


def _scale_compensated_surface(
    surface: FiniteStateLikelihoodSurface,
    *,
    rate_matrix: numpy.ndarray,
    scale_factor: float,
) -> FiniteStateLikelihoodSurface:
    return FiniteStateLikelihoodSurface(
        model_name=surface.model_name,
        state_labels=surface.state_labels,
        root_prior=surface.root_prior,
        root_prior_source=surface.root_prior_source,
        observation_policy=surface.observation_policy,
        transition_matrix_for_branch_length=lambda branch_length: (
            transition_probability_matrix(
                rate_matrix * scale_factor,
                branch_length / scale_factor,
            )
        ),
        leaf_likelihood_vector_for_state=surface.leaf_likelihood_vector_for_state,
        normalization_policy=surface.normalization_policy,
        parameter_values=surface.parameter_values,
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
    state_index = {state: index for index, state in enumerate(PROTEIN_STATE_ORDER)}
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
    state_index = {state: index for index, state in enumerate(PROTEIN_STATE_ORDER)}
    rate_matrix = numpy.full(
        (len(PROTEIN_STATE_ORDER), len(PROTEIN_STATE_ORDER)),
        0.02,
        dtype=float,
    )
    numpy.fill_diagonal(rate_matrix, 0.0)
    for (left_state, right_state), rate in boosted_pairs.items():
        left_index = state_index[left_state]
        right_index = state_index[right_state]
        rate_matrix[left_index, right_index] = rate
        rate_matrix[right_index, left_index] = rate
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    return rate_matrix
