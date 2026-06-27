from __future__ import annotations

import math
from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
    compress_alignment_site_patterns_from_records,
    compute_ctmc_expected_substitution_rate,
    evaluate_f81_tree_likelihood,
    evaluate_finite_state_tree_log_likelihood,
    evaluate_gtr_tree_likelihood,
    evaluate_hky85_tree_likelihood,
    evaluate_jc69_tree_likelihood,
    evaluate_k80_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.f81 import f81_rate_matrix
from bijux_phylogenetics.phylo.likelihood.finite_state import FiniteStateLikelihoodSurface
from bijux_phylogenetics.phylo.likelihood.gtr import gtr_rate_matrix
from bijux_phylogenetics.phylo.likelihood.hky85 import hky85_rate_matrix
from bijux_phylogenetics.phylo.likelihood.jc69 import jc69_rate_matrix
from bijux_phylogenetics.phylo.likelihood.k80 import k80_rate_matrix
from bijux_phylogenetics.phylo.likelihood.nucleotide_finite_state import (
    build_selected_nucleotide_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix

FIXTURES = Path(__file__).parent / "fixtures"
_NUCLEOTIDE_BASE_FREQUENCIES = {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}
_GTR_EXCHANGEABILITIES = {
    "AC": 1.0,
    "AG": 4.5,
    "AT": 0.8,
    "CG": 1.6,
    "CT": 2.4,
    "GT": 3.1,
}


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


@pytest.mark.parametrize(
    ("model_name", "tree_fixture", "alignment_fixture", "kwargs", "report_evaluator", "rate_matrix"),
    [
        (
            "jc69",
            "jc69_likelihood_tree_2_taxa.nwk",
            "jc69_likelihood_alignment_2_taxa.fasta",
            {},
            evaluate_jc69_tree_likelihood,
            jc69_rate_matrix(),
        ),
        (
            "k80",
            "k80_likelihood_tree_2_taxa.nwk",
            "k80_likelihood_alignment_2_taxa.fasta",
            {"kappa": 4.0},
            evaluate_k80_tree_likelihood,
            k80_rate_matrix(4.0),
        ),
        (
            "f81",
            "f81_likelihood_tree_2_taxa.nwk",
            "f81_likelihood_alignment_2_taxa.fasta",
            {"base_frequencies": _NUCLEOTIDE_BASE_FREQUENCIES},
            evaluate_f81_tree_likelihood,
            f81_rate_matrix(_NUCLEOTIDE_BASE_FREQUENCIES),
        ),
        (
            "hky85",
            "hky85_likelihood_tree_2_taxa.nwk",
            "hky85_likelihood_alignment_2_taxa.fasta",
            {"base_frequencies": _NUCLEOTIDE_BASE_FREQUENCIES, "kappa": 4.0},
            evaluate_hky85_tree_likelihood,
            hky85_rate_matrix(_NUCLEOTIDE_BASE_FREQUENCIES, kappa=4.0),
        ),
        (
            "gtr",
            "gtr_likelihood_tree_2_taxa.nwk",
            "gtr_likelihood_alignment_2_taxa.fasta",
            {
                "base_frequencies": _NUCLEOTIDE_BASE_FREQUENCIES,
                "exchangeabilities": _GTR_EXCHANGEABILITIES,
            },
            evaluate_gtr_tree_likelihood,
            gtr_rate_matrix(
                _GTR_EXCHANGEABILITIES,
                base_frequencies=_NUCLEOTIDE_BASE_FREQUENCIES,
            ),
        ),
    ],
)
def test_nucleotide_ctmc_families_declare_expected_rate_one_and_preserve_scale_compensated_likelihood(
    model_name: str,
    tree_fixture: str,
    alignment_fixture: str,
    kwargs: dict[str, object],
    report_evaluator,
    rate_matrix: numpy.ndarray,
) -> None:
    tree = load_tree(fixture("trees", tree_fixture))
    records = load_fasta_alignment(fixture("alignments", alignment_fixture))
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    surface = build_selected_nucleotide_finite_state_surface(
        records,
        model_name=model_name,
        owner_name="nucleotide normalization policy test",
        **kwargs,
    )
    report = report_evaluator(tree, records, **kwargs)

    assert surface.normalization_policy == EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY
    assert report.normalization_policy == EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY
    assert math.isclose(
        compute_ctmc_expected_substitution_rate(
            rate_matrix,
            _surface_stationary_distribution(surface),
            state_labels=surface.state_labels,
        ),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    scaled_log_likelihood = evaluate_finite_state_tree_log_likelihood(
        tree,
        compressed_patterns,
        surface=_scale_compensated_surface(
            surface,
            rate_matrix=rate_matrix,
            scale_factor=7.0,
        ),
    )

    assert math.isclose(
        report.log_likelihood,
        scaled_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def _surface_stationary_distribution(
    surface: FiniteStateLikelihoodSurface,
) -> numpy.ndarray:
    if surface.model_name in {"JC69", "K80"}:
        return surface.root_prior
    return surface.root_prior


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
