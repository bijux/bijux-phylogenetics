from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.phylo.likelihood import (
    optimize_empirical_protein_branch_lengths_from_alignment,
    optimize_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
    optimize_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment,
)
import bijux_phylogenetics.phylo.likelihood.fixed_topology_nucleotide_joint_optimization as joint_optimization

pytestmark = pytest.mark.slow

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


@pytest.mark.parametrize(
    ("report_builder", "expected_classification"),
    [
        (
            lambda: optimize_empirical_protein_branch_lengths_from_alignment(
                fixture(
                    "trees", "empirical_protein_branch_optimization_start_tree_2_taxa.nwk"
                ),
                fixture(
                    "alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"
                ),
                rate_matrix=_compact_polar_rate_matrix(),
                likelihood_model="fixed-rate",
                root_prior=_biased_root_prior(),
                matrix_label="compact-polar",
                max_coordinate_passes=4,
            ),
            "failed",
        ),
        (
            lambda: optimize_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment(
                fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
                fixture(
                    "alignments",
                    "empirical_protein_invariant_mixture_alignment_2_taxa.fasta",
                ),
                rate_matrix=_compact_polar_rate_matrix(),
                root_prior=_biased_root_prior(),
                matrix_label="compact-polar",
                initial_invariant_proportion=0.1,
            ),
            "stable",
        ),
        (
            lambda: optimize_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
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
                initial_invariant_proportion=0.31,
                lower_invariant_proportion_bound=0.3,
                upper_invariant_proportion_bound=0.5,
            ),
            "boundary",
        ),
        (
            lambda: joint_optimization.optimize_fixed_topology_nucleotide_branches_and_model_from_alignment(
                fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"),
                fixture(
                    "alignments", "gtr_exchangeability_optimization_alignment_2_taxa.fasta"
                ),
                model_name="gtr",
                max_joint_passes=4,
                max_model_coordinate_passes=12,
            ),
            None,
        ),
    ],
)
def test_workflow_reports_serialize_fit_assessment(
    report_builder,
    expected_classification: str | None,
) -> None:
    report = report_builder()
    payload = asdict(report)
    fit_assessment = payload["fit_assessment"]

    assert fit_assessment["summary"]
    assert fit_assessment["boundary_warning_count"] >= 0
    assert isinstance(fit_assessment["warning_kinds"], list)
    if expected_classification is not None:
        assert fit_assessment["classification"] == expected_classification
    else:
        assert fit_assessment["classification"] in {"failed", "weak", "boundary"}


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
