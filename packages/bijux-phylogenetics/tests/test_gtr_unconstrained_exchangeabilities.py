from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_gtr_tree_likelihood,
    evaluate_gtr_tree_likelihood_from_unconstrained_exchangeabilities,
)
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    parameterize_dna_exchangeability_simplex,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_gtr_likelihood_from_unconstrained_exchangeabilities_matches_constrained_surface() -> (
    None
):
    tree = load_tree(fixture("trees", "gtr_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "gtr_likelihood_alignment_2_taxa.fasta")
    )
    constrained_exchangeabilities = {
        "AC": 1.0,
        "AG": 2.0,
        "AT": 0.5,
        "CG": 1.5,
        "CT": 1.75,
        "GT": 1.25,
    }
    parameterization = parameterize_dna_exchangeability_simplex(
        constrained_exchangeabilities
    )

    constrained_report = evaluate_gtr_tree_likelihood(
        tree,
        records,
        exchangeabilities=constrained_exchangeabilities,
        base_frequencies=[0.4, 0.1, 0.2, 0.3],
    )
    unconstrained_report = (
        evaluate_gtr_tree_likelihood_from_unconstrained_exchangeabilities(
            tree,
            records,
            unconstrained_exchangeabilities=parameterization.unconstrained_values,
            base_frequencies=[0.4, 0.1, 0.2, 0.3],
        )
    )

    assert math.isclose(
        unconstrained_report.log_likelihood,
        constrained_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        unconstrained_report.exchangeability_ac,
        constrained_report.exchangeability_ac,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        unconstrained_report.exchangeability_ag,
        constrained_report.exchangeability_ag,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        unconstrained_report.exchangeability_gt,
        constrained_report.exchangeability_gt,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
