from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_f81_tree_likelihood,
    evaluate_gtr_tree_likelihood,
    evaluate_hky85_tree_likelihood,
    optimize_gtr_exchangeabilities_from_alignment,
    optimize_hky85_kappa,
)

pytestmark = pytest.mark.slow

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_gtr_exchangeability_optimization_recovers_dominant_pattern_over_hky() -> None:
    tree = load_tree(
        fixture("trees", "gtr_exchangeability_optimization_tree_2_taxa.nwk")
    )
    records = load_fasta_alignment(
        fixture("alignments", "gtr_exchangeability_optimization_alignment_2_taxa.fasta")
    )

    f81_report = evaluate_f81_tree_likelihood(tree, records)
    hky_optimization = optimize_hky85_kappa(tree, records)
    hky_report = evaluate_hky85_tree_likelihood(
        tree,
        records,
        kappa=hky_optimization.optimized_kappa,
    )
    gtr_optimization = optimize_gtr_exchangeabilities_from_alignment(
        fixture("trees", "gtr_exchangeability_optimization_tree_2_taxa.nwk"),
        fixture(
            "alignments", "gtr_exchangeability_optimization_alignment_2_taxa.fasta"
        ),
    )
    gtr_report = evaluate_gtr_tree_likelihood(
        tree,
        records,
        exchangeabilities={
            "AC": gtr_optimization.exchangeability_ac,
            "AG": gtr_optimization.exchangeability_ag,
            "AT": gtr_optimization.exchangeability_at,
            "CG": gtr_optimization.exchangeability_cg,
            "CT": gtr_optimization.exchangeability_ct,
            "GT": gtr_optimization.exchangeability_gt,
        },
    )

    assert gtr_optimization.tree_newick == "(A:0.15,B:0.15);"
    assert gtr_optimization.site_count == 256
    assert gtr_optimization.pattern_count == 15
    assert gtr_optimization.base_frequency_source == "estimated"
    assert gtr_optimization.exchangeability_anchor == "AC=1"
    assert math.isclose(
        gtr_optimization.exchangeability_ac,
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert gtr_optimization.exchangeability_ag > gtr_optimization.exchangeability_ct
    assert gtr_optimization.exchangeability_ag > gtr_optimization.exchangeability_gt
    assert gtr_optimization.exchangeability_ag > gtr_optimization.exchangeability_cg
    assert gtr_optimization.exchangeability_ag > gtr_optimization.exchangeability_at
    assert math.isclose(
        gtr_optimization.exchangeability_ag,
        10.436093496916191,
        rel_tol=0.0,
        abs_tol=5e-6,
    )
    assert math.isclose(
        gtr_optimization.exchangeability_at,
        0.624578992788826,
        rel_tol=0.0,
        abs_tol=5e-6,
    )
    assert math.isclose(
        gtr_optimization.exchangeability_cg,
        2.351009315272331,
        rel_tol=0.0,
        abs_tol=5e-6,
    )
    assert math.isclose(
        gtr_optimization.exchangeability_ct,
        2.066395062525222,
        rel_tol=0.0,
        abs_tol=5e-6,
    )
    assert math.isclose(
        gtr_optimization.exchangeability_gt,
        5.880958916954807,
        rel_tol=0.0,
        abs_tol=5e-6,
    )
    assert gtr_optimization.parameter_count == 8
    assert math.isclose(
        gtr_optimization.initial_log_likelihood,
        f81_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert gtr_optimization.optimized_log_likelihood > hky_report.log_likelihood
    assert gtr_optimization.optimized_aic < hky_report.aic
    assert gtr_optimization.optimization_pass_count >= 1
    assert gtr_optimization.function_evaluation_count > 1
    assert gtr_optimization.optimization_pass_count == 24
    assert gtr_optimization.converged is True
    assert math.isclose(
        gtr_report.log_likelihood,
        gtr_optimization.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        gtr_report.aic,
        gtr_optimization.optimized_aic,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
