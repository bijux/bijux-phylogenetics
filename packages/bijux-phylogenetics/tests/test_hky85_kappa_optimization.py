from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_f81_tree_likelihood,
    evaluate_hky85_tree_likelihood,
    evaluate_jc69_tree_likelihood,
    evaluate_k80_tree_likelihood,
    optimize_hky85_kappa_from_alignment,
    optimize_k80_kappa,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_hky85_kappa_optimization_prefers_hky_fixture_over_jc69_f81_and_k80() -> None:
    tree = load_tree(fixture("trees", "hky85_kappa_optimization_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "hky85_kappa_optimization_alignment_2_taxa.fasta")
    )

    jc69_report = evaluate_jc69_tree_likelihood(tree, records)
    f81_report = evaluate_f81_tree_likelihood(tree, records)
    k80_optimization = optimize_k80_kappa(tree, records)
    k80_report = evaluate_k80_tree_likelihood(
        tree,
        records,
        kappa=k80_optimization.optimized_kappa,
    )
    hky85_optimization = optimize_hky85_kappa_from_alignment(
        fixture("trees", "hky85_kappa_optimization_tree_2_taxa.nwk"),
        fixture("alignments", "hky85_kappa_optimization_alignment_2_taxa.fasta"),
        initial_kappa=1.0,
    )
    hky85_report = evaluate_hky85_tree_likelihood(
        tree,
        records,
        kappa=hky85_optimization.optimized_kappa,
    )

    jc69_aic = -2.0 * jc69_report.log_likelihood
    k80_aic = (-2.0 * k80_report.log_likelihood) + 2.0

    assert hky85_optimization.tree_newick == "(A:0.15,B:0.15);"
    assert hky85_optimization.site_count == 128
    assert hky85_optimization.pattern_count == 16
    assert hky85_optimization.base_frequency_source == "estimated"
    assert hky85_optimization.base_frequency_a > 0.4
    assert hky85_optimization.base_frequency_t > 0.25
    assert hky85_optimization.base_frequency_c < 0.12
    assert hky85_optimization.base_frequency_g < 0.16
    assert hky85_optimization.parameter_count == 4
    assert hky85_optimization.initial_kappa == 1.0
    assert hky85_optimization.optimized_kappa > hky85_optimization.initial_kappa
    assert math.isclose(
        hky85_optimization.optimized_kappa,
        5.55724924899235,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert hky85_optimization.optimized_log_likelihood > jc69_report.log_likelihood
    assert hky85_optimization.optimized_log_likelihood > f81_report.log_likelihood
    assert hky85_optimization.optimized_log_likelihood > k80_report.log_likelihood
    assert hky85_optimization.optimized_aic < jc69_aic
    assert hky85_optimization.optimized_aic < f81_report.aic
    assert hky85_optimization.optimized_aic < k80_aic
    assert hky85_optimization.function_evaluation_count > 1
    assert hky85_optimization.converged is True
    assert math.isclose(
        hky85_report.log_likelihood,
        hky85_optimization.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_hky85_kappa_optimization_reevaluates_to_same_aic() -> None:
    tree = load_tree(fixture("trees", "hky85_kappa_optimization_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "hky85_kappa_optimization_alignment_2_taxa.fasta")
    )

    optimization = optimize_hky85_kappa_from_alignment(
        fixture("trees", "hky85_kappa_optimization_tree_2_taxa.nwk"),
        fixture("alignments", "hky85_kappa_optimization_alignment_2_taxa.fasta"),
        initial_kappa=1.0,
    )
    reevaluated = evaluate_hky85_tree_likelihood(
        tree,
        records,
        kappa=optimization.optimized_kappa,
    )

    assert math.isclose(
        reevaluated.base_frequency_a,
        optimization.base_frequency_a,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        reevaluated.base_frequency_c,
        optimization.base_frequency_c,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        reevaluated.base_frequency_g,
        optimization.base_frequency_g,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        reevaluated.base_frequency_t,
        optimization.base_frequency_t,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        reevaluated.log_likelihood,
        optimization.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        reevaluated.aic,
        optimization.optimized_aic,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
