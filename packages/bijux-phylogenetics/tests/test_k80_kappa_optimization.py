from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_jc69_tree_likelihood,
    evaluate_k80_tree_likelihood,
    optimize_k80_kappa_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_k80_kappa_optimization_improves_transition_biased_fixture() -> None:
    tree_path = fixture("trees", "k80_kappa_optimization_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments",
        "k80_kappa_optimization_alignment_2_taxa.fasta",
    )

    report = optimize_k80_kappa_from_alignment(
        tree_path,
        alignment_path,
        initial_kappa=1.0,
    )

    assert report.tree_newick == "(A:0.15,B:0.15);"
    assert report.site_count == 13
    assert report.pattern_count == 9
    assert report.initial_kappa == 1.0
    assert report.optimized_kappa > report.initial_kappa
    assert math.isclose(
        report.optimized_kappa,
        9.17375541937712,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert report.optimized_log_likelihood > report.initial_log_likelihood
    assert report.function_evaluation_count > 1
    assert report.converged is True


def test_k80_kappa_optimization_beats_jc69_on_transition_biased_fixture() -> None:
    tree = load_tree(fixture("trees", "k80_kappa_optimization_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "k80_kappa_optimization_alignment_2_taxa.fasta")
    )
    jc69_report = evaluate_jc69_tree_likelihood(tree, records)
    optimized_report = optimize_k80_kappa_from_alignment(
        fixture("trees", "k80_kappa_optimization_tree_2_taxa.nwk"),
        fixture("alignments", "k80_kappa_optimization_alignment_2_taxa.fasta"),
        initial_kappa=1.0,
    )
    reevaluated_k80_report = evaluate_k80_tree_likelihood(
        tree,
        records,
        kappa=optimized_report.optimized_kappa,
    )

    assert math.isclose(
        optimized_report.initial_log_likelihood,
        jc69_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert reevaluated_k80_report.log_likelihood > jc69_report.log_likelihood
    assert math.isclose(
        reevaluated_k80_report.log_likelihood,
        optimized_report.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
