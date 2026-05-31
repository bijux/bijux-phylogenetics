from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_jc69_tree_likelihood,
    optimize_jc69_branch_lengths_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_jc69_branch_optimization_improves_two_tip_likelihood() -> None:
    tree_path = fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments",
        "jc69_branch_optimization_alignment_2_taxa.fasta",
    )

    report = optimize_jc69_branch_lengths_from_alignment(tree_path, alignment_path)
    optimized_tree = loads_newick(report.optimized_tree_newick)
    optimized_branch_sum = sum(
        float(child.branch_length or 0.0)
        for _parent, child in optimized_tree.iter_edges()
    )
    expected_branch_sum = -(3.0 / 4.0) * math.log(1.0 - ((4.0 / 3.0) * 0.25))

    assert report.initial_tree_newick == "(A:0.8,B:0.8);"
    assert report.branch_count == 2
    assert report.site_count == 8
    assert report.pattern_count == 6
    assert report.optimized_log_likelihood > report.initial_log_likelihood
    assert report.function_evaluation_count > 1
    assert report.converged is True
    assert any(step.accepted for step in report.steps)
    assert math.isclose(
        optimized_branch_sum,
        expected_branch_sum,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    reevaluated_report = evaluate_jc69_tree_likelihood(
        optimized_tree,
        load_fasta_alignment(alignment_path),
    )
    assert math.isclose(
        reevaluated_report.log_likelihood,
        report.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_jc69_branch_optimization_beats_starting_tree_directly() -> None:
    starting_tree = load_tree(
        fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk")
    )
    records = load_fasta_alignment(
        fixture("alignments", "jc69_branch_optimization_alignment_2_taxa.fasta")
    )
    starting_report = evaluate_jc69_tree_likelihood(starting_tree, records)
    optimization_report = optimize_jc69_branch_lengths_from_alignment(
        fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"),
        fixture("alignments", "jc69_branch_optimization_alignment_2_taxa.fasta"),
    )

    assert optimization_report.initial_log_likelihood == starting_report.log_likelihood
    assert optimization_report.optimized_log_likelihood > starting_report.log_likelihood
