from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_f81_tree_likelihood,
    evaluate_gtr_tree_likelihood,
    evaluate_hky85_tree_likelihood,
    evaluate_jc69_tree_likelihood,
    evaluate_k80_tree_likelihood,
)
import bijux_phylogenetics.phylo.likelihood.fixed_topology_branch_lengths as branch_optimization
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


@pytest.mark.parametrize(
    ("model_name", "alignment_name", "optimization_kwargs", "evaluator"),
    [
        (
            "jc69",
            "jc69_branch_optimization_alignment_2_taxa.fasta",
            {},
            evaluate_jc69_tree_likelihood,
        ),
        (
            "k80",
            "k80_likelihood_alignment_2_taxa.fasta",
            {"kappa": 4.0},
            evaluate_k80_tree_likelihood,
        ),
        (
            "f81",
            "f81_likelihood_alignment_2_taxa.fasta",
            {
                "base_frequencies": {
                    "A": 0.4,
                    "C": 0.1,
                    "G": 0.2,
                    "T": 0.3,
                }
            },
            evaluate_f81_tree_likelihood,
        ),
        (
            "hky85",
            "hky85_likelihood_alignment_2_taxa.fasta",
            {
                "kappa": 4.0,
                "base_frequencies": {
                    "A": 0.4,
                    "C": 0.1,
                    "G": 0.2,
                    "T": 0.3,
                },
            },
            evaluate_hky85_tree_likelihood,
        ),
        (
            "gtr",
            "gtr_likelihood_alignment_2_taxa.fasta",
            {
                "exchangeabilities": {
                    "AC": 1.0,
                    "AG": 4.5,
                    "AT": 0.8,
                    "CG": 1.6,
                    "CT": 2.4,
                    "GT": 3.1,
                },
                "base_frequencies": {
                    "A": 0.4,
                    "C": 0.1,
                    "G": 0.2,
                    "T": 0.3,
                },
            },
            evaluate_gtr_tree_likelihood,
        ),
    ],
)
def test_fixed_topology_nucleotide_branch_optimization_improves_selected_model(
    model_name: str,
    alignment_name: str,
    optimization_kwargs: dict[str, object],
    evaluator,
) -> None:
    tree_path = fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk")
    alignment_path = fixture("alignments", alignment_name)

    report = branch_optimization.optimize_fixed_topology_nucleotide_branch_lengths_from_alignment(
        tree_path,
        alignment_path,
        model_name=model_name,
        **optimization_kwargs,
    )

    optimized_tree = loads_newick(report.optimized_tree_newick)
    reevaluated_report = evaluator(
        optimized_tree,
        load_fasta_alignment(alignment_path),
        **optimization_kwargs,
    )

    assert report.model_name == model_name.upper()
    assert report.taxa == ["A", "B"]
    assert report.branch_count == 2
    assert report.initial_tree_newick == "(A:0.8,B:0.8);"
    assert report.parameter_count == len(report.fixed_parameter_values)
    assert report.optimized_log_likelihood > report.initial_log_likelihood
    assert report.function_evaluation_count > 1
    assert all(row.optimized_branch_length >= 0.0 for row in report.branches)
    assert any(
        not math.isclose(
            row.initial_branch_length,
            row.optimized_branch_length,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for row in report.branches
    )
    assert math.isclose(
        reevaluated_report.log_likelihood,
        report.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_fixed_topology_nucleotide_branch_optimization_rejects_starting_lengths_outside_bounds_before_evaluation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tree = load_tree(fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta")
    )
    for _parent, child in tree.iter_edges():
        child.branch_length = 5.5

    def fail_if_called(*args: object, **kwargs: object) -> object:
        raise AssertionError(
            "fixed-topology nucleotide branch objective should not run"
        )

    monkeypatch.setattr(
        branch_optimization,
        "evaluate_selected_nucleotide_log_likelihood_from_patterns",
        fail_if_called,
    )

    with pytest.raises(
        InvalidBranchLengthError,
        match=(
            "fixed-topology nucleotide branch optimization requires every starting "
            "branch length to lie within the declared bounds"
        ),
    ):
        branch_optimization.optimize_fixed_topology_nucleotide_branch_lengths(
            tree,
            records,
            model_name="jc69",
            lower_branch_length_bound=0.0,
            upper_branch_length_bound=5.0,
        )
