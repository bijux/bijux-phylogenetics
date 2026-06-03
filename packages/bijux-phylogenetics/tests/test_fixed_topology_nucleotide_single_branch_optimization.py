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
def test_fixed_topology_nucleotide_single_branch_optimization_keeps_other_branches_fixed(
    model_name: str,
    alignment_name: str,
    optimization_kwargs: dict[str, object],
    evaluator,
) -> None:
    tree = load_tree(fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"))
    alignment_path = fixture("alignments", alignment_name)
    branch_ids_by_name = {
        child.name or "": child.node_id or "" for _parent, child in tree.iter_edges()
    }
    selected_branch_id = branch_ids_by_name["A"]
    unchanged_branch_id = branch_ids_by_name["B"]

    report = branch_optimization.optimize_fixed_topology_nucleotide_single_branch_length_from_alignment(
        fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"),
        alignment_path,
        model_name=model_name,
        branch_id=selected_branch_id,
        **optimization_kwargs,
    )

    optimized_tree = loads_newick(report.optimized_tree_newick)
    optimized_lengths = {
        child.node_id or "": float(child.branch_length or 0.0)
        for _parent, child in optimized_tree.iter_edges()
    }
    starting_lengths = {
        child.node_id or "": float(child.branch_length or 0.0)
        for _parent, child in tree.iter_edges()
    }
    reevaluated_report = evaluator(
        optimized_tree,
        load_fasta_alignment(alignment_path),
        **optimization_kwargs,
    )

    assert report.model_name == model_name.upper()
    assert report.branch_count == 2
    assert report.selected_branch.branch_id == selected_branch_id
    assert report.selected_branch.child_name == "A"
    assert report.selected_branch.descendant_taxa == ["A"]
    assert report.unchanged_branch_count == 1
    assert report.unchanged_branch_ids == [unchanged_branch_id]
    assert report.optimized_log_likelihood >= report.initial_log_likelihood
    assert math.isclose(
        optimized_lengths[unchanged_branch_id],
        starting_lengths[unchanged_branch_id],
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        report.selected_branch.initial_branch_length,
        report.selected_branch.optimized_branch_length,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.selected_branch.initial_branch_length,
        starting_lengths[selected_branch_id],
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.selected_branch.optimized_branch_length,
        optimized_lengths[selected_branch_id],
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        reevaluated_report.log_likelihood,
        report.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_fixed_topology_nucleotide_single_branch_optimization_rejects_missing_branch_before_evaluation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tree = load_tree(fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta")
    )

    def fail_if_called(*args: object, **kwargs: object) -> object:
        raise AssertionError("single-branch likelihood evaluator should not run")

    monkeypatch.setattr(
        branch_optimization,
        "evaluate_selected_nucleotide_log_likelihood_from_patterns",
        fail_if_called,
    )

    with pytest.raises(
        ValueError, match="tree does not contain branch_id 'missing-branch'"
    ):
        branch_optimization.optimize_fixed_topology_nucleotide_single_branch_length(
            tree,
            records,
            model_name="jc69",
            branch_id="missing-branch",
        )
