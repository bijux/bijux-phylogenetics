from __future__ import annotations

import math
from pathlib import Path

import pytest

import bijux_phylogenetics.phylo.likelihood.fixed_topology_nucleotide_joint_optimization as joint_optimization

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_fixed_topology_nucleotide_joint_optimization_rejects_models_without_joint_model_updates() -> (
    None
):
    with pytest.raises(
        ValueError,
        match="fixed-topology nucleotide joint optimization model must be one of",
    ):
        joint_optimization.validate_fixed_topology_nucleotide_joint_optimization_model(
            "jc69"
        )


def test_fixed_topology_k80_joint_optimization_records_branch_and_model_updates() -> (
    None
):
    report = joint_optimization.optimize_fixed_topology_nucleotide_branches_and_model_from_alignment(
        fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"),
        fixture("alignments", "k80_kappa_optimization_alignment_2_taxa.fasta"),
        model_name="k80",
        initial_kappa=1.0,
        upper_branch_length_bound=5.0,
        max_joint_passes=4,
    )

    assert report.model_name == "K80"
    assert report.initial_tree_newick == "(A:0.8,B:0.8);"
    assert report.optimized_log_likelihood > report.initial_log_likelihood
    assert report.joint_optimization_pass_count >= 1
    assert report.convergence_reason in {
        "joint-schedule-converged",
        "max-joint-passes-exhausted",
    }
    assert [row.update_kind for row in report.update_rows[:2]] == [
        "branch-lengths",
        "substitution-parameters",
    ]
    assert report.update_rows[0].optimized_branch_count == 2
    assert report.update_rows[0].updated_parameter_names == []
    assert report.update_rows[1].optimized_branch_count == 0
    assert report.update_rows[1].updated_parameter_names == ["kappa"]
    assert report.parameter_rows[0].parameter_name == "kappa"
    assert (
        report.parameter_rows[0].optimized_value
        > report.parameter_rows[0].initial_value
    )
    assert any(
        not math.isclose(
            row.initial_branch_length,
            row.optimized_branch_length,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for row in report.branch_rows
    )


def test_fixed_topology_gtr_joint_optimization_reports_free_exchangeability_updates() -> (
    None
):
    report = joint_optimization.optimize_fixed_topology_nucleotide_branches_and_model_from_alignment(
        fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"),
        fixture(
            "alignments", "gtr_exchangeability_optimization_alignment_2_taxa.fasta"
        ),
        model_name="gtr",
        max_joint_passes=4,
        max_model_coordinate_passes=12,
    )

    model_rows = [
        row
        for row in report.update_rows
        if row.update_kind == "substitution-parameters"
    ]
    assert report.model_name == "GTR"
    assert report.fixed_parameter_values == {"AC": 1.0}
    assert report.parameter_count == 8
    assert report.base_frequency_source == "estimated"
    assert report.optimized_log_likelihood > report.initial_log_likelihood
    assert model_rows
    assert model_rows[0].updated_parameter_names == ["AG", "AT", "CG", "CT", "GT"]
    assert [row.parameter_name for row in report.parameter_rows] == [
        "AG",
        "AT",
        "CG",
        "CT",
        "GT",
    ]
    assert any(row.optimized_value > row.initial_value for row in report.parameter_rows)
