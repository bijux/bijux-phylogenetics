from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.phylo.likelihood.fixed_topology_nucleotide_joint_optimization as joint_optimization

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_fixed_topology_nucleotide_joint_restart_policy_validation_rejects_unknown_name() -> (
    None
):
    with pytest.raises(
        ValueError,
        match="fixed-topology nucleotide joint restart policy must be one of",
    ):
        joint_optimization.validate_fixed_topology_nucleotide_joint_restart_policy(
            "retry-everything"
        )


def test_fixed_topology_nucleotide_joint_restart_none_keeps_single_attempt() -> None:
    report = joint_optimization.optimize_fixed_topology_nucleotide_branches_and_model_with_restarts_from_alignment(
        fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"),
        fixture("alignments", "k80_kappa_optimization_alignment_2_taxa.fasta"),
        model_name="k80",
        restart_policy="none",
        max_restart_count=3,
        initial_kappa=1.0,
        max_joint_passes=4,
    )

    assert report.restart_policy == "none"
    assert report.attempt_count == 1
    assert report.selected_attempt_index == 1
    assert report.attempt_rows[0].trigger_reason == "initial-attempt"
    assert report.attempt_rows[0].selected_best is True
    assert report.selected_report.model_name == "K80"


def test_fixed_topology_gtr_joint_restart_policy_records_attempts_and_selects_best_available_solution() -> (
    None
):
    report = joint_optimization.optimize_fixed_topology_nucleotide_branches_and_model_with_restarts_from_alignment(
        fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"),
        fixture("alignments", "gtr_exchangeability_optimization_alignment_2_taxa.fasta"),
        model_name="gtr",
        restart_policy="restart-on-nonconvergence-or-boundary",
        max_restart_count=2,
        max_joint_passes=3,
        max_model_coordinate_passes=6,
    )

    assert report.model_name == "GTR"
    assert report.attempt_count == 3
    assert report.selected_solution_reason == "best-available-attempt"
    assert [row.attempt_index for row in report.attempt_rows] == [1, 2, 3]
    assert report.attempt_rows[0].trigger_reason == "initial-attempt"
    assert report.attempt_rows[1].trigger_reason == "restart-after-nonconvergence"
    assert report.attempt_rows[2].trigger_reason == "restart-after-boundary-warning"
    assert sum(1 for row in report.attempt_rows if row.selected_best) == 1
    assert report.attempt_rows[report.selected_attempt_index - 1].selected_best is True
    assert report.selected_report.optimized_log_likelihood == report.attempt_rows[
        report.selected_attempt_index - 1
    ].optimized_log_likelihood
    assert all(row.boundary_warning_count == 0 for row in report.attempt_rows)
    assert report.attempt_rows[1].branch_boundary_count == 1
    assert report.attempt_rows[2].branch_boundary_count == 1
