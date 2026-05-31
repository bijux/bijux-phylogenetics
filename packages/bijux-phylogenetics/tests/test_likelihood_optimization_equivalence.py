from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
import bijux_phylogenetics.phylo.likelihood.fixed_topology_branch_lengths as branch_optimization
import bijux_phylogenetics.phylo.likelihood.fixed_topology_nucleotide_joint_optimization as joint_optimization
import bijux_phylogenetics.phylo.likelihood.optimization_likelihood_equivalence as optimization_equivalence
import bijux_phylogenetics.phylo.likelihood.substitution_parameters as substitution_optimization

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_fixed_topology_branch_optimization_equivalence_replays_custom_root_prior() -> (
    None
):
    alignment_path = fixture("alignments", "k80_likelihood_alignment_2_taxa.fasta")
    report = branch_optimization.optimize_fixed_topology_nucleotide_branch_lengths_from_alignment(
        fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"),
        alignment_path,
        model_name="k80",
        kappa=4.0,
        root_prior_policy="provided",
        root_prior={"A": 0.7, "C": 0.1, "G": 0.1, "T": 0.1},
    )

    equivalence_report = (
        optimization_equivalence.check_nucleotide_likelihood_optimization_equivalence_from_alignment(
            report,
            alignment_path,
        )
    )

    assert report.root_prior_source == "provided"
    assert report.root_prior_values == [0.7, 0.1, 0.1, 0.1]
    assert (
        equivalence_report.optimization_surface
        == "fixed-topology-nucleotide-branch-length-optimization"
    )
    assert equivalence_report.model_name == "K80"
    assert equivalence_report.parameter_values == {"kappa": 4.0}
    assert equivalence_report.root_prior_source == "provided"
    assert equivalence_report.equivalent is True
    assert math.isclose(
        equivalence_report.independently_rescored_log_likelihood,
        report.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_fixed_topology_single_branch_optimization_equivalence_preserves_selected_branch_result() -> (
    None
):
    tree = load_tree(fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"))
    branch_ids_by_name = {
        child.name or "": child.node_id or ""
        for _parent, child in tree.iter_edges()
    }
    alignment_path = fixture("alignments", "f81_likelihood_alignment_2_taxa.fasta")
    report = branch_optimization.optimize_fixed_topology_nucleotide_single_branch_length_from_alignment(
        fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"),
        alignment_path,
        model_name="f81",
        branch_id=branch_ids_by_name["A"],
        base_frequencies={"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
        root_prior_policy="provided",
        root_prior={"A": 0.55, "C": 0.15, "G": 0.1, "T": 0.2},
    )

    equivalence_report = (
        optimization_equivalence.check_nucleotide_likelihood_optimization_equivalence_from_alignment(
            report,
            alignment_path,
        )
    )

    assert report.root_prior_source == "provided"
    assert report.root_prior_values == [0.55, 0.15, 0.1, 0.2]
    assert (
        equivalence_report.optimization_surface
        == "fixed-topology-nucleotide-single-branch-optimization"
    )
    assert equivalence_report.model_name == "F81"
    assert equivalence_report.parameter_values == {
        "A": 0.4,
        "C": 0.1,
        "G": 0.2,
        "T": 0.3,
    }
    assert equivalence_report.root_prior_source == "provided"
    assert equivalence_report.equivalent is True
    assert math.isclose(
        equivalence_report.independently_rescored_log_likelihood,
        report.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_substitution_parameter_optimization_equivalence_replays_optimized_model_state() -> (
    None
):
    alignment_path = fixture("alignments", "gtr_exchangeability_optimization_alignment_2_taxa.fasta")
    report = substitution_optimization.optimize_nucleotide_substitution_parameters_from_alignment(
        fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"),
        alignment_path,
        model_name="gtr",
        max_coordinate_passes=12,
    )

    equivalence_report = (
        optimization_equivalence.check_nucleotide_likelihood_optimization_equivalence_from_alignment(
            report,
            alignment_path,
        )
    )

    assert (
        equivalence_report.optimization_surface
        == "nucleotide-substitution-parameter-optimization"
    )
    assert equivalence_report.model_name == "GTR"
    assert equivalence_report.parameter_values["AC"] == 1.0
    assert equivalence_report.parameter_values["AG"] == report.parameter_rows[0].optimized_value
    assert equivalence_report.equivalent is True
    assert math.isclose(
        equivalence_report.independently_rescored_log_likelihood,
        report.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_joint_optimization_equivalence_replays_final_branch_and_model_state() -> None:
    alignment_path = fixture("alignments", "gtr_exchangeability_optimization_alignment_2_taxa.fasta")
    report = joint_optimization.optimize_fixed_topology_nucleotide_branches_and_model_from_alignment(
        fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"),
        alignment_path,
        model_name="gtr",
        max_joint_passes=4,
        max_model_coordinate_passes=12,
    )

    equivalence_report = (
        optimization_equivalence.check_nucleotide_likelihood_optimization_equivalence_from_alignment(
            report,
            alignment_path,
        )
    )

    assert equivalence_report.optimization_surface == "fixed-topology-nucleotide-joint-optimization"
    assert equivalence_report.model_name == "GTR"
    assert equivalence_report.parameter_values["AC"] == 1.0
    assert equivalence_report.root_prior_source == report.root_prior_source
    assert equivalence_report.equivalent is True
    assert math.isclose(
        equivalence_report.independently_rescored_log_likelihood,
        report.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_joint_restart_equivalence_replays_selected_attempt_result() -> None:
    alignment_path = fixture("alignments", "gtr_exchangeability_optimization_alignment_2_taxa.fasta")
    report = joint_optimization.optimize_fixed_topology_nucleotide_branches_and_model_with_restarts_from_alignment(
        fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"),
        alignment_path,
        model_name="gtr",
        restart_policy="restart-on-nonconvergence-or-boundary",
        max_restart_count=2,
        max_joint_passes=1,
        max_model_coordinate_passes=1,
        upper_branch_length_bound=1.0,
    )

    equivalence_report = (
        optimization_equivalence.check_nucleotide_likelihood_optimization_equivalence_from_alignment(
            report,
            alignment_path,
        )
    )

    assert (
        equivalence_report.optimization_surface
        == "fixed-topology-nucleotide-joint-optimization-restarts"
    )
    assert equivalence_report.model_name == report.selected_report.model_name
    assert equivalence_report.equivalent is True
    assert math.isclose(
        equivalence_report.independently_rescored_log_likelihood,
        report.selected_report.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_optimization_equivalence_reports_mismatched_stored_likelihood() -> None:
    alignment_path = fixture("alignments", "k80_likelihood_alignment_2_taxa.fasta")
    report = branch_optimization.optimize_fixed_topology_nucleotide_branch_lengths_from_alignment(
        fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"),
        alignment_path,
        model_name="k80",
        kappa=4.0,
    )
    report.optimized_log_likelihood += 0.5

    equivalence_report = (
        optimization_equivalence.check_nucleotide_likelihood_optimization_equivalence_from_alignment(
            report,
            alignment_path,
        )
    )

    assert equivalence_report.equivalent is False
    assert math.isclose(equivalence_report.absolute_difference, 0.5, rel_tol=0.0, abs_tol=1e-12)


def test_optimization_equivalence_tolerance_validation_rejects_negative_values() -> None:
    try:
        optimization_equivalence.validate_nucleotide_likelihood_optimization_equivalence_tolerances(
            absolute_tolerance=-1e-9,
            relative_tolerance=1e-12,
        )
    except ValueError as error:
        assert "absolute_tolerance" in str(error)
    else:
        raise AssertionError("negative absolute tolerance should fail")
