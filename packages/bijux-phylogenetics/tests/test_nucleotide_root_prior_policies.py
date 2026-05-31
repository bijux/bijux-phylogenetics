from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_f81_tree_likelihood,
    evaluate_nucleotide_marginal_ancestral_probabilities_from_alignment,
    evaluate_nucleotide_site_log_likelihoods_from_alignment,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_models import (
    resolve_selected_nucleotide_likelihood_specification,
)

FIXTURES = Path(__file__).parent / "fixtures"

_F81_BASE_FREQUENCIES = {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}
_PROVIDED_ROOT_PRIOR = {"A": 0.7, "C": 0.1, "G": 0.1, "T": 0.1}


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_selected_nucleotide_specification_supports_all_root_prior_policies() -> None:
    records = load_fasta_alignment(
        fixture("alignments", "f81_likelihood_alignment_2_taxa.fasta")
    )

    equal_specification = resolve_selected_nucleotide_likelihood_specification(
        records,
        model_name="f81",
        owner_name="F81 selected likelihood policy test",
        base_frequencies=_F81_BASE_FREQUENCIES,
        root_prior_policy="equal",
    )
    empirical_specification = resolve_selected_nucleotide_likelihood_specification(
        records,
        model_name="f81",
        owner_name="F81 selected likelihood policy test",
        base_frequencies=_F81_BASE_FREQUENCIES,
        root_prior_policy="empirical",
    )
    stationary_specification = resolve_selected_nucleotide_likelihood_specification(
        records,
        model_name="f81",
        owner_name="F81 selected likelihood policy test",
        base_frequencies=_F81_BASE_FREQUENCIES,
        root_prior_policy="stationary",
    )
    fixed_specification = resolve_selected_nucleotide_likelihood_specification(
        records,
        model_name="f81",
        owner_name="F81 selected likelihood policy test",
        base_frequencies=_F81_BASE_FREQUENCIES,
        root_prior_policy="fixed-state",
        fixed_root_state="A",
    )
    provided_specification = resolve_selected_nucleotide_likelihood_specification(
        records,
        model_name="f81",
        owner_name="F81 selected likelihood policy test",
        base_frequencies=_F81_BASE_FREQUENCIES,
        root_prior_policy="provided",
        root_prior=_PROVIDED_ROOT_PRIOR,
    )

    assert equal_specification.root_prior_source == "equal"
    assert empirical_specification.root_prior_source == "empirical"
    assert stationary_specification.root_prior_source == "stationary"
    assert fixed_specification.root_prior_source == "fixed-state"
    assert provided_specification.root_prior_source == "provided"

    assert numpy.allclose(
        equal_specification.root_prior,
        numpy.array([0.25, 0.25, 0.25, 0.25], dtype=float),
    )
    assert numpy.allclose(
        empirical_specification.root_prior,
        numpy.array([0.375, 0.125, 0.125, 0.375], dtype=float),
    )
    assert numpy.allclose(
        stationary_specification.root_prior,
        numpy.array([0.4, 0.1, 0.2, 0.3], dtype=float),
    )
    assert numpy.allclose(
        fixed_specification.root_prior,
        numpy.array([1.0, 0.0, 0.0, 0.0], dtype=float),
    )
    assert numpy.allclose(
        provided_specification.root_prior,
        numpy.array([0.7, 0.1, 0.1, 0.1], dtype=float),
    )


def test_f81_fixed_tree_likelihood_changes_with_root_prior_policy() -> None:
    tree = load_tree(fixture("trees", "f81_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "f81_likelihood_alignment_2_taxa.fasta")
    )

    default_report = evaluate_f81_tree_likelihood(
        tree,
        records,
        base_frequencies=_F81_BASE_FREQUENCIES,
    )
    equal_report = evaluate_f81_tree_likelihood(
        tree,
        records,
        base_frequencies=_F81_BASE_FREQUENCIES,
        root_prior_policy="equal",
    )
    empirical_report = evaluate_f81_tree_likelihood(
        tree,
        records,
        base_frequencies=_F81_BASE_FREQUENCIES,
        root_prior_policy="empirical",
    )
    stationary_report = evaluate_f81_tree_likelihood(
        tree,
        records,
        base_frequencies=_F81_BASE_FREQUENCIES,
        root_prior_policy="stationary",
    )
    fixed_report = evaluate_f81_tree_likelihood(
        tree,
        records,
        base_frequencies=_F81_BASE_FREQUENCIES,
        root_prior_policy="fixed-state",
        fixed_root_state="A",
    )
    provided_report = evaluate_f81_tree_likelihood(
        tree,
        records,
        base_frequencies=_F81_BASE_FREQUENCIES,
        root_prior_policy="provided",
        root_prior=_PROVIDED_ROOT_PRIOR,
    )

    assert math.isclose(
        default_report.log_likelihood,
        stationary_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        equal_report.log_likelihood,
        stationary_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        empirical_report.log_likelihood,
        stationary_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        fixed_report.log_likelihood,
        stationary_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        provided_report.log_likelihood,
        stationary_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_site_log_likelihood_export_matches_direct_f81_total_for_provided_root_prior() -> (
    None
):
    tree_path = fixture("trees", "f81_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture("alignments", "f81_likelihood_alignment_2_taxa.fasta")
    direct_report = evaluate_f81_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        base_frequencies=_F81_BASE_FREQUENCIES,
        root_prior_policy="provided",
        root_prior=_PROVIDED_ROOT_PRIOR,
    )
    site_report = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        tree_path,
        alignment_path,
        model_name="f81",
        base_frequencies=_F81_BASE_FREQUENCIES,
        root_prior_policy="provided",
        root_prior=_PROVIDED_ROOT_PRIOR,
    )

    assert math.isclose(
        site_report.log_likelihood,
        direct_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        sum(row.log_likelihood for row in site_report.site_log_likelihoods),
        site_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_marginal_ancestral_probabilities_respect_fixed_root_state() -> None:
    report = evaluate_nucleotide_marginal_ancestral_probabilities_from_alignment(
        fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta"),
        model_name="jc69",
        root_prior_policy="fixed-state",
        fixed_root_state="A",
    )

    root_node_id = report.posterior_rows[0].node_id
    site_one_rows = [
        row
        for row in report.posterior_rows
        if row.node_id == root_node_id and row.site_position == 1
    ]
    posterior_by_state = {row.state: row.posterior_probability for row in site_one_rows}

    assert posterior_by_state["A"] == 1.0
    assert posterior_by_state["C"] == 0.0
    assert posterior_by_state["G"] == 0.0
    assert posterior_by_state["T"] == 0.0
