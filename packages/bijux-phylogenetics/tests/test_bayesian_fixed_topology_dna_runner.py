from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    build_fixed_topology_dna_model_definition,
    build_fixed_topology_dna_proposal_schedule,
    run_fixed_topology_dna_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
)
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_fixed_topology_dna_runner_samples_hky85_branch_lengths_and_model_parameters() -> (
    None
):
    model_definition = build_fixed_topology_dna_model_definition(
        substitution_model_name="HKY85",
        branch_length_prior=build_exponential_branch_length_prior(rate=4.0),
        substitution_parameter_prior_bundle=build_substitution_parameter_prior_bundle(
            kappa_prior=build_exponential_positive_substitution_parameter_prior(
                rate=1.5
            ),
            base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
                expected_component_names=("A", "C", "G", "T"),
                concentration_parameters={"A": 2.0, "C": 2.0, "G": 2.0, "T": 2.0},
            ),
        ),
    )
    proposal_schedule = build_fixed_topology_dna_proposal_schedule(
        model_definition=model_definition,
        branch_length_move_weight=1.0,
        branch_length_log_scale_standard_deviation=0.25,
        kappa_move_weight=1.0,
        kappa_log_scale_standard_deviation=0.35,
        base_frequency_move_weight=1.0,
        base_frequency_coordinate_standard_deviation=0.45,
    )

    report = run_fixed_topology_dna_metropolis_hastings(
        tree=load_tree(fixture("trees", "hky85_likelihood_tree_2_taxa.nwk")),
        records=load_fasta_alignment(
            fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta")
        ),
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        iteration_count=9,
        sample_every=1,
        seed=0,
    )

    assert report.chain_report.accepted_count >= 1
    assert len(report.posterior_rows) == len(report.chain_report.sampled_states)
    assert all(row.substitution_model_name == "HKY85" for row in report.posterior_rows)
    assert len({row.topology_id for row in report.posterior_rows}) == 1
    assert all(
        "branch-lengths" in row.prior_component_log_priors
        for row in report.posterior_rows
    )
    assert all(
        "substitution:kappa" in row.prior_component_log_priors
        for row in report.posterior_rows
    )
    assert all(
        "substitution:base-frequencies" in row.prior_component_log_priors
        for row in report.posterior_rows
    )
    assert all(row.branch_lengths for row in report.posterior_rows)
    assert all("kappa" in row.scalar_parameters for row in report.posterior_rows)
    assert all(
        "base-frequencies" in row.vector_parameters for row in report.posterior_rows
    )
    assert any(
        changed_field == "scalar_parameters.kappa"
        for step_row in report.chain_report.step_rows
        for changed_field in step_row.proposal_changed_fields
    )
    assert any(
        changed_field.startswith("tree.branch_length:")
        for step_row in report.chain_report.step_rows
        for changed_field in step_row.proposal_changed_fields
    )
    assert any(
        changed_field.startswith("vector_parameters.base-frequencies.")
        for step_row in report.chain_report.step_rows
        for changed_field in step_row.proposal_changed_fields
    )
    assert all(
        state.model_parameters.categorical_parameters["substitution-model"] == "HKY85"
        for state in report.chain_report.sampled_states
    )
