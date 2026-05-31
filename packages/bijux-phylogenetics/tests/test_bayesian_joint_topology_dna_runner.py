from __future__ import annotations

from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    build_fixed_topology_dna_model_definition,
    build_fixed_topology_dna_proposal_schedule,
)
from bijux_phylogenetics.bayesian.joint_topology_dna import (
    build_joint_topology_dna_model_definition,
    build_joint_topology_dna_proposal_schedule,
    run_joint_topology_dna_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_exponential_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    build_uniform_rooted_tree_topology_prior,
)
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


def test_joint_topology_dna_runner_samples_multiple_topologies_and_model_parameters() -> (
    None
):
    sequence_model_definition = build_fixed_topology_dna_model_definition(
        substitution_model_name="K80",
        branch_length_prior=build_exponential_branch_length_prior(rate=4.0),
        substitution_parameter_prior_bundle=build_substitution_parameter_prior_bundle(
            kappa_prior=build_exponential_positive_substitution_parameter_prior(
                rate=1.5
            )
        ),
    )
    joint_model_definition = build_joint_topology_dna_model_definition(
        sequence_model_definition=sequence_model_definition,
        topology_prior=build_uniform_rooted_tree_topology_prior(["A", "B", "C", "D"]),
    )
    joint_proposal_schedule = build_joint_topology_dna_proposal_schedule(
        sequence_proposal_schedule=build_fixed_topology_dna_proposal_schedule(
            model_definition=sequence_model_definition,
            branch_length_move_weight=1.0,
            branch_length_log_scale_standard_deviation=0.25,
            kappa_move_weight=1.0,
            kappa_log_scale_standard_deviation=0.35,
        ),
        nni_move_weight=1.0,
    )

    report = run_joint_topology_dna_metropolis_hastings(
        tree=_build_ambiguous_start_tree(),
        records=_build_ambiguous_alignment_records(),
        model_definition=joint_model_definition,
        proposal_schedule=joint_proposal_schedule,
        iteration_count=20,
        sample_every=1,
        seed=0,
    )

    assert report.chain_report.accepted_count >= 1
    assert report.distinct_topology_count > 1
    assert len(report.distinct_topology_ids) == report.distinct_topology_count
    assert len(report.posterior_rows) == len(report.chain_report.sampled_states)
    assert all(row.substitution_model_name == "K80" for row in report.posterior_rows)
    assert all(
        "tree-topology" in row.prior_component_log_priors
        for row in report.posterior_rows
    )
    assert all(
        "branch-lengths" in row.prior_component_log_priors
        for row in report.posterior_rows
    )
    assert all(
        "substitution:kappa" in row.prior_component_log_priors
        for row in report.posterior_rows
    )
    assert all(row.branch_lengths for row in report.posterior_rows)
    assert all("kappa" in row.scalar_parameters for row in report.posterior_rows)
    assert all(
        state.model_parameters.categorical_parameters["substitution-model"] == "K80"
        for state in report.chain_report.sampled_states
    )
    assert any(
        changed_field.startswith("tree.topology:nni:")
        for step_row in report.chain_report.step_rows
        for changed_field in step_row.proposal_changed_fields
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


def _build_ambiguous_start_tree() -> PhyloTree:
    tree = PhyloTree.from_newick("(((A:0.1,B:0.1):0.1,C:0.1):0.1,D:0.1);")
    tree.rooted = True
    return tree


def _build_ambiguous_alignment_records() -> list[AlignmentRecord]:
    return [
        AlignmentRecord(identifier="A", sequence="A"),
        AlignmentRecord(identifier="B", sequence="C"),
        AlignmentRecord(identifier="C", sequence="G"),
        AlignmentRecord(identifier="D", sequence="T"),
    ]
