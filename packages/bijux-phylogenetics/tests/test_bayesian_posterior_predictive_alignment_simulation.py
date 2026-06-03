from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    FixedTopologyDnaRunReport,
    build_fixed_topology_dna_model_definition,
    build_fixed_topology_dna_proposal_schedule,
)
from bijux_phylogenetics.bayesian.fixed_topology_partitioned_dna import (
    build_fixed_topology_partitioned_dna_model_definition,
    build_fixed_topology_partitioned_dna_proposal_schedule,
    run_fixed_topology_partitioned_dna_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.joint_topology_dna import (
    build_joint_topology_dna_model_definition,
    build_joint_topology_dna_proposal_schedule,
    run_joint_topology_dna_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import MetropolisHastingsRunReport
from bijux_phylogenetics.bayesian.partition_model_priors import (
    build_partition_model_prior_bundle,
    build_partition_parameter_linkage_plan,
    build_partition_substitution_model_definition,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    build_posterior_predictive_simulation_definition,
    simulate_fixed_topology_dna_posterior_predictive,
    simulate_fixed_topology_partitioned_dna_posterior_predictive,
    simulate_joint_topology_dna_posterior_predictive,
)
from bijux_phylogenetics.bayesian.state import (
    build_bayesian_model_parameter_state,
    build_bayesian_phylogenetic_state,
    build_bayesian_prior_component_state,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    build_uniform_rooted_tree_topology_prior,
)
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import parse_locus_partitions
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_posterior_predictive_definition_rejects_unsupported_policy() -> None:
    with pytest.raises(
        PhylogeneticsError, match="supported posterior-sample selection"
    ):
        build_posterior_predictive_simulation_definition(
            replicate_count=2,
            sample_selection_policy="without-replacement",
        )


def test_fixed_topology_dna_posterior_predictive_uses_selected_posterior_sample_parameters() -> (
    None
):
    run_report = _build_manual_fixed_topology_dna_run_report()
    observed_records = [
        AlignmentRecord(identifier="A", sequence="A" * 400),
        AlignmentRecord(identifier="B", sequence="A" * 400),
    ]

    at_heavy = simulate_fixed_topology_dna_posterior_predictive(
        run_report=run_report,
        records=observed_records,
        definition=build_posterior_predictive_simulation_definition(
            replicate_count=1,
            seed=1,
        ),
    )
    gc_heavy = simulate_fixed_topology_dna_posterior_predictive(
        run_report=run_report,
        records=observed_records,
        definition=build_posterior_predictive_simulation_definition(
            replicate_count=1,
            seed=0,
        ),
    )

    assert at_heavy.replicates[0].posterior_sample_index == 0
    assert gc_heavy.replicates[0].posterior_sample_index == 1

    at_gc_fraction = next(
        row.value
        for row in at_heavy.replicate_statistic_rows
        if row.statistic_name == "gc-fraction"
    )
    gc_gc_fraction = next(
        row.value
        for row in gc_heavy.replicate_statistic_rows
        if row.statistic_name == "gc-fraction"
    )

    assert at_gc_fraction < 0.1
    assert gc_gc_fraction > 0.9


def test_joint_topology_dna_posterior_predictive_emits_observed_and_replicate_statistics() -> (
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
    run_report = run_joint_topology_dna_metropolis_hastings(
        tree=_build_ambiguous_start_tree(),
        records=_build_ambiguous_alignment_records(),
        model_definition=build_joint_topology_dna_model_definition(
            sequence_model_definition=sequence_model_definition,
            topology_prior=build_uniform_rooted_tree_topology_prior(
                ["A", "B", "C", "D"]
            ),
        ),
        proposal_schedule=build_joint_topology_dna_proposal_schedule(
            sequence_proposal_schedule=build_fixed_topology_dna_proposal_schedule(
                model_definition=sequence_model_definition,
                branch_length_move_weight=1.0,
                branch_length_log_scale_standard_deviation=0.25,
                kappa_move_weight=1.0,
                kappa_log_scale_standard_deviation=0.35,
            ),
            nni_move_weight=1.0,
        ),
        iteration_count=12,
        sample_every=1,
        seed=0,
    )

    report = simulate_joint_topology_dna_posterior_predictive(
        run_report=run_report,
        records=_build_ambiguous_alignment_records(),
        definition=build_posterior_predictive_simulation_definition(
            replicate_count=3,
            seed=0,
        ),
    )

    assert report.model_name == "K80"
    assert len(report.observed_statistic_rows) == 2
    assert len(report.replicates) == 3
    assert len(report.replicate_statistic_rows) == 6
    assert len(report.statistic_summary_rows) == 2
    assert all(
        replicate.posterior_sample_index < len(run_report.chain_report.sampled_states)
        for replicate in report.replicates
    )


def test_partitioned_dna_posterior_predictive_emits_full_alignment_replicates() -> None:
    locus_partitions = parse_locus_partitions(
        fixture("metadata", "partitioned_dna_partitions_4_taxa.txt")
    )
    model_definition = build_fixed_topology_partitioned_dna_model_definition(
        locus_partitions=locus_partitions,
        partition_prior_bundle=build_partition_model_prior_bundle(
            partition_models=tuple(
                build_partition_substitution_model_definition(
                    partition_name=partition.name,
                    model_name="K80",
                )
                for partition in locus_partitions
            ),
            linkage_plan=build_partition_parameter_linkage_plan(
                partition_names=tuple(partition.name for partition in locus_partitions),
                linkage_policies={"kappa": "linked"},
            ),
            substitution_prior_bundle=build_substitution_parameter_prior_bundle(
                kappa_prior=build_exponential_positive_substitution_parameter_prior(
                    rate=1.2
                )
            ),
        ),
        branch_length_prior=build_exponential_branch_length_prior(rate=4.0),
    )
    proposal_schedule = build_fixed_topology_partitioned_dna_proposal_schedule(
        model_definition=model_definition,
        branch_length_move_weight=1.0,
        branch_length_log_scale_standard_deviation=0.2,
        kappa_move_weight=1.0,
        kappa_log_scale_standard_deviation=0.3,
        linkage_move_weight=2.0,
    )
    observed_records = load_fasta_alignment(
        fixture("alignments", "partitioned_dna_alignment_4_taxa.fasta")
    )
    run_report = run_fixed_topology_partitioned_dna_metropolis_hastings(
        tree=load_tree(fixture("trees", "partitioned_dna_tree_4_taxa.nwk")),
        records=observed_records,
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        iteration_count=12,
        sample_every=1,
        seed=4,
    )

    report = simulate_fixed_topology_partitioned_dna_posterior_predictive(
        run_report=run_report,
        records=observed_records,
        definition=build_posterior_predictive_simulation_definition(
            replicate_count=2,
            seed=0,
        ),
    )

    assert report.model_name == "partitioned-dna"
    assert len(report.replicates) == 2
    assert len(report.replicate_statistic_rows) == 4
    assert all(
        len(replicate.records[0].sequence) == len(observed_records[0].sequence)
        for replicate in report.replicates
    )


def _build_manual_fixed_topology_dna_run_report() -> FixedTopologyDnaRunReport:
    model_definition = build_fixed_topology_dna_model_definition(
        substitution_model_name="HKY85",
        branch_length_prior=build_exponential_branch_length_prior(rate=3.0),
        substitution_parameter_prior_bundle=build_substitution_parameter_prior_bundle(
            kappa_prior=build_exponential_positive_substitution_parameter_prior(
                rate=1.0
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
        branch_length_log_scale_standard_deviation=0.2,
        kappa_move_weight=1.0,
        kappa_log_scale_standard_deviation=0.2,
        base_frequency_move_weight=1.0,
        base_frequency_coordinate_standard_deviation=0.2,
    )
    tree = PhyloTree.from_newick("(A:0.01,B:0.01);")
    tree.rooted = True
    first_state = build_bayesian_phylogenetic_state(
        tree=tree,
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters={"substitution-model": "HKY85"},
            scalar_parameters={"kappa": 2.0},
            vector_parameters={
                "base-frequencies": {
                    "A": 0.49,
                    "C": 0.01,
                    "G": 0.01,
                    "T": 0.49,
                }
            },
        ),
        prior_components=[
            build_bayesian_prior_component_state(
                component_name="branch-lengths",
                family="exponential",
                log_prior=-0.2,
            ),
            build_bayesian_prior_component_state(
                component_name="substitution:kappa",
                family="exponential",
                log_prior=-0.1,
            ),
            build_bayesian_prior_component_state(
                component_name="substitution:base-frequencies",
                family="dirichlet",
                log_prior=-0.1,
            ),
        ],
        log_likelihood=-1.0,
    )
    second_state = build_bayesian_phylogenetic_state(
        tree=tree,
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters={"substitution-model": "HKY85"},
            scalar_parameters={"kappa": 2.0},
            vector_parameters={
                "base-frequencies": {
                    "A": 0.01,
                    "C": 0.49,
                    "G": 0.49,
                    "T": 0.01,
                }
            },
        ),
        prior_components=[
            build_bayesian_prior_component_state(
                component_name="branch-lengths",
                family="exponential",
                log_prior=-0.2,
            ),
            build_bayesian_prior_component_state(
                component_name="substitution:kappa",
                family="exponential",
                log_prior=-0.1,
            ),
            build_bayesian_prior_component_state(
                component_name="substitution:base-frequencies",
                family="dirichlet",
                log_prior=-0.1,
            ),
        ],
        log_likelihood=-1.0,
    )
    return FixedTopologyDnaRunReport(
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        observation_policy="reject",
        chain_report=MetropolisHastingsRunReport(
            iteration_count=1,
            sample_every=1,
            seed=0,
            accepted_count=1,
            rejected_count=0,
            acceptance_rate=1.0,
            initial_state=first_state,
            final_state=second_state,
            sampled_states=[first_state, second_state],
            step_rows=[],
        ),
        posterior_rows=[],
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
