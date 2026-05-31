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
from bijux_phylogenetics.bayesian.metropolis_hastings import MetropolisHastingsRunReport
from bijux_phylogenetics.bayesian.partition_model_priors import (
    build_partition_model_prior_bundle,
    build_partition_parameter_linkage_plan,
    build_partition_substitution_model_definition,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    build_posterior_missing_nucleotide_definition,
    summarize_fixed_topology_dna_posterior_missing_states,
    summarize_fixed_topology_partitioned_dna_posterior_missing_states,
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
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import parse_locus_partitions
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_posterior_missing_nucleotide_definition_rejects_dna_state_as_missing_symbol() -> (
    None
):
    with pytest.raises(PhylogeneticsError, match="distinct from A, C, G, and T"):
        build_posterior_missing_nucleotide_definition(
            records=[AlignmentRecord(identifier="A", sequence="AN")],
            missing_state_symbols=("A", "N"),
        )


def test_fixed_topology_dna_missing_state_summary_recovers_masked_known_state() -> None:
    report = summarize_fixed_topology_dna_posterior_missing_states(
        run_report=_build_manual_fixed_topology_dna_run_report(),
        definition=build_posterior_missing_nucleotide_definition(
            records=[
                AlignmentRecord(identifier="A", sequence="A"),
                AlignmentRecord(identifier="B", sequence="N"),
            ],
            posterior_probability_threshold=0.5,
        ),
    )

    assert report.sample_count == 1
    assert report.masked_site_count == 1
    assert report.sampled_substitution_models == ["HKY85"]
    assert report.warnings == []

    summary_row = report.site_summary_rows[0]
    assert summary_row.taxon == "B"
    assert summary_row.site_position == 1
    assert summary_row.observed_symbol == "N"
    assert summary_row.consensus_state == "A"
    assert summary_row.exported_state == "A"
    assert summary_row.max_posterior_probability > 0.5
    assert summary_row.max_posterior_probability < 1.0
    assert report.sequence_records == [
        type(report.sequence_records[0])(
            identifier="B",
            masked_site_count=1,
            sequence="A",
        )
    ]

    probability_by_state = {
        row.state: row.posterior_probability for row in report.state_probability_rows
    }
    assert probability_by_state["A"] > probability_by_state["C"]
    assert probability_by_state["A"] > probability_by_state["G"]
    assert probability_by_state["A"] > probability_by_state["T"]


def test_partitioned_dna_missing_state_summary_preserves_global_site_positions() -> (
    None
):
    locus_partitions = parse_locus_partitions(
        fixture("metadata", "partitioned_dna_partitions_4_taxa.txt")
    )
    observed_records = load_fasta_alignment(
        fixture("alignments", "partitioned_dna_alignment_4_taxa.fasta")
    )
    masked_records = [
        AlignmentRecord(
            identifier=record.identifier,
            sequence=record.sequence[:2]
            + ("N" if record.identifier == "A" else record.sequence[2])
            + record.sequence[3:],
        )
        for record in observed_records
    ]
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
    run_report = run_fixed_topology_partitioned_dna_metropolis_hastings(
        tree=load_tree(fixture("trees", "partitioned_dna_tree_4_taxa.nwk")),
        records=observed_records,
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        iteration_count=8,
        sample_every=1,
        seed=4,
    )

    report = summarize_fixed_topology_partitioned_dna_posterior_missing_states(
        run_report=run_report,
        definition=build_posterior_missing_nucleotide_definition(
            records=masked_records,
            locus_partitions=locus_partitions,
            partition_models=tuple(
                build_partition_substitution_model_definition(
                    partition_name=partition.name,
                    model_name="K80",
                )
                for partition in locus_partitions
            ),
        ),
    )

    assert report.masked_site_count == 1
    assert report.site_summary_rows[0].taxon == "A"
    assert report.site_summary_rows[0].site_position == 3
    assert len(report.sequence_records[0].sequence) == len(observed_records[0].sequence)


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
    tree = PhyloTree.from_newick("(A:0.05,B:0.05);")
    tree.rooted = True
    sampled_state = build_bayesian_phylogenetic_state(
        tree=tree,
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters={"substitution-model": "HKY85"},
            scalar_parameters={"kappa": 2.0},
            vector_parameters={
                "base-frequencies": {
                    "A": 0.55,
                    "C": 0.15,
                    "G": 0.15,
                    "T": 0.15,
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
            initial_state=sampled_state,
            final_state=sampled_state,
            sampled_states=[sampled_state],
            step_rows=[],
        ),
        posterior_rows=[],
    )
