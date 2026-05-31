from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
)
from bijux_phylogenetics.bayesian.fixed_topology_partitioned_dna import (
    build_fixed_topology_partitioned_dna_model_definition,
    build_fixed_topology_partitioned_dna_proposal_schedule,
    run_fixed_topology_partitioned_dna_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    build_partition_model_prior_bundle,
    build_partition_parameter_linkage_plan,
    build_partition_substitution_model_definition,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_exponential_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
)
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.partitions import parse_locus_partitions

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_fixed_topology_partitioned_dna_runner_samples_partition_likelihoods_and_linkage_states() -> (
    None
):
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

    report = run_fixed_topology_partitioned_dna_metropolis_hastings(
        tree=load_tree(fixture("trees", "partitioned_dna_tree_4_taxa.nwk")),
        records=load_fasta_alignment(
            fixture("alignments", "partitioned_dna_alignment_4_taxa.fasta")
        ),
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        iteration_count=20,
        sample_every=1,
        seed=4,
    )

    assert report.chain_report.accepted_count >= 1
    assert len(report.posterior_rows) == len(report.chain_report.sampled_states)
    assert all(len(row.partition_rows) == 2 for row in report.posterior_rows)
    assert all(
        "branch-lengths" in row.prior_component_log_priors
        for row in report.posterior_rows
    )
    assert all(
        math.isclose(
            sum(partition_row.log_likelihood for partition_row in row.partition_rows),
            row.log_likelihood,
            rel_tol=0.0,
            abs_tol=1e-9,
        )
        for row in report.posterior_rows
    )
    assert all(
        [partition_row.partition_name for partition_row in row.partition_rows]
        == ["gene_alpha", "gene_beta"]
        for row in report.posterior_rows
    )
    assert all(
        all(
            "kappa" in partition_row.scalar_parameters
            for partition_row in row.partition_rows
        )
        for row in report.posterior_rows
    )
    assert any(
        changed_field.startswith("categorical_parameters.partition-linkage:kappa:")
        for step_row in report.chain_report.step_rows
        for changed_field in step_row.proposal_changed_fields
    )
    linkage_layouts = {
        tuple(
            partition_row.linkage_groups["kappa"]
            for partition_row in row.partition_rows
        )
        for row in report.posterior_rows
    }
    assert ("kappa-shared", "kappa-shared") in linkage_layouts
    assert ("gene_alpha", "gene_beta") in linkage_layouts
