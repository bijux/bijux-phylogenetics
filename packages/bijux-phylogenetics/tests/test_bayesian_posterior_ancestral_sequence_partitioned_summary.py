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
from bijux_phylogenetics.bayesian.posterior_ancestral_sequences import (
    build_posterior_ancestral_sequence_definition,
    summarize_nucleotide_posterior_ancestral_sequences,
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


def test_posterior_ancestral_sequence_summary_preserves_partitioned_global_site_coordinates() -> (
    None
):
    locus_partitions = parse_locus_partitions(
        fixture("metadata", "partitioned_dna_partitions_4_taxa.txt")
    )
    partition_models = tuple(
        build_partition_substitution_model_definition(
            partition_name=partition.name,
            model_name="K80",
        )
        for partition in locus_partitions
    )
    model_definition = build_fixed_topology_partitioned_dna_model_definition(
        locus_partitions=locus_partitions,
        partition_prior_bundle=build_partition_model_prior_bundle(
            partition_models=partition_models,
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
    records = load_fasta_alignment(
        fixture("alignments", "partitioned_dna_alignment_4_taxa.fasta")
    )
    run_report = run_fixed_topology_partitioned_dna_metropolis_hastings(
        tree=load_tree(fixture("trees", "partitioned_dna_tree_4_taxa.nwk")),
        records=records,
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        iteration_count=20,
        sample_every=1,
        seed=4,
    )

    summary = summarize_nucleotide_posterior_ancestral_sequences(
        run_report.chain_report.sampled_states,
        definition=build_posterior_ancestral_sequence_definition(
            records=records,
            posterior_probability_threshold=0.5,
            minimum_clade_posterior_support=0.5,
            locus_partitions=locus_partitions,
            partition_models=partition_models,
        ),
    )

    assert summary.sampled_substitution_models == [
        "partitioned-dna[gene_alpha=K80,gene_beta=K80]"
    ]
    assert summary.distinct_topology_count == 1
    assert summary.sequence_records
    assert all(
        len(record.sequence) == len(records[0].sequence)
        for record in summary.sequence_records
    )
    assert all(
        1 <= row.site_position <= len(records[0].sequence)
        for row in summary.site_summary_rows
    )
    assert any(row.site_position == 1 for row in summary.site_summary_rows)
    assert any(
        row.site_position == len(records[0].sequence)
        for row in summary.site_summary_rows
    )
    _assert_site_probabilities_sum_to_clade_support(summary)


def _assert_site_probabilities_sum_to_clade_support(summary) -> None:
    grouped_probabilities: dict[tuple[str, int], list[float]] = {}
    clade_support_by_key: dict[tuple[str, int], float] = {}
    for row in summary.state_probability_rows:
        key = (row.clade_id, row.site_position)
        grouped_probabilities.setdefault(key, []).append(
            row.marginal_posterior_probability
        )
        clade_support_by_key.setdefault(key, row.clade_posterior_probability)
    for key, probabilities in grouped_probabilities.items():
        assert math.isclose(
            sum(probabilities),
            clade_support_by_key[key],
            rel_tol=0.0,
            abs_tol=1e-9,
        )
