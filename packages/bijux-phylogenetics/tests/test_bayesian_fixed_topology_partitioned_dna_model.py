from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
    build_fixed_branch_length_prior,
)
from bijux_phylogenetics.bayesian.fixed_topology_partitioned_dna import (
    FixedTopologyPartitionedDnaModelDefinition,
    FixedTopologyPartitionedDnaProposalSchedule,
    build_fixed_topology_partitioned_dna_model_definition,
    build_fixed_topology_partitioned_dna_proposal_schedule,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    build_partition_model_prior_bundle,
    build_partition_parameter_linkage_plan,
    build_partition_substitution_model_definition,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
)
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    LocusSegment,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def _mixed_partition_prior_bundle():
    return build_substitution_parameter_prior_bundle(
        kappa_prior=build_exponential_positive_substitution_parameter_prior(rate=1.5),
        exchangeability_prior=build_dirichlet_simplex_substitution_parameter_prior(
            expected_component_names=("AC", "AG", "AT", "CG", "CT", "GT"),
            concentration_parameters={
                "AC": 2.0,
                "AG": 2.0,
                "AT": 2.0,
                "CG": 2.0,
                "CT": 2.0,
                "GT": 2.0,
            },
        ),
        base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
            expected_component_names=("A", "C", "G", "T"),
            concentration_parameters={"A": 2.0, "C": 2.0, "G": 2.0, "T": 2.0},
        ),
    )


def _locus_partition(name: str, start: int, end: int) -> LocusPartition:
    return LocusPartition(
        name=name,
        segments=(LocusSegment(start=start, end=end),),
        total_sites=(end - start) + 1,
        data_type="DNA",
    )


def test_build_fixed_topology_partitioned_dna_model_definition_tracks_targets_and_linkage() -> (
    None
):
    model_definition = build_fixed_topology_partitioned_dna_model_definition(
        locus_partitions=(
            _locus_partition("gene_alpha", 1, 6),
            _locus_partition("gene_beta", 7, 12),
        ),
        partition_prior_bundle=build_partition_model_prior_bundle(
            partition_models=(
                build_partition_substitution_model_definition(
                    partition_name="gene_alpha",
                    model_name="HKY85",
                ),
                build_partition_substitution_model_definition(
                    partition_name="gene_beta",
                    model_name="GTR",
                ),
            ),
            linkage_plan=build_partition_parameter_linkage_plan(
                partition_names=("gene_alpha", "gene_beta"),
                linkage_policies={
                    "kappa": "unlinked",
                    "base-frequencies": "linked",
                    "exchangeabilities": "unlinked",
                },
            ),
            substitution_prior_bundle=_mixed_partition_prior_bundle(),
        ),
        branch_length_prior=build_exponential_branch_length_prior(rate=4.0),
    )

    assert isinstance(model_definition, FixedTopologyPartitionedDnaModelDefinition)
    assert tuple(partition.name for partition in model_definition.locus_partitions) == (
        "gene_alpha",
        "gene_beta",
    )
    assert tuple(model.model_name for model in model_definition.partition_models) == (
        "HKY85",
        "GTR",
    )
    assert model_definition.active_parameter_targets == (
        "base-frequencies",
        "exchangeabilities",
        "kappa",
    )
    assert model_definition.linkage_eligible_targets == ("base-frequencies",)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        (
            {
                "locus_partitions": (
                    _locus_partition("gene_alpha", 1, 6),
                    _locus_partition("gene_beta", 7, 12),
                ),
                "partition_prior_bundle": build_partition_model_prior_bundle(
                    partition_models=(
                        build_partition_substitution_model_definition(
                            partition_name="gene_alpha",
                            model_name="HKY85+G",
                        ),
                        build_partition_substitution_model_definition(
                            partition_name="gene_beta",
                            model_name="HKY85",
                        ),
                    ),
                    linkage_plan=build_partition_parameter_linkage_plan(
                        partition_names=("gene_alpha", "gene_beta"),
                        linkage_policies={
                            "kappa": "linked",
                            "base-frequencies": "linked",
                            "gamma-alpha": "linked",
                        },
                    ),
                    substitution_prior_bundle=build_substitution_parameter_prior_bundle(
                        kappa_prior=build_exponential_positive_substitution_parameter_prior(
                            rate=1.5
                        ),
                        base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
                            expected_component_names=("A", "C", "G", "T"),
                            concentration_parameters={
                                "A": 2.0,
                                "C": 2.0,
                                "G": 2.0,
                                "T": 2.0,
                            },
                        ),
                        gamma_alpha_prior=build_exponential_positive_substitution_parameter_prior(
                            rate=0.8
                        ),
                    ),
                ),
                "branch_length_prior": build_exponential_branch_length_prior(rate=4.0),
            },
            "without gamma-rate or invariant-site modifiers",
        ),
        (
            {
                "locus_partitions": (
                    _locus_partition("gene_alpha", 1, 6),
                    _locus_partition("gene_beta", 7, 12),
                ),
                "partition_prior_bundle": build_partition_model_prior_bundle(
                    partition_models=(
                        build_partition_substitution_model_definition(
                            partition_name="gene_alpha",
                            model_name="HKY85",
                        ),
                        build_partition_substitution_model_definition(
                            partition_name="gene_beta",
                            model_name="HKY85",
                        ),
                    ),
                    linkage_plan=build_partition_parameter_linkage_plan(
                        partition_names=("gene_alpha", "gene_beta"),
                        linkage_policies={
                            "kappa": "linked",
                            "base-frequencies": "linked",
                        },
                    ),
                    substitution_prior_bundle=build_substitution_parameter_prior_bundle(
                        kappa_prior=build_exponential_positive_substitution_parameter_prior(
                            rate=1.5
                        ),
                        base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
                            expected_component_names=("A", "C", "G", "T"),
                            concentration_parameters={
                                "A": 2.0,
                                "C": 2.0,
                                "G": 2.0,
                                "T": 2.0,
                            },
                        ),
                    ),
                ),
                "branch_length_prior": build_fixed_branch_length_prior(fixed_value=0.1),
            },
            "does not support fixed branch-length priors",
        ),
    ],
)
def test_build_fixed_topology_partitioned_dna_model_definition_rejects_unsupported_inputs(
    kwargs: dict[str, object],
    match: str,
) -> None:
    with pytest.raises(PhylogeneticsError, match=match):
        build_fixed_topology_partitioned_dna_model_definition(**kwargs)


def test_build_fixed_topology_partitioned_dna_proposal_schedule_requires_active_moves() -> (
    None
):
    model_definition = build_fixed_topology_partitioned_dna_model_definition(
        locus_partitions=(
            _locus_partition("gene_alpha", 1, 6),
            _locus_partition("gene_beta", 7, 12),
        ),
        partition_prior_bundle=build_partition_model_prior_bundle(
            partition_models=(
                build_partition_substitution_model_definition(
                    partition_name="gene_alpha",
                    model_name="HKY85",
                ),
                build_partition_substitution_model_definition(
                    partition_name="gene_beta",
                    model_name="HKY85",
                ),
            ),
            linkage_plan=build_partition_parameter_linkage_plan(
                partition_names=("gene_alpha", "gene_beta"),
                linkage_policies={
                    "kappa": "linked",
                    "base-frequencies": "linked",
                },
            ),
            substitution_prior_bundle=build_substitution_parameter_prior_bundle(
                kappa_prior=build_exponential_positive_substitution_parameter_prior(
                    rate=1.5
                ),
                base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
                    expected_component_names=("A", "C", "G", "T"),
                    concentration_parameters={
                        "A": 2.0,
                        "C": 2.0,
                        "G": 2.0,
                        "T": 2.0,
                    },
                ),
            ),
        ),
        branch_length_prior=build_exponential_branch_length_prior(rate=4.0),
    )

    schedule = build_fixed_topology_partitioned_dna_proposal_schedule(
        model_definition=model_definition,
        branch_length_move_weight=1.0,
        branch_length_log_scale_standard_deviation=0.25,
        kappa_move_weight=1.0,
        kappa_log_scale_standard_deviation=0.3,
        base_frequency_move_weight=1.0,
        base_frequency_coordinate_standard_deviation=0.35,
        linkage_move_weight=0.5,
    )

    assert isinstance(schedule, FixedTopologyPartitionedDnaProposalSchedule)
    assert schedule.kappa_move_weight == 1.0
    assert schedule.linkage_move_weight == 0.5

    with pytest.raises(PhylogeneticsError, match="positive move weight"):
        build_fixed_topology_partitioned_dna_proposal_schedule(
            model_definition=model_definition,
            branch_length_move_weight=1.0,
            branch_length_log_scale_standard_deviation=0.25,
            kappa_move_weight=0.0,
            kappa_log_scale_standard_deviation=0.3,
            base_frequency_move_weight=1.0,
            base_frequency_coordinate_standard_deviation=0.35,
        )
