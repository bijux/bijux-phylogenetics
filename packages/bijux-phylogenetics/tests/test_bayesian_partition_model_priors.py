from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.partition_model_priors import (
    PartitionSubstitutionParameterState,
    build_partition_model_prior_bundle,
    build_partition_parameter_linkage_plan,
    build_partition_substitution_model_definition,
    evaluate_partition_model_log_prior,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_beta_probability_substitution_parameter_prior,
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
)
from bijux_phylogenetics.phylo.alignment.partitions import parse_locus_partitions
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def _partition_names() -> tuple[str, ...]:
    return tuple(
        partition.name
        for partition in parse_locus_partitions(
            fixture("example_multilocus_partitions.txt")
        )
    )


def _gtr_gamma_prior_bundle():
    return build_substitution_parameter_prior_bundle(
        exchangeability_prior=build_dirichlet_simplex_substitution_parameter_prior(
            expected_component_names=("AC", "AG", "AT", "CG", "CT", "GT"),
            concentration_parameters=(2.0, 3.0, 4.0, 5.0, 6.0, 7.0),
        ),
        base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
            expected_component_names=("A", "C", "G", "T"),
            concentration_parameters=(4.0, 3.0, 2.0, 5.0),
        ),
        gamma_alpha_prior=build_exponential_positive_substitution_parameter_prior(
            rate=0.75
        ),
    )


def test_partition_model_priors_distinguish_linked_and_unlinked_parameter_terms() -> (
    None
):
    partition_models = tuple(
        build_partition_substitution_model_definition(
            partition_name=partition_name,
            model_name="GTR+G",
        )
        for partition_name in _partition_names()
    )
    linked_bundle = build_partition_model_prior_bundle(
        partition_models=partition_models,
        linkage_plan=build_partition_parameter_linkage_plan(
            partition_names=_partition_names(),
            linkage_policies={
                "exchangeabilities": "linked",
                "base-frequencies": "linked",
                "gamma-alpha": "linked",
            },
        ),
        substitution_prior_bundle=_gtr_gamma_prior_bundle(),
    )
    unlinked_bundle = build_partition_model_prior_bundle(
        partition_models=partition_models,
        linkage_plan=build_partition_parameter_linkage_plan(
            partition_names=_partition_names(),
            linkage_policies={
                "exchangeabilities": "unlinked",
                "base-frequencies": "unlinked",
                "gamma-alpha": "unlinked",
            },
        ),
        substitution_prior_bundle=_gtr_gamma_prior_bundle(),
    )

    linked_report = evaluate_partition_model_log_prior(
        prior_bundle=linked_bundle,
        partition_parameter_states=(
            PartitionSubstitutionParameterState(
                partition_name="gene_alpha",
                exchangeabilities=(1.0, 2.0, 1.0, 3.0, 2.0, 1.0),
                base_frequencies=(0.3, 0.2, 0.1, 0.4),
                gamma_alpha=0.7,
            ),
            PartitionSubstitutionParameterState(
                partition_name="gene_beta",
                exchangeabilities=(1.0, 2.0, 1.0, 3.0, 2.0, 1.0),
                base_frequencies=(0.3, 0.2, 0.1, 0.4),
                gamma_alpha=0.7,
            ),
            PartitionSubstitutionParameterState(
                partition_name="gene_gamma",
                exchangeabilities=(1.0, 2.0, 1.0, 3.0, 2.0, 1.0),
                base_frequencies=(0.3, 0.2, 0.1, 0.4),
                gamma_alpha=0.7,
            ),
        ),
    )
    unlinked_report = evaluate_partition_model_log_prior(
        prior_bundle=unlinked_bundle,
        partition_parameter_states=(
            PartitionSubstitutionParameterState(
                partition_name="gene_alpha",
                exchangeabilities=(1.0, 2.0, 1.0, 3.0, 2.0, 1.0),
                base_frequencies=(0.3, 0.2, 0.1, 0.4),
                gamma_alpha=0.7,
            ),
            PartitionSubstitutionParameterState(
                partition_name="gene_beta",
                exchangeabilities=(2.0, 1.0, 1.0, 2.0, 3.0, 1.0),
                base_frequencies=(0.25, 0.15, 0.25, 0.35),
                gamma_alpha=1.2,
            ),
            PartitionSubstitutionParameterState(
                partition_name="gene_gamma",
                exchangeabilities=(1.0, 1.0, 2.0, 1.0, 2.0, 3.0),
                base_frequencies=(0.2, 0.3, 0.2, 0.3),
                gamma_alpha=0.9,
            ),
        ),
    )

    assert linked_report.partition_count == 3
    assert linked_report.parameter_count == 3
    assert all(
        row.partition_names == ("gene_alpha", "gene_beta", "gene_gamma")
        for row in linked_report.rows
    )
    assert unlinked_report.parameter_count == 9
    assert len(unlinked_report.rows) == 9
    assert not math.isclose(
        linked_report.total_log_prior,
        unlinked_report.total_log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_partition_model_priors_respect_per_partition_model_targets() -> None:
    prior_bundle = build_partition_model_prior_bundle(
        partition_models=(
            build_partition_substitution_model_definition(
                partition_name="gene_alpha",
                model_name="GTR+G+I",
            ),
            build_partition_substitution_model_definition(
                partition_name="gene_beta",
                model_name="HKY85+G",
            ),
            build_partition_substitution_model_definition(
                partition_name="gene_gamma",
                model_name="JC69",
            ),
        ),
        linkage_plan=build_partition_parameter_linkage_plan(
            partition_names=("gene_alpha", "gene_beta", "gene_gamma"),
            linkage_policies={
                "kappa": "unlinked",
                "exchangeabilities": "unlinked",
                "base-frequencies": "unlinked",
                "gamma-alpha": "unlinked",
                "invariant-proportion": "unlinked",
            },
        ),
        substitution_prior_bundle=build_substitution_parameter_prior_bundle(
            kappa_prior=build_exponential_positive_substitution_parameter_prior(
                rate=0.4
            ),
            exchangeability_prior=build_dirichlet_simplex_substitution_parameter_prior(
                expected_component_names=("AC", "AG", "AT", "CG", "CT", "GT"),
                concentration_parameters=(2.0, 2.0, 2.0, 2.0, 2.0, 2.0),
            ),
            base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
                expected_component_names=("A", "C", "G", "T"),
                concentration_parameters=(3.0, 2.0, 2.0, 3.0),
            ),
            gamma_alpha_prior=build_exponential_positive_substitution_parameter_prior(
                rate=0.9
            ),
            invariant_proportion_prior=build_beta_probability_substitution_parameter_prior(
                alpha=2.0,
                beta=5.0,
            ),
        ),
    )

    report = evaluate_partition_model_log_prior(
        prior_bundle=prior_bundle,
        partition_parameter_states=(
            PartitionSubstitutionParameterState(
                partition_name="gene_alpha",
                exchangeabilities=(1.0, 2.0, 1.0, 3.0, 2.0, 1.0),
                base_frequencies=(0.3, 0.2, 0.1, 0.4),
                gamma_alpha=0.7,
                invariant_proportion=0.1,
            ),
            PartitionSubstitutionParameterState(
                partition_name="gene_beta",
                kappa=2.5,
                base_frequencies=(0.25, 0.15, 0.25, 0.35),
                gamma_alpha=1.2,
            ),
            PartitionSubstitutionParameterState(partition_name="gene_gamma"),
        ),
    )

    target_counts = {
        target_name: sum(row.target_name == target_name for row in report.rows)
        for target_name in {
            "kappa",
            "exchangeabilities",
            "base-frequencies",
            "gamma-alpha",
            "invariant-proportion",
        }
    }

    assert report.partition_count == 3
    assert report.parameter_count == 7
    assert target_counts == {
        "kappa": 1,
        "exchangeabilities": 1,
        "base-frequencies": 2,
        "gamma-alpha": 2,
        "invariant-proportion": 1,
    }
    assert {
        tuple(row.partition_names)
        for row in report.rows
        if row.target_name == "base-frequencies"
    } == {
        ("gene_alpha",),
        ("gene_beta",),
    }
    assert all("gene_gamma" not in row.partition_names for row in report.rows)


def test_partition_model_priors_reject_mismatched_linked_parameter_values() -> None:
    prior_bundle = build_partition_model_prior_bundle(
        partition_models=(
            build_partition_substitution_model_definition(
                partition_name="gene_alpha",
                model_name="GTR+G",
            ),
            build_partition_substitution_model_definition(
                partition_name="gene_beta",
                model_name="GTR+G",
            ),
        ),
        linkage_plan=build_partition_parameter_linkage_plan(
            partition_names=("gene_alpha", "gene_beta"),
            linkage_policies={"base-frequencies": "linked"},
        ),
        substitution_prior_bundle=build_substitution_parameter_prior_bundle(
            base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
                expected_component_names=("A", "C", "G", "T"),
                concentration_parameters=(3.0, 2.0, 2.0, 3.0),
            )
        ),
    )

    with pytest.raises(
        PhylogeneticsError,
        match="linked partition parameter group received mismatched realized values",
    ):
        evaluate_partition_model_log_prior(
            prior_bundle=prior_bundle,
            partition_parameter_states=(
                PartitionSubstitutionParameterState(
                    partition_name="gene_alpha",
                    base_frequencies=(0.3, 0.2, 0.1, 0.4),
                ),
                PartitionSubstitutionParameterState(
                    partition_name="gene_beta",
                    base_frequencies=(0.25, 0.15, 0.25, 0.35),
                ),
            ),
        )


def test_partition_model_prior_bundle_rejects_unused_prior_targets() -> None:
    with pytest.raises(
        PhylogeneticsError,
        match="configures substitution priors that no partition model can use",
    ):
        build_partition_model_prior_bundle(
            partition_models=(
                build_partition_substitution_model_definition(
                    partition_name="gene_gamma",
                    model_name="JC69",
                ),
            ),
            linkage_plan=build_partition_parameter_linkage_plan(
                partition_names=("gene_gamma",)
            ),
            substitution_prior_bundle=build_substitution_parameter_prior_bundle(
                kappa_prior=build_exponential_positive_substitution_parameter_prior(
                    rate=0.5
                )
            ),
        )
