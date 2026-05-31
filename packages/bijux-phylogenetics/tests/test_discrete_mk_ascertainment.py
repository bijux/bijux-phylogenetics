from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import load_discrete_dataset
from bijux_phylogenetics.comparative.discrete_mk import (
    compare_discrete_mk_model_ranking,
    fit_discrete_mk_model,
    fit_discrete_mk_model_from_dataset,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_discrete_mk_lewis_ascertainment_changes_variable_only_likelihood() -> None:
    uncorrected = fit_discrete_mk_model(
        fixture("trees", "example_tree.nwk"),
        fixture("metadata", "example_traits_discrete_mk_variable_only_four_taxa.tsv"),
        trait="state",
        taxon_column="taxon",
        model="equal-rates",
    )
    corrected = fit_discrete_mk_model(
        fixture("trees", "example_tree.nwk"),
        fixture("metadata", "example_traits_discrete_mk_variable_only_four_taxa.tsv"),
        trait="state",
        taxon_column="taxon",
        model="equal-rates",
        ascertainment_policy="lewis-variable-only",
    )

    assert uncorrected.ascertainment_policy == "none"
    assert uncorrected.ascertainment_conditioning_log_probability is None
    assert uncorrected.invariant_pattern_log_probability is None
    assert corrected.ascertainment_policy == "lewis-variable-only"
    assert corrected.ascertainment_conditioning_log_probability is not None
    assert corrected.invariant_pattern_log_probability is not None
    assert corrected.likelihood_constant_policy.endswith(
        "declared-lewis-mk-ascertainment-policy"
    )
    assert corrected.log_likelihood > uncorrected.log_likelihood
    assert not math.isclose(
        corrected.log_likelihood,
        uncorrected.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_discrete_mk_lewis_ascertainment_dataset_surface_matches_path_surface() -> None:
    dataset = load_discrete_dataset(
        fixture("trees", "example_tree.nwk"),
        fixture("metadata", "example_traits_discrete_mk_variable_only_four_taxa.tsv"),
        trait="state",
        taxon_column="taxon",
    )

    path_report = fit_discrete_mk_model(
        dataset.tree_path,
        dataset.traits_path,
        trait=dataset.trait,
        taxon_column=dataset.taxon_column,
        model="equal-rates",
        ascertainment_policy="lewis-variable-only",
    )
    dataset_report = fit_discrete_mk_model_from_dataset(
        dataset,
        model="equal-rates",
        ascertainment_policy="lewis-variable-only",
    )

    assert dataset_report.log_likelihood == path_report.log_likelihood
    assert (
        dataset_report.ascertainment_conditioning_log_probability
        == path_report.ascertainment_conditioning_log_probability
    )
    assert (
        dataset_report.invariant_pattern_log_probability
        == path_report.invariant_pattern_log_probability
    )


def test_discrete_mk_model_ranking_carries_ascertainment_policy() -> None:
    report = compare_discrete_mk_model_ranking(
        fixture("trees", "example_tree.nwk"),
        fixture("metadata", "example_traits_discrete_mk_variable_only_four_taxa.tsv"),
        trait="state",
        taxon_column="taxon",
        ascertainment_policy="lewis-variable-only",
    )

    assert report.ascertainment_policy == "lewis-variable-only"
    assert all(
        row.likelihood_constant_policy == report.likelihood_constant_policy
        for row in report.rows
        if row.comparable
    )
