from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.ancestral.discrete import DiscreteTransitionRateRow
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    build_exponential_discrete_trait_rate_prior,
    build_gamma_discrete_trait_rate_prior,
    build_lognormal_discrete_trait_rate_prior,
    evaluate_discrete_trait_rate_log_prior,
)
from bijux_phylogenetics.comparative.discrete_mk import fit_discrete_mk_model
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_phytools_comparative_fixture,
)
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


def _gamma_log_density(value: float, *, shape: float, scale: float) -> float:
    return (
        ((shape - 1.0) * math.log(value))
        - (value / scale)
        - math.lgamma(shape)
        - (shape * math.log(scale))
    )


def _lognormal_log_density(
    value: float,
    *,
    log_mean: float,
    log_standard_deviation: float,
) -> float:
    return (
        -math.log(value)
        - math.log(log_standard_deviation)
        - (0.5 * math.log(2.0 * math.pi))
        - (((math.log(value) - log_mean) ** 2) / (2.0 * (log_standard_deviation**2)))
    )


def test_equal_rates_discrete_trait_prior_collapses_to_one_shared_parameter() -> None:
    fit_report = fit_discrete_mk_model(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="region_state",
        taxon_column="taxon",
        model="equal-rates",
    )
    prior_model = build_exponential_discrete_trait_rate_prior(rate=0.75)

    report = evaluate_discrete_trait_rate_log_prior(
        model=fit_report.model,
        transition_rate_rows=fit_report.transition_rate_rows,
        prior_model=prior_model,
    )

    shared_rate = next(
        row.rate for row in fit_report.transition_rate_rows if row.transition_allowed
    )
    expected_log_prior = math.log(0.75) - (0.75 * shared_rate)

    assert report.model == "equal-rates"
    assert report.parameter_count == 1
    assert report.rows[0].parameter_name == "shared-rate"
    assert len(report.rows[0].transition_pairs) == 6
    assert math.isclose(
        report.total_log_prior,
        expected_log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


@pytest.mark.slow
def test_symmetric_discrete_trait_prior_scores_each_bidirectional_pair_once() -> None:
    fit_report = fit_discrete_mk_model(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="region_state",
        taxon_column="taxon",
        model="symmetric",
    )
    prior_model = build_gamma_discrete_trait_rate_prior(shape=2.5, scale=0.75)

    report = evaluate_discrete_trait_rate_log_prior(
        model=fit_report.model,
        transition_rate_rows=fit_report.transition_rate_rows,
        prior_model=prior_model,
    )

    expected_total_log_prior = math.fsum(
        _gamma_log_density(row.rate_value, shape=2.5, scale=0.75) for row in report.rows
    )

    assert report.parameter_count == 3
    assert {row.parameter_name for row in report.rows} == {
        "north<->south",
        "north<->west",
        "south<->west",
    }
    assert all(len(row.transition_pairs) == 2 for row in report.rows)
    assert math.isclose(
        report.total_log_prior,
        expected_total_log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_ard_discrete_trait_prior_scores_each_directional_rate() -> None:
    fixture_entry = get_shared_phytools_comparative_fixture(
        "phytools_discrete_ard_binary_twenty_four_taxa"
    )
    fit_report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="all-rates-different",
    )
    prior_model = build_lognormal_discrete_trait_rate_prior(
        log_mean=math.log(0.5),
        log_standard_deviation=0.8,
    )

    report = evaluate_discrete_trait_rate_log_prior(
        model=fit_report.model,
        transition_rate_rows=fit_report.transition_rate_rows,
        prior_model=prior_model,
    )

    expected_total_log_prior = math.fsum(
        _lognormal_log_density(
            row.rate_value,
            log_mean=math.log(0.5),
            log_standard_deviation=0.8,
        )
        for row in report.rows
    )

    assert report.parameter_count == 2
    assert {row.parameter_name for row in report.rows} == {"0->1", "1->0"}
    assert all(len(row.transition_pairs) == 1 for row in report.rows)
    assert math.isclose(
        report.total_log_prior,
        expected_total_log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_discrete_trait_rate_prior_changes_combined_mk_posterior_score() -> None:
    fit_report = fit_discrete_mk_model(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="region_state",
        taxon_column="taxon",
        model="equal-rates",
    )
    first_prior = build_exponential_discrete_trait_rate_prior(rate=0.25)
    second_prior = build_exponential_discrete_trait_rate_prior(rate=1.25)

    first_report = evaluate_discrete_trait_rate_log_prior(
        model=fit_report.model,
        transition_rate_rows=fit_report.transition_rate_rows,
        prior_model=first_prior,
    )
    second_report = evaluate_discrete_trait_rate_log_prior(
        model=fit_report.model,
        transition_rate_rows=fit_report.transition_rate_rows,
        prior_model=second_prior,
    )

    first_posterior_score = fit_report.log_likelihood + first_report.total_log_prior
    second_posterior_score = fit_report.log_likelihood + second_report.total_log_prior

    assert not math.isclose(
        first_report.total_log_prior,
        second_report.total_log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        first_posterior_score,
        second_posterior_score,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_symmetric_discrete_trait_prior_rejects_mismatched_reverse_rates() -> None:
    prior_model = build_exponential_discrete_trait_rate_prior(rate=0.5)

    with pytest.raises(
        PhylogeneticsError,
        match="symmetric prior evaluation requires matched forward and reverse transition rates",
    ):
        evaluate_discrete_trait_rate_log_prior(
            model="symmetric",
            transition_rate_rows=[
                DiscreteTransitionRateRow(
                    source_state="A",
                    target_state="B",
                    transition_allowed=True,
                    step_distance=1,
                    rate=0.25,
                ),
                DiscreteTransitionRateRow(
                    source_state="B",
                    target_state="A",
                    transition_allowed=True,
                    step_distance=1,
                    rate=0.5,
                ),
            ],
            prior_model=prior_model,
        )
