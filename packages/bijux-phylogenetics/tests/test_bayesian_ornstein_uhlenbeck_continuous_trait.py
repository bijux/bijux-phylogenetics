from __future__ import annotations

import csv
import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.continuous_trait_location_priors import (
    build_fixed_continuous_trait_location_prior,
    build_normal_continuous_trait_location_prior,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    build_exponential_continuous_trait_scalar_prior,
    build_fixed_continuous_trait_scalar_prior,
)
from bijux_phylogenetics.bayesian.ornstein_uhlenbeck_continuous_trait import (
    build_ornstein_uhlenbeck_continuous_trait_model_definition,
    build_ornstein_uhlenbeck_continuous_trait_proposal_schedule,
    run_ornstein_uhlenbeck_continuous_trait_metropolis_hastings,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_ou_continuous_trait_model_requires_sampled_alpha_prior() -> None:
    with pytest.raises(PhylogeneticsError, match="non-fixed alpha prior"):
        build_ornstein_uhlenbeck_continuous_trait_model_definition(
            alpha_prior=build_fixed_continuous_trait_scalar_prior(fixed_value=0.5),
            optimum_prior=build_normal_continuous_trait_location_prior(
                mean=0.0,
                standard_deviation=2.0,
            ),
            sigma_squared_prior=build_exponential_continuous_trait_scalar_prior(
                rate=1.0
            ),
        )


def test_ou_continuous_trait_model_requires_sampled_optimum_prior() -> None:
    with pytest.raises(PhylogeneticsError, match="non-fixed optimum prior"):
        build_ornstein_uhlenbeck_continuous_trait_model_definition(
            alpha_prior=build_exponential_continuous_trait_scalar_prior(rate=1.0),
            optimum_prior=build_fixed_continuous_trait_location_prior(fixed_value=0.0),
            sigma_squared_prior=build_exponential_continuous_trait_scalar_prior(
                rate=1.0
            ),
        )


def test_ou_continuous_trait_model_requires_sampled_sigma_squared_prior() -> None:
    with pytest.raises(PhylogeneticsError, match="non-fixed sigma-squared prior"):
        build_ornstein_uhlenbeck_continuous_trait_model_definition(
            alpha_prior=build_exponential_continuous_trait_scalar_prior(rate=1.0),
            optimum_prior=build_normal_continuous_trait_location_prior(
                mean=0.0,
                standard_deviation=2.0,
            ),
            sigma_squared_prior=build_fixed_continuous_trait_scalar_prior(
                fixed_value=0.5
            ),
        )


def test_ou_continuous_trait_runner_emits_sampled_parameter_summaries() -> None:
    model_definition = _build_model_definition()
    report = run_ornstein_uhlenbeck_continuous_trait_metropolis_hastings(
        tree=_load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk"),
        tip_values={"A": 0.2, "B": 0.5, "C": 1.0, "D": 1.4},
        model_definition=model_definition,
        proposal_schedule=build_ornstein_uhlenbeck_continuous_trait_proposal_schedule(
            model_definition=model_definition,
            alpha_move_weight=1.0,
            alpha_log_scale_standard_deviation=0.4,
            optimum_move_weight=1.0,
            optimum_standard_deviation=0.25,
            sigma_squared_move_weight=1.0,
            sigma_squared_log_scale_standard_deviation=0.4,
        ),
        iteration_count=24,
        sample_every=1,
        seed=7,
    )

    assert report.taxa == ["A", "B", "C", "D"]
    assert len(report.posterior_rows) == 25
    assert all(math.isfinite(row.log_likelihood) for row in report.posterior_rows)
    assert all(math.isfinite(row.posterior_log_score) for row in report.posterior_rows)
    assert all(row.alpha > 0.0 for row in report.posterior_rows)
    assert all(row.sigma_squared > 0.0 for row in report.posterior_rows)
    assert {
        step_row.proposal_changed_fields for step_row in report.chain_report.step_rows
    } == {
        ("scalar_parameters.alpha",),
        ("scalar_parameters.optimum",),
        ("scalar_parameters.sigma-squared",),
    }

    summary_by_parameter = {
        summary.parameter_name: summary for summary in report.parameter_summaries
    }
    assert set(summary_by_parameter) == {"alpha", "optimum", "sigma-squared"}
    for summary in summary_by_parameter.values():
        assert summary.sample_count == len(report.posterior_rows)
        assert math.isfinite(summary.posterior_mean)
        assert summary.hpd_95_lower <= summary.hpd_95_upper
    warning_kinds = [warning.kind for warning in report.identifiability_warnings]
    assert "small_sample_size" in warning_kinds


def test_ou_continuous_trait_weak_fixture_reports_broad_and_boundary_warnings() -> None:
    model_definition = build_ornstein_uhlenbeck_continuous_trait_model_definition(
        alpha_prior=build_exponential_continuous_trait_scalar_prior(rate=0.2),
        optimum_prior=build_normal_continuous_trait_location_prior(
            mean=0.0,
            standard_deviation=5.0,
        ),
        sigma_squared_prior=build_exponential_continuous_trait_scalar_prior(rate=0.2),
    )
    report = run_ornstein_uhlenbeck_continuous_trait_metropolis_hastings(
        tree=_load_rooted_tree_fixture("example_tree_six_taxa.nwk"),
        tip_values=_load_tip_values_from_tsv(
            name="example_traits_comparative_multiple.tsv",
            trait_name="response_growth",
        ),
        model_definition=model_definition,
        proposal_schedule=build_ornstein_uhlenbeck_continuous_trait_proposal_schedule(
            model_definition=model_definition,
            alpha_move_weight=2.0,
            alpha_log_scale_standard_deviation=0.9,
            optimum_move_weight=1.0,
            optimum_standard_deviation=0.5,
            sigma_squared_move_weight=1.0,
            sigma_squared_log_scale_standard_deviation=0.7,
        ),
        iteration_count=80,
        sample_every=1,
        seed=9,
    )

    warning_kinds = [warning.kind for warning in report.identifiability_warnings]

    assert "boundary_alpha" in warning_kinds
    assert "broad_alpha_posterior" in warning_kinds
    alpha_summary = next(
        summary
        for summary in report.parameter_summaries
        if summary.parameter_name == "alpha"
    )
    assert alpha_summary.hpd_95_lower < 0.1
    assert (alpha_summary.hpd_95_upper - alpha_summary.hpd_95_lower) > 1.0


def _build_model_definition():
    return build_ornstein_uhlenbeck_continuous_trait_model_definition(
        alpha_prior=build_exponential_continuous_trait_scalar_prior(rate=1.0),
        optimum_prior=build_normal_continuous_trait_location_prior(
            mean=0.0,
            standard_deviation=2.0,
        ),
        sigma_squared_prior=build_exponential_continuous_trait_scalar_prior(rate=1.0),
    )


def _load_tip_values_from_tsv(*, name: str, trait_name: str) -> dict[str, float]:
    with fixture("metadata", name).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {
            row["taxon"]: float(row[trait_name])
            for row in reader
            if row["taxon"] and row[trait_name]
        }


def _load_rooted_tree_fixture(name: str):
    tree = load_tree(fixture("trees", name))
    tree.rooted = True
    return tree
