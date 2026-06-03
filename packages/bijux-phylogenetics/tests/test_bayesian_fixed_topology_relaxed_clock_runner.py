from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.clock_model_priors import (
    build_exponential_clock_model_scalar_prior,
)
from bijux_phylogenetics.bayesian.fixed_topology_relaxed_clock import (
    build_fixed_topology_relaxed_clock_model_definition,
    build_fixed_topology_relaxed_clock_proposal_schedule,
    run_fixed_topology_relaxed_clock_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    build_crown_conditioned_yule_tree_prior,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_fixed_topology_relaxed_clock_runner_emits_branch_rate_and_node_age_summaries() -> (
    None
):
    model_definition = build_fixed_topology_relaxed_clock_model_definition(
        rate_policy="independent",
        time_tree_prior=build_crown_conditioned_yule_tree_prior(speciation_rate=0.4),
        mean_clock_rate_prior=build_exponential_clock_model_scalar_prior(rate=3.0),
        log_standard_deviation_prior=build_exponential_clock_model_scalar_prior(
            rate=4.0
        ),
    )
    proposal_schedule = build_fixed_topology_relaxed_clock_proposal_schedule(
        model_definition=model_definition,
        mean_clock_rate_move_weight=1.0,
        mean_clock_rate_log_scale_standard_deviation=0.12,
        log_standard_deviation_move_weight=1.0,
        log_standard_deviation_log_scale_standard_deviation=0.1,
        node_height_move_weight=2.0,
        node_height_slide_standard_deviation=0.04,
        tree_height_move_weight=1.0,
        tree_height_log_scale_standard_deviation=0.08,
    )

    report = run_fixed_topology_relaxed_clock_metropolis_hastings(
        substitution_tree=_load_rooted_tree_fixture(
            "relaxed_rate_summary_substitution_tree_4_taxa.nwk"
        ),
        initial_dated_tree=_load_rooted_tree_fixture(
            "relaxed_rate_summary_dated_tree_4_taxa.nwk"
        ),
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        iteration_count=30,
        sample_every=1,
        seed=9,
    )

    assert len(report.posterior_rows) == 31
    assert len(report.branch_rate_summaries) == 6
    assert len(report.node_age_summaries) == 3
    assert all(
        row.topology_id == report.posterior_rows[0].topology_id
        for row in report.posterior_rows
    )
    assert all(math.isfinite(row.log_likelihood) for row in report.posterior_rows)
    assert all(math.isfinite(row.posterior_log_score) for row in report.posterior_rows)
    assert all(row.mean_clock_rate > 0.0 for row in report.posterior_rows)
    assert all(row.log_standard_deviation > 0.0 for row in report.posterior_rows)
    changed_fields = {
        step_row.proposal_changed_fields for step_row in report.chain_report.step_rows
    }
    assert ("scalar_parameters.mean-clock-rate",) in changed_fields
    assert ("scalar_parameters.log-standard-deviation",) in changed_fields
    assert ("tree.branch_lengths",) in changed_fields
    assert any(
        proposal_changed_fields[0].startswith("tree.node_height:")
        for proposal_changed_fields in changed_fields
    )
    assert any(
        summary.hpd_95_lower < summary.hpd_95_upper
        for summary in report.branch_rate_summaries
    )
    root_summary = next(
        summary for summary in report.node_age_summaries if summary.node_kind == "root"
    )
    assert root_summary.sample_count == len(report.posterior_rows)
    assert root_summary.hpd_95_lower <= root_summary.posterior_mean
    assert root_summary.posterior_mean <= root_summary.hpd_95_upper
    assert any(
        summary.node_kind == "internal" and summary.hpd_95_lower < summary.hpd_95_upper
        for summary in report.node_age_summaries
    )


def test_fixed_topology_relaxed_clock_runner_rejects_unrooted_substitution_tree() -> (
    None
):
    model_definition = build_fixed_topology_relaxed_clock_model_definition(
        rate_policy="independent",
        time_tree_prior=build_crown_conditioned_yule_tree_prior(speciation_rate=0.4),
        mean_clock_rate_prior=build_exponential_clock_model_scalar_prior(rate=3.0),
        log_standard_deviation_prior=build_exponential_clock_model_scalar_prior(
            rate=4.0
        ),
    )
    proposal_schedule = build_fixed_topology_relaxed_clock_proposal_schedule(
        model_definition=model_definition,
        mean_clock_rate_move_weight=1.0,
        mean_clock_rate_log_scale_standard_deviation=0.12,
        log_standard_deviation_move_weight=1.0,
        log_standard_deviation_log_scale_standard_deviation=0.1,
        node_height_move_weight=2.0,
        node_height_slide_standard_deviation=0.04,
        tree_height_move_weight=1.0,
        tree_height_log_scale_standard_deviation=0.08,
    )
    unrooted_substitution_tree = _load_rooted_tree_fixture(
        "relaxed_rate_summary_substitution_tree_4_taxa.nwk"
    )
    unrooted_substitution_tree.rooted = False

    with pytest.raises(PhylogeneticsError, match="rooted substitution_tree"):
        run_fixed_topology_relaxed_clock_metropolis_hastings(
            substitution_tree=unrooted_substitution_tree,
            initial_dated_tree=_load_rooted_tree_fixture(
                "relaxed_rate_summary_dated_tree_4_taxa.nwk"
            ),
            model_definition=model_definition,
            proposal_schedule=proposal_schedule,
            iteration_count=4,
            sample_every=1,
            seed=0,
        )


def _load_rooted_tree_fixture(name: str):
    tree = load_tree(fixture("trees", name))
    tree.rooted = True
    return tree
