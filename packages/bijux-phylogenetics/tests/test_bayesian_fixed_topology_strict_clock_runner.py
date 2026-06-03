from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.calibration_priors import (
    load_calibration_prior_definitions,
)
from bijux_phylogenetics.bayesian.clock_model_priors import (
    build_exponential_clock_model_scalar_prior,
)
from bijux_phylogenetics.bayesian.fixed_topology_strict_clock import (
    build_fixed_topology_strict_clock_model_definition,
    build_fixed_topology_strict_clock_proposal_schedule,
    run_fixed_topology_strict_clock_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    build_crown_conditioned_yule_tree_prior,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_fixed_topology_strict_clock_runner_emits_clock_rate_and_node_age_summaries() -> (
    None
):
    model_definition = build_fixed_topology_strict_clock_model_definition(
        time_tree_prior=build_crown_conditioned_yule_tree_prior(speciation_rate=0.4),
        global_clock_rate_prior=build_exponential_clock_model_scalar_prior(rate=1.5),
        calibration_priors=load_calibration_prior_definitions(
            fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
            fixture("metadata", "strict_clock_calibrations_root_recent_4_taxa.tsv"),
        ),
        initial_global_clock_rate=0.5,
    )
    proposal_schedule = build_fixed_topology_strict_clock_proposal_schedule(
        model_definition=model_definition,
        global_clock_rate_move_weight=1.0,
        global_clock_rate_log_scale_standard_deviation=0.12,
    )

    report = run_fixed_topology_strict_clock_metropolis_hastings(
        substitution_tree=_build_strict_clock_substitution_tree_fixture(clock_rate=0.5),
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        iteration_count=40,
        sample_every=1,
        seed=7,
    )

    assert len(report.posterior_rows) == 41
    assert report.clock_rate_summary.sample_count == len(report.posterior_rows)
    assert (
        report.clock_rate_summary.hpd_95_lower < report.clock_rate_summary.hpd_95_upper
    )
    assert len(report.node_age_summaries) == 3
    assert all(
        row.topology_id == report.posterior_rows[0].topology_id
        for row in report.posterior_rows
    )
    assert all(math.isfinite(row.total_log_prior) for row in report.posterior_rows)
    assert all(math.isfinite(row.log_likelihood) for row in report.posterior_rows)
    assert all(row.global_clock_rate > 0.0 for row in report.posterior_rows)
    assert {
        step_row.proposal_changed_fields for step_row in report.chain_report.step_rows
    } == {("scalar_parameters.global-clock-rate", "tree.branch_lengths")}
    root_summary = next(
        summary for summary in report.node_age_summaries if summary.node_kind == "root"
    )
    assert root_summary.sample_count == len(report.posterior_rows)
    assert root_summary.hpd_95_lower <= root_summary.posterior_mean
    assert root_summary.posterior_mean <= root_summary.hpd_95_upper


def test_fixed_topology_strict_clock_runner_clock_rate_posterior_changes_under_calibration() -> (
    None
):
    substitution_tree = _build_strict_clock_substitution_tree_fixture(clock_rate=0.5)
    recent_model_definition = build_fixed_topology_strict_clock_model_definition(
        time_tree_prior=build_crown_conditioned_yule_tree_prior(speciation_rate=0.4),
        global_clock_rate_prior=build_exponential_clock_model_scalar_prior(rate=1.5),
        calibration_priors=load_calibration_prior_definitions(
            fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
            fixture("metadata", "strict_clock_calibrations_root_recent_4_taxa.tsv"),
        ),
        initial_global_clock_rate=0.5,
    )
    deep_model_definition = build_fixed_topology_strict_clock_model_definition(
        time_tree_prior=build_crown_conditioned_yule_tree_prior(speciation_rate=0.4),
        global_clock_rate_prior=build_exponential_clock_model_scalar_prior(rate=1.5),
        calibration_priors=load_calibration_prior_definitions(
            fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
            fixture("metadata", "strict_clock_calibrations_root_deep_4_taxa.tsv"),
        ),
    )
    proposal_schedule = build_fixed_topology_strict_clock_proposal_schedule(
        model_definition=recent_model_definition,
        global_clock_rate_move_weight=1.0,
        global_clock_rate_log_scale_standard_deviation=0.12,
    )

    recent_report = run_fixed_topology_strict_clock_metropolis_hastings(
        substitution_tree=substitution_tree,
        model_definition=recent_model_definition,
        proposal_schedule=proposal_schedule,
        iteration_count=60,
        sample_every=1,
        seed=11,
    )
    deep_report = run_fixed_topology_strict_clock_metropolis_hastings(
        substitution_tree=substitution_tree,
        model_definition=deep_model_definition,
        proposal_schedule=proposal_schedule,
        iteration_count=60,
        sample_every=1,
        seed=11,
    )

    assert (
        deep_report.clock_rate_summary.posterior_mean
        < recent_report.clock_rate_summary.posterior_mean
    )
    assert (
        deep_report.posterior_rows[-1].root_age
        > recent_report.posterior_rows[-1].root_age
    )


def test_fixed_topology_strict_clock_runner_seeds_rate_from_root_calibration_center() -> (
    None
):
    model_definition = build_fixed_topology_strict_clock_model_definition(
        time_tree_prior=build_crown_conditioned_yule_tree_prior(speciation_rate=0.4),
        global_clock_rate_prior=build_exponential_clock_model_scalar_prior(rate=1.5),
        calibration_priors=load_calibration_prior_definitions(
            fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
            fixture("metadata", "strict_clock_calibrations_root_deep_4_taxa.tsv"),
        ),
    )
    proposal_schedule = build_fixed_topology_strict_clock_proposal_schedule(
        model_definition=model_definition,
        global_clock_rate_move_weight=1.0,
        global_clock_rate_log_scale_standard_deviation=0.12,
    )

    report = run_fixed_topology_strict_clock_metropolis_hastings(
        substitution_tree=_build_strict_clock_substitution_tree_fixture(clock_rate=0.5),
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        iteration_count=4,
        sample_every=1,
        seed=3,
    )

    assert math.isclose(
        report.posterior_rows[0].global_clock_rate,
        0.25,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.posterior_rows[0].root_age,
        6.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_fixed_topology_strict_clock_runner_rejects_unrooted_substitution_tree() -> (
    None
):
    model_definition = build_fixed_topology_strict_clock_model_definition(
        time_tree_prior=build_crown_conditioned_yule_tree_prior(speciation_rate=0.4),
        global_clock_rate_prior=build_exponential_clock_model_scalar_prior(rate=1.5),
        initial_global_clock_rate=0.5,
    )
    proposal_schedule = build_fixed_topology_strict_clock_proposal_schedule(
        model_definition=model_definition,
        global_clock_rate_move_weight=1.0,
        global_clock_rate_log_scale_standard_deviation=0.12,
    )
    unrooted_substitution_tree = _build_strict_clock_substitution_tree_fixture(
        clock_rate=0.5
    )
    unrooted_substitution_tree.rooted = False

    with pytest.raises(PhylogeneticsError, match="rooted substitution_tree"):
        run_fixed_topology_strict_clock_metropolis_hastings(
            substitution_tree=unrooted_substitution_tree,
            model_definition=model_definition,
            proposal_schedule=proposal_schedule,
            iteration_count=4,
            sample_every=1,
            seed=0,
        )


def _build_strict_clock_substitution_tree_fixture(*, clock_rate: float) -> PhyloTree:
    dated_tree = load_tree(fixture("trees", "strict_clock_time_tree_4_taxa.nwk"))
    dated_tree.rooted = True
    substitution_tree = dated_tree.copy()
    substitution_tree.rooted = True
    for _parent, child in substitution_tree.iter_edges():
        child.branch_length = float(child.branch_length or 0.0) * clock_rate
    substitution_tree.refresh()
    return substitution_tree
