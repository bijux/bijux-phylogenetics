from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.calibration_priors import CalibrationPriorDefinition
from bijux_phylogenetics.bayesian.clock_model_priors import (
    build_exponential_clock_model_scalar_prior,
    build_fixed_clock_model_scalar_prior,
)
from bijux_phylogenetics.bayesian.fixed_topology_relaxed_clock import (
    FIXED_TOPOLOGY_RELAXED_CLOCK_MODELS,
    FixedTopologyRelaxedClockModelDefinition,
    FixedTopologyRelaxedClockProposalSchedule,
    build_fixed_topology_relaxed_clock_model_definition,
    build_fixed_topology_relaxed_clock_proposal_schedule,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    build_crown_conditioned_yule_tree_prior,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_fixed_topology_relaxed_clock_model_definition_normalizes_rate_policy() -> None:
    model_definition = build_fixed_topology_relaxed_clock_model_definition(
        rate_policy="Independent",
        time_tree_prior=build_crown_conditioned_yule_tree_prior(speciation_rate=0.5),
        mean_clock_rate_prior=build_exponential_clock_model_scalar_prior(rate=2.0),
        log_standard_deviation_prior=build_exponential_clock_model_scalar_prior(
            rate=3.0
        ),
    )

    assert FIXED_TOPOLOGY_RELAXED_CLOCK_MODELS == ("relaxed-lognormal",)
    assert isinstance(model_definition, FixedTopologyRelaxedClockModelDefinition)
    assert model_definition.rate_policy == "independent"
    assert model_definition.calibration_priors == []


def test_fixed_topology_relaxed_clock_model_definition_requires_sampled_clock_hyperpriors() -> (
    None
):
    with pytest.raises(PhylogeneticsError, match="non-fixed mean_clock_rate_prior"):
        build_fixed_topology_relaxed_clock_model_definition(
            rate_policy="independent",
            time_tree_prior=build_crown_conditioned_yule_tree_prior(
                speciation_rate=0.5
            ),
            mean_clock_rate_prior=build_fixed_clock_model_scalar_prior(fixed_value=0.2),
            log_standard_deviation_prior=build_exponential_clock_model_scalar_prior(
                rate=3.0
            ),
        )

    with pytest.raises(
        PhylogeneticsError,
        match="non-fixed log_standard_deviation_prior",
    ):
        build_fixed_topology_relaxed_clock_model_definition(
            rate_policy="independent",
            time_tree_prior=build_crown_conditioned_yule_tree_prior(
                speciation_rate=0.5
            ),
            mean_clock_rate_prior=build_exponential_clock_model_scalar_prior(rate=2.0),
            log_standard_deviation_prior=build_fixed_clock_model_scalar_prior(
                fixed_value=0.5
            ),
        )


def test_fixed_topology_relaxed_clock_model_definition_requires_unique_calibrations() -> (
    None
):
    duplicate_calibration = CalibrationPriorDefinition(
        calibration_id="cal-root",
        requested_distribution="fixed",
        family="fixed",
        target_kind="clade",
        target_label="root",
        descendant_taxa=["A", "B", "C", "D"],
        node_id="root",
        node_kind="root",
        minimum_age=2.0,
        maximum_age=2.0,
        translated=False,
        translation_note=None,
        fixed_age=2.0,
        fixed_tolerance=1e-12,
    )

    with pytest.raises(PhylogeneticsError, match="calibration ids to be unique"):
        build_fixed_topology_relaxed_clock_model_definition(
            rate_policy="independent",
            time_tree_prior=build_crown_conditioned_yule_tree_prior(
                speciation_rate=0.5
            ),
            mean_clock_rate_prior=build_exponential_clock_model_scalar_prior(rate=2.0),
            log_standard_deviation_prior=build_exponential_clock_model_scalar_prior(
                rate=3.0
            ),
            calibration_priors=[duplicate_calibration, duplicate_calibration],
        )


def test_fixed_topology_relaxed_clock_proposal_schedule_requires_positive_moves() -> (
    None
):
    model_definition = build_fixed_topology_relaxed_clock_model_definition(
        rate_policy="autocorrelated",
        time_tree_prior=build_crown_conditioned_yule_tree_prior(speciation_rate=0.5),
        mean_clock_rate_prior=build_exponential_clock_model_scalar_prior(rate=2.0),
        log_standard_deviation_prior=build_exponential_clock_model_scalar_prior(
            rate=3.0
        ),
    )

    proposal_schedule = build_fixed_topology_relaxed_clock_proposal_schedule(
        model_definition=model_definition,
        mean_clock_rate_move_weight=1.0,
        mean_clock_rate_log_scale_standard_deviation=0.2,
        log_standard_deviation_move_weight=1.0,
        log_standard_deviation_log_scale_standard_deviation=0.25,
        node_height_move_weight=2.0,
        node_height_slide_standard_deviation=0.1,
        tree_height_move_weight=1.0,
        tree_height_log_scale_standard_deviation=0.15,
    )

    assert isinstance(proposal_schedule, FixedTopologyRelaxedClockProposalSchedule)
    assert proposal_schedule.node_height_move_weight == 2.0
    assert proposal_schedule.tree_height_log_scale_standard_deviation == 0.15

    with pytest.raises(PhylogeneticsError, match="mean_clock_rate_move_weight"):
        build_fixed_topology_relaxed_clock_proposal_schedule(
            model_definition=model_definition,
            mean_clock_rate_move_weight=0.0,
            mean_clock_rate_log_scale_standard_deviation=0.2,
            log_standard_deviation_move_weight=1.0,
            log_standard_deviation_log_scale_standard_deviation=0.25,
            node_height_move_weight=2.0,
            node_height_slide_standard_deviation=0.1,
            tree_height_move_weight=1.0,
            tree_height_log_scale_standard_deviation=0.15,
        )


def test_fixed_topology_relaxed_clock_model_definition_rejects_unknown_rate_policy() -> (
    None
):
    with pytest.raises(PhylogeneticsError, match="supported relaxed-clock rate policy"):
        build_fixed_topology_relaxed_clock_model_definition(
            rate_policy="unsupported",
            time_tree_prior=build_crown_conditioned_yule_tree_prior(
                speciation_rate=0.5
            ),
            mean_clock_rate_prior=build_exponential_clock_model_scalar_prior(rate=2.0),
            log_standard_deviation_prior=build_exponential_clock_model_scalar_prior(
                rate=3.0
            ),
        )
