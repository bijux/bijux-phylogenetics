from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.calibration_priors import (
    CalibrationPriorDefinition,
    load_calibration_prior_definitions,
)
from bijux_phylogenetics.bayesian.clock_model_priors import (
    build_exponential_clock_model_scalar_prior,
    build_fixed_clock_model_scalar_prior,
)
from bijux_phylogenetics.bayesian.fixed_topology_strict_clock import (
    build_fixed_topology_strict_clock_model_definition,
    build_fixed_topology_strict_clock_proposal_schedule,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    build_crown_conditioned_yule_tree_prior,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_fixed_topology_strict_clock_model_definition_preserves_priors_and_tuning() -> (
    None
):
    calibration_priors = load_calibration_prior_definitions(
        fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
        fixture("metadata", "strict_clock_calibrations_root_recent_4_taxa.tsv"),
    )

    model_definition = build_fixed_topology_strict_clock_model_definition(
        time_tree_prior=build_crown_conditioned_yule_tree_prior(speciation_rate=0.5),
        global_clock_rate_prior=build_exponential_clock_model_scalar_prior(rate=2.0),
        calibration_priors=calibration_priors,
        initial_global_clock_rate=0.6,
        branch_length_tolerance=1e-10,
    )
    proposal_schedule = build_fixed_topology_strict_clock_proposal_schedule(
        model_definition=model_definition,
        global_clock_rate_move_weight=1.0,
        global_clock_rate_log_scale_standard_deviation=0.15,
    )

    assert model_definition.time_tree_prior.family == "crown-conditioned-yule"
    assert model_definition.global_clock_rate_prior.family == "exponential"
    assert model_definition.calibration_priors == calibration_priors
    assert model_definition.initial_global_clock_rate == 0.6
    assert model_definition.branch_length_tolerance == 1e-10
    assert proposal_schedule.global_clock_rate_move_weight == 1.0
    assert proposal_schedule.global_clock_rate_log_scale_standard_deviation == 0.15


def test_fixed_topology_strict_clock_model_definition_requires_sampled_clock_rate() -> (
    None
):
    with pytest.raises(PhylogeneticsError, match="non-fixed global_clock_rate_prior"):
        build_fixed_topology_strict_clock_model_definition(
            time_tree_prior=build_crown_conditioned_yule_tree_prior(
                speciation_rate=0.5
            ),
            global_clock_rate_prior=build_fixed_clock_model_scalar_prior(
                fixed_value=0.5
            ),
        )


def test_fixed_topology_strict_clock_model_definition_requires_unique_calibrations() -> (
    None
):
    duplicate_calibration = CalibrationPriorDefinition(
        calibration_id="cal-root",
        requested_distribution="uniform",
        family="uniform",
        target_kind="clade",
        target_label="A|B|C|D",
        descendant_taxa=["A", "B", "C", "D"],
        node_id="root",
        node_kind="root",
        minimum_age=2.8,
        maximum_age=3.2,
        translated=False,
        translation_note=None,
    )

    with pytest.raises(PhylogeneticsError, match="calibration ids to be unique"):
        build_fixed_topology_strict_clock_model_definition(
            time_tree_prior=build_crown_conditioned_yule_tree_prior(
                speciation_rate=0.5
            ),
            global_clock_rate_prior=build_exponential_clock_model_scalar_prior(
                rate=2.0
            ),
            calibration_priors=[duplicate_calibration, duplicate_calibration],
        )


def test_fixed_topology_strict_clock_proposal_schedule_requires_model_definition() -> (
    None
):
    with pytest.raises(
        PhylogeneticsError,
        match="requires one FixedTopologyStrictClockModelDefinition",
    ):
        build_fixed_topology_strict_clock_proposal_schedule(
            model_definition="not-a-model-definition",
            global_clock_rate_move_weight=1.0,
            global_clock_rate_log_scale_standard_deviation=0.15,
        )
