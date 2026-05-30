from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.bayesian.calibration_priors import CalibrationPriorDefinition
from bijux_phylogenetics.bayesian.clock_model_priors import ClockModelScalarPriorModel
from bijux_phylogenetics.bayesian.metropolis_hastings import MetropolisHastingsRunReport
from bijux_phylogenetics.bayesian.time_tree_priors import (
    BirthDeathTreePriorModel,
    ConstantPopulationCoalescentPriorModel,
    SkylineCoalescentPriorModel,
    YuleTreePriorModel,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXED_TOPOLOGY_RELAXED_CLOCK_MODELS = ("relaxed-lognormal",)
_RELAXED_CLOCK_RATE_POLICIES = ("independent", "autocorrelated")


@dataclass(frozen=True, slots=True)
class FixedTopologyRelaxedClockModelDefinition:
    """One validated fixed-topology relaxed-clock posterior model definition."""

    rate_policy: str
    time_tree_prior: (
        YuleTreePriorModel
        | BirthDeathTreePriorModel
        | ConstantPopulationCoalescentPriorModel
        | SkylineCoalescentPriorModel
    )
    mean_clock_rate_prior: ClockModelScalarPriorModel
    log_standard_deviation_prior: ClockModelScalarPriorModel
    calibration_priors: list[CalibrationPriorDefinition]
    initial_mean_clock_rate: float | None = None
    initial_log_standard_deviation: float | None = None


@dataclass(frozen=True, slots=True)
class FixedTopologyRelaxedClockProposalSchedule:
    """One validated proposal schedule for relaxed-clock dating sampling."""

    mean_clock_rate_move_weight: float
    mean_clock_rate_log_scale_standard_deviation: float
    log_standard_deviation_move_weight: float
    log_standard_deviation_log_scale_standard_deviation: float
    node_height_move_weight: float
    node_height_slide_standard_deviation: float
    tree_height_move_weight: float
    tree_height_log_scale_standard_deviation: float


@dataclass(frozen=True, slots=True)
class FixedTopologyRelaxedClockPosteriorRow:
    """One sampled posterior row from a fixed-topology relaxed-clock chain."""

    sample_index: int
    iteration_index: int
    topology_id: str
    model_name: str
    rate_policy: str
    total_log_prior: float
    log_likelihood: float
    posterior_log_score: float
    prior_component_log_priors: dict[str, float]
    mean_clock_rate: float
    log_standard_deviation: float
    root_age: float
    tree_newick: str


@dataclass(frozen=True, slots=True)
class FixedTopologyRelaxedClockBranchRateSummary:
    """One posterior summary for one implied relaxed-clock branch rate."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    sample_count: int
    posterior_mean: float
    hpd_95_lower: float
    hpd_95_upper: float
    minimum_rate: float
    maximum_rate: float


@dataclass(frozen=True, slots=True)
class FixedTopologyRelaxedClockNodeAgeSummary:
    """One posterior summary for one sampled internal-node age."""

    node_id: str
    node_kind: str
    child_name: str | None
    descendant_taxa: list[str]
    sample_count: int
    posterior_mean: float
    hpd_95_lower: float
    hpd_95_upper: float
    minimum_age: float
    maximum_age: float


@dataclass(frozen=True, slots=True)
class FixedTopologyRelaxedClockRunReport:
    """One completed fixed-topology relaxed-clock posterior run."""

    model_definition: FixedTopologyRelaxedClockModelDefinition
    proposal_schedule: FixedTopologyRelaxedClockProposalSchedule
    chain_report: MetropolisHastingsRunReport
    posterior_rows: list[FixedTopologyRelaxedClockPosteriorRow]
    branch_rate_summaries: list[FixedTopologyRelaxedClockBranchRateSummary]
    node_age_summaries: list[FixedTopologyRelaxedClockNodeAgeSummary]


def build_fixed_topology_relaxed_clock_model_definition(
    *,
    rate_policy: str,
    time_tree_prior: (
        YuleTreePriorModel
        | BirthDeathTreePriorModel
        | ConstantPopulationCoalescentPriorModel
        | SkylineCoalescentPriorModel
    ),
    mean_clock_rate_prior: ClockModelScalarPriorModel,
    log_standard_deviation_prior: ClockModelScalarPriorModel,
    calibration_priors: list[CalibrationPriorDefinition] | None = None,
    initial_mean_clock_rate: float | None = None,
    initial_log_standard_deviation: float | None = None,
) -> FixedTopologyRelaxedClockModelDefinition:
    """Build one validated fixed-topology relaxed-clock posterior model."""
    validated_rate_policy = _validate_rate_policy(rate_policy)
    validated_time_tree_prior = _validate_time_tree_prior_model(time_tree_prior)
    validated_mean_clock_rate_prior = _validate_clock_model_scalar_prior(
        prior_model=mean_clock_rate_prior,
        field_name="mean_clock_rate_prior",
        owner_name="fixed-topology relaxed-clock posterior model",
    )
    if validated_mean_clock_rate_prior.family == "fixed":
        raise PhylogeneticsError(
            "fixed-topology relaxed-clock posterior model requires one non-fixed mean_clock_rate_prior so mean clock rate remains sampled",
            code="fixed_topology_relaxed_clock_mean_clock_rate_prior_fixed",
        )
    validated_log_standard_deviation_prior = _validate_clock_model_scalar_prior(
        prior_model=log_standard_deviation_prior,
        field_name="log_standard_deviation_prior",
        owner_name="fixed-topology relaxed-clock posterior model",
    )
    if validated_log_standard_deviation_prior.family == "fixed":
        raise PhylogeneticsError(
            "fixed-topology relaxed-clock posterior model requires one non-fixed log_standard_deviation_prior so relaxed-clock variance remains sampled",
            code="fixed_topology_relaxed_clock_log_standard_deviation_prior_fixed",
        )
    return FixedTopologyRelaxedClockModelDefinition(
        rate_policy=validated_rate_policy,
        time_tree_prior=validated_time_tree_prior,
        mean_clock_rate_prior=validated_mean_clock_rate_prior,
        log_standard_deviation_prior=validated_log_standard_deviation_prior,
        calibration_priors=_validate_calibration_priors(calibration_priors or []),
        initial_mean_clock_rate=_validate_optional_positive_finite_float(
            value=initial_mean_clock_rate,
            field_name="initial_mean_clock_rate",
            owner_name="fixed-topology relaxed-clock posterior model",
        ),
        initial_log_standard_deviation=_validate_optional_positive_finite_float(
            value=initial_log_standard_deviation,
            field_name="initial_log_standard_deviation",
            owner_name="fixed-topology relaxed-clock posterior model",
        ),
    )


def build_fixed_topology_relaxed_clock_proposal_schedule(
    *,
    model_definition: FixedTopologyRelaxedClockModelDefinition,
    mean_clock_rate_move_weight: float,
    mean_clock_rate_log_scale_standard_deviation: float,
    log_standard_deviation_move_weight: float,
    log_standard_deviation_log_scale_standard_deviation: float,
    node_height_move_weight: float,
    node_height_slide_standard_deviation: float,
    tree_height_move_weight: float,
    tree_height_log_scale_standard_deviation: float,
) -> FixedTopologyRelaxedClockProposalSchedule:
    """Build one validated proposal schedule for relaxed-clock dating."""
    if not isinstance(model_definition, FixedTopologyRelaxedClockModelDefinition):
        raise PhylogeneticsError(
            "fixed-topology relaxed-clock proposal schedule requires one FixedTopologyRelaxedClockModelDefinition",
            code="fixed_topology_relaxed_clock_proposal_schedule_model_definition_type_invalid",
        )
    return FixedTopologyRelaxedClockProposalSchedule(
        mean_clock_rate_move_weight=_validate_positive_finite_float(
            value=mean_clock_rate_move_weight,
            field_name="mean_clock_rate_move_weight",
            owner_name="fixed-topology relaxed-clock proposal schedule",
        ),
        mean_clock_rate_log_scale_standard_deviation=_validate_positive_finite_float(
            value=mean_clock_rate_log_scale_standard_deviation,
            field_name="mean_clock_rate_log_scale_standard_deviation",
            owner_name="fixed-topology relaxed-clock proposal schedule",
        ),
        log_standard_deviation_move_weight=_validate_positive_finite_float(
            value=log_standard_deviation_move_weight,
            field_name="log_standard_deviation_move_weight",
            owner_name="fixed-topology relaxed-clock proposal schedule",
        ),
        log_standard_deviation_log_scale_standard_deviation=_validate_positive_finite_float(
            value=log_standard_deviation_log_scale_standard_deviation,
            field_name="log_standard_deviation_log_scale_standard_deviation",
            owner_name="fixed-topology relaxed-clock proposal schedule",
        ),
        node_height_move_weight=_validate_positive_finite_float(
            value=node_height_move_weight,
            field_name="node_height_move_weight",
            owner_name="fixed-topology relaxed-clock proposal schedule",
        ),
        node_height_slide_standard_deviation=_validate_positive_finite_float(
            value=node_height_slide_standard_deviation,
            field_name="node_height_slide_standard_deviation",
            owner_name="fixed-topology relaxed-clock proposal schedule",
        ),
        tree_height_move_weight=_validate_positive_finite_float(
            value=tree_height_move_weight,
            field_name="tree_height_move_weight",
            owner_name="fixed-topology relaxed-clock proposal schedule",
        ),
        tree_height_log_scale_standard_deviation=_validate_positive_finite_float(
            value=tree_height_log_scale_standard_deviation,
            field_name="tree_height_log_scale_standard_deviation",
            owner_name="fixed-topology relaxed-clock proposal schedule",
        ),
    )


def _validate_rate_policy(rate_policy: str) -> str:
    normalized_rate_policy = rate_policy.strip().casefold()
    if normalized_rate_policy not in _RELAXED_CLOCK_RATE_POLICIES:
        raise PhylogeneticsError(
            "fixed-topology relaxed-clock posterior model requires one supported relaxed-clock rate policy",
            code="fixed_topology_relaxed_clock_rate_policy_invalid",
            details={
                "rate_policy": rate_policy,
                "allowed_rate_policies": list(_RELAXED_CLOCK_RATE_POLICIES),
            },
        )
    return normalized_rate_policy


def _validate_time_tree_prior_model(
    time_tree_prior: object,
) -> (
    YuleTreePriorModel
    | BirthDeathTreePriorModel
    | ConstantPopulationCoalescentPriorModel
    | SkylineCoalescentPriorModel
):
    if isinstance(
        time_tree_prior,
        (
            YuleTreePriorModel,
            BirthDeathTreePriorModel,
            ConstantPopulationCoalescentPriorModel,
            SkylineCoalescentPriorModel,
        ),
    ):
        return time_tree_prior
    raise PhylogeneticsError(
        "fixed-topology relaxed-clock posterior model requires one supported rooted time-tree prior model",
        code="fixed_topology_relaxed_clock_time_tree_prior_type_invalid",
    )


def _validate_clock_model_scalar_prior(
    *,
    prior_model: object,
    field_name: str,
    owner_name: str,
) -> ClockModelScalarPriorModel:
    if not isinstance(prior_model, ClockModelScalarPriorModel):
        raise PhylogeneticsError(
            f"{owner_name} requires one ClockModelScalarPriorModel for {field_name}",
            code="fixed_topology_relaxed_clock_clock_prior_type_invalid",
            details={"field_name": field_name},
        )
    return prior_model


def _validate_calibration_priors(
    calibration_priors: list[CalibrationPriorDefinition],
) -> list[CalibrationPriorDefinition]:
    seen_calibration_ids: set[str] = set()
    validated_calibration_priors: list[CalibrationPriorDefinition] = []
    for calibration_prior in calibration_priors:
        if not isinstance(calibration_prior, CalibrationPriorDefinition):
            raise PhylogeneticsError(
                "fixed-topology relaxed-clock posterior model requires every calibration prior to be one CalibrationPriorDefinition",
                code="fixed_topology_relaxed_clock_calibration_prior_type_invalid",
            )
        if calibration_prior.calibration_id in seen_calibration_ids:
            raise PhylogeneticsError(
                "fixed-topology relaxed-clock posterior model requires calibration ids to be unique",
                code="fixed_topology_relaxed_clock_calibration_prior_duplicate",
                details={"calibration_id": calibration_prior.calibration_id},
            )
        seen_calibration_ids.add(calibration_prior.calibration_id)
        validated_calibration_priors.append(calibration_prior)
    return validated_calibration_priors


def _validate_positive_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    if not math.isfinite(value) or value <= 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be positive and finite",
            code="fixed_topology_relaxed_clock_positive_float_invalid",
            details={field_name: value},
        )
    return float(format(value, ".15g"))


def _validate_optional_positive_finite_float(
    *,
    value: float | None,
    field_name: str,
    owner_name: str,
) -> float | None:
    if value is None:
        return None
    return _validate_positive_finite_float(
        value=value,
        field_name=field_name,
        owner_name=owner_name,
    )
