from __future__ import annotations

from dataclasses import dataclass
import math
from statistics import mean

from bijux_phylogenetics.bayesian.calibration_priors import (
    CalibrationPriorDefinition,
    evaluate_calibration_tree_log_prior,
)
from bijux_phylogenetics.bayesian.clock_model_priors import (
    ClockModelScalarPriorModel,
    evaluate_clock_model_scalar_log_prior,
)
from bijux_phylogenetics.bayesian.clock_models import (
    build_relaxed_lognormal_clock_model,
    evaluate_relaxed_lognormal_clock_tree_log_prior,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsRunReport,
    propose_clock_rate_move,
    propose_global_tree_height_scaling_move,
    propose_node_height_sliding_move,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.posterior_sets.diagnostics import (
    highest_posterior_density_interval,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
    build_bayesian_prior_component_state,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    BirthDeathTreePriorModel,
    ConstantPopulationCoalescentPriorModel,
    SkylineCoalescentPriorModel,
    YuleTreePriorModel,
    evaluate_birth_death_tree_log_prior,
    evaluate_constant_population_coalescent_tree_log_prior,
    evaluate_skyline_coalescent_tree_log_prior,
    evaluate_yule_tree_log_prior,
)
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    APE_ULTRAMETRIC_TOLERANCE,
    summarize_ultrametric_tip_depths,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXED_TOPOLOGY_RELAXED_CLOCK_MODELS = ("relaxed-lognormal",)
_RELAXED_CLOCK_RATE_POLICIES = ("independent", "autocorrelated")
_MODEL_NAME = "relaxed-lognormal"
_CLOCK_MODEL_PARAMETER_NAME = "clock-model"
_CLOCK_RATE_POLICY_PARAMETER_NAME = "clock-rate-policy"
_MEAN_CLOCK_RATE_PARAMETER_NAME = "mean-clock-rate"
_LOG_STANDARD_DEVIATION_PARAMETER_NAME = "log-standard-deviation"


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


def run_fixed_topology_relaxed_clock_metropolis_hastings(
    *,
    substitution_tree: PhyloTree,
    initial_dated_tree: PhyloTree,
    model_definition: FixedTopologyRelaxedClockModelDefinition,
    proposal_schedule: FixedTopologyRelaxedClockProposalSchedule,
    iteration_count: int,
    sample_every: int = 1,
    seed: int = 0,
) -> FixedTopologyRelaxedClockRunReport:
    """Run one fixed-topology relaxed-clock posterior sampler."""
    if not isinstance(substitution_tree, PhyloTree):
        raise PhylogeneticsError(
            "fixed-topology relaxed-clock posterior runner requires one PhyloTree substitution_tree",
            code="fixed_topology_relaxed_clock_substitution_tree_type_invalid",
        )
    if not isinstance(initial_dated_tree, PhyloTree):
        raise PhylogeneticsError(
            "fixed-topology relaxed-clock posterior runner requires one PhyloTree initial_dated_tree",
            code="fixed_topology_relaxed_clock_dated_tree_type_invalid",
        )
    if not isinstance(model_definition, FixedTopologyRelaxedClockModelDefinition):
        raise PhylogeneticsError(
            "fixed-topology relaxed-clock posterior runner requires one FixedTopologyRelaxedClockModelDefinition",
            code="fixed_topology_relaxed_clock_model_definition_type_invalid",
        )
    if not isinstance(proposal_schedule, FixedTopologyRelaxedClockProposalSchedule):
        raise PhylogeneticsError(
            "fixed-topology relaxed-clock posterior runner requires one FixedTopologyRelaxedClockProposalSchedule",
            code="fixed_topology_relaxed_clock_proposal_schedule_type_invalid",
        )
    fixed_substitution_tree = substitution_tree.copy()
    fixed_substitution_tree.rooted = substitution_tree.rooted
    fixed_substitution_tree.refresh()
    fixed_dated_tree = initial_dated_tree.copy()
    fixed_dated_tree.rooted = initial_dated_tree.rooted
    fixed_dated_tree.refresh()
    if fixed_substitution_tree.rooted is not True:
        raise PhylogeneticsError(
            "fixed-topology relaxed-clock posterior runner requires one rooted substitution_tree",
            code="fixed_topology_relaxed_clock_substitution_tree_rooting_invalid",
        )
    if fixed_dated_tree.rooted is not True:
        raise PhylogeneticsError(
            "fixed-topology relaxed-clock posterior runner requires one rooted initial_dated_tree",
            code="fixed_topology_relaxed_clock_dated_tree_rooting_invalid",
        )
    initial_model_parameters = _build_initial_model_parameters(
        substitution_tree=fixed_substitution_tree,
        initial_dated_tree=fixed_dated_tree,
        model_definition=model_definition,
    )
    initial_state = score_bayesian_phylogenetic_state(
        tree=fixed_dated_tree,
        model_parameters=initial_model_parameters,
        update_prior_components=lambda state: (
            _build_fixed_topology_relaxed_clock_prior_components(
                state=state,
                model_definition=model_definition,
                fixed_topology_id=None,
            )
        ),
        update_log_likelihood=lambda state: (
            _evaluate_fixed_topology_relaxed_clock_log_likelihood(
                state=state,
                substitution_tree=fixed_substitution_tree,
                model_definition=model_definition,
                fixed_topology_id=None,
            )
        ),
    )
    fixed_topology_id = initial_state.tree.topology_id
    chain_report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=lambda current_state, rng: (
            _propose_fixed_topology_relaxed_clock_state(
                current_state=current_state,
                rng=rng,
                proposal_schedule=proposal_schedule,
            )
        ),
        update_prior_components=lambda state: (
            _build_fixed_topology_relaxed_clock_prior_components(
                state=state,
                model_definition=model_definition,
                fixed_topology_id=fixed_topology_id,
            )
        ),
        update_log_likelihood=lambda state: (
            _evaluate_fixed_topology_relaxed_clock_log_likelihood(
                state=state,
                substitution_tree=fixed_substitution_tree,
                model_definition=model_definition,
                fixed_topology_id=fixed_topology_id,
            )
        ),
        iteration_count=iteration_count,
        sample_every=sample_every,
        seed=seed,
    )
    posterior_rows = _build_fixed_topology_relaxed_clock_posterior_rows(
        chain_report=chain_report,
        fixed_topology_id=fixed_topology_id,
        rate_policy=model_definition.rate_policy,
    )
    branch_rate_summaries = _build_fixed_topology_relaxed_clock_branch_rate_summaries(
        chain_report=chain_report,
        substitution_tree=fixed_substitution_tree,
        model_definition=model_definition,
        fixed_topology_id=fixed_topology_id,
    )
    node_age_summaries = _build_fixed_topology_relaxed_clock_node_age_summaries(
        chain_report=chain_report,
        fixed_topology_id=fixed_topology_id,
    )
    return FixedTopologyRelaxedClockRunReport(
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        chain_report=chain_report,
        posterior_rows=posterior_rows,
        branch_rate_summaries=branch_rate_summaries,
        node_age_summaries=node_age_summaries,
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


def _build_initial_model_parameters(
    *,
    substitution_tree: PhyloTree,
    initial_dated_tree: PhyloTree,
    model_definition: FixedTopologyRelaxedClockModelDefinition,
):
    initial_mean_clock_rate = (
        model_definition.initial_mean_clock_rate
        if model_definition.initial_mean_clock_rate is not None
        else _estimate_initial_mean_clock_rate(
            substitution_tree=substitution_tree,
            dated_tree=initial_dated_tree,
            rate_policy=model_definition.rate_policy,
        )
    )
    initial_log_standard_deviation = (
        model_definition.initial_log_standard_deviation
        if model_definition.initial_log_standard_deviation is not None
        else _estimate_initial_log_standard_deviation(
            substitution_tree=substitution_tree,
            dated_tree=initial_dated_tree,
            rate_policy=model_definition.rate_policy,
            mean_clock_rate=initial_mean_clock_rate,
        )
    )
    return build_bayesian_model_parameter_state(
        categorical_parameters={
            _CLOCK_MODEL_PARAMETER_NAME: _MODEL_NAME,
            _CLOCK_RATE_POLICY_PARAMETER_NAME: model_definition.rate_policy,
        },
        scalar_parameters={
            _MEAN_CLOCK_RATE_PARAMETER_NAME: initial_mean_clock_rate,
            _LOG_STANDARD_DEVIATION_PARAMETER_NAME: initial_log_standard_deviation,
        },
    )


def _estimate_initial_mean_clock_rate(
    *,
    substitution_tree: PhyloTree,
    dated_tree: PhyloTree,
    rate_policy: str,
) -> float:
    branch_rows = _evaluate_relaxed_clock_branch_rows(
        substitution_tree=substitution_tree,
        dated_tree=dated_tree,
        rate_policy=rate_policy,
        mean_clock_rate=1.0,
        log_standard_deviation=1.0,
    )
    return float(format(mean(row.branch_rate for row in branch_rows), ".15g"))


def _estimate_initial_log_standard_deviation(
    *,
    substitution_tree: PhyloTree,
    dated_tree: PhyloTree,
    rate_policy: str,
    mean_clock_rate: float,
) -> float:
    branch_rows = _evaluate_relaxed_clock_branch_rows(
        substitution_tree=substitution_tree,
        dated_tree=dated_tree,
        rate_policy=rate_policy,
        mean_clock_rate=mean_clock_rate,
        log_standard_deviation=1.0,
    )
    log_rates = [math.log(row.branch_rate) for row in branch_rows]
    average_log_rate = mean(log_rates)
    variance = mean(
        (log_rate - average_log_rate) * (log_rate - average_log_rate)
        for log_rate in log_rates
    )
    return float(format(max(math.sqrt(variance), 0.05), ".15g"))


def _build_fixed_topology_relaxed_clock_prior_components(
    *,
    state: BayesianPhylogeneticState,
    model_definition: FixedTopologyRelaxedClockModelDefinition,
    fixed_topology_id: str | None,
) -> list[BayesianPriorComponentState]:
    mean_clock_rate, log_standard_deviation = (
        _require_fixed_topology_relaxed_clock_state_consistency(
            state=state,
            model_definition=model_definition,
            fixed_topology_id=fixed_topology_id,
        )
    )
    dated_tree = state.tree.to_tree()
    time_tree_prior_report = _evaluate_time_tree_log_prior(
        tree=dated_tree,
        prior_model=model_definition.time_tree_prior,
    )
    prior_components = [
        build_bayesian_prior_component_state(
            component_name="time-tree",
            family=time_tree_prior_report.family,
            log_prior=time_tree_prior_report.log_prior,
            parameter_values=_time_tree_prior_parameter_values(
                model_definition.time_tree_prior
            ),
        ),
        build_bayesian_prior_component_state(
            component_name="clock:mean-clock-rate",
            family=model_definition.mean_clock_rate_prior.family,
            log_prior=evaluate_clock_model_scalar_log_prior(
                parameter_value=mean_clock_rate,
                prior_model=model_definition.mean_clock_rate_prior,
                parameter_name=_MEAN_CLOCK_RATE_PARAMETER_NAME,
            ),
            parameter_values=model_definition.mean_clock_rate_prior.parameter_values(),
        ),
        build_bayesian_prior_component_state(
            component_name="clock:log-standard-deviation",
            family=model_definition.log_standard_deviation_prior.family,
            log_prior=evaluate_clock_model_scalar_log_prior(
                parameter_value=log_standard_deviation,
                prior_model=model_definition.log_standard_deviation_prior,
                parameter_name=_LOG_STANDARD_DEVIATION_PARAMETER_NAME,
            ),
            parameter_values=model_definition.log_standard_deviation_prior.parameter_values(),
        ),
    ]
    if model_definition.calibration_priors:
        calibration_prior_report = evaluate_calibration_tree_log_prior(
            dated_tree,
            model_definition.calibration_priors,
        )
        prior_components.extend(
            build_bayesian_prior_component_state(
                component_name=f"calibration:{row.calibration_id}",
                family=row.family,
                log_prior=row.log_prior_contribution,
                parameter_values=row.parameter_values,
            )
            for row in calibration_prior_report.calibration_rows
        )
    return prior_components


def _evaluate_fixed_topology_relaxed_clock_log_likelihood(
    *,
    state: BayesianPhylogeneticState,
    substitution_tree: PhyloTree,
    model_definition: FixedTopologyRelaxedClockModelDefinition,
    fixed_topology_id: str | None,
) -> float:
    mean_clock_rate, log_standard_deviation = (
        _require_fixed_topology_relaxed_clock_state_consistency(
            state=state,
            model_definition=model_definition,
            fixed_topology_id=fixed_topology_id,
        )
    )
    report = evaluate_relaxed_lognormal_clock_tree_log_prior(
        substitution_tree,
        state.tree.to_tree(),
        build_relaxed_lognormal_clock_model(
            rate_policy=model_definition.rate_policy,
            mean_clock_rate=mean_clock_rate,
            log_standard_deviation=log_standard_deviation,
        ),
    )
    return report.total_log_prior


def _propose_fixed_topology_relaxed_clock_state(
    *,
    current_state: BayesianPhylogeneticState,
    rng,
    proposal_schedule: FixedTopologyRelaxedClockProposalSchedule,
):
    weighted_moves = [
        (
            proposal_schedule.mean_clock_rate_move_weight,
            lambda: propose_clock_rate_move(
                current_state,
                rng,
                log_scale_standard_deviation=(
                    proposal_schedule.mean_clock_rate_log_scale_standard_deviation
                ),
                parameter_name=_MEAN_CLOCK_RATE_PARAMETER_NAME,
            ),
        ),
        (
            proposal_schedule.log_standard_deviation_move_weight,
            lambda: propose_clock_rate_move(
                current_state,
                rng,
                log_scale_standard_deviation=(
                    proposal_schedule.log_standard_deviation_log_scale_standard_deviation
                ),
                parameter_name=_LOG_STANDARD_DEVIATION_PARAMETER_NAME,
            ),
        ),
        (
            proposal_schedule.node_height_move_weight,
            lambda: propose_node_height_sliding_move(
                current_state,
                rng,
                height_slide_standard_deviation=(
                    proposal_schedule.node_height_slide_standard_deviation
                ),
            ),
        ),
        (
            proposal_schedule.tree_height_move_weight,
            lambda: propose_global_tree_height_scaling_move(
                current_state,
                rng,
                log_scale_standard_deviation=(
                    proposal_schedule.tree_height_log_scale_standard_deviation
                ),
            ),
        ),
    ]
    total_weight = math.fsum(weight for weight, _move in weighted_moves)
    move_threshold = rng.random() * total_weight
    cumulative_weight = 0.0
    for weight, move in weighted_moves:
        cumulative_weight += weight
        if move_threshold <= cumulative_weight:
            return move()
    return weighted_moves[-1][1]()


def _build_fixed_topology_relaxed_clock_posterior_rows(
    *,
    chain_report: MetropolisHastingsRunReport,
    fixed_topology_id: str,
    rate_policy: str,
) -> list[FixedTopologyRelaxedClockPosteriorRow]:
    posterior_rows: list[FixedTopologyRelaxedClockPosteriorRow] = []
    for sample_index, state in enumerate(chain_report.sampled_states):
        if state.tree.topology_id != fixed_topology_id:
            raise PhylogeneticsError(
                "fixed-topology relaxed-clock posterior trace detected one topology change in sampled states",
                code="fixed_topology_relaxed_clock_trace_topology_changed",
                details={
                    "expected_topology_id": fixed_topology_id,
                    "observed_topology_id": state.tree.topology_id,
                    "sample_index": sample_index,
                },
            )
        tree = state.tree.to_tree()
        posterior_rows.append(
            FixedTopologyRelaxedClockPosteriorRow(
                sample_index=sample_index,
                iteration_index=sample_index * chain_report.sample_every,
                topology_id=state.tree.topology_id,
                model_name=_MODEL_NAME,
                rate_policy=rate_policy,
                total_log_prior=state.total_log_prior,
                log_likelihood=state.log_likelihood,
                posterior_log_score=state.posterior_log_score,
                prior_component_log_priors={
                    component.component_name: component.log_prior
                    for component in state.prior_components
                },
                mean_clock_rate=state.model_parameters.scalar_parameters[
                    _MEAN_CLOCK_RATE_PARAMETER_NAME
                ],
                log_standard_deviation=state.model_parameters.scalar_parameters[
                    _LOG_STANDARD_DEVIATION_PARAMETER_NAME
                ],
                root_age=_compute_root_age(tree),
                tree_newick=tree.to_newick(),
            )
        )
    return posterior_rows


def _build_fixed_topology_relaxed_clock_branch_rate_summaries(
    *,
    chain_report: MetropolisHastingsRunReport,
    substitution_tree: PhyloTree,
    model_definition: FixedTopologyRelaxedClockModelDefinition,
    fixed_topology_id: str,
) -> list[FixedTopologyRelaxedClockBranchRateSummary]:
    branch_metadata_by_id: dict[str, tuple[str | None, list[str]]] = {}
    branch_rate_series_by_id: dict[str, list[float]] = {}
    for state in chain_report.sampled_states:
        if state.tree.topology_id != fixed_topology_id:
            raise PhylogeneticsError(
                "fixed-topology relaxed-clock branch-rate summary requires a fixed topology across sampled states",
                code="fixed_topology_relaxed_clock_branch_rate_summary_topology_changed",
            )
        branch_rows = _evaluate_relaxed_clock_branch_rows(
            substitution_tree=substitution_tree,
            dated_tree=state.tree.to_tree(),
            rate_policy=model_definition.rate_policy,
            mean_clock_rate=state.model_parameters.scalar_parameters[
                _MEAN_CLOCK_RATE_PARAMETER_NAME
            ],
            log_standard_deviation=state.model_parameters.scalar_parameters[
                _LOG_STANDARD_DEVIATION_PARAMETER_NAME
            ],
        )
        for branch_row in branch_rows:
            branch_metadata_by_id[branch_row.branch_id] = (
                branch_row.child_name,
                list(branch_row.descendant_taxa),
            )
            branch_rate_series_by_id.setdefault(branch_row.branch_id, []).append(
                branch_row.branch_rate
            )
    return [
        FixedTopologyRelaxedClockBranchRateSummary(
            branch_id=branch_id,
            child_name=branch_metadata_by_id[branch_id][0],
            descendant_taxa=branch_metadata_by_id[branch_id][1],
            sample_count=len(rates),
            posterior_mean=float(format(mean(rates), ".15g")),
            hpd_95_lower=float(
                format(highest_posterior_density_interval(rates)[0], ".15g")
            ),
            hpd_95_upper=float(
                format(highest_posterior_density_interval(rates)[1], ".15g")
            ),
            minimum_rate=float(format(min(rates), ".15g")),
            maximum_rate=float(format(max(rates), ".15g")),
        )
        for branch_id, rates in sorted(
            branch_rate_series_by_id.items(),
            key=lambda item: (
                len(branch_metadata_by_id[item[0]][1]),
                branch_metadata_by_id[item[0]][1],
            ),
        )
    ]


def _build_fixed_topology_relaxed_clock_node_age_summaries(
    *,
    chain_report: MetropolisHastingsRunReport,
    fixed_topology_id: str,
) -> list[FixedTopologyRelaxedClockNodeAgeSummary]:
    node_metadata_by_id: dict[str, tuple[str, str | None, list[str]]] = {}
    node_age_series_by_id: dict[str, list[float]] = {}
    for state in chain_report.sampled_states:
        if state.tree.topology_id != fixed_topology_id:
            raise PhylogeneticsError(
                "fixed-topology relaxed-clock node-age summary requires a fixed topology across sampled states",
                code="fixed_topology_relaxed_clock_node_age_summary_topology_changed",
            )
        tree = state.tree.to_tree()
        node_age_by_id = _compute_internal_node_ages(tree)
        for node in tree.iter_nodes(order="preorder"):
            if node.is_leaf():
                continue
            if node.node_id is None:
                raise PhylogeneticsError(
                    "fixed-topology relaxed-clock node-age summary requires stable node ids",
                    code="fixed_topology_relaxed_clock_node_age_summary_node_id_missing",
                )
            node_metadata_by_id[node.node_id] = (
                "root" if node is tree.root else "internal",
                node.name,
                list(node.descendant_taxa),
            )
            node_age_series_by_id.setdefault(node.node_id, []).append(
                node_age_by_id[node.node_id]
            )
    return [
        FixedTopologyRelaxedClockNodeAgeSummary(
            node_id=node_id,
            node_kind=node_metadata_by_id[node_id][0],
            child_name=node_metadata_by_id[node_id][1],
            descendant_taxa=node_metadata_by_id[node_id][2],
            sample_count=len(ages),
            posterior_mean=float(format(mean(ages), ".15g")),
            hpd_95_lower=float(
                format(highest_posterior_density_interval(ages)[0], ".15g")
            ),
            hpd_95_upper=float(
                format(highest_posterior_density_interval(ages)[1], ".15g")
            ),
            minimum_age=float(format(min(ages), ".15g")),
            maximum_age=float(format(max(ages), ".15g")),
        )
        for node_id, ages in sorted(
            node_age_series_by_id.items(),
            key=lambda item: (
                0 if node_metadata_by_id[item[0]][0] == "root" else 1,
                len(node_metadata_by_id[item[0]][2]),
                node_metadata_by_id[item[0]][2],
            ),
        )
    ]


def _require_fixed_topology_relaxed_clock_state_consistency(
    *,
    state: BayesianPhylogeneticState,
    model_definition: FixedTopologyRelaxedClockModelDefinition,
    fixed_topology_id: str | None,
) -> tuple[float, float]:
    model_name = state.model_parameters.categorical_parameters.get(
        _CLOCK_MODEL_PARAMETER_NAME
    )
    if model_name != _MODEL_NAME:
        raise PhylogeneticsError(
            "fixed-topology relaxed-clock posterior model requires every sampled state to preserve the configured clock-model label",
            code="fixed_topology_relaxed_clock_state_model_label_invalid",
            details={
                "observed_model_name": model_name,
                "expected_model_name": _MODEL_NAME,
            },
        )
    rate_policy = state.model_parameters.categorical_parameters.get(
        _CLOCK_RATE_POLICY_PARAMETER_NAME
    )
    if rate_policy != model_definition.rate_policy:
        raise PhylogeneticsError(
            "fixed-topology relaxed-clock posterior model requires every sampled state to preserve the configured clock-rate-policy label",
            code="fixed_topology_relaxed_clock_state_rate_policy_label_invalid",
            details={
                "observed_rate_policy": rate_policy,
                "expected_rate_policy": model_definition.rate_policy,
            },
        )
    if fixed_topology_id is not None and state.tree.topology_id != fixed_topology_id:
        raise PhylogeneticsError(
            "fixed-topology relaxed-clock posterior model requires every sampled state to preserve the initial rooted topology",
            code="fixed_topology_relaxed_clock_state_topology_changed",
            details={
                "expected_topology_id": fixed_topology_id,
                "observed_topology_id": state.tree.topology_id,
            },
        )
    try:
        mean_clock_rate = state.model_parameters.scalar_parameters[
            _MEAN_CLOCK_RATE_PARAMETER_NAME
        ]
        log_standard_deviation = state.model_parameters.scalar_parameters[
            _LOG_STANDARD_DEVIATION_PARAMETER_NAME
        ]
    except KeyError as error:
        raise PhylogeneticsError(
            "fixed-topology relaxed-clock posterior model requires every sampled state to include both relaxed-clock scalar parameters",
            code="fixed_topology_relaxed_clock_state_scalar_parameter_missing",
            details={"parameter_name": str(error)},
        ) from error
    return mean_clock_rate, log_standard_deviation


def _evaluate_time_tree_log_prior(
    *,
    tree: PhyloTree,
    prior_model: (
        YuleTreePriorModel
        | BirthDeathTreePriorModel
        | ConstantPopulationCoalescentPriorModel
        | SkylineCoalescentPriorModel
    ),
):
    if isinstance(prior_model, YuleTreePriorModel):
        return evaluate_yule_tree_log_prior(tree, prior_model)
    if isinstance(prior_model, BirthDeathTreePriorModel):
        return evaluate_birth_death_tree_log_prior(tree, prior_model)
    if isinstance(prior_model, ConstantPopulationCoalescentPriorModel):
        return evaluate_constant_population_coalescent_tree_log_prior(
            tree,
            prior_model,
        )
    if isinstance(prior_model, SkylineCoalescentPriorModel):
        return evaluate_skyline_coalescent_tree_log_prior(tree, prior_model)
    raise AssertionError("unsupported time-tree prior model reached evaluation")


def _time_tree_prior_parameter_values(
    prior_model: (
        YuleTreePriorModel
        | BirthDeathTreePriorModel
        | ConstantPopulationCoalescentPriorModel
        | SkylineCoalescentPriorModel
    ),
) -> dict[str, float]:
    if isinstance(prior_model, YuleTreePriorModel):
        return {"speciation_rate": prior_model.speciation_rate}
    if isinstance(prior_model, BirthDeathTreePriorModel):
        return {
            "speciation_rate": prior_model.speciation_rate,
            "extinction_rate": prior_model.extinction_rate,
            "sampling_fraction": prior_model.sampling_fraction,
        }
    if isinstance(prior_model, ConstantPopulationCoalescentPriorModel):
        return {"effective_population_size": prior_model.effective_population_size}
    if isinstance(prior_model, SkylineCoalescentPriorModel):
        parameter_values: dict[str, float] = {
            "epoch_count": float(len(prior_model.epochs))
        }
        for index, epoch in enumerate(prior_model.epochs, start=1):
            parameter_values[f"epoch_{index}_younger_boundary_age"] = (
                epoch.younger_boundary_age
            )
            if epoch.older_boundary_age is not None:
                parameter_values[f"epoch_{index}_older_boundary_age"] = (
                    epoch.older_boundary_age
                )
            parameter_values[f"epoch_{index}_effective_population_size"] = (
                epoch.effective_population_size
            )
        return parameter_values
    raise AssertionError("unsupported time-tree prior model reached parameter export")


def _evaluate_relaxed_clock_branch_rows(
    *,
    substitution_tree: PhyloTree,
    dated_tree: PhyloTree,
    rate_policy: str,
    mean_clock_rate: float,
    log_standard_deviation: float,
):
    report = evaluate_relaxed_lognormal_clock_tree_log_prior(
        substitution_tree,
        dated_tree,
        build_relaxed_lognormal_clock_model(
            rate_policy=rate_policy,
            mean_clock_rate=mean_clock_rate,
            log_standard_deviation=log_standard_deviation,
        ),
    )
    return report.branch_rows


def _compute_internal_node_ages(tree: PhyloTree) -> dict[str, float]:
    root_age = _compute_root_age(tree)
    depth_by_node_id: dict[str, float] = {}

    def visit(node: TreeNode, current_depth: float) -> None:
        if node.node_id is None:
            raise PhylogeneticsError(
                "fixed-topology relaxed-clock posterior summaries require stable node ids",
                code="fixed_topology_relaxed_clock_summary_node_id_missing",
            )
        depth_by_node_id[node.node_id] = current_depth
        for child in node.children:
            visit(child, current_depth + float(child.branch_length or 0.0))

    visit(tree.root, 0.0)
    return {
        node_id: float(format(root_age - depth, ".15g"))
        for node_id, depth in depth_by_node_id.items()
        if depth_by_node_id[node_id] <= root_age + 1e-12
    }


def _compute_root_age(tree: PhyloTree) -> float:
    tip_depth_by_label = _tip_depth_by_label(tree)
    ultrametric_summary = summarize_ultrametric_tip_depths(
        tip_depth_by_label,
        tolerance=APE_ULTRAMETRIC_TOLERANCE,
    )
    if ultrametric_summary.ultrametric is not True:
        raise PhylogeneticsError(
            "fixed-topology relaxed-clock posterior summaries require one ultrametric dated tree",
            code="fixed_topology_relaxed_clock_summary_tree_not_ultrametric",
            details={
                "minimum_tip_depth": ultrametric_summary.minimum_tip_depth,
                "maximum_tip_depth": ultrametric_summary.maximum_tip_depth,
            },
        )
    return float(format(ultrametric_summary.root_age, ".15g"))


def _tip_depth_by_label(tree: PhyloTree) -> dict[str, float]:
    tip_depth_by_label: dict[str, float] = {}

    def visit(node: TreeNode, current_depth: float) -> None:
        if node.is_leaf():
            if node.name is None:
                raise PhylogeneticsError(
                    "fixed-topology relaxed-clock posterior summaries require every tip to have a name",
                    code="fixed_topology_relaxed_clock_summary_tip_name_missing",
                )
            tip_depth_by_label[node.name] = current_depth
            return
        for child in node.children:
            visit(child, current_depth + float(child.branch_length or 0.0))

    visit(tree.root, 0.0)
    if not tip_depth_by_label:
        raise PhylogeneticsError(
            "fixed-topology relaxed-clock posterior summaries require at least one named tip",
            code="fixed_topology_relaxed_clock_summary_tip_missing",
        )
    return tip_depth_by_label
