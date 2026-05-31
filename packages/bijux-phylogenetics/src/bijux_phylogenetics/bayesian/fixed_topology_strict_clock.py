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
    build_strict_clock_rate_model,
    evaluate_strict_clock_tree_log_prior,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsProposal,
    MetropolisHastingsRunReport,
    build_metropolis_hastings_proposal,
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

FIXED_TOPOLOGY_STRICT_CLOCK_MODELS = ("strict-clock",)
_MODEL_NAME = "strict-clock"
_CLOCK_MODEL_PARAMETER_NAME = "clock-model"
_GLOBAL_CLOCK_RATE_PARAMETER_NAME = "global-clock-rate"


@dataclass(frozen=True, slots=True)
class FixedTopologyStrictClockModelDefinition:
    """One validated fixed-topology strict-clock posterior model definition."""

    time_tree_prior: (
        YuleTreePriorModel
        | BirthDeathTreePriorModel
        | ConstantPopulationCoalescentPriorModel
        | SkylineCoalescentPriorModel
    )
    global_clock_rate_prior: ClockModelScalarPriorModel
    calibration_priors: list[CalibrationPriorDefinition]
    initial_global_clock_rate: float | None = None
    branch_length_tolerance: float = 1e-12


@dataclass(frozen=True, slots=True)
class FixedTopologyStrictClockProposalSchedule:
    """One validated proposal schedule for strict-clock dating sampling."""

    global_clock_rate_move_weight: float
    global_clock_rate_log_scale_standard_deviation: float


@dataclass(frozen=True, slots=True)
class FixedTopologyStrictClockPosteriorRow:
    """One sampled posterior row from a fixed-topology strict-clock chain."""

    sample_index: int
    iteration_index: int
    topology_id: str
    model_name: str
    total_log_prior: float
    log_likelihood: float
    posterior_log_score: float
    prior_component_log_priors: dict[str, float]
    global_clock_rate: float
    root_age: float
    tree_newick: str


@dataclass(frozen=True, slots=True)
class FixedTopologyStrictClockRateSummary:
    """One posterior summary for the sampled global strict-clock rate."""

    sample_count: int
    posterior_mean: float
    hpd_95_lower: float
    hpd_95_upper: float
    minimum_rate: float
    maximum_rate: float


@dataclass(frozen=True, slots=True)
class FixedTopologyStrictClockNodeAgeSummary:
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
class FixedTopologyStrictClockRunReport:
    """One completed fixed-topology strict-clock posterior run."""

    model_definition: FixedTopologyStrictClockModelDefinition
    proposal_schedule: FixedTopologyStrictClockProposalSchedule
    chain_report: MetropolisHastingsRunReport
    posterior_rows: list[FixedTopologyStrictClockPosteriorRow]
    clock_rate_summary: FixedTopologyStrictClockRateSummary
    node_age_summaries: list[FixedTopologyStrictClockNodeAgeSummary]


def build_fixed_topology_strict_clock_model_definition(
    *,
    time_tree_prior: (
        YuleTreePriorModel
        | BirthDeathTreePriorModel
        | ConstantPopulationCoalescentPriorModel
        | SkylineCoalescentPriorModel
    ),
    global_clock_rate_prior: ClockModelScalarPriorModel,
    calibration_priors: list[CalibrationPriorDefinition] | None = None,
    initial_global_clock_rate: float | None = None,
    branch_length_tolerance: float = 1e-12,
) -> FixedTopologyStrictClockModelDefinition:
    """Build one validated fixed-topology strict-clock posterior model."""
    validated_time_tree_prior = _validate_time_tree_prior_model(time_tree_prior)
    validated_global_clock_rate_prior = _validate_clock_model_scalar_prior(
        prior_model=global_clock_rate_prior,
        field_name="global_clock_rate_prior",
        owner_name="fixed-topology strict-clock posterior model",
    )
    if validated_global_clock_rate_prior.family == "fixed":
        raise PhylogeneticsError(
            "fixed-topology strict-clock posterior model requires one non-fixed global_clock_rate_prior so clock rate remains sampled",
            code="fixed_topology_strict_clock_global_clock_rate_prior_fixed",
        )
    return FixedTopologyStrictClockModelDefinition(
        time_tree_prior=validated_time_tree_prior,
        global_clock_rate_prior=validated_global_clock_rate_prior,
        calibration_priors=_validate_calibration_priors(calibration_priors or []),
        initial_global_clock_rate=_validate_optional_positive_finite_float(
            value=initial_global_clock_rate,
            field_name="initial_global_clock_rate",
            owner_name="fixed-topology strict-clock posterior model",
        ),
        branch_length_tolerance=_validate_nonnegative_finite_float(
            value=branch_length_tolerance,
            field_name="branch_length_tolerance",
            owner_name="fixed-topology strict-clock posterior model",
        ),
    )


def build_fixed_topology_strict_clock_proposal_schedule(
    *,
    model_definition: FixedTopologyStrictClockModelDefinition,
    global_clock_rate_move_weight: float,
    global_clock_rate_log_scale_standard_deviation: float,
) -> FixedTopologyStrictClockProposalSchedule:
    """Build one validated proposal schedule for strict-clock dating."""
    if not isinstance(model_definition, FixedTopologyStrictClockModelDefinition):
        raise PhylogeneticsError(
            "fixed-topology strict-clock proposal schedule requires one FixedTopologyStrictClockModelDefinition",
            code="fixed_topology_strict_clock_proposal_schedule_model_definition_type_invalid",
        )
    return FixedTopologyStrictClockProposalSchedule(
        global_clock_rate_move_weight=_validate_positive_finite_float(
            value=global_clock_rate_move_weight,
            field_name="global_clock_rate_move_weight",
            owner_name="fixed-topology strict-clock proposal schedule",
        ),
        global_clock_rate_log_scale_standard_deviation=_validate_positive_finite_float(
            value=global_clock_rate_log_scale_standard_deviation,
            field_name="global_clock_rate_log_scale_standard_deviation",
            owner_name="fixed-topology strict-clock proposal schedule",
        ),
    )


def run_fixed_topology_strict_clock_metropolis_hastings(
    *,
    substitution_tree: PhyloTree,
    model_definition: FixedTopologyStrictClockModelDefinition,
    proposal_schedule: FixedTopologyStrictClockProposalSchedule,
    iteration_count: int,
    sample_every: int = 1,
    seed: int = 0,
) -> FixedTopologyStrictClockRunReport:
    """Run one fixed-topology strict-clock posterior sampler."""
    if not isinstance(substitution_tree, PhyloTree):
        raise PhylogeneticsError(
            "fixed-topology strict-clock posterior runner requires one PhyloTree substitution_tree",
            code="fixed_topology_strict_clock_substitution_tree_type_invalid",
        )
    if not isinstance(model_definition, FixedTopologyStrictClockModelDefinition):
        raise PhylogeneticsError(
            "fixed-topology strict-clock posterior runner requires one FixedTopologyStrictClockModelDefinition",
            code="fixed_topology_strict_clock_model_definition_type_invalid",
        )
    if not isinstance(proposal_schedule, FixedTopologyStrictClockProposalSchedule):
        raise PhylogeneticsError(
            "fixed-topology strict-clock posterior runner requires one FixedTopologyStrictClockProposalSchedule",
            code="fixed_topology_strict_clock_proposal_schedule_type_invalid",
        )
    fixed_substitution_tree = substitution_tree.copy()
    fixed_substitution_tree.rooted = substitution_tree.rooted
    fixed_substitution_tree.refresh()
    if fixed_substitution_tree.rooted is not True:
        raise PhylogeneticsError(
            "fixed-topology strict-clock posterior runner requires one rooted substitution_tree",
            code="fixed_topology_strict_clock_substitution_tree_rooting_invalid",
        )
    initial_global_clock_rate = _resolve_initial_global_clock_rate(
        substitution_tree=fixed_substitution_tree,
        model_definition=model_definition,
    )
    initial_tree = _derive_dated_tree_from_strict_clock_rate(
        substitution_tree=fixed_substitution_tree,
        global_clock_rate=initial_global_clock_rate,
    )
    initial_model_parameters = build_bayesian_model_parameter_state(
        categorical_parameters={_CLOCK_MODEL_PARAMETER_NAME: _MODEL_NAME},
        scalar_parameters={
            _GLOBAL_CLOCK_RATE_PARAMETER_NAME: initial_global_clock_rate,
        },
    )
    initial_state = score_bayesian_phylogenetic_state(
        tree=initial_tree,
        model_parameters=initial_model_parameters,
        update_prior_components=lambda state: (
            _build_fixed_topology_strict_clock_prior_components(
                state=state,
                model_definition=model_definition,
                fixed_topology_id=None,
            )
        ),
        update_log_likelihood=lambda state: (
            _evaluate_fixed_topology_strict_clock_log_likelihood(
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
            _propose_fixed_topology_strict_clock_state(
                current_state=current_state,
                rng=rng,
                substitution_tree=fixed_substitution_tree,
                model_definition=model_definition,
                proposal_schedule=proposal_schedule,
            )
        ),
        update_prior_components=lambda state: (
            _build_fixed_topology_strict_clock_prior_components(
                state=state,
                model_definition=model_definition,
                fixed_topology_id=fixed_topology_id,
            )
        ),
        update_log_likelihood=lambda state: (
            _evaluate_fixed_topology_strict_clock_log_likelihood(
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
    posterior_rows = _build_fixed_topology_strict_clock_posterior_rows(
        chain_report=chain_report,
        fixed_topology_id=fixed_topology_id,
    )
    return FixedTopologyStrictClockRunReport(
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        chain_report=chain_report,
        posterior_rows=posterior_rows,
        clock_rate_summary=_build_fixed_topology_strict_clock_rate_summary(
            chain_report=chain_report,
        ),
        node_age_summaries=_build_fixed_topology_strict_clock_node_age_summaries(
            chain_report=chain_report,
            fixed_topology_id=fixed_topology_id,
        ),
    )


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
        "fixed-topology strict-clock posterior model requires one supported rooted time-tree prior model",
        code="fixed_topology_strict_clock_time_tree_prior_type_invalid",
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
            code="fixed_topology_strict_clock_clock_prior_type_invalid",
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
                "fixed-topology strict-clock posterior model requires every calibration prior to be one CalibrationPriorDefinition",
                code="fixed_topology_strict_clock_calibration_prior_type_invalid",
            )
        if calibration_prior.calibration_id in seen_calibration_ids:
            raise PhylogeneticsError(
                "fixed-topology strict-clock posterior model requires calibration ids to be unique",
                code="fixed_topology_strict_clock_calibration_prior_duplicate",
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
            code="fixed_topology_strict_clock_positive_float_invalid",
            details={field_name: value},
        )
    return float(format(value, ".15g"))


def _validate_nonnegative_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    if not math.isfinite(value) or value < 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be non-negative and finite",
            code="fixed_topology_strict_clock_nonnegative_float_invalid",
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


def _resolve_initial_global_clock_rate(
    *,
    substitution_tree: PhyloTree,
    model_definition: FixedTopologyStrictClockModelDefinition,
) -> float:
    if model_definition.initial_global_clock_rate is not None:
        return model_definition.initial_global_clock_rate
    substitution_root_depth = _compute_root_age(substitution_tree)
    calibration_age_centers = [
        calibration_age_center
        for calibration_age_center in (
            _representative_calibration_age(prior_definition)
            for prior_definition in model_definition.calibration_priors
            if prior_definition.node_kind == "root"
        )
        if calibration_age_center is not None and calibration_age_center > 0.0
    ]
    if calibration_age_centers:
        return float(
            format(substitution_root_depth / mean(calibration_age_centers), ".15g")
        )
    return 1.0


def _representative_calibration_age(
    prior_definition: CalibrationPriorDefinition,
) -> float | None:
    if prior_definition.fixed_age is not None:
        return prior_definition.fixed_age
    if (
        prior_definition.minimum_age is not None
        and prior_definition.maximum_age is not None
    ):
        return (prior_definition.minimum_age + prior_definition.maximum_age) / 2.0
    if prior_definition.mean_age is not None:
        return prior_definition.mean_age
    if (
        prior_definition.offset_age is not None
        and prior_definition.exponential_mean is not None
    ):
        return prior_definition.offset_age + prior_definition.exponential_mean
    return None


def _build_fixed_topology_strict_clock_prior_components(
    *,
    state: BayesianPhylogeneticState,
    model_definition: FixedTopologyStrictClockModelDefinition,
    fixed_topology_id: str | None,
) -> list[BayesianPriorComponentState]:
    global_clock_rate = _require_fixed_topology_strict_clock_state_consistency(
        state=state,
        model_definition=model_definition,
        fixed_topology_id=fixed_topology_id,
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
            component_name="clock:global-clock-rate",
            family=model_definition.global_clock_rate_prior.family,
            log_prior=evaluate_clock_model_scalar_log_prior(
                parameter_value=global_clock_rate,
                prior_model=model_definition.global_clock_rate_prior,
                parameter_name=_GLOBAL_CLOCK_RATE_PARAMETER_NAME,
            ),
            parameter_values=model_definition.global_clock_rate_prior.parameter_values(),
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


def _evaluate_fixed_topology_strict_clock_log_likelihood(
    *,
    state: BayesianPhylogeneticState,
    substitution_tree: PhyloTree,
    model_definition: FixedTopologyStrictClockModelDefinition,
    fixed_topology_id: str | None,
) -> float:
    global_clock_rate = _require_fixed_topology_strict_clock_state_consistency(
        state=state,
        model_definition=model_definition,
        fixed_topology_id=fixed_topology_id,
    )
    report = evaluate_strict_clock_tree_log_prior(
        substitution_tree,
        state.tree.to_tree(),
        build_strict_clock_rate_model(
            global_clock_rate=global_clock_rate,
            branch_length_tolerance=model_definition.branch_length_tolerance,
        ),
    )
    return report.total_log_prior


def _propose_fixed_topology_strict_clock_state(
    *,
    current_state: BayesianPhylogeneticState,
    rng,
    substitution_tree: PhyloTree,
    model_definition: FixedTopologyStrictClockModelDefinition,
    proposal_schedule: FixedTopologyStrictClockProposalSchedule,
) -> MetropolisHastingsProposal:
    current_clock_rate = current_state.model_parameters.scalar_parameters.get(
        _GLOBAL_CLOCK_RATE_PARAMETER_NAME
    )
    if current_clock_rate is None:
        return build_metropolis_hastings_proposal(
            changed_fields=(
                "scalar_parameters.global-clock-rate",
                "tree.branch_lengths",
            ),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=(
                "fixed-topology strict-clock proposal requires one "
                "'global-clock-rate' scalar parameter"
            ),
        )
    try:
        validated_current_clock_rate = _validate_positive_finite_float(
            value=current_clock_rate,
            field_name=_GLOBAL_CLOCK_RATE_PARAMETER_NAME,
            owner_name="fixed-topology strict-clock proposal",
        )
        scale_factor = math.exp(
            rng.gauss(
                0.0,
                proposal_schedule.global_clock_rate_log_scale_standard_deviation,
            )
        )
        validated_proposed_clock_rate = _validate_positive_finite_float(
            value=validated_current_clock_rate * scale_factor,
            field_name=_GLOBAL_CLOCK_RATE_PARAMETER_NAME,
            owner_name="fixed-topology strict-clock proposal",
        )
    except (OverflowError, PhylogeneticsError) as error:
        return build_metropolis_hastings_proposal(
            changed_fields=(
                "scalar_parameters.global-clock-rate",
                "tree.branch_lengths",
            ),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=str(error),
        )
    proposed_tree = _derive_dated_tree_from_strict_clock_rate(
        substitution_tree=substitution_tree,
        global_clock_rate=validated_proposed_clock_rate,
    )
    invalid_reason = _invalid_prior_support_reason(
        tree=proposed_tree,
        model_definition=model_definition,
    )
    if invalid_reason is not None:
        return build_metropolis_hastings_proposal(
            changed_fields=(
                "scalar_parameters.global-clock-rate",
                "tree.branch_lengths",
            ),
            log_forward_density=0.0,
            log_reverse_density=0.0,
            is_valid=False,
            invalid_reason=invalid_reason,
        )
    proposed_scalar_parameters = dict(current_state.model_parameters.scalar_parameters)
    proposed_scalar_parameters[_GLOBAL_CLOCK_RATE_PARAMETER_NAME] = (
        validated_proposed_clock_rate
    )
    log_forward_density = _lognormal_scaling_density(
        current_value=validated_current_clock_rate,
        proposed_value=validated_proposed_clock_rate,
        log_scale_standard_deviation=(
            proposal_schedule.global_clock_rate_log_scale_standard_deviation
        ),
    )
    log_reverse_density = _lognormal_scaling_density(
        current_value=validated_proposed_clock_rate,
        proposed_value=validated_current_clock_rate,
        log_scale_standard_deviation=(
            proposal_schedule.global_clock_rate_log_scale_standard_deviation
        ),
    )
    return build_metropolis_hastings_proposal(
        changed_fields=(
            "scalar_parameters.global-clock-rate",
            "tree.branch_lengths",
        ),
        log_forward_density=log_forward_density,
        log_reverse_density=log_reverse_density,
        is_valid=True,
        proposed_tree=proposed_tree,
        proposed_model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters=current_state.model_parameters.categorical_parameters,
            scalar_parameters=proposed_scalar_parameters,
            vector_parameters=current_state.model_parameters.vector_parameters,
        ),
    )


def _build_fixed_topology_strict_clock_posterior_rows(
    *,
    chain_report: MetropolisHastingsRunReport,
    fixed_topology_id: str,
) -> list[FixedTopologyStrictClockPosteriorRow]:
    posterior_rows: list[FixedTopologyStrictClockPosteriorRow] = []
    for sample_index, state in enumerate(chain_report.sampled_states):
        if state.tree.topology_id != fixed_topology_id:
            raise PhylogeneticsError(
                "fixed-topology strict-clock posterior trace detected one topology change in sampled states",
                code="fixed_topology_strict_clock_trace_topology_changed",
            )
        tree = state.tree.to_tree()
        posterior_rows.append(
            FixedTopologyStrictClockPosteriorRow(
                sample_index=sample_index,
                iteration_index=sample_index * chain_report.sample_every,
                topology_id=state.tree.topology_id,
                model_name=_MODEL_NAME,
                total_log_prior=state.total_log_prior,
                log_likelihood=state.log_likelihood,
                posterior_log_score=state.posterior_log_score,
                prior_component_log_priors={
                    component.component_name: component.log_prior
                    for component in state.prior_components
                },
                global_clock_rate=state.model_parameters.scalar_parameters[
                    _GLOBAL_CLOCK_RATE_PARAMETER_NAME
                ],
                root_age=_compute_root_age(tree),
                tree_newick=tree.to_newick(),
            )
        )
    return posterior_rows


def _build_fixed_topology_strict_clock_rate_summary(
    *,
    chain_report: MetropolisHastingsRunReport,
) -> FixedTopologyStrictClockRateSummary:
    rates = [
        state.model_parameters.scalar_parameters[_GLOBAL_CLOCK_RATE_PARAMETER_NAME]
        for state in chain_report.sampled_states
    ]
    interval = highest_posterior_density_interval(rates)
    return FixedTopologyStrictClockRateSummary(
        sample_count=len(rates),
        posterior_mean=float(format(mean(rates), ".15g")),
        hpd_95_lower=float(format(interval[0], ".15g")),
        hpd_95_upper=float(format(interval[1], ".15g")),
        minimum_rate=float(format(min(rates), ".15g")),
        maximum_rate=float(format(max(rates), ".15g")),
    )


def _build_fixed_topology_strict_clock_node_age_summaries(
    *,
    chain_report: MetropolisHastingsRunReport,
    fixed_topology_id: str,
) -> list[FixedTopologyStrictClockNodeAgeSummary]:
    node_metadata_by_id: dict[str, tuple[str, str | None, list[str]]] = {}
    node_age_series_by_id: dict[str, list[float]] = {}
    for state in chain_report.sampled_states:
        if state.tree.topology_id != fixed_topology_id:
            raise PhylogeneticsError(
                "fixed-topology strict-clock node-age summary requires a fixed topology across sampled states",
                code="fixed_topology_strict_clock_node_age_summary_topology_changed",
            )
        tree = state.tree.to_tree()
        node_age_by_id = _compute_internal_node_ages(tree)
        for node in tree.iter_nodes(order="preorder"):
            if node.is_leaf():
                continue
            if node.node_id is None:
                raise PhylogeneticsError(
                    "fixed-topology strict-clock node-age summary requires stable node ids",
                    code="fixed_topology_strict_clock_node_age_summary_node_id_missing",
                )
            node_metadata_by_id[node.node_id] = (
                "root" if node is tree.root else "internal",
                node.name,
                list(node.descendant_taxa),
            )
            node_age_series_by_id.setdefault(node.node_id, []).append(
                node_age_by_id[node.node_id]
            )
    summaries: list[FixedTopologyStrictClockNodeAgeSummary] = []
    for node_id, ages in sorted(
        node_age_series_by_id.items(),
        key=lambda item: (
            0 if node_metadata_by_id[item[0]][0] == "root" else 1,
            len(node_metadata_by_id[item[0]][2]),
            node_metadata_by_id[item[0]][2],
        ),
    ):
        interval = highest_posterior_density_interval(ages)
        summaries.append(
            FixedTopologyStrictClockNodeAgeSummary(
                node_id=node_id,
                node_kind=node_metadata_by_id[node_id][0],
                child_name=node_metadata_by_id[node_id][1],
                descendant_taxa=node_metadata_by_id[node_id][2],
                sample_count=len(ages),
                posterior_mean=float(format(mean(ages), ".15g")),
                hpd_95_lower=float(format(interval[0], ".15g")),
                hpd_95_upper=float(format(interval[1], ".15g")),
                minimum_age=float(format(min(ages), ".15g")),
                maximum_age=float(format(max(ages), ".15g")),
            )
        )
    return summaries


def _require_fixed_topology_strict_clock_state_consistency(
    *,
    state: BayesianPhylogeneticState,
    model_definition: FixedTopologyStrictClockModelDefinition,
    fixed_topology_id: str | None,
) -> float:
    model_name = state.model_parameters.categorical_parameters.get(
        _CLOCK_MODEL_PARAMETER_NAME
    )
    if model_name != _MODEL_NAME:
        raise PhylogeneticsError(
            "fixed-topology strict-clock posterior model requires every sampled state to preserve the configured clock-model label",
            code="fixed_topology_strict_clock_state_model_label_invalid",
            details={
                "observed_model_name": model_name,
                "expected_model_name": _MODEL_NAME,
            },
        )
    if fixed_topology_id is not None and state.tree.topology_id != fixed_topology_id:
        raise PhylogeneticsError(
            "fixed-topology strict-clock posterior model requires every sampled state to preserve the initial rooted topology",
            code="fixed_topology_strict_clock_state_topology_changed",
            details={
                "expected_topology_id": fixed_topology_id,
                "observed_topology_id": state.tree.topology_id,
            },
        )
    try:
        return state.model_parameters.scalar_parameters[
            _GLOBAL_CLOCK_RATE_PARAMETER_NAME
        ]
    except KeyError as error:
        raise PhylogeneticsError(
            "fixed-topology strict-clock posterior model requires every sampled state to include one global-clock-rate scalar parameter",
            code="fixed_topology_strict_clock_state_scalar_parameter_missing",
            details={"parameter_name": str(error)},
        ) from error


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


def _invalid_prior_support_reason(
    *,
    tree: PhyloTree,
    model_definition: FixedTopologyStrictClockModelDefinition,
) -> str | None:
    time_tree_prior_report = _evaluate_time_tree_log_prior(
        tree=tree,
        prior_model=model_definition.time_tree_prior,
    )
    if not math.isfinite(time_tree_prior_report.log_prior):
        return "strict-clock proposal moved dated tree outside time-tree prior support"
    if not model_definition.calibration_priors:
        return None
    calibration_prior_report = evaluate_calibration_tree_log_prior(
        tree,
        model_definition.calibration_priors,
    )
    if any(
        not math.isfinite(row.log_prior_contribution)
        for row in calibration_prior_report.calibration_rows
    ):
        return (
            "strict-clock proposal moved dated tree outside calibration prior support"
        )
    return None


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


def _derive_dated_tree_from_strict_clock_rate(
    *,
    substitution_tree: PhyloTree,
    global_clock_rate: float,
) -> PhyloTree:
    validated_global_clock_rate = _validate_positive_finite_float(
        value=global_clock_rate,
        field_name=_GLOBAL_CLOCK_RATE_PARAMETER_NAME,
        owner_name="fixed-topology strict-clock tree derivation",
    )
    substitution_root_age = _compute_root_age(substitution_tree)
    if substitution_root_age <= 0.0:
        raise PhylogeneticsError(
            "fixed-topology strict-clock posterior runner requires one substitution_tree with positive ultrametric depth",
            code="fixed_topology_strict_clock_substitution_tree_depth_invalid",
        )
    dated_tree = substitution_tree.copy()
    dated_tree.rooted = substitution_tree.rooted
    for _parent, child in dated_tree.iter_edges():
        if child.branch_length is None:
            raise PhylogeneticsError(
                "fixed-topology strict-clock tree derivation requires explicit substitution branch lengths on every edge",
                code="fixed_topology_strict_clock_substitution_branch_length_missing",
            )
        child.branch_length = float(
            format(child.branch_length / validated_global_clock_rate, ".15g")
        )
    dated_tree.refresh()
    return dated_tree


def _compute_internal_node_ages(tree: PhyloTree) -> dict[str, float]:
    root_age = _compute_root_age(tree)
    depth_by_node_id: dict[str, float] = {}

    def visit(node: TreeNode, current_depth: float) -> None:
        if node.node_id is None:
            raise PhylogeneticsError(
                "fixed-topology strict-clock posterior summaries require stable node ids",
                code="fixed_topology_strict_clock_summary_node_id_missing",
            )
        depth_by_node_id[node.node_id] = current_depth
        for child in node.children:
            visit(child, current_depth + float(child.branch_length or 0.0))

    visit(tree.root, 0.0)
    return {
        node_id: float(format(root_age - depth, ".15g"))
        for node_id, depth in depth_by_node_id.items()
        if depth <= root_age + 1e-12
    }


def _compute_root_age(tree: PhyloTree) -> float:
    tip_depth_by_label = _tip_depth_by_label(tree)
    ultrametric_summary = summarize_ultrametric_tip_depths(
        tip_depth_by_label,
        tolerance=APE_ULTRAMETRIC_TOLERANCE,
    )
    if ultrametric_summary.ultrametric is not True:
        raise PhylogeneticsError(
            "fixed-topology strict-clock posterior requires one ultrametric tree",
            code="fixed_topology_strict_clock_tree_not_ultrametric",
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
                    "fixed-topology strict-clock posterior requires every tip to have a name",
                    code="fixed_topology_strict_clock_tip_name_missing",
                )
            tip_depth_by_label[node.name] = current_depth
            return
        for child in node.children:
            if child.branch_length is None:
                raise PhylogeneticsError(
                    "fixed-topology strict-clock posterior requires explicit branch lengths on every edge",
                    code="fixed_topology_strict_clock_branch_length_missing",
                )
            visit(child, current_depth + child.branch_length)

    visit(tree.root, 0.0)
    if not tip_depth_by_label:
        raise PhylogeneticsError(
            "fixed-topology strict-clock posterior requires at least one named tip",
            code="fixed_topology_strict_clock_tip_missing",
        )
    return tip_depth_by_label


def _lognormal_scaling_density(
    *,
    current_value: float,
    proposed_value: float,
    log_scale_standard_deviation: float,
) -> float:
    log_ratio = math.log(proposed_value / current_value)
    variance = log_scale_standard_deviation * log_scale_standard_deviation
    return float(
        format(
            (-math.log(proposed_value * log_scale_standard_deviation))
            - 0.5 * math.log(2.0 * math.pi)
            - (log_ratio * log_ratio) / (2.0 * variance),
            ".15g",
        )
    )
