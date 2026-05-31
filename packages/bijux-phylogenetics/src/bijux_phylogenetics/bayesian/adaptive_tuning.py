from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import math
import random

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsProposal,
    MetropolisHastingsRunReport,
    MetropolisHastingsStepRow,
    build_metropolis_hastings_proposal,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

_ADAPTIVE_TUNING_ACTIONS = ("decrease", "frozen", "hold", "increase")


@dataclass(frozen=True, slots=True)
class AdaptiveTuningController:
    """One validated burn-in tuning policy for a scalar proposal scale."""

    proposal_name: str
    scale_parameter_name: str
    initial_scale: float
    target_acceptance_rate: float
    burnin_iteration_count: int
    adaptation_window_size: int
    decrease_factor: float
    increase_factor: float
    minimum_scale: float
    maximum_scale: float


@dataclass(frozen=True, slots=True)
class AdaptiveTuningWindowRow:
    """One window-level adaptive tuning summary."""

    window_index: int
    window_start_iteration: int
    window_end_iteration: int
    within_burnin: bool
    attempted_count: int
    accepted_count: int
    acceptance_rate: float
    target_acceptance_rate: float
    scale_before_window: float
    scale_after_window: float
    action: str


@dataclass(frozen=True, slots=True)
class AdaptiveTuningReport:
    """One completed adaptive tuning report for a Metropolis-Hastings run."""

    proposal_name: str
    scale_parameter_name: str
    initial_scale: float
    final_scale: float
    target_acceptance_rate: float
    burnin_iteration_count: int
    adaptation_window_size: int
    freeze_iteration_index: int
    burnin_sample_count: int
    retained_sample_count: int
    window_rows: list[AdaptiveTuningWindowRow]


@dataclass(frozen=True, slots=True)
class AdaptiveMetropolisHastingsRunReport:
    """One completed adaptive Metropolis-Hastings run with a tuning report."""

    chain_report: MetropolisHastingsRunReport
    burnin_iteration_count: int
    freeze_iteration_index: int
    burnin_sample_count: int
    retained_sample_count: int
    retained_sampled_states: list[BayesianPhylogeneticState]
    tuning_report: AdaptiveTuningReport


AdaptiveTuningProposal = Callable[
    [BayesianPhylogeneticState, random.Random, float],
    MetropolisHastingsProposal,
]
AdaptivePriorUpdate = Callable[
    [BayesianPhylogeneticState],
    list[BayesianPriorComponentState],
]
AdaptiveLikelihoodUpdate = Callable[[BayesianPhylogeneticState], float]


def build_adaptive_tuning_controller(
    *,
    proposal_name: str,
    scale_parameter_name: str,
    initial_scale: float,
    target_acceptance_rate: float,
    burnin_iteration_count: int,
    adaptation_window_size: int,
    decrease_factor: float = 0.5,
    increase_factor: float = 2.0,
    minimum_scale: float = 1e-6,
    maximum_scale: float = 1e6,
) -> AdaptiveTuningController:
    """Build one validated controller for burn-in-only proposal adaptation."""
    validated_proposal_name = _validate_nonblank_name(
        value=proposal_name,
        field_name="proposal_name",
        owner_name="adaptive tuning controller",
    )
    validated_scale_parameter_name = _validate_nonblank_name(
        value=scale_parameter_name,
        field_name="scale_parameter_name",
        owner_name="adaptive tuning controller",
    )
    validated_initial_scale = _validate_positive_finite_float(
        value=initial_scale,
        field_name="initial_scale",
        owner_name="adaptive tuning controller",
    )
    validated_target_acceptance_rate = _validate_probability_rate(
        value=target_acceptance_rate,
        field_name="target_acceptance_rate",
        owner_name="adaptive tuning controller",
    )
    validated_burnin_iteration_count = _validate_nonnegative_integer(
        value=burnin_iteration_count,
        field_name="burnin_iteration_count",
        owner_name="adaptive tuning controller",
    )
    validated_adaptation_window_size = _validate_positive_integer(
        value=adaptation_window_size,
        field_name="adaptation_window_size",
        owner_name="adaptive tuning controller",
    )
    validated_decrease_factor = _validate_open_interval_float(
        value=decrease_factor,
        field_name="decrease_factor",
        owner_name="adaptive tuning controller",
        lower_bound=0.0,
        upper_bound=1.0,
    )
    validated_increase_factor = _validate_greater_than_float(
        value=increase_factor,
        field_name="increase_factor",
        owner_name="adaptive tuning controller",
        lower_bound=1.0,
    )
    validated_minimum_scale = _validate_positive_finite_float(
        value=minimum_scale,
        field_name="minimum_scale",
        owner_name="adaptive tuning controller",
    )
    validated_maximum_scale = _validate_positive_finite_float(
        value=maximum_scale,
        field_name="maximum_scale",
        owner_name="adaptive tuning controller",
    )
    if validated_minimum_scale > validated_maximum_scale:
        raise PhylogeneticsError(
            "adaptive tuning controller requires 'minimum_scale' to be less than or equal to 'maximum_scale'",
            code="adaptive_tuning_scale_bounds_invalid",
        )
    if (
        not validated_minimum_scale
        <= validated_initial_scale
        <= validated_maximum_scale
    ):
        raise PhylogeneticsError(
            "adaptive tuning controller requires 'initial_scale' to lie within the configured scale bounds",
            code="adaptive_tuning_initial_scale_out_of_bounds",
            details={
                "initial_scale": validated_initial_scale,
                "minimum_scale": validated_minimum_scale,
                "maximum_scale": validated_maximum_scale,
            },
        )
    return AdaptiveTuningController(
        proposal_name=validated_proposal_name,
        scale_parameter_name=validated_scale_parameter_name,
        initial_scale=validated_initial_scale,
        target_acceptance_rate=validated_target_acceptance_rate,
        burnin_iteration_count=validated_burnin_iteration_count,
        adaptation_window_size=validated_adaptation_window_size,
        decrease_factor=validated_decrease_factor,
        increase_factor=validated_increase_factor,
        minimum_scale=validated_minimum_scale,
        maximum_scale=validated_maximum_scale,
    )


def build_adaptive_tuning_window_row(
    *,
    window_index: int,
    window_start_iteration: int,
    window_end_iteration: int,
    within_burnin: bool,
    attempted_count: int,
    accepted_count: int,
    target_acceptance_rate: float,
    scale_before_window: float,
    scale_after_window: float,
    action: str,
) -> AdaptiveTuningWindowRow:
    """Build one validated adaptive tuning window summary row."""
    validated_window_index = _validate_positive_integer(
        value=window_index,
        field_name="window_index",
        owner_name="adaptive tuning window row",
    )
    validated_window_start_iteration = _validate_positive_integer(
        value=window_start_iteration,
        field_name="window_start_iteration",
        owner_name="adaptive tuning window row",
    )
    validated_window_end_iteration = _validate_positive_integer(
        value=window_end_iteration,
        field_name="window_end_iteration",
        owner_name="adaptive tuning window row",
    )
    if validated_window_start_iteration > validated_window_end_iteration:
        raise PhylogeneticsError(
            "adaptive tuning window row requires 'window_start_iteration' to be less than or equal to 'window_end_iteration'",
            code="adaptive_tuning_window_iteration_range_invalid",
        )
    validated_attempted_count = _validate_positive_integer(
        value=attempted_count,
        field_name="attempted_count",
        owner_name="adaptive tuning window row",
    )
    validated_accepted_count = _validate_nonnegative_integer(
        value=accepted_count,
        field_name="accepted_count",
        owner_name="adaptive tuning window row",
    )
    if validated_accepted_count > validated_attempted_count:
        raise PhylogeneticsError(
            "adaptive tuning window row requires 'accepted_count' to be less than or equal to 'attempted_count'",
            code="adaptive_tuning_window_acceptance_count_invalid",
        )
    validated_target_acceptance_rate = _validate_probability_rate(
        value=target_acceptance_rate,
        field_name="target_acceptance_rate",
        owner_name="adaptive tuning window row",
    )
    validated_scale_before_window = _validate_positive_finite_float(
        value=scale_before_window,
        field_name="scale_before_window",
        owner_name="adaptive tuning window row",
    )
    validated_scale_after_window = _validate_positive_finite_float(
        value=scale_after_window,
        field_name="scale_after_window",
        owner_name="adaptive tuning window row",
    )
    validated_action = _validate_adaptive_tuning_action(action)
    return AdaptiveTuningWindowRow(
        window_index=validated_window_index,
        window_start_iteration=validated_window_start_iteration,
        window_end_iteration=validated_window_end_iteration,
        within_burnin=bool(within_burnin),
        attempted_count=validated_attempted_count,
        accepted_count=validated_accepted_count,
        acceptance_rate=validated_accepted_count / validated_attempted_count,
        target_acceptance_rate=validated_target_acceptance_rate,
        scale_before_window=validated_scale_before_window,
        scale_after_window=validated_scale_after_window,
        action=validated_action,
    )


def build_adaptive_tuning_report(
    *,
    controller: AdaptiveTuningController,
    freeze_iteration_index: int,
    burnin_sample_count: int,
    retained_sample_count: int,
    window_rows: list[AdaptiveTuningWindowRow],
) -> AdaptiveTuningReport:
    """Build one validated adaptive tuning report."""
    validated_freeze_iteration_index = _validate_positive_integer(
        value=freeze_iteration_index,
        field_name="freeze_iteration_index",
        owner_name="adaptive tuning report",
    )
    expected_freeze_iteration_index = controller.burnin_iteration_count + 1
    if validated_freeze_iteration_index != expected_freeze_iteration_index:
        raise PhylogeneticsError(
            "adaptive tuning report requires one freeze iteration immediately after burn-in",
            code="adaptive_tuning_report_freeze_iteration_invalid",
            details={
                "freeze_iteration_index": validated_freeze_iteration_index,
                "expected_freeze_iteration_index": expected_freeze_iteration_index,
            },
        )
    validated_burnin_sample_count = _validate_nonnegative_integer(
        value=burnin_sample_count,
        field_name="burnin_sample_count",
        owner_name="adaptive tuning report",
    )
    validated_retained_sample_count = _validate_nonnegative_integer(
        value=retained_sample_count,
        field_name="retained_sample_count",
        owner_name="adaptive tuning report",
    )
    if not window_rows:
        raise PhylogeneticsError(
            "adaptive tuning report requires at least one window row",
            code="adaptive_tuning_report_window_rows_empty",
        )
    _validate_adaptive_tuning_window_rows(
        controller=controller,
        window_rows=window_rows,
    )
    return AdaptiveTuningReport(
        proposal_name=controller.proposal_name,
        scale_parameter_name=controller.scale_parameter_name,
        initial_scale=controller.initial_scale,
        final_scale=window_rows[-1].scale_after_window,
        target_acceptance_rate=controller.target_acceptance_rate,
        burnin_iteration_count=controller.burnin_iteration_count,
        adaptation_window_size=controller.adaptation_window_size,
        freeze_iteration_index=validated_freeze_iteration_index,
        burnin_sample_count=validated_burnin_sample_count,
        retained_sample_count=validated_retained_sample_count,
        window_rows=list(window_rows),
    )


def run_adaptive_tuned_metropolis_hastings_sampler(
    *,
    initial_state: BayesianPhylogeneticState,
    propose_state: AdaptiveTuningProposal,
    tuning_controller: AdaptiveTuningController,
    update_prior_components: AdaptivePriorUpdate,
    update_log_likelihood: AdaptiveLikelihoodUpdate,
    iteration_count: int,
    sample_every: int = 1,
    seed: int = 0,
) -> AdaptiveMetropolisHastingsRunReport:
    """Run one burn-in-adaptive Metropolis-Hastings chain with a frozen post-burn-in scale."""
    validated_iteration_count = _validate_positive_integer(
        value=iteration_count,
        field_name="iteration_count",
        owner_name="adaptive metropolis-hastings sampler",
    )
    validated_sample_every = _validate_positive_integer(
        value=sample_every,
        field_name="sample_every",
        owner_name="adaptive metropolis-hastings sampler",
    )
    validated_seed = _validate_integer_seed(seed)
    if tuning_controller.burnin_iteration_count >= validated_iteration_count:
        raise PhylogeneticsError(
            "adaptive metropolis-hastings sampler requires 'iteration_count' to exceed the tuning-controller burn-in so retained posterior samples begin after tuning freezes",
            code="adaptive_metropolis_hastings_burnin_not_before_retention",
            details={
                "iteration_count": validated_iteration_count,
                "burnin_iteration_count": tuning_controller.burnin_iteration_count,
            },
        )
    rng = random.Random(validated_seed)  # nosec B311
    current_state = initial_state
    sampled_states: list[BayesianPhylogeneticState] = [current_state]
    retained_sampled_states: list[BayesianPhylogeneticState] = []
    step_rows: list[MetropolisHastingsStepRow] = []
    tuning_window_rows: list[AdaptiveTuningWindowRow] = []
    accepted_count = 0
    current_scale = tuning_controller.initial_scale
    scale_before_window = current_scale
    window_start_iteration = 1
    window_attempted_count = 0
    window_accepted_count = 0

    for iteration_index in range(1, validated_iteration_count + 1):
        previous_state = current_state
        raw_proposal = propose_state(current_state, rng, current_scale)
        proposal = _validate_metropolis_hastings_proposal_instance(raw_proposal)
        proposed_state: BayesianPhylogeneticState | None = None
        log_acceptance_ratio: float | None = None
        accepted = False
        if proposal.is_valid:
            proposed_state = score_bayesian_phylogenetic_state(
                tree=proposal.proposed_tree,
                model_parameters=proposal.proposed_model_parameters,
                update_prior_components=update_prior_components,
                update_log_likelihood=update_log_likelihood,
            )
            log_acceptance_ratio = (
                proposed_state.posterior_log_score
                - current_state.posterior_log_score
                + proposal.log_hastings_ratio
            )
            accepted = _accept_metropolis_hastings_proposal(
                log_acceptance_ratio=log_acceptance_ratio,
                rng=rng,
            )
            if accepted:
                current_state = proposed_state
                accepted_count += 1
        window_attempted_count += 1
        if accepted:
            window_accepted_count += 1
        if iteration_index % validated_sample_every == 0:
            sampled_states.append(current_state)
            if iteration_index > tuning_controller.burnin_iteration_count:
                retained_sampled_states.append(current_state)
        step_rows.append(
            MetropolisHastingsStepRow(
                iteration_index=iteration_index,
                proposal_changed_fields=proposal.changed_fields,
                proposal_valid=proposal.is_valid,
                proposal_invalid_reason=proposal.invalid_reason,
                log_forward_density=proposal.log_forward_density,
                log_reverse_density=proposal.log_reverse_density,
                accepted=accepted,
                log_hastings_ratio=proposal.log_hastings_ratio,
                current_posterior_log_score=previous_state.posterior_log_score,
                proposed_posterior_log_score=(
                    proposed_state.posterior_log_score
                    if proposed_state is not None
                    else None
                ),
                log_acceptance_ratio=log_acceptance_ratio,
                recorded_posterior_log_score=current_state.posterior_log_score,
            )
        )
        if _should_close_tuning_window(
            iteration_index=iteration_index,
            iteration_count=validated_iteration_count,
            burnin_iteration_count=tuning_controller.burnin_iteration_count,
            adaptation_window_size=tuning_controller.adaptation_window_size,
            window_attempted_count=window_attempted_count,
        ):
            current_scale, action = _resolve_next_adaptive_scale(
                controller=tuning_controller,
                window_end_iteration=iteration_index,
                current_scale=current_scale,
                accepted_count=window_accepted_count,
                attempted_count=window_attempted_count,
            )
            tuning_window_rows.append(
                build_adaptive_tuning_window_row(
                    window_index=len(tuning_window_rows) + 1,
                    window_start_iteration=window_start_iteration,
                    window_end_iteration=iteration_index,
                    within_burnin=(
                        iteration_index <= tuning_controller.burnin_iteration_count
                    ),
                    attempted_count=window_attempted_count,
                    accepted_count=window_accepted_count,
                    target_acceptance_rate=tuning_controller.target_acceptance_rate,
                    scale_before_window=scale_before_window,
                    scale_after_window=current_scale,
                    action=action,
                )
            )
            scale_before_window = current_scale
            window_start_iteration = iteration_index + 1
            window_attempted_count = 0
            window_accepted_count = 0

    chain_report = MetropolisHastingsRunReport(
        iteration_count=validated_iteration_count,
        sample_every=validated_sample_every,
        seed=validated_seed,
        accepted_count=accepted_count,
        rejected_count=validated_iteration_count - accepted_count,
        acceptance_rate=accepted_count / validated_iteration_count,
        initial_state=initial_state,
        final_state=current_state,
        sampled_states=sampled_states,
        step_rows=step_rows,
    )
    tuning_report = build_adaptive_tuning_report(
        controller=tuning_controller,
        freeze_iteration_index=tuning_controller.burnin_iteration_count + 1,
        burnin_sample_count=len(sampled_states) - len(retained_sampled_states),
        retained_sample_count=len(retained_sampled_states),
        window_rows=tuning_window_rows,
    )
    return AdaptiveMetropolisHastingsRunReport(
        chain_report=chain_report,
        burnin_iteration_count=tuning_controller.burnin_iteration_count,
        freeze_iteration_index=tuning_report.freeze_iteration_index,
        burnin_sample_count=tuning_report.burnin_sample_count,
        retained_sample_count=tuning_report.retained_sample_count,
        retained_sampled_states=retained_sampled_states,
        tuning_report=tuning_report,
    )


def _validate_nonblank_name(
    *,
    value: str,
    field_name: str,
    owner_name: str,
) -> str:
    if not isinstance(value, str) or not value.strip():
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one nonblank name",
            code="adaptive_tuning_name_invalid",
            details={"field_name": field_name},
        )
    return value.strip()


def _validate_integer_seed(seed: int) -> int:
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise PhylogeneticsError(
            "adaptive metropolis-hastings sampler requires 'seed' to be one integer",
            code="adaptive_tuning_seed_invalid",
        )
    return seed


def _validate_adaptive_tuning_action(action: str) -> str:
    validated_action = _validate_nonblank_name(
        value=action,
        field_name="action",
        owner_name="adaptive tuning window row",
    )
    if validated_action not in _ADAPTIVE_TUNING_ACTIONS:
        raise PhylogeneticsError(
            "adaptive tuning window row requires one supported action",
            code="adaptive_tuning_action_invalid",
            details={
                "action": action,
                "allowed_actions": list(_ADAPTIVE_TUNING_ACTIONS),
            },
        )
    return validated_action


def _validate_positive_integer(
    *,
    value: int,
    field_name: str,
    owner_name: str,
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one integer",
            code="adaptive_tuning_integer_required",
            details={"field_name": field_name},
        )
    if value <= 0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be positive",
            code="adaptive_tuning_positive_integer_required",
            details={"field_name": field_name},
        )
    return value


def _validate_nonnegative_integer(
    *,
    value: int,
    field_name: str,
    owner_name: str,
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one integer",
            code="adaptive_tuning_integer_required",
            details={"field_name": field_name},
        )
    if value < 0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be nonnegative",
            code="adaptive_tuning_nonnegative_integer_required",
            details={"field_name": field_name},
        )
    return value


def _validate_positive_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be numeric",
            code="adaptive_tuning_float_required",
            details={"field_name": field_name},
        )
    normalized_value = float(value)
    if not math.isfinite(normalized_value):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be finite",
            code="adaptive_tuning_finite_float_required",
            details={"field_name": field_name},
        )
    if normalized_value <= 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be strictly positive",
            code="adaptive_tuning_positive_float_required",
            details={"field_name": field_name},
        )
    return normalized_value


def _validate_probability_rate(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be numeric",
            code="adaptive_tuning_float_required",
            details={"field_name": field_name},
        )
    normalized_value = float(value)
    if not math.isfinite(normalized_value):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be finite",
            code="adaptive_tuning_finite_float_required",
            details={"field_name": field_name},
        )
    if not 0.0 < normalized_value < 1.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to lie strictly between 0 and 1",
            code="adaptive_tuning_probability_rate_invalid",
            details={"field_name": field_name},
        )
    return normalized_value


def _validate_open_interval_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
    lower_bound: float,
    upper_bound: float,
) -> float:
    normalized_value = _validate_positive_finite_float(
        value=value,
        field_name=field_name,
        owner_name=owner_name,
    )
    if not lower_bound < normalized_value < upper_bound:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to lie strictly between {lower_bound} and {upper_bound}",
            code="adaptive_tuning_interval_float_invalid",
            details={
                "field_name": field_name,
                "lower_bound": lower_bound,
                "upper_bound": upper_bound,
            },
        )
    return normalized_value


def _validate_greater_than_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
    lower_bound: float,
) -> float:
    normalized_value = _validate_positive_finite_float(
        value=value,
        field_name=field_name,
        owner_name=owner_name,
    )
    if normalized_value <= lower_bound:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be greater than {lower_bound}",
            code="adaptive_tuning_lower_bound_invalid",
            details={"field_name": field_name, "lower_bound": lower_bound},
        )
    return normalized_value


def _validate_adaptive_tuning_window_rows(
    *,
    controller: AdaptiveTuningController,
    window_rows: list[AdaptiveTuningWindowRow],
) -> None:
    previous_window_end_iteration = 0
    adaptation_ended = False
    for expected_window_index, window_row in enumerate(window_rows, start=1):
        if window_row.window_index != expected_window_index:
            raise PhylogeneticsError(
                "adaptive tuning report requires consecutive window indexes",
                code="adaptive_tuning_window_index_invalid",
                details={
                    "expected_window_index": expected_window_index,
                    "observed_window_index": window_row.window_index,
                },
            )
        if window_row.window_start_iteration != previous_window_end_iteration + 1:
            raise PhylogeneticsError(
                "adaptive tuning report requires contiguous window iteration ranges",
                code="adaptive_tuning_window_contiguity_invalid",
                details={
                    "expected_window_start_iteration": previous_window_end_iteration
                    + 1,
                    "observed_window_start_iteration": window_row.window_start_iteration,
                },
            )
        if not math.isclose(
            window_row.target_acceptance_rate,
            controller.target_acceptance_rate,
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise PhylogeneticsError(
                "adaptive tuning report requires every window row to use the controller target acceptance rate",
                code="adaptive_tuning_window_target_rate_invalid",
            )
        if window_row.within_burnin:
            if adaptation_ended:
                raise PhylogeneticsError(
                    "adaptive tuning report cannot resume burn-in windows after tuning has frozen",
                    code="adaptive_tuning_window_phase_invalid",
                )
        else:
            adaptation_ended = True
            if window_row.action != "frozen":
                raise PhylogeneticsError(
                    "adaptive tuning report requires post-burn-in windows to record the frozen action",
                    code="adaptive_tuning_window_frozen_action_missing",
                )
            if not math.isclose(
                window_row.scale_before_window,
                window_row.scale_after_window,
                rel_tol=0.0,
                abs_tol=1e-12,
            ):
                raise PhylogeneticsError(
                    "adaptive tuning report requires post-burn-in windows to preserve one frozen scale",
                    code="adaptive_tuning_window_frozen_scale_invalid",
                )
        previous_window_end_iteration = window_row.window_end_iteration


def _validate_metropolis_hastings_proposal_instance(
    proposal: MetropolisHastingsProposal,
) -> MetropolisHastingsProposal:
    if not isinstance(proposal, MetropolisHastingsProposal):
        raise PhylogeneticsError(
            "adaptive metropolis-hastings sampler requires the proposal callback to return one MetropolisHastingsProposal",
            code="adaptive_tuning_proposal_type_invalid",
        )
    return build_metropolis_hastings_proposal(
        changed_fields=proposal.changed_fields,
        log_forward_density=proposal.log_forward_density,
        log_reverse_density=proposal.log_reverse_density,
        is_valid=proposal.is_valid,
        invalid_reason=proposal.invalid_reason,
        proposed_tree=proposal.proposed_tree,
        proposed_model_parameters=proposal.proposed_model_parameters,
    )


def _accept_metropolis_hastings_proposal(
    *,
    log_acceptance_ratio: float,
    rng: random.Random,
) -> bool:
    if log_acceptance_ratio >= 0.0:
        return True
    return math.log(rng.random()) <= log_acceptance_ratio


def _should_close_tuning_window(
    *,
    iteration_index: int,
    iteration_count: int,
    burnin_iteration_count: int,
    adaptation_window_size: int,
    window_attempted_count: int,
) -> bool:
    return window_attempted_count >= adaptation_window_size or iteration_index in (
        burnin_iteration_count,
        iteration_count,
    )


def _resolve_next_adaptive_scale(
    *,
    controller: AdaptiveTuningController,
    window_end_iteration: int,
    current_scale: float,
    accepted_count: int,
    attempted_count: int,
) -> tuple[float, str]:
    acceptance_rate = accepted_count / attempted_count
    if window_end_iteration > controller.burnin_iteration_count:
        return current_scale, "frozen"
    if acceptance_rate < controller.target_acceptance_rate:
        return (
            max(controller.minimum_scale, current_scale * controller.decrease_factor),
            "decrease",
        )
    if acceptance_rate > controller.target_acceptance_rate:
        return (
            min(controller.maximum_scale, current_scale * controller.increase_factor),
            "increase",
        )
    return current_scale, "hold"
