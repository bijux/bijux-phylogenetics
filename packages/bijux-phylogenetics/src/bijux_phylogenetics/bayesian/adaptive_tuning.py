from __future__ import annotations

from dataclasses import dataclass
import math

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
    if not validated_minimum_scale <= validated_initial_scale <= validated_maximum_scale:
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
