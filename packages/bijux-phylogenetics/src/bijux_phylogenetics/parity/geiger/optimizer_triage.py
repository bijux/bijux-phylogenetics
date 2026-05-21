from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .runner import GeigerParityObservation, GeigerParityReport

_OBJECTIVE_ONLY_PARAMETERS = {"log_likelihood", "aic", "aicc"}


@dataclass(frozen=True, slots=True)
class GeigerOptimizerTriageRow:
    """One governed optimizer-disagreement triage row for a geiger parity case."""

    case_id: str
    function_name: str
    model_name: str
    status: str
    mismatch_reason: str | None
    mismatch_type: str
    parameter_surface_comparable: bool
    same_likelihood_within_tolerance: bool | None
    same_parameter_surface_within_tolerance: bool | None
    reference_log_likelihood: float | None
    bijux_log_likelihood: float | None
    objective_delta: float | None
    reference_parameter_name: str | None
    bijux_parameter_name: str | None
    reference_parameter_value: float | None
    bijux_parameter_value: float | None
    parameter_delta: float | None
    reference_boundary_detected: bool
    bijux_boundary_detected: bool
    boundary_solution_detected: bool
    reference_trace_row_count: int | None
    bijux_trace_row_count: int | None
    reference_local_optimum_count: int | None
    bijux_local_optimum_count: int | None
    reference_optimizer_trace: list[dict[str, object]] | None
    bijux_optimizer_trace: list[dict[str, object]] | None


def build_geiger_optimizer_triage_rows(
    observations: list[GeigerParityObservation],
) -> list[GeigerOptimizerTriageRow]:
    """Classify governed parity observations by optimizer-disagreement type."""

    return [_build_optimizer_triage_row(observation) for observation in observations]


def write_geiger_optimizer_triage_table(
    path: Path,
    report: GeigerParityReport,
) -> Path:
    """Write one row per governed geiger optimizer-triage observation."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_id",
                "function_name",
                "model_name",
                "status",
                "mismatch_reason",
                "mismatch_type",
                "parameter_surface_comparable",
                "same_likelihood_within_tolerance",
                "same_parameter_surface_within_tolerance",
                "reference_log_likelihood",
                "bijux_log_likelihood",
                "objective_delta",
                "reference_parameter_name",
                "bijux_parameter_name",
                "reference_parameter_value",
                "bijux_parameter_value",
                "parameter_delta",
                "reference_boundary_detected",
                "bijux_boundary_detected",
                "boundary_solution_detected",
                "reference_trace_row_count",
                "bijux_trace_row_count",
                "reference_local_optimum_count",
                "bijux_local_optimum_count",
                "reference_optimizer_trace",
                "bijux_optimizer_trace",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.optimizer_triage_rows:
            payload = asdict(row)
            payload["reference_optimizer_trace"] = json.dumps(
                row.reference_optimizer_trace,
                sort_keys=True,
            )
            payload["bijux_optimizer_trace"] = json.dumps(
                row.bijux_optimizer_trace,
                sort_keys=True,
            )
            writer.writerow(payload)
    return path


def _build_optimizer_triage_row(
    observation: GeigerParityObservation,
) -> GeigerOptimizerTriageRow:
    if observation.status == "skipped":
        return GeigerOptimizerTriageRow(
            case_id=observation.case_id,
            function_name=observation.function_name,
            model_name=observation.model_name,
            status=observation.status,
            mismatch_reason=observation.mismatch_reason,
            mismatch_type="reference_environment_unavailable",
            parameter_surface_comparable=False,
            same_likelihood_within_tolerance=None,
            same_parameter_surface_within_tolerance=None,
            reference_log_likelihood=None,
            bijux_log_likelihood=None,
            objective_delta=None,
            reference_parameter_name=None,
            bijux_parameter_name=None,
            reference_parameter_value=None,
            bijux_parameter_value=None,
            parameter_delta=None,
            reference_boundary_detected=False,
            bijux_boundary_detected=False,
            boundary_solution_detected=False,
            reference_trace_row_count=None,
            bijux_trace_row_count=None,
            reference_local_optimum_count=None,
            bijux_local_optimum_count=None,
            reference_optimizer_trace=None,
            bijux_optimizer_trace=None,
        )

    if observation.reference_error is not None:
        return _error_row(observation, mismatch_type="reference_execution_failure")
    if observation.bijux_error is not None:
        return _error_row(observation, mismatch_type="bijux_execution_failure")
    if observation.reference_summary is None or observation.bijux_summary is None:
        return _error_row(observation, mismatch_type="comparison_not_available")

    reference_log_likelihood = _optional_float(
        observation.reference_summary.get("log_likelihood")
    )
    bijux_log_likelihood = _optional_float(
        observation.bijux_summary.get("log_likelihood")
    )
    same_likelihood = _same_numeric_value(
        reference_log_likelihood,
        bijux_log_likelihood,
        tolerance=observation.tolerance,
    )
    objective_delta = (
        None
        if reference_log_likelihood is None or bijux_log_likelihood is None
        else abs(reference_log_likelihood - bijux_log_likelihood)
    )

    reference_parameter_name = _optional_text(
        observation.reference_summary.get("parameter_name")
    )
    bijux_parameter_name = _optional_text(
        observation.bijux_summary.get("parameter_name")
    )
    reference_parameter_value = _optional_float(
        observation.reference_summary.get("parameter_value")
    )
    bijux_parameter_value = _optional_float(
        observation.bijux_summary.get("parameter_value")
    )
    canonical_reference = _canonical_primary_parameter(
        observation.model_name,
        reference_parameter_name,
        reference_parameter_value,
    )
    canonical_bijux = _canonical_primary_parameter(
        observation.model_name,
        bijux_parameter_name,
        bijux_parameter_value,
    )
    parameter_delta = _parameter_delta(
        canonical_reference[1],
        canonical_bijux[1],
    )

    parameter_surface_comparable, same_parameter_surface = _parameter_surface_relation(
        observation,
    )
    reference_boundary_detected = _boundary_detected(observation.reference_summary)
    bijux_boundary_detected = _boundary_detected(observation.bijux_summary)
    boundary_solution_detected = reference_boundary_detected or bijux_boundary_detected

    reference_trace = _optimizer_trace(observation.reference_summary)
    bijux_trace = _optimizer_trace(observation.bijux_summary)
    reference_local_optimum_count = _local_optimum_count(
        reference_trace,
        model_name=observation.model_name,
        parameter_name=reference_parameter_name,
    )
    bijux_local_optimum_count = _local_optimum_count(
        bijux_trace,
        model_name=observation.model_name,
        parameter_name=bijux_parameter_name,
    )

    mismatch_type = _mismatch_type(
        same_likelihood=same_likelihood,
        parameter_surface_comparable=parameter_surface_comparable,
        same_parameter_surface=same_parameter_surface,
        boundary_solution_detected=boundary_solution_detected,
    )

    return GeigerOptimizerTriageRow(
        case_id=observation.case_id,
        function_name=observation.function_name,
        model_name=observation.model_name,
        status=observation.status,
        mismatch_reason=observation.mismatch_reason,
        mismatch_type=mismatch_type,
        parameter_surface_comparable=parameter_surface_comparable,
        same_likelihood_within_tolerance=same_likelihood,
        same_parameter_surface_within_tolerance=same_parameter_surface,
        reference_log_likelihood=reference_log_likelihood,
        bijux_log_likelihood=bijux_log_likelihood,
        objective_delta=objective_delta,
        reference_parameter_name=canonical_reference[0],
        bijux_parameter_name=canonical_bijux[0],
        reference_parameter_value=canonical_reference[1],
        bijux_parameter_value=canonical_bijux[1],
        parameter_delta=parameter_delta,
        reference_boundary_detected=reference_boundary_detected,
        bijux_boundary_detected=bijux_boundary_detected,
        boundary_solution_detected=boundary_solution_detected,
        reference_trace_row_count=None
        if reference_trace is None
        else len(reference_trace),
        bijux_trace_row_count=None if bijux_trace is None else len(bijux_trace),
        reference_local_optimum_count=reference_local_optimum_count,
        bijux_local_optimum_count=bijux_local_optimum_count,
        reference_optimizer_trace=reference_trace,
        bijux_optimizer_trace=bijux_trace,
    )


def _error_row(
    observation: GeigerParityObservation,
    *,
    mismatch_type: str,
) -> GeigerOptimizerTriageRow:
    return GeigerOptimizerTriageRow(
        case_id=observation.case_id,
        function_name=observation.function_name,
        model_name=observation.model_name,
        status=observation.status,
        mismatch_reason=observation.mismatch_reason,
        mismatch_type=mismatch_type,
        parameter_surface_comparable=False,
        same_likelihood_within_tolerance=None,
        same_parameter_surface_within_tolerance=None,
        reference_log_likelihood=None,
        bijux_log_likelihood=None,
        objective_delta=None,
        reference_parameter_name=None,
        bijux_parameter_name=None,
        reference_parameter_value=None,
        bijux_parameter_value=None,
        parameter_delta=None,
        reference_boundary_detected=False,
        bijux_boundary_detected=False,
        boundary_solution_detected=False,
        reference_trace_row_count=None,
        bijux_trace_row_count=None,
        reference_local_optimum_count=None,
        bijux_local_optimum_count=None,
        reference_optimizer_trace=None,
        bijux_optimizer_trace=None,
    )


def _optional_float(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        numeric = float(value)
        return None if math.isnan(numeric) else numeric
    return None


def _optional_text(value: object) -> str | None:
    return value if isinstance(value, str) and value else None


def _same_numeric_value(
    left: float | None,
    right: float | None,
    *,
    tolerance: float,
) -> bool | None:
    if left is None or right is None:
        return None
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def _canonical_primary_parameter(
    model_name: str,
    parameter_name: str | None,
    parameter_value: float | None,
) -> tuple[str | None, float | None]:
    if parameter_name is None or parameter_value is None:
        return parameter_name, parameter_value
    if model_name == "EB":
        if parameter_name == "a":
            return "rate_change", -parameter_value
        if parameter_name == "rate_change":
            return "rate_change", parameter_value
    return parameter_name, parameter_value


def _parameter_delta(
    left: float | None,
    right: float | None,
) -> float | None:
    if left is None or right is None:
        return None
    return abs(left - right)


def _parameter_surface_relation(
    observation: GeigerParityObservation,
) -> tuple[bool, bool | None]:
    reference_rows = observation.reference_rows
    bijux_rows = observation.bijux_rows
    if not reference_rows or not bijux_rows:
        return False, None
    if "model" in reference_rows[0] or "model" in bijux_rows[0]:
        return False, None

    reference_surface = _parameter_surface_map(
        rows=reference_rows,
        model_name=observation.model_name,
    )
    bijux_surface = _parameter_surface_map(
        rows=bijux_rows,
        model_name=observation.model_name,
    )
    if reference_surface is None or bijux_surface is None:
        return False, None
    if set(reference_surface) != set(bijux_surface):
        return True, False
    for key, reference_value in reference_surface.items():
        bijux_value = bijux_surface[key]
        if isinstance(reference_value, tuple) and isinstance(bijux_value, tuple):
            if reference_value[:2] != bijux_value[:2]:
                return True, False
            if not math.isclose(
                reference_value[2],
                bijux_value[2],
                rel_tol=observation.tolerance,
                abs_tol=observation.tolerance,
            ):
                return True, False
            continue
        if not math.isclose(
            float(reference_value),
            float(bijux_value),
            rel_tol=observation.tolerance,
            abs_tol=observation.tolerance,
        ):
            return True, False
    return True, True


def _parameter_surface_map(
    *,
    rows: list[dict[str, object]],
    model_name: str,
) -> dict[str, float | tuple[bool, int, float]] | None:
    if not rows:
        return None
    first_row = rows[0]
    if "parameter" in first_row:
        surface: dict[str, float] = {}
        for row in rows:
            parameter_name = _optional_text(row.get("parameter"))
            value = _optional_float(row.get("value"))
            if parameter_name is None or value is None:
                continue
            if parameter_name in _OBJECTIVE_ONLY_PARAMETERS:
                continue
            canonical_name, canonical_value = _canonical_primary_parameter(
                model_name,
                parameter_name,
                value,
            )
            if canonical_name is None or canonical_value is None:
                continue
            surface[canonical_name] = canonical_value
        return surface
    if "source_state" in first_row and "target_state" in first_row:
        surface = {}
        for row in rows:
            source_state = _optional_text(row.get("source_state"))
            target_state = _optional_text(row.get("target_state"))
            rate = _optional_float(row.get("rate"))
            if source_state is None or target_state is None or rate is None:
                continue
            identifier = f"{source_state}->{target_state}"
            surface[identifier] = (
                bool(row.get("transition_allowed", False)),
                int(row.get("step_distance", 0)),
                rate,
            )
        return surface
    return None


def _boundary_detected(summary: dict[str, object]) -> bool:
    if summary.get("hit_lower_parameter_boundary") is True:
        return True
    if summary.get("hit_upper_parameter_boundary") is True:
        return True
    optimizer_result = summary.get("optimizer_result")
    if not isinstance(optimizer_result, dict):
        return False
    return bool(
        optimizer_result.get("hit_lower_parameter_bound")
        or optimizer_result.get("hit_upper_parameter_bound")
    )


def _optimizer_trace(summary: dict[str, object]) -> list[dict[str, object]] | None:
    optimizer_result = summary.get("optimizer_result")
    if not isinstance(optimizer_result, dict):
        return None
    trace = optimizer_result.get("attempt_rows") or optimizer_result.get("profile_rows")
    if not isinstance(trace, list):
        return None
    return [dict(row) for row in trace if isinstance(row, dict)]


def _local_optimum_count(
    trace: list[dict[str, object]] | None,
    *,
    model_name: str,
    parameter_name: str | None,
) -> int | None:
    if not trace:
        return None
    canonical_parameter_name = _canonical_primary_parameter(
        model_name,
        parameter_name,
        0.0 if parameter_name is not None else None,
    )[0]
    parameter_key = _trace_parameter_key(trace, canonical_parameter_name)
    if parameter_key is None:
        return None
    ordered_points = []
    for row in trace:
        parameter_value = _optional_float(row.get(parameter_key))
        log_likelihood = _optional_float(row.get("log_likelihood"))
        if parameter_value is None or log_likelihood is None:
            continue
        if canonical_parameter_name == "rate_change" and parameter_key == "a":
            parameter_value = -parameter_value
        ordered_points.append((parameter_value, log_likelihood))
    if not ordered_points:
        return None
    profile = sorted(
        {
            round(parameter, 12): (parameter, log_likelihood)
            for parameter, log_likelihood in ordered_points
        }.values(),
        key=lambda item: item[0],
    )
    if len(profile) == 1:
        return 1
    tolerance = 1e-9
    local_optimum_count = 0
    for index, (_, log_likelihood) in enumerate(profile):
        left = None if index == 0 else profile[index - 1][1]
        right = None if index == len(profile) - 1 else profile[index + 1][1]
        if left is None:
            if log_likelihood >= right - tolerance:
                local_optimum_count += 1
            continue
        if right is None:
            if log_likelihood >= left - tolerance:
                local_optimum_count += 1
            continue
        if log_likelihood >= left - tolerance and log_likelihood >= right - tolerance:
            local_optimum_count += 1
    return local_optimum_count


def _trace_parameter_key(
    trace: list[dict[str, object]],
    parameter_name: str | None,
) -> str | None:
    preferred_keys = []
    if parameter_name is not None:
        preferred_keys.append(parameter_name)
    preferred_keys.extend(
        [
            "parameter_value",
            "lambda",
            "kappa",
            "delta",
            "alpha",
            "a",
            "rate_change",
        ]
    )
    for key in preferred_keys:
        if any(_optional_float(row.get(key)) is not None for row in trace):
            return key
    return None


def _mismatch_type(
    *,
    same_likelihood: bool | None,
    parameter_surface_comparable: bool,
    same_parameter_surface: bool | None,
    boundary_solution_detected: bool,
) -> str:
    if same_likelihood is None:
        return "comparison_not_available"
    if not parameter_surface_comparable or same_parameter_surface is None:
        if boundary_solution_detected:
            return "boundary_solution_review"
        return "parameter_surface_not_applicable"
    if same_likelihood and same_parameter_surface:
        return "no_algorithm_mismatch"
    if same_likelihood and not same_parameter_surface:
        return "same_likelihood_different_parameters"
    if (not same_likelihood) and same_parameter_surface:
        return "different_likelihood_same_parameters"
    if boundary_solution_detected:
        return "boundary_solution_review"
    return "different_likelihood_different_parameters"
