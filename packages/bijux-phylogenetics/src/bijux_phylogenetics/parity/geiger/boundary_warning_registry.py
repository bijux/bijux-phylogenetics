from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import TYPE_CHECKING

from .registry import GeigerParityCase, list_geiger_parity_cases

if TYPE_CHECKING:
    from .runner import GeigerParityObservation, GeigerParityReport


@dataclass(frozen=True, slots=True)
class GeigerBoundaryWarningRow:
    """One governed record of a Bijux-versus-geiger boundary-warning interpretation."""

    case_id: str
    function_name: str
    model_name: str
    status: str
    scope: str
    affected_parameter: str | None
    lower_bound: float | None
    upper_bound: float | None
    reference_parameter_value: float | None
    bijux_parameter_value: float | None
    reference_hit_lower_boundary: bool | None
    reference_hit_upper_boundary: bool | None
    bijux_hit_lower_boundary: bool | None
    bijux_hit_upper_boundary: bool | None
    reference_near_lower_boundary: bool | None
    reference_near_upper_boundary: bool | None
    bijux_near_lower_boundary: bool | None
    bijux_near_upper_boundary: bool | None
    reference_flat_likelihood_near_boundary: bool | None
    bijux_flat_likelihood_near_boundary: bool | None
    reference_boundary_warning_kinds: list[str]
    bijux_boundary_warning_kinds: list[str]
    reference_boundary_dominates_interpretation: bool | None
    bijux_boundary_dominates_interpretation: bool | None
    reference_stable_conclusion_supported: bool | None
    bijux_stable_conclusion_supported: bool | None
    affected_parameter_match: bool | None
    stable_conclusion_supported_match: bool | None
    boundary_evidence: str


@dataclass(frozen=True, slots=True)
class _BoundaryAssessment:
    affected_parameter: str
    lower_bound: float
    upper_bound: float
    parameter_value: float
    hit_lower_boundary: bool
    hit_upper_boundary: bool
    near_lower_boundary: bool
    near_upper_boundary: bool
    flat_likelihood_near_boundary: bool
    boundary_warning_kinds: list[str]
    boundary_dominates_interpretation: bool
    stable_conclusion_supported: bool


def build_geiger_boundary_warning_rows(
    observations: list[GeigerParityObservation],
) -> list[GeigerBoundaryWarningRow]:
    """Build the governed geiger boundary-warning registry rows."""

    case_by_id = {case.case_id: case for case in list_geiger_parity_cases()}
    return [
        _build_boundary_warning_row(observation, case_by_id[observation.case_id])
        for observation in observations
    ]


def write_geiger_boundary_warning_table(
    path: Path,
    report: GeigerParityReport,
) -> Path:
    """Write the governed geiger boundary-warning registry as TSV."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_id",
                "function_name",
                "model_name",
                "status",
                "scope",
                "affected_parameter",
                "lower_bound",
                "upper_bound",
                "reference_parameter_value",
                "bijux_parameter_value",
                "reference_hit_lower_boundary",
                "reference_hit_upper_boundary",
                "bijux_hit_lower_boundary",
                "bijux_hit_upper_boundary",
                "reference_near_lower_boundary",
                "reference_near_upper_boundary",
                "bijux_near_lower_boundary",
                "bijux_near_upper_boundary",
                "reference_flat_likelihood_near_boundary",
                "bijux_flat_likelihood_near_boundary",
                "reference_boundary_warning_kinds",
                "bijux_boundary_warning_kinds",
                "reference_boundary_dominates_interpretation",
                "bijux_boundary_dominates_interpretation",
                "reference_stable_conclusion_supported",
                "bijux_stable_conclusion_supported",
                "affected_parameter_match",
                "stable_conclusion_supported_match",
                "boundary_evidence",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.boundary_warning_rows:
            payload = asdict(row)
            payload["reference_boundary_warning_kinds"] = json.dumps(
                row.reference_boundary_warning_kinds
            )
            payload["bijux_boundary_warning_kinds"] = json.dumps(
                row.bijux_boundary_warning_kinds
            )
            writer.writerow(payload)
    return path


def _build_boundary_warning_row(
    observation: GeigerParityObservation,
    case: GeigerParityCase,
) -> GeigerBoundaryWarningRow:
    if case.operation != "fit-continuous" or case.model_name in {"BM", "white"}:
        return GeigerBoundaryWarningRow(
            case_id=observation.case_id,
            function_name=observation.function_name,
            model_name=observation.model_name,
            status=observation.status,
            scope="outside-current-fitcontinuous-boundary-warning-scope",
            affected_parameter=None,
            lower_bound=None,
            upper_bound=None,
            reference_parameter_value=None,
            bijux_parameter_value=None,
            reference_hit_lower_boundary=None,
            reference_hit_upper_boundary=None,
            bijux_hit_lower_boundary=None,
            bijux_hit_upper_boundary=None,
            reference_near_lower_boundary=None,
            reference_near_upper_boundary=None,
            bijux_near_lower_boundary=None,
            bijux_near_upper_boundary=None,
            reference_flat_likelihood_near_boundary=None,
            bijux_flat_likelihood_near_boundary=None,
            reference_boundary_warning_kinds=[],
            bijux_boundary_warning_kinds=[],
            reference_boundary_dominates_interpretation=None,
            bijux_boundary_dominates_interpretation=None,
            reference_stable_conclusion_supported=None,
            bijux_stable_conclusion_supported=None,
            affected_parameter_match=None,
            stable_conclusion_supported_match=None,
            boundary_evidence=(
                "This governed boundary-warning registry tranche covers only "
                "parameterized geiger::fitContinuous surfaces where lower, upper, "
                "near-boundary, and flat-likelihood review are meaningful."
            ),
        )
    reference_assessment = _summary_boundary_assessment(
        case, observation.reference_summary
    )
    bijux_assessment = _summary_boundary_assessment(case, observation.bijux_summary)
    return GeigerBoundaryWarningRow(
        case_id=observation.case_id,
        function_name=observation.function_name,
        model_name=observation.model_name,
        status=observation.status,
        scope="fitcontinuous-parameter-boundary-review",
        affected_parameter=(
            None if bijux_assessment is None else bijux_assessment.affected_parameter
        ),
        lower_bound=(
            None if bijux_assessment is None else bijux_assessment.lower_bound
        ),
        upper_bound=(
            None if bijux_assessment is None else bijux_assessment.upper_bound
        ),
        reference_parameter_value=(
            None
            if reference_assessment is None
            else reference_assessment.parameter_value
        ),
        bijux_parameter_value=(
            None if bijux_assessment is None else bijux_assessment.parameter_value
        ),
        reference_hit_lower_boundary=(
            None
            if reference_assessment is None
            else reference_assessment.hit_lower_boundary
        ),
        reference_hit_upper_boundary=(
            None
            if reference_assessment is None
            else reference_assessment.hit_upper_boundary
        ),
        bijux_hit_lower_boundary=(
            None if bijux_assessment is None else bijux_assessment.hit_lower_boundary
        ),
        bijux_hit_upper_boundary=(
            None if bijux_assessment is None else bijux_assessment.hit_upper_boundary
        ),
        reference_near_lower_boundary=(
            None
            if reference_assessment is None
            else reference_assessment.near_lower_boundary
        ),
        reference_near_upper_boundary=(
            None
            if reference_assessment is None
            else reference_assessment.near_upper_boundary
        ),
        bijux_near_lower_boundary=(
            None if bijux_assessment is None else bijux_assessment.near_lower_boundary
        ),
        bijux_near_upper_boundary=(
            None if bijux_assessment is None else bijux_assessment.near_upper_boundary
        ),
        reference_flat_likelihood_near_boundary=(
            None
            if reference_assessment is None
            else reference_assessment.flat_likelihood_near_boundary
        ),
        bijux_flat_likelihood_near_boundary=(
            None
            if bijux_assessment is None
            else bijux_assessment.flat_likelihood_near_boundary
        ),
        reference_boundary_warning_kinds=(
            []
            if reference_assessment is None
            else reference_assessment.boundary_warning_kinds
        ),
        bijux_boundary_warning_kinds=(
            [] if bijux_assessment is None else bijux_assessment.boundary_warning_kinds
        ),
        reference_boundary_dominates_interpretation=(
            None
            if reference_assessment is None
            else reference_assessment.boundary_dominates_interpretation
        ),
        bijux_boundary_dominates_interpretation=(
            None
            if bijux_assessment is None
            else bijux_assessment.boundary_dominates_interpretation
        ),
        reference_stable_conclusion_supported=(
            None
            if reference_assessment is None
            else reference_assessment.stable_conclusion_supported
        ),
        bijux_stable_conclusion_supported=(
            None
            if bijux_assessment is None
            else bijux_assessment.stable_conclusion_supported
        ),
        affected_parameter_match=_optional_str(
            None
            if reference_assessment is None
            else reference_assessment.affected_parameter
        )
        == _optional_str(
            None if bijux_assessment is None else bijux_assessment.affected_parameter
        ),
        stable_conclusion_supported_match=(
            None
            if reference_assessment is None or bijux_assessment is None
            else reference_assessment.stable_conclusion_supported
            == bijux_assessment.stable_conclusion_supported
        ),
        boundary_evidence=_boundary_evidence(
            case,
            reference_assessment=reference_assessment,
            bijux_assessment=bijux_assessment,
        ),
    )


def _summary_boundary_assessment(
    case: GeigerParityCase,
    summary: dict[str, object] | None,
) -> _BoundaryAssessment | None:
    if summary is None:
        return None
    parameter_name = _optional_str(summary.get("parameter_name"))
    parameter_value = _optional_float(summary.get("parameter_value"))
    bounds = _parameter_bounds(case)
    if parameter_name is None or parameter_value is None or bounds is None:
        return None
    lower_bound, upper_bound = bounds
    near_boundary_tolerance = max((upper_bound - lower_bound) / 20.0, 1e-6)
    hit_lower_boundary = bool(summary.get("hit_lower_parameter_boundary"))
    hit_upper_boundary = bool(summary.get("hit_upper_parameter_boundary"))
    near_lower_boundary = hit_lower_boundary or (
        parameter_value <= lower_bound + near_boundary_tolerance
    )
    near_upper_boundary = hit_upper_boundary or (
        parameter_value >= upper_bound - near_boundary_tolerance
    )
    boundary_warning_kinds = _boundary_warning_kinds(
        case,
        summary=summary,
        near_lower_boundary=near_lower_boundary,
        near_upper_boundary=near_upper_boundary,
    )
    flat_likelihood_near_boundary = any(
        kind in {"flat_likelihood", "flat_likelihood_profile"}
        for kind in boundary_warning_kinds
    )
    boundary_dominates_interpretation = (
        near_lower_boundary or near_upper_boundary
    ) and (
        flat_likelihood_near_boundary
        or any(
            kind.startswith("boundary_") or kind in _boundary_limit_warning_kinds(case)
            for kind in boundary_warning_kinds
        )
    )
    return _BoundaryAssessment(
        affected_parameter=parameter_name,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        parameter_value=parameter_value,
        hit_lower_boundary=hit_lower_boundary,
        hit_upper_boundary=hit_upper_boundary,
        near_lower_boundary=near_lower_boundary,
        near_upper_boundary=near_upper_boundary,
        flat_likelihood_near_boundary=flat_likelihood_near_boundary,
        boundary_warning_kinds=boundary_warning_kinds,
        boundary_dominates_interpretation=boundary_dominates_interpretation,
        stable_conclusion_supported=not boundary_dominates_interpretation,
    )


def _boundary_warning_kinds(
    case: GeigerParityCase,
    *,
    summary: dict[str, object],
    near_lower_boundary: bool,
    near_upper_boundary: bool,
) -> list[str]:
    explicit = _string_list(summary.get("identifiability_warning_kinds"))
    if explicit:
        return [kind for kind in explicit if _is_boundary_warning_kind(case, kind)]
    kinds: list[str] = []
    parameter_name = _optional_str(summary.get("parameter_name"))
    if parameter_name is None:
        return kinds
    if bool(summary.get("hit_lower_parameter_boundary")) or bool(
        summary.get("hit_upper_parameter_boundary")
    ):
        kinds.append(f"boundary_{parameter_name}")
    if _flat_likelihood_from_summary(summary):
        kinds.append(
            "flat_likelihood_profile" if case.model_name == "EB" else "flat_likelihood"
        )
    if case.model_name == "OU" and near_lower_boundary:
        kinds.append("weak_pull_to_optimum")
    if case.model_name == "EB" and near_lower_boundary:
        kinds.append("brownian_like_rate_change")
    if case.model_name == "lambda" and near_lower_boundary:
        kinds.append("weak_phylogenetic_signal")
    if case.model_name == "lambda" and near_upper_boundary:
        kinds.append("brownian_limit")
    if case.model_name == "kappa" and near_lower_boundary:
        kinds.append("punctuational_limit")
    if case.model_name == "kappa" and near_upper_boundary:
        kinds.append("upper_search_limit")
    if case.model_name == "delta" and near_lower_boundary:
        kinds.append("early_change_limit")
    if case.model_name == "delta" and near_upper_boundary:
        kinds.append("late_change_limit")
    return kinds


def _flat_likelihood_from_summary(summary: dict[str, object]) -> bool:
    optimizer_result = summary.get("optimizer_result")
    if not isinstance(optimizer_result, dict):
        return False
    profile_rows = optimizer_result.get("profile_rows")
    if isinstance(profile_rows, list):
        log_likelihoods = sorted(
            (
                value
                for row in profile_rows
                if isinstance(row, dict)
                for value in [_optional_float(row.get("log_likelihood"))]
                if value is not None
            ),
            reverse=True,
        )
        return (
            len(log_likelihoods) > 1 and (log_likelihoods[0] - log_likelihoods[1]) < 0.5
        )
    attempt_rows = optimizer_result.get("attempt_rows")
    if not isinstance(attempt_rows, list):
        return False
    log_likelihoods = sorted(
        (
            value
            for row in attempt_rows
            if isinstance(row, dict)
            for value in [_attempt_row_log_likelihood(row)]
            if value is not None and math.isfinite(value)
        ),
        reverse=True,
    )
    return len(log_likelihoods) > 1 and (log_likelihoods[0] - log_likelihoods[1]) < 0.5


def _attempt_row_log_likelihood(row: dict[str, object]) -> float | None:
    named = _optional_float(row.get("log_likelihood"))
    if named is not None:
        return named
    numeric_columns = sorted(
        (key for key in row if key.isdigit()),
        key=lambda value: int(value),
    )
    if len(numeric_columns) < 3:
        return None
    return _optional_float(row.get(numeric_columns[-2]))


def _parameter_bounds(case: GeigerParityCase) -> tuple[float, float] | None:
    if case.model_name == "lambda":
        return (0.0, 1.0) if case.lambda_bounds is None else case.lambda_bounds
    if case.model_name == "kappa":
        return (0.0, 1.0) if case.kappa_bounds is None else case.kappa_bounds
    if case.model_name == "delta":
        return (0.0, 3.0) if case.delta_bounds is None else case.delta_bounds
    if case.model_name == "OU":
        return (0.0, 10.0) if case.ou_bounds is None else case.ou_bounds
    if case.model_name == "EB":
        return (
            (0.0, 50.0) if case.early_burst_bounds is None else case.early_burst_bounds
        )
    return None


def _is_boundary_warning_kind(case: GeigerParityCase, kind: str) -> bool:
    return kind.startswith("boundary_") or kind in {
        "flat_likelihood",
        "flat_likelihood_profile",
        *_boundary_limit_warning_kinds(case),
    }


def _boundary_limit_warning_kinds(case: GeigerParityCase) -> set[str]:
    if case.model_name == "OU":
        return {"weak_pull_to_optimum"}
    if case.model_name == "EB":
        return {"brownian_like_rate_change"}
    if case.model_name == "lambda":
        return {"weak_phylogenetic_signal", "brownian_limit"}
    if case.model_name == "kappa":
        return {"punctuational_limit", "upper_search_limit"}
    if case.model_name == "delta":
        return {"early_change_limit", "late_change_limit"}
    return set()


def _boundary_evidence(
    case: GeigerParityCase,
    *,
    reference_assessment: _BoundaryAssessment | None,
    bijux_assessment: _BoundaryAssessment | None,
) -> str:
    if reference_assessment is None or bijux_assessment is None:
        return (
            "Boundary-warning evidence is unavailable because one side did not retain "
            "a parameterized fitContinuous summary for this governed case."
        )
    return (
        f"{case.model_name} boundary review compares lower={reference_assessment.lower_bound:g}, "
        f"upper={reference_assessment.upper_bound:g}, reference warning kinds "
        f"{reference_assessment.boundary_warning_kinds}, and bijux warning kinds "
        f"{bijux_assessment.boundary_warning_kinds}."
    )


def _optional_float(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        text = value.strip()
        if not text or text.upper() == "NA":
            return None
        try:
            return float(text)
        except ValueError:
            return None
    return None


def _optional_str(value: object) -> str | None:
    return value if isinstance(value, str) and value else None


def _string_list(value: object) -> list[str]:
    if not isinstance(value, list):
        return []
    return [item for item in value if isinstance(item, str) and item]
