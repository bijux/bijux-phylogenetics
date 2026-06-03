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
class GeigerParameterizationRegistryRow:
    """One governed record of a Bijux-versus-geiger parameterization contract."""

    case_id: str
    function_name: str
    model_name: str
    reference_surface_contract: str
    bijux_surface_contract: str
    status: str
    canonical_parameter_name: str | None
    reference_parameter_name: str | None
    bijux_parameter_name: str | None
    reference_parameter_value: float | None
    bijux_parameter_value: float | None
    converted_reference_parameter_value: float | None
    converted_bijux_parameter_value: float | None
    parameter_match_after_conversion: bool | None
    parameter_conversion_rule: str
    reference_parameter_bounds: dict[str, float] | None
    bijux_parameter_bounds: dict[str, float] | None
    converted_reference_parameter_bounds: dict[str, float] | None
    parameter_bounds_match_after_conversion: bool | None
    bounds_conversion_rule: str
    root_state_parameterization_reference: str
    root_state_parameterization_bijux: str
    root_state_match_within_tolerance: bool | None
    variance_parameterization_reference: str
    variance_parameterization_bijux: str
    variance_match_within_tolerance: bool | None
    likelihood_constants_policy_reference: str
    likelihood_constants_policy_bijux: str
    likelihood_constants_comparison_policy: str
    log_likelihood_match_within_tolerance: bool | None
    expected_divergence: bool
    expected_divergence_kind: str | None
    expected_divergence_evidence: str | None


def build_geiger_parameterization_registry_rows(
    observations: list[GeigerParityObservation],
) -> list[GeigerParameterizationRegistryRow]:
    """Build the governed parameterization registry for geiger live parity."""

    case_by_id = {case.case_id: case for case in list_geiger_parity_cases()}
    return [
        _build_parameterization_registry_row(
            observation, case_by_id[observation.case_id]
        )
        for observation in observations
    ]


def write_geiger_parameterization_registry_table(
    path: Path,
    report: GeigerParityReport,
) -> Path:
    """Write the governed geiger parameterization registry as TSV."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_id",
                "function_name",
                "model_name",
                "reference_surface_contract",
                "bijux_surface_contract",
                "status",
                "canonical_parameter_name",
                "reference_parameter_name",
                "bijux_parameter_name",
                "reference_parameter_value",
                "bijux_parameter_value",
                "converted_reference_parameter_value",
                "converted_bijux_parameter_value",
                "parameter_match_after_conversion",
                "parameter_conversion_rule",
                "reference_parameter_bounds",
                "bijux_parameter_bounds",
                "converted_reference_parameter_bounds",
                "parameter_bounds_match_after_conversion",
                "bounds_conversion_rule",
                "root_state_parameterization_reference",
                "root_state_parameterization_bijux",
                "root_state_match_within_tolerance",
                "variance_parameterization_reference",
                "variance_parameterization_bijux",
                "variance_match_within_tolerance",
                "likelihood_constants_policy_reference",
                "likelihood_constants_policy_bijux",
                "likelihood_constants_comparison_policy",
                "log_likelihood_match_within_tolerance",
                "expected_divergence",
                "expected_divergence_kind",
                "expected_divergence_evidence",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.parameterization_registry_rows:
            payload = asdict(row)
            payload["reference_parameter_bounds"] = json.dumps(
                row.reference_parameter_bounds,
                sort_keys=True,
            )
            payload["bijux_parameter_bounds"] = json.dumps(
                row.bijux_parameter_bounds,
                sort_keys=True,
            )
            payload["converted_reference_parameter_bounds"] = json.dumps(
                row.converted_reference_parameter_bounds,
                sort_keys=True,
            )
            writer.writerow(payload)
    return path


def _build_parameterization_registry_row(
    observation: GeigerParityObservation,
    case: GeigerParityCase,
) -> GeigerParameterizationRegistryRow:
    canonical_parameter_name = _canonical_parameter_name(case, observation)
    canonical_reference_parameter_value = _optional_float(
        None
        if observation.reference_summary is None
        else observation.reference_summary.get("parameter_value")
    )
    canonical_bijux_parameter_value = _optional_float(
        None
        if observation.bijux_summary is None
        else observation.bijux_summary.get("parameter_value")
    )
    reference_parameter_name, reference_parameter_value = _reference_parameter_surface(
        case,
        canonical_parameter_name=canonical_parameter_name,
        canonical_parameter_value=canonical_reference_parameter_value,
    )
    bijux_parameter_name, bijux_parameter_value = _bijux_parameter_surface(
        case,
        canonical_parameter_name=canonical_parameter_name,
        canonical_parameter_value=canonical_bijux_parameter_value,
    )
    bijux_bounds = _bijux_parameter_bounds(case)
    reference_bounds = _reference_parameter_bounds(case)
    converted_reference_bounds = _converted_reference_parameter_bounds(
        case, reference_bounds
    )
    parameter_match_after_conversion = _same_numeric_value(
        canonical_reference_parameter_value,
        canonical_bijux_parameter_value,
        tolerance=observation.tolerance,
    )
    parameter_bounds_match_after_conversion = _same_bounds(
        converted_reference_bounds,
        bijux_bounds,
        tolerance=observation.tolerance,
    )
    root_state_match = _same_numeric_value(
        _optional_float(
            None
            if observation.reference_summary is None
            else observation.reference_summary.get("root_state")
        ),
        _optional_float(
            None
            if observation.bijux_summary is None
            else observation.bijux_summary.get("root_state")
        ),
        tolerance=observation.tolerance,
    )
    variance_match = _same_numeric_value(
        _optional_float(
            None
            if observation.reference_summary is None
            else observation.reference_summary.get("rate")
        ),
        _optional_float(
            None
            if observation.bijux_summary is None
            else observation.bijux_summary.get("rate")
        ),
        tolerance=observation.tolerance,
    )
    log_likelihood_match = _same_numeric_value(
        _optional_float(
            None
            if observation.reference_summary is None
            else observation.reference_summary.get("log_likelihood")
        ),
        _optional_float(
            None
            if observation.bijux_summary is None
            else observation.bijux_summary.get("log_likelihood")
        ),
        tolerance=observation.tolerance,
    )
    expected_divergence_kind = _expected_divergence_kind(case)
    expected_divergence_evidence = _expected_divergence_evidence(
        case,
        reference_parameter_name=reference_parameter_name,
        bijux_parameter_name=bijux_parameter_name,
        reference_bounds=reference_bounds,
        bijux_bounds=bijux_bounds,
        converted_reference_bounds=converted_reference_bounds,
    )
    return GeigerParameterizationRegistryRow(
        case_id=observation.case_id,
        function_name=observation.function_name,
        model_name=observation.model_name,
        reference_surface_contract=observation.function_name,
        bijux_surface_contract=_bijux_surface_contract(case),
        status=observation.status,
        canonical_parameter_name=canonical_parameter_name,
        reference_parameter_name=reference_parameter_name,
        bijux_parameter_name=bijux_parameter_name,
        reference_parameter_value=reference_parameter_value,
        bijux_parameter_value=bijux_parameter_value,
        converted_reference_parameter_value=canonical_reference_parameter_value,
        converted_bijux_parameter_value=canonical_bijux_parameter_value,
        parameter_match_after_conversion=parameter_match_after_conversion,
        parameter_conversion_rule=_parameter_conversion_rule(case),
        reference_parameter_bounds=reference_bounds,
        bijux_parameter_bounds=bijux_bounds,
        converted_reference_parameter_bounds=converted_reference_bounds,
        parameter_bounds_match_after_conversion=parameter_bounds_match_after_conversion,
        bounds_conversion_rule=_bounds_conversion_rule(case),
        root_state_parameterization_reference=_root_state_parameterization_reference(
            case
        ),
        root_state_parameterization_bijux=_root_state_parameterization_bijux(case),
        root_state_match_within_tolerance=root_state_match,
        variance_parameterization_reference=_variance_parameterization_reference(case),
        variance_parameterization_bijux=_variance_parameterization_bijux(case),
        variance_match_within_tolerance=variance_match,
        likelihood_constants_policy_reference=_likelihood_constants_policy_reference(
            case
        ),
        likelihood_constants_policy_bijux=_likelihood_constants_policy_bijux(case),
        likelihood_constants_comparison_policy=_likelihood_constants_comparison_policy(
            case
        ),
        log_likelihood_match_within_tolerance=log_likelihood_match,
        expected_divergence=expected_divergence_kind is not None,
        expected_divergence_kind=expected_divergence_kind,
        expected_divergence_evidence=expected_divergence_evidence,
    )


def _canonical_parameter_name(
    case: GeigerParityCase,
    observation: GeigerParityObservation,
) -> str | None:
    if observation.bijux_summary is not None:
        value = observation.bijux_summary.get("parameter_name")
        if isinstance(value, str) and value:
            return value
    if observation.reference_summary is not None:
        value = observation.reference_summary.get("parameter_name")
        if isinstance(value, str) and value:
            return value
    if case.operation == "fit-continuous":
        if case.model_name == "OU":
            return "alpha"
        if case.model_name == "EB":
            return "rate_change"
        if case.model_name in {"lambda", "kappa", "delta"}:
            return case.model_name
    if case.operation == "fit-discrete-mk" and case.discrete_transform_name == "EB":
        return "a"
    if case.operation == "fit-discrete-mk" and case.discrete_transform_name in {
        "lambda",
        "kappa",
        "delta",
    }:
        return case.discrete_transform_name
    return None


def _reference_parameter_surface(
    case: GeigerParityCase,
    *,
    canonical_parameter_name: str | None,
    canonical_parameter_value: float | None,
) -> tuple[str | None, float | None]:
    if canonical_parameter_name is None or canonical_parameter_value is None:
        return None, None
    if case.operation == "fit-continuous" and case.model_name == "EB":
        return "a", -canonical_parameter_value
    return canonical_parameter_name, canonical_parameter_value


def _bijux_parameter_surface(
    case: GeigerParityCase,
    *,
    canonical_parameter_name: str | None,
    canonical_parameter_value: float | None,
) -> tuple[str | None, float | None]:
    if canonical_parameter_name is None or canonical_parameter_value is None:
        return None, None
    if case.operation == "fit-continuous" and case.model_name == "EB":
        return "rate_change", canonical_parameter_value
    return canonical_parameter_name, canonical_parameter_value


def _bijux_parameter_bounds(
    case: GeigerParityCase,
) -> dict[str, float] | None:
    bounds = _case_bounds(case)
    if bounds is None:
        return None
    return {"lower": bounds[0], "upper": bounds[1]}


def _reference_parameter_bounds(
    case: GeigerParityCase,
) -> dict[str, float] | None:
    bounds = _case_bounds(case)
    if bounds is None:
        return None
    if case.operation == "fit-continuous" and case.model_name == "EB":
        lower_public, upper_public = bounds
        upper_a = -1e-6 if lower_public <= 0.0 else -lower_public
        return {"lower": -upper_public, "upper": upper_a}
    return {"lower": bounds[0], "upper": bounds[1]}


def _converted_reference_parameter_bounds(
    case: GeigerParityCase,
    reference_bounds: dict[str, float] | None,
) -> dict[str, float] | None:
    if reference_bounds is None:
        return None
    if case.operation == "fit-continuous" and case.model_name == "EB":
        return {
            "lower": -reference_bounds["upper"],
            "upper": -reference_bounds["lower"],
        }
    return dict(reference_bounds)


def _case_bounds(case: GeigerParityCase) -> tuple[float, float] | None:
    if case.operation == "fit-continuous":
        if case.model_name == "lambda":
            return case.lambda_bounds
        if case.model_name == "kappa":
            return case.kappa_bounds
        if case.model_name == "delta":
            return case.delta_bounds
        if case.model_name == "OU":
            return case.ou_bounds
        if case.model_name == "EB":
            return case.early_burst_bounds
        return None
    if case.operation != "fit-discrete-mk":
        return None
    if case.discrete_transform_name == "lambda":
        return case.lambda_bounds
    if case.discrete_transform_name == "kappa":
        return case.kappa_bounds
    if case.discrete_transform_name == "delta":
        return case.delta_bounds
    if case.discrete_transform_name == "EB":
        return case.early_burst_bounds
    return None


def _parameter_conversion_rule(case: GeigerParityCase) -> str:
    if case.operation == "fit-continuous" and case.model_name == "EB":
        return (
            "geiger raw early-burst parameter a is sign-flipped into the Bijux public "
            "rate_change surface: rate_change = -a"
        )
    return "no parameter conversion required"


def _bounds_conversion_rule(case: GeigerParityCase) -> str:
    if case.operation == "fit-continuous" and case.model_name == "EB":
        return (
            "geiger raw early-burst bounds are recorded on a and converted into the "
            "Bijux public rate_change surface by sign flip and interval reversal"
        )
    return "no bounds conversion required"


def _root_state_parameterization_reference(case: GeigerParityCase) -> str:
    if case.operation != "fit-continuous":
        return "not-applicable-for-discrete-and-model-comparison-surfaces"
    if case.model_name == "model-comparison":
        return "not-applicable-for-model-comparison-summary"
    return "geiger fitContinuous opt$z0 exposed in parity summaries as root_state"


def _root_state_parameterization_bijux(case: GeigerParityCase) -> str:
    if case.operation != "fit-continuous":
        return "not-applicable-for-discrete-and-model-comparison-surfaces"
    if case.model_name == "model-comparison":
        return "not-applicable-for-model-comparison-summary"
    return "Bijux intercept-only theta exposed in parity summaries as root_state"


def _variance_parameterization_reference(case: GeigerParityCase) -> str:
    if case.operation != "fit-continuous":
        return "not-applicable-for-discrete-and-model-comparison-surfaces"
    if case.model_name == "model-comparison":
        return "not-applicable-for-model-comparison-summary"
    return "geiger fitContinuous opt$sigsq exposed in parity summaries as rate"


def _variance_parameterization_bijux(case: GeigerParityCase) -> str:
    if case.operation != "fit-continuous":
        return "not-applicable-for-discrete-and-model-comparison-surfaces"
    if case.model_name == "model-comparison":
        return "not-applicable-for-model-comparison-summary"
    return "Bijux intercept-only sigma_squared exposed in parity summaries as rate"


def _likelihood_constants_policy_reference(case: GeigerParityCase) -> str:
    if case.operation == "compare-fitcontinuous-models":
        return "model-comparison summaries do not expose one case-level log-likelihood"
    return (
        "raw geiger summary log_likelihood is compared directly without a post-hoc "
        "likelihood-constant conversion in this registry tranche"
    )


def _likelihood_constants_policy_bijux(case: GeigerParityCase) -> str:
    if case.operation == "compare-fitcontinuous-models":
        return "model-comparison summaries do not expose one case-level log-likelihood"
    return (
        "raw Bijux summary log_likelihood is compared directly without a post-hoc "
        "likelihood-constant conversion in this registry tranche"
    )


def _likelihood_constants_comparison_policy(case: GeigerParityCase) -> str:
    if case.operation == "compare-fitcontinuous-models":
        return "not-applicable-for-model-comparison-summary"
    return "direct-summary-log-likelihood-comparison"


def _expected_divergence_kind(case: GeigerParityCase) -> str | None:
    if case.operation == "fit-continuous" and case.model_name == "EB":
        return "continuous-early-burst-sign-and-bound-convention"
    return None


def _expected_divergence_evidence(
    case: GeigerParityCase,
    *,
    reference_parameter_name: str | None,
    bijux_parameter_name: str | None,
    reference_bounds: dict[str, float] | None,
    bijux_bounds: dict[str, float] | None,
    converted_reference_bounds: dict[str, float] | None,
) -> str | None:
    if case.operation == "fit-continuous" and case.model_name == "EB":
        return (
            "The governed live geiger runner exposes the raw early-burst surface on "
            f"{reference_parameter_name}, while Bijux exposes {bijux_parameter_name}; "
            f"the raw geiger bounds {reference_bounds} convert to "
            f"{converted_reference_bounds}, whereas the Bijux public bounds are "
            f"{bijux_bounds}. This keeps the sign flip explicit and also preserves "
            "the near-zero exclusion on the geiger raw a surface."
        )
    return None


def _bijux_surface_contract(case: GeigerParityCase) -> str:
    if case.operation == "fit-continuous":
        return f"{case.python_function_name}(mode='{case.python_mode}')"
    if case.operation == "compare-fitcontinuous-models":
        return f"{case.python_function_name}(modes=governed-fitcontinuous-order)"
    if case.discrete_transform_name is None:
        return f"{case.python_function_name}(model='{case.python_mode}')"
    return (
        f"{case.python_function_name}(model='{case.python_mode}', "
        f"transform='{case.discrete_transform_name}')"
    )


def _optional_float(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        numeric = float(value)
        return None if math.isnan(numeric) else numeric
    return None


def _same_numeric_value(
    left: float | None,
    right: float | None,
    *,
    tolerance: float,
) -> bool | None:
    if left is None or right is None:
        return None
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def _same_bounds(
    left: dict[str, float] | None,
    right: dict[str, float] | None,
    *,
    tolerance: float,
) -> bool | None:
    if left is None or right is None:
        return None if left is None and right is None else False
    bounds_tolerance = min(tolerance, 1e-9)
    return math.isclose(
        left["lower"],
        right["lower"],
        rel_tol=bounds_tolerance,
        abs_tol=bounds_tolerance,
    ) and math.isclose(
        left["upper"],
        right["upper"],
        rel_tol=bounds_tolerance,
        abs_tol=bounds_tolerance,
    )
