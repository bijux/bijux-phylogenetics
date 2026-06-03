from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
import math
from pathlib import Path
from typing import TYPE_CHECKING

from bijux_phylogenetics.comparative.evolutionary_modes import (
    FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD,
    FITCONTINUOUS_MODEL_CONFIDENCE_WEIGHT_BASIS,
)

from .registry import GeigerParityCase, list_geiger_parity_cases

if TYPE_CHECKING:
    from .runner import GeigerParityObservation, GeigerParityReport


@dataclass(frozen=True, slots=True)
class GeigerModelConfidenceRow:
    """One governed fitContinuous model-confidence comparison row."""

    case_id: str
    function_name: str
    model_name: str
    status: str
    candidate_model: str
    weight_basis: str
    delta_threshold: float
    reference_best_model: str | None
    bijux_best_model: str | None
    reference_rank: int | None
    bijux_rank: int | None
    reference_comparable: bool
    bijux_comparable: bool
    reference_delta_aic: float | None
    bijux_delta_aic: float | None
    reference_delta_aicc: float | None
    bijux_delta_aicc: float | None
    reference_akaike_weight: float | None
    bijux_akaike_weight: float | None
    akaike_weight_match_within_tolerance: bool | None
    reference_within_delta_aic_threshold: bool | None
    bijux_within_delta_aic_threshold: bool | None
    within_delta_aic_threshold_match: bool | None
    reference_within_delta_aicc_threshold: bool | None
    bijux_within_delta_aicc_threshold: bool | None
    within_delta_aicc_threshold_match: bool | None
    reference_selected: bool
    bijux_selected: bool
    selected_match: bool
    reference_uncertainty_class: str
    bijux_uncertainty_class: str
    uncertainty_class_match: bool
    reference_uncertainty_language: str
    bijux_uncertainty_language: str
    uncertainty_language_match: bool
    confidence_evidence: str


def build_geiger_model_confidence_rows(
    observations: list[GeigerParityObservation],
) -> list[GeigerModelConfidenceRow]:
    """Build the governed fitContinuous model-confidence registry rows."""

    case_by_id = {case.case_id: case for case in list_geiger_parity_cases()}
    built_rows: list[GeigerModelConfidenceRow] = []
    for observation in observations:
        case = case_by_id[observation.case_id]
        if case.operation != "compare-fitcontinuous-models":
            continue
        built_rows.extend(_build_model_confidence_rows(observation, case))
    return built_rows


def write_geiger_model_confidence_table(
    path: Path,
    report: GeigerParityReport,
) -> Path:
    """Write the governed fitContinuous model-confidence table as TSV."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_id",
                "function_name",
                "model_name",
                "status",
                "candidate_model",
                "weight_basis",
                "delta_threshold",
                "reference_best_model",
                "bijux_best_model",
                "reference_rank",
                "bijux_rank",
                "reference_comparable",
                "bijux_comparable",
                "reference_delta_aic",
                "bijux_delta_aic",
                "reference_delta_aicc",
                "bijux_delta_aicc",
                "reference_akaike_weight",
                "bijux_akaike_weight",
                "akaike_weight_match_within_tolerance",
                "reference_within_delta_aic_threshold",
                "bijux_within_delta_aic_threshold",
                "within_delta_aic_threshold_match",
                "reference_within_delta_aicc_threshold",
                "bijux_within_delta_aicc_threshold",
                "within_delta_aicc_threshold_match",
                "reference_selected",
                "bijux_selected",
                "selected_match",
                "reference_uncertainty_class",
                "bijux_uncertainty_class",
                "uncertainty_class_match",
                "reference_uncertainty_language",
                "bijux_uncertainty_language",
                "uncertainty_language_match",
                "confidence_evidence",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.model_confidence_rows:
            writer.writerow(asdict(row))
    return path


def _build_model_confidence_rows(
    observation: GeigerParityObservation,
    case: GeigerParityCase,
) -> list[GeigerModelConfidenceRow]:
    reference_rows = _rows_by_model(observation.reference_rows)
    bijux_rows = _rows_by_model(observation.bijux_rows)
    model_order = list(dict.fromkeys([*reference_rows.keys(), *bijux_rows.keys()]))
    reference_confidence = _confidence_surface(reference_rows)
    bijux_confidence = _confidence_surface(bijux_rows)
    built_rows: list[GeigerModelConfidenceRow] = []
    for candidate_model in model_order:
        reference_row = reference_rows.get(candidate_model)
        bijux_row = bijux_rows.get(candidate_model)
        reference_weight = reference_confidence.weights.get(candidate_model)
        bijux_weight = bijux_confidence.weights.get(candidate_model)
        reference_within_aic = reference_confidence.within_delta_aic.get(
            candidate_model
        )
        bijux_within_aic = bijux_confidence.within_delta_aic.get(candidate_model)
        reference_within_aicc = reference_confidence.within_delta_aicc.get(
            candidate_model
        )
        bijux_within_aicc = bijux_confidence.within_delta_aicc.get(candidate_model)
        built_rows.append(
            GeigerModelConfidenceRow(
                case_id=observation.case_id,
                function_name=observation.function_name,
                model_name=observation.model_name,
                status=observation.status,
                candidate_model=candidate_model,
                weight_basis=FITCONTINUOUS_MODEL_CONFIDENCE_WEIGHT_BASIS,
                delta_threshold=FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD,
                reference_best_model=reference_confidence.best_model,
                bijux_best_model=bijux_confidence.best_model,
                reference_rank=_optional_int(
                    None if reference_row is None else reference_row.get("rank")
                ),
                bijux_rank=_optional_int(
                    None if bijux_row is None else bijux_row.get("rank")
                ),
                reference_comparable=_row_comparable(reference_row),
                bijux_comparable=_row_comparable(bijux_row),
                reference_delta_aic=_optional_float(
                    None if reference_row is None else reference_row.get("delta_aic")
                ),
                bijux_delta_aic=_optional_float(
                    None if bijux_row is None else bijux_row.get("delta_aic")
                ),
                reference_delta_aicc=_optional_float(
                    None if reference_row is None else reference_row.get("delta_aicc")
                ),
                bijux_delta_aicc=_optional_float(
                    None if bijux_row is None else bijux_row.get("delta_aicc")
                ),
                reference_akaike_weight=reference_weight,
                bijux_akaike_weight=bijux_weight,
                akaike_weight_match_within_tolerance=_match_optional_numeric(
                    reference_weight,
                    bijux_weight,
                    tolerance=observation.tolerance,
                ),
                reference_within_delta_aic_threshold=reference_within_aic,
                bijux_within_delta_aic_threshold=bijux_within_aic,
                within_delta_aic_threshold_match=_match_optional_bool(
                    reference_within_aic,
                    bijux_within_aic,
                ),
                reference_within_delta_aicc_threshold=reference_within_aicc,
                bijux_within_delta_aicc_threshold=bijux_within_aicc,
                within_delta_aicc_threshold_match=_match_optional_bool(
                    reference_within_aicc,
                    bijux_within_aicc,
                ),
                reference_selected=_row_selected(reference_row),
                bijux_selected=_row_selected(bijux_row),
                selected_match=_row_selected(reference_row) == _row_selected(bijux_row),
                reference_uncertainty_class=reference_confidence.uncertainty_class,
                bijux_uncertainty_class=bijux_confidence.uncertainty_class,
                uncertainty_class_match=(
                    reference_confidence.uncertainty_class
                    == bijux_confidence.uncertainty_class
                ),
                reference_uncertainty_language=reference_confidence.uncertainty_language,
                bijux_uncertainty_language=bijux_confidence.uncertainty_language,
                uncertainty_language_match=(
                    reference_confidence.uncertainty_class
                    == bijux_confidence.uncertainty_class
                ),
                confidence_evidence=(
                    "Akaike weights are derived from comparable delta AICc rows only; "
                    "non-comparable candidates retain blank weights while threshold "
                    "review remains governed by delta AIC and delta AICc against the "
                    f"shared {FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD:.1f}-unit policy."
                ),
            )
        )
    return built_rows


@dataclass(frozen=True, slots=True)
class _ConfidenceSurface:
    best_model: str | None
    weights: dict[str, float | None]
    within_delta_aic: dict[str, bool | None]
    within_delta_aicc: dict[str, bool | None]
    uncertainty_class: str
    uncertainty_language: str


def _confidence_surface(
    rows_by_model: dict[str, dict[str, object]],
) -> _ConfidenceSurface:
    comparable_rows = [
        row
        for row in rows_by_model.values()
        if _row_comparable(row)
        and _optional_float(row.get("delta_aicc")) is not None
        and _optional_float(row.get("delta_aic")) is not None
    ]
    if not comparable_rows:
        return _ConfidenceSurface(
            best_model=None,
            weights=dict.fromkeys(rows_by_model),
            within_delta_aic=dict.fromkeys(rows_by_model),
            within_delta_aicc=dict.fromkeys(rows_by_model),
            uncertainty_class="unresolved",
            uncertainty_language=(
                "model confidence is unresolved because no comparable candidate retained "
                "one finite AICc surface for Akaike-weight review"
            ),
        )
    raw_weights = {}
    for row in comparable_rows:
        delta_aicc = _optional_float(row.get("delta_aicc"))
        if delta_aicc is None:
            continue
        raw_weights[row["model"]] = math.exp(-0.5 * delta_aicc)
    weight_total = sum(raw_weights.values())
    weights = {
        model: (raw_weight / weight_total if weight_total else 0.0)
        for model, raw_weight in raw_weights.items()
    }
    weight_map = dict.fromkeys(rows_by_model)
    weight_map.update(weights)
    within_delta_aic = dict.fromkeys(rows_by_model)
    within_delta_aicc = dict.fromkeys(rows_by_model)
    for row in comparable_rows:
        model = _row_model(row)
        delta_aic = _optional_float(row.get("delta_aic"))
        delta_aicc = _optional_float(row.get("delta_aicc"))
        within_delta_aic[model] = (
            None
            if delta_aic is None
            else delta_aic <= FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD
        )
        within_delta_aicc[model] = (
            None
            if delta_aicc is None
            else delta_aicc <= FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD
        )
    best_row = min(
        comparable_rows,
        key=lambda row: (
            _optional_float(row.get("delta_aicc")),
            _optional_float(row.get("delta_aic")),
            _row_model(row),
        ),
    )
    best_model = _row_model(best_row)
    best_weight = weight_map.get(best_model)
    nearby_models = [
        model
        for model, within_threshold in within_delta_aicc.items()
        if model != best_model and within_threshold is True
    ]
    if nearby_models:
        uncertainty_class = "limited"
        tied_models = ", ".join([best_model, *nearby_models])
        uncertainty_language = (
            "model confidence is limited because "
            f"{tied_models} remain within "
            f"{FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD:.1f} AICc units of the "
            f"selected surface; {best_model} carries Akaike weight {best_weight:.3f}"
        )
    elif best_weight is not None and best_weight < 0.9:
        uncertainty_class = "moderate"
        uncertainty_language = (
            f"model confidence is moderate because {best_model} carries Akaike "
            f"weight {best_weight:.3f} even though no runner-up remains within "
            f"{FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD:.1f} AICc units"
        )
    else:
        uncertainty_class = "strong"
        uncertainty_language = (
            f"model confidence is strong because {best_model} carries Akaike "
            f"weight {0.0 if best_weight is None else best_weight:.3f} and no "
            "runner-up remains within "
            f"{FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD:.1f} AICc units"
        )
    return _ConfidenceSurface(
        best_model=best_model,
        weights=weight_map,
        within_delta_aic=within_delta_aic,
        within_delta_aicc=within_delta_aicc,
        uncertainty_class=uncertainty_class,
        uncertainty_language=uncertainty_language,
    )


def _rows_by_model(
    rows: list[dict[str, object]] | None,
) -> dict[str, dict[str, object]]:
    if rows is None:
        return {}
    return {
        model: row
        for row in rows
        if isinstance((model := row.get("model")), str) and model
    }


def _row_model(row: dict[str, object]) -> str:
    model = row.get("model")
    if not isinstance(model, str) or not model:
        raise ValueError("model-confidence rows require one non-empty model name")
    return model


def _row_comparable(row: dict[str, object] | None) -> bool:
    return bool(row is not None and row.get("comparable") is True)


def _row_selected(row: dict[str, object] | None) -> bool:
    return bool(row is not None and row.get("selected") is True)


def _match_optional_numeric(
    reference_value: float | None,
    bijux_value: float | None,
    *,
    tolerance: float,
) -> bool | None:
    if reference_value is None or bijux_value is None:
        return None if reference_value is None and bijux_value is None else False
    return math.isclose(
        reference_value,
        bijux_value,
        rel_tol=tolerance,
        abs_tol=tolerance,
    )


def _match_optional_bool(
    reference_value: bool | None,
    bijux_value: bool | None,
) -> bool | None:
    if reference_value is None or bijux_value is None:
        return None if reference_value is None and bijux_value is None else False
    return reference_value == bijux_value


def _optional_float(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _optional_int(value: object) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return None
