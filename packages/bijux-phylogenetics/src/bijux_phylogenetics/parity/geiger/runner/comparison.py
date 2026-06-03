from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from bijux_phylogenetics.parity.geiger.registry import GeigerParityCase


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_rows_table(path: Path) -> list[dict[str, object]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows: list[dict[str, object]] = []
        for row in reader:
            parsed: dict[str, object] = {}
            for key, value in row.items():
                if value is None or value == "":
                    parsed[key] = ""
                    continue
                if key in {
                    "parameter",
                    "model",
                    "comparability_note",
                    "source_state",
                    "target_state",
                }:
                    parsed[key] = value
                    continue
                if value.lower() in {"true", "false"}:
                    parsed[key] = value.lower() == "true"
                    continue
                try:
                    parsed[key] = int(value)
                    continue
                except ValueError:
                    pass
                try:
                    parsed[key] = float(value)
                    continue
                except ValueError:
                    parsed[key] = value
            rows.append(parsed)
        return rows


def optional_payload_string(payload: dict[str, object], key: str) -> str | None:
    value = payload.get(key)
    return value if isinstance(value, str) else None


def field_tolerance(case: GeigerParityCase, key: str) -> float:
    if case.field_tolerances and key in case.field_tolerances:
        return case.field_tolerances[key]
    return case.tolerance


def row_field_tolerance(case: GeigerParityCase, key: str) -> float:
    if case.row_field_tolerances and key in case.row_field_tolerances:
        return case.row_field_tolerances[key]
    return case.tolerance


def isclose(left: object, right: object, *, tolerance: float) -> bool:
    if isinstance(left, bool) or isinstance(right, bool):
        return left == right
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        return math.isclose(
            float(left),
            float(right),
            rel_tol=tolerance,
            abs_tol=tolerance,
        )
    return left == right


def mismatch_reason(
    case: GeigerParityCase,
    *,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
) -> str | None:
    if reference_summary is None or bijux_summary is None:
        return "summary_missing"
    for key in case.comparison_fields:
        if key not in reference_summary:
            return f"reference_summary_missing:{key}"
        if key not in bijux_summary:
            return f"bijux_summary_missing:{key}"
        if not isclose(
            reference_summary[key],
            bijux_summary[key],
            tolerance=field_tolerance(case, key),
        ):
            return f"summary_field_mismatch:{key}"
    return None


def row_mismatch_reason(
    case: GeigerParityCase,
    *,
    reference_rows: list[dict[str, object]] | None,
    bijux_rows: list[dict[str, object]] | None,
) -> str | None:
    if case.row_comparison_policy == "summary-only":
        return None
    if reference_rows is None or bijux_rows is None:
        return "rows_missing"
    normalized_reference_rows = normalized_rows(case, reference_rows)
    normalized_bijux_rows = normalized_rows(case, bijux_rows)
    if len(normalized_reference_rows) != len(normalized_bijux_rows):
        return "row_count_mismatch"
    for reference_row, bijux_row in zip(
        normalized_reference_rows,
        normalized_bijux_rows,
        strict=True,
    ):
        reference_id = row_identifier(reference_row)
        bijux_id = row_identifier(bijux_row)
        if reference_id != bijux_id:
            return "row_identifier_mismatch"
        if set(reference_row) != set(bijux_row):
            return "row_field_set_mismatch"
        for key in reference_row:
            if not isclose(
                reference_row.get(key),
                bijux_row.get(key),
                tolerance=row_field_tolerance(case, key),
            ):
                return f"row_field_mismatch:{reference_id}:{key}"
    return None


def row_identifier(row: dict[str, object]) -> object:
    identifier = row.get("parameter", row.get("model"))
    if identifier is not None:
        return identifier
    source_state = row.get("source_state")
    target_state = row.get("target_state")
    if source_state is not None and target_state is not None:
        return f"{source_state}->{target_state}"
    return None


def normalized_rows(
    case: GeigerParityCase,
    rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    normalized = [dict(row) for row in rows]
    if case.operation != "compare-fitcontinuous-models":
        if case.operation == "fit-discrete-mk":
            normalized.sort(
                key=lambda row: (
                    str(row.get("source_state", "")),
                    str(row.get("target_state", "")),
                )
            )
        return normalized
    for row in normalized:
        row.pop("rank", None)
    normalized.sort(key=lambda row: str(row.get("model", "")))
    return normalized


def parameter_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for parameter in ("root_state", "rate", "log_likelihood", "aic", "aicc"):
        value = summary.get(parameter)
        if value in {None, ""}:
            continue
        rows.append({"parameter": parameter, "value": value})
    parameter_name = summary.get("parameter_name")
    parameter_value = summary.get("parameter_value")
    if (
        isinstance(parameter_name, str)
        and parameter_name
        and parameter_value
        not in {
            None,
            "",
        }
    ):
        rows.append({"parameter": parameter_name, "value": parameter_value})
    return rows


def comparison_rows(report) -> list[dict[str, object]]:
    return [
        {
            "model": row.model,
            "rank": "" if row.rank is None else row.rank,
            "parameter_count": row.parameter_count,
            "log_likelihood": row.log_likelihood,
            "aic": row.aic,
            "aicc": row.aicc,
            "delta_aic": row.delta_aic,
            "delta_aicc": row.delta_aicc,
            "selected": row.selected,
            "comparable": row.comparable,
            "likelihood_constant_policy": (
                ""
                if row.likelihood_constant_policy is None
                else row.likelihood_constant_policy
            ),
            "comparability_note": (
                "" if row.comparability_note is None else row.comparability_note
            ),
        }
        for row in report.rows
    ]
