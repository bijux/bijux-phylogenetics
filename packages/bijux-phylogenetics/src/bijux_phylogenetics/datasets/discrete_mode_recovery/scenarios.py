from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.comparative.discrete_mode_recovery import (
    DiscreteModeRecoveryScenario,
)
from bijux_phylogenetics.simulation import DiscreteHistoryRateRow

from .models import DiscreteModeRecoveryPanelDataset


def load_discrete_mode_recovery_panel_scenarios(
    dataset: DiscreteModeRecoveryPanelDataset,
) -> list[DiscreteModeRecoveryScenario]:
    return _load_scenarios(dataset.simulation_cases_path, dataset.dataset_root)


def count_discrete_mode_recovery_panel_scenarios(
    path: Path,
    dataset_root: Path,
) -> int:
    return len(_load_scenarios(path, dataset_root))


def _load_scenarios(
    path: Path,
    dataset_root: Path,
) -> list[DiscreteModeRecoveryScenario]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            DiscreteModeRecoveryScenario(
                case_id=row["case_id"],
                label=row["label"],
                generating_model=row["generating_model"],
                expected_selected_model=(
                    None
                    if not row["expected_selected_model"]
                    else row["expected_selected_model"]
                ),
                states=_split_items(row["states"]),
                rate_rows=_parse_rate_rows(row["rate_rows"]),
                root_state=row["root_state"],
                seed=int(row["seed"]),
                tree_path=dataset_root / row["tree_file"],
                transform=_optional_string(row.get("transform", "")),
                transform_parameter_value=_optional_float(
                    row.get("transform_parameter_value", "")
                ),
                rate_tolerance=_optional_float(row["rate_tolerance"]),
                parameter_tolerances=_parse_parameter_tolerances(
                    row.get("parameter_tolerances", "")
                ),
                lambda_bounds=_parse_bounds(
                    row.get("lambda_bounds", ""),
                    default=(0.0, 1.0),
                ),
                kappa_bounds=_parse_bounds(
                    row.get("kappa_bounds", ""),
                    default=(0.0, 1.0),
                ),
                delta_bounds=_parse_bounds(
                    row.get("delta_bounds", ""),
                    default=(0.006737947, 3.0),
                ),
                early_burst_bounds=_parse_bounds(
                    row.get("early_burst_bounds", ""),
                    default=(-10.0, 10.0),
                ),
                expected_overparameterized=(
                    row["expected_overparameterized"] == "true"
                ),
                expected_warning_kinds=_split_items(row["expected_warning_kinds"]),
                notes=row["notes"],
            )
            for row in reader
        ]


def _parse_rate_rows(value: str) -> list[DiscreteHistoryRateRow]:
    rows: list[DiscreteHistoryRateRow] = []
    for item in (entry for entry in value.split(";") if entry):
        pair, rate_text = item.split("=")
        source_state, target_state = pair.split(">")
        rows.append(
            DiscreteHistoryRateRow(
                source_state=source_state,
                target_state=target_state,
                rate=float(rate_text),
            )
        )
    return rows


def _optional_float(value: str) -> float | None:
    return None if not value else float(value)


def _optional_string(value: str) -> str | None:
    return value if value else None


def _parse_parameter_tolerances(value: str) -> dict[str, float]:
    tolerances: dict[str, float] = {}
    for item in _split_items(value):
        parameter, tolerance = item.split("=")
        tolerances[parameter] = float(tolerance)
    return tolerances


def _parse_bounds(
    value: str,
    *,
    default: tuple[float, float],
) -> tuple[float, float]:
    if not value:
        return default
    left, right = value.split(":")
    return (float(left), float(right))


def _split_items(value: str) -> list[str]:
    if not value:
        return []
    return [item for item in value.split(",") if item]
