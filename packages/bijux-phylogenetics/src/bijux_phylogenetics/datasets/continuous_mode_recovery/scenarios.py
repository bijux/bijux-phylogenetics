from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.comparative.continuous_mode_recovery import (
    ContinuousModeRecoveryScenario,
)

from .models import ContinuousModeRecoveryPanelDataset


def load_continuous_mode_recovery_panel_scenarios(
    dataset: ContinuousModeRecoveryPanelDataset,
) -> list[ContinuousModeRecoveryScenario]:
    return _load_scenarios(dataset.simulation_cases_path, dataset.dataset_root)


def count_continuous_mode_recovery_panel_scenarios(
    path: Path,
    dataset_root: Path,
) -> int:
    return len(_load_scenarios(path, dataset_root))


def _load_scenarios(
    path: Path,
    dataset_root: Path,
) -> list[ContinuousModeRecoveryScenario]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            ContinuousModeRecoveryScenario(
                case_id=row["case_id"],
                label=row["label"],
                generating_model=row["generating_model"],
                expected_selected_model=(
                    None
                    if not row["expected_selected_model"]
                    else row["expected_selected_model"]
                ),
                root_state=float(row["root_state"]),
                sigma=float(row["sigma"]),
                seed=int(row["seed"]),
                tree_path=dataset_root / row["tree_file"],
                alpha=_optional_float(row["alpha"]),
                theta=_optional_float(row["theta"]),
                rate_change=_optional_float(row["rate_change"]),
                lambda_value=_optional_float(row["lambda_value"]),
                kappa=_optional_float(row["kappa"]),
                delta=_optional_float(row["delta"]),
                candidate_modes=tuple(
                    item for item in row["candidate_modes"].split(",") if item
                ),
                lambda_bounds=(float(row["lambda_lower"]), float(row["lambda_upper"])),
                kappa_bounds=(float(row["kappa_lower"]), float(row["kappa_upper"])),
                delta_bounds=(float(row["delta_lower"]), float(row["delta_upper"])),
                ou_bounds=(float(row["ou_lower"]), float(row["ou_upper"])),
                early_burst_bounds=(
                    float(row["early_burst_lower"]),
                    float(row["early_burst_upper"]),
                ),
                parameter_tolerances=_build_parameter_tolerances(row),
                expected_warning_kinds=_split_items(row["expected_warning_kinds"]),
                notes=row["notes"],
            )
            for row in reader
        ]


def _build_parameter_tolerances(row: dict[str, str]) -> dict[str, float]:
    tolerances: dict[str, float] = {}
    for field_name, parameter_name in (
        ("sigma_squared_tolerance", "sigma_squared"),
        ("alpha_tolerance", "alpha"),
        ("theta_tolerance", "theta"),
        ("rate_change_tolerance", "rate_change"),
        ("lambda_tolerance", "lambda"),
        ("kappa_tolerance", "kappa"),
        ("delta_tolerance", "delta"),
    ):
        if row[field_name]:
            tolerances[parameter_name] = float(row[field_name])
    return tolerances


def _optional_float(value: str) -> float | None:
    return None if not value else float(value)


def _split_items(value: str) -> list[str]:
    if not value:
        return []
    return [item for item in value.split(",") if item]
