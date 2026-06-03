from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

FITCONTINUOUS_REFERENCE_MODEL_ORDER = (
    "BM",
    "white",
    "lambda",
    "kappa",
    "delta",
    "OU",
    "EB",
)


@dataclass(frozen=True, slots=True)
class GeigerParityCase:
    """One governed live `geiger` parity case."""

    case_id: str
    fixture_id: str
    function_name: str
    python_function_name: str
    operation: str
    model_name: str
    python_mode: str
    input_fixtures: tuple[Path, ...]
    tolerance: float
    trait_name: str
    taxon_column: str | None = None
    optimizer_settings: dict[str, object] | None = None
    candidate_model_names: tuple[str, ...] | None = None
    reference_control: dict[str, object] | None = None
    coarse_grid_point_count: int | None = None
    fine_grid_point_count: int | None = None
    initial_parameter_value: float | None = None
    comparison_fields: tuple[str, ...] = ()
    row_comparison_policy: str = "full"
    lambda_bounds: tuple[float, float] | None = None
    discrete_transform_name: str | None = None
    kappa_bounds: tuple[float, float] | None = None
    delta_bounds: tuple[float, float] | None = None
    ou_bounds: tuple[float, float] | None = None
    early_burst_bounds: tuple[float, float] | None = None
    field_tolerances: dict[str, float] | None = None
    row_field_tolerances: dict[str, float] | None = None
