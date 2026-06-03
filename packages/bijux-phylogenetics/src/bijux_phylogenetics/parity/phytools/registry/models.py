from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.simulation import DiscreteHistoryRateRow


@dataclass(frozen=True, slots=True)
class PhytoolsParityCase:
    """One governed live `phytools` parity case."""

    case_id: str
    fixture_id: str
    function_name: str
    python_function_name: str
    operation: str
    input_fixtures: tuple[Path, ...]
    tolerance: float
    trait_name: str
    taxon_column: str | None = None
    discrete_model: str | None = None
    root_prior_mode: str = "equal"
    permutation_count: int | None = None
    permutation_seed: int | None = None
    stochastic_map_replicate_count: int | None = None
    stochastic_map_seed: int | None = None
    density_resolution: int | None = None
    focal_state: str | None = None
    simulation_states: tuple[str, ...] | None = None
    simulation_rate_rows: tuple[DiscreteHistoryRateRow, ...] | None = None
    simulation_root_state: str | None = None
    simulation_root_state_probabilities: dict[str, float] | None = None
    simulation_replicate_count: int | None = None
    simulation_seed: int | None = None
    continuous_root_state: float | None = None
    continuous_sigma_squared: float | None = None
    continuous_replicate_count: int | None = None
    continuous_seed: int | None = None
    continuous_trait_names: tuple[str, ...] | None = None
    continuous_root_states: tuple[float, ...] | None = None
    continuous_covariance_matrix: tuple[tuple[float, ...], ...] | None = None
    comparative_formula: str | None = None
    comparative_predictors: tuple[str, ...] | None = None
    comparative_lambda_value: float | None = None
    field_tolerances: dict[str, float] | None = None
    row_field_tolerances: dict[str, float] | None = None
    compare_rows: bool = True
