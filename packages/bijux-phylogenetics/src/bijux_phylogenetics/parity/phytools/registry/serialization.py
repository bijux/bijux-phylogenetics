from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from .models import PhytoolsParityCase


def write_case_file(path: Path, case: PhytoolsParityCase) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "case_id": case.case_id,
        "fixture_id": case.fixture_id,
        "function_name": case.function_name,
        "operation": case.operation,
        "input_fixtures": [str(path) for path in case.input_fixtures],
        "trait_name": case.trait_name,
        "taxon_column": case.taxon_column,
        "discrete_model": case.discrete_model,
        "root_prior_mode": case.root_prior_mode,
        "tolerance": case.tolerance,
        "permutation_count": case.permutation_count,
        "permutation_seed": case.permutation_seed,
        "stochastic_map_replicate_count": case.stochastic_map_replicate_count,
        "stochastic_map_seed": case.stochastic_map_seed,
        "density_resolution": case.density_resolution,
        "focal_state": case.focal_state,
        "simulation_states": case.simulation_states,
        "simulation_rate_rows": (
            None
            if case.simulation_rate_rows is None
            else [asdict(row) for row in case.simulation_rate_rows]
        ),
        "simulation_root_state": case.simulation_root_state,
        "simulation_root_state_probabilities": case.simulation_root_state_probabilities,
        "simulation_replicate_count": case.simulation_replicate_count,
        "simulation_seed": case.simulation_seed,
        "continuous_root_state": case.continuous_root_state,
        "continuous_sigma_squared": case.continuous_sigma_squared,
        "continuous_replicate_count": case.continuous_replicate_count,
        "continuous_seed": case.continuous_seed,
        "continuous_trait_names": case.continuous_trait_names,
        "continuous_root_states": case.continuous_root_states,
        "continuous_covariance_matrix": case.continuous_covariance_matrix,
        "comparative_formula": case.comparative_formula,
        "comparative_predictors": case.comparative_predictors,
        "comparative_lambda_value": case.comparative_lambda_value,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path
