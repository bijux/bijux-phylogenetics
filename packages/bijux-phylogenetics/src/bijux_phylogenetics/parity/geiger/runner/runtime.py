from __future__ import annotations

from importlib import metadata
import json
import os
from pathlib import Path
import shutil

# Parity helpers invoke repository-owned reference commands under governed paths.
import subprocess  # nosec B404

from bijux_phylogenetics.parity.geiger.registry import (
    GeigerParityCase,
    list_geiger_parity_cases,
)


def repository_root() -> Path:
    return Path(__file__).resolve().parents[7]


def geiger_runner_path() -> Path:
    return (
        Path(__file__).resolve().parents[3]
        / "resources"
        / "reference"
        / "geiger_parity_runner.R"
    )


def failure_root() -> Path:
    return repository_root() / "artifacts" / "geiger-parity-failures"


def reference_environment() -> dict[str, str]:
    environment = dict(os.environ)
    r_library = repository_root() / "artifacts" / "r-lib"
    if "R_LIBS_USER" not in environment and r_library.is_dir():
        environment["R_LIBS_USER"] = str(r_library)
    return environment


def bijux_version() -> str:
    try:
        return metadata.version("bijux-phylogenetics")
    except metadata.PackageNotFoundError:
        return "0.1.0"


def bijux_commit() -> str | None:
    git_executable = shutil.which("git")
    if git_executable is None:
        return None
    result = subprocess.run(  # nosec B603
        [git_executable, "rev-parse", "--short", "HEAD"],
        capture_output=True,
        check=False,
        cwd=repository_root(),
        text=True,
    )
    if result.returncode != 0:
        return None
    commit = result.stdout.strip()
    return commit or None


def selected_cases(case_ids: list[str] | None) -> list[GeigerParityCase]:
    registry = {case.case_id: case for case in list_geiger_parity_cases()}
    if case_ids is None:
        return list(registry.values())
    missing = [case_id for case_id in case_ids if case_id not in registry]
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"unknown geiger parity case id(s): {missing_text}")
    return [registry[case_id] for case_id in case_ids]


def write_case_file(path: Path, case: GeigerParityCase) -> Path:
    payload = {
        "case_id": case.case_id,
        "fixture_id": case.fixture_id,
        "function_name": case.function_name,
        "python_function_name": case.python_function_name,
        "operation": case.operation,
        "model_name": case.model_name,
        "python_mode": case.python_mode,
        "input_fixtures": [str(item) for item in case.input_fixtures],
        "tolerance": case.tolerance,
        "trait_name": case.trait_name,
        "taxon_column": case.taxon_column,
        "optimizer_settings": case.optimizer_settings,
        "candidate_model_names": None
        if case.candidate_model_names is None
        else list(case.candidate_model_names),
        "reference_control": case.reference_control,
        "coarse_grid_point_count": case.coarse_grid_point_count,
        "fine_grid_point_count": case.fine_grid_point_count,
        "initial_parameter_value": case.initial_parameter_value,
        "comparison_fields": list(case.comparison_fields),
        "row_comparison_policy": case.row_comparison_policy,
        "lambda_bounds": None
        if case.lambda_bounds is None
        else list(case.lambda_bounds),
        "discrete_transform_name": case.discrete_transform_name,
        "kappa_bounds": None if case.kappa_bounds is None else list(case.kappa_bounds),
        "delta_bounds": None if case.delta_bounds is None else list(case.delta_bounds),
        "ou_bounds": None if case.ou_bounds is None else list(case.ou_bounds),
        "early_burst_bounds": None
        if case.early_burst_bounds is None
        else list(case.early_burst_bounds),
        "field_tolerances": case.field_tolerances,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path
