from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess  # nosec B404
import sys
import textwrap
import time
import tracemalloc

from bijux_phylogenetics.comparative.evolutionary_modes import (
    fit_continuous_evolutionary_mode,
)

from .case_definitions import ContinuousCaseDefinition


def measure(callback):
    tracemalloc.start()
    started = time.perf_counter()
    result = callback()
    runtime_seconds = time.perf_counter() - started
    _, peak_memory_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result, round(runtime_seconds, 15), peak_memory_bytes


def measure_continuous_fit(
    *,
    case_definition: ContinuousCaseDefinition,
    tree_path: Path,
    traits_path: Path,
):
    if case_definition.timeout_seconds is None:
        return measure(
            lambda: fit_continuous_evolutionary_mode(
                tree_path,
                traits_path,
                trait="value",
                mode=case_definition.fit_mode,
                search_controls=case_definition.search_controls,
                lambda_bounds=(0.0, 1.0),
            )
        )
    payload = run_continuous_fit_subprocess(
        case_definition=case_definition,
        tree_path=tree_path,
        traits_path=traits_path,
    )
    return (
        continuous_report_from_payload(
            tree_path=tree_path,
            traits_path=traits_path,
            payload=payload,
        ),
        float(payload["runtime_seconds"]),
        int(payload["peak_memory_bytes"]),
    )


def run_continuous_fit_subprocess(
    *,
    case_definition: ContinuousCaseDefinition,
    tree_path: Path,
    traits_path: Path,
) -> dict[str, object]:
    command = [
        sys.executable,
        "-c",
        textwrap.dedent(
            """
            import json
            import tracemalloc
            import time
            from pathlib import Path

            from bijux_phylogenetics.comparative.evolutionary_modes import (
                fit_continuous_evolutionary_mode,
            )

            tree_path = Path(__import__("sys").argv[1])
            traits_path = Path(__import__("sys").argv[2])
            mode = __import__("sys").argv[3]
            coarse_grid_point_count = __import__("sys").argv[4]
            fine_grid_point_count = __import__("sys").argv[5]

            search_controls = None
            if coarse_grid_point_count != "none":
                from bijux_phylogenetics.comparative.evolutionary_modes import (
                    ContinuousModeSearchControls,
                )

                search_controls = ContinuousModeSearchControls(
                    coarse_grid_point_count=int(coarse_grid_point_count),
                    fine_grid_point_count=int(fine_grid_point_count),
                )

            tracemalloc.start()
            started = time.perf_counter()
            report = fit_continuous_evolutionary_mode(
                tree_path,
                traits_path,
                trait="value",
                mode=mode,
                search_controls=search_controls,
                lambda_bounds=(0.0, 1.0),
            )
            runtime_seconds = time.perf_counter() - started
            _, peak_memory_bytes = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            diagnostics = report.optimizer_diagnostics
            boundary = report.boundary_assessment
            payload = {
                "runtime_seconds": runtime_seconds,
                "peak_memory_bytes": peak_memory_bytes,
                "taxon_count": report.taxon_count,
                "mode": report.mode,
                "parameter_name": report.parameter_name,
                "parameter_value": report.parameter_value,
                "rate": report.rate,
                "log_likelihood": report.log_likelihood,
                "aic": report.aic,
                "aicc": report.aicc,
                "optimizer_name": None if diagnostics is None else diagnostics.optimizer_name,
                "optimizer_step_count": None if diagnostics is None else diagnostics.function_evaluation_count,
                "converged": None if diagnostics is None else diagnostics.converged,
                "hit_lower_parameter_boundary": None if diagnostics is None else diagnostics.hit_lower_boundary,
                "hit_upper_parameter_boundary": None if diagnostics is None else diagnostics.hit_upper_boundary,
                "stable_conclusion_supported": None if boundary is None else boundary.stable_conclusion_supported,
            }
            print(json.dumps(payload))
            """
        ),
        str(tree_path),
        str(traits_path),
        case_definition.fit_mode,
        (
            "none"
            if case_definition.search_controls is None
            else str(case_definition.search_controls.coarse_grid_point_count)
        ),
        (
            "none"
            if case_definition.search_controls is None
            else str(case_definition.search_controls.fine_grid_point_count)
        ),
    ]
    result = subprocess.run(  # nosec B603
        command,
        capture_output=True,
        check=False,
        text=True,
        timeout=case_definition.timeout_seconds,
        env=dict(os.environ),
    )
    if result.returncode != 0:
        raise RuntimeError(
            result.stderr.strip() or result.stdout.strip() or "subprocess fit failed"
        )
    return json.loads(result.stdout)


def continuous_report_from_payload(
    *,
    tree_path: Path,
    traits_path: Path,
    payload: dict[str, object],
):
    from types import SimpleNamespace

    optimizer_name = optional_string(payload.get("optimizer_name"))
    optimizer_step_count = optional_int(payload.get("optimizer_step_count"))
    converged = payload.get("converged")
    hit_lower_parameter_boundary = payload.get("hit_lower_parameter_boundary")
    hit_upper_parameter_boundary = payload.get("hit_upper_parameter_boundary")
    stable_conclusion_supported = payload.get("stable_conclusion_supported")
    return SimpleNamespace(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_count=optional_int(payload.get("taxon_count")),
        mode=optional_string(payload.get("mode")),
        parameter_name=optional_string(payload.get("parameter_name")),
        parameter_value=optional_float(payload.get("parameter_value")),
        rate=optional_float(payload.get("rate")),
        log_likelihood=optional_float(payload.get("log_likelihood")),
        aic=optional_float(payload.get("aic")),
        aicc=optional_float(payload.get("aicc")),
        optimizer_diagnostics=(
            None
            if optimizer_name is None and optimizer_step_count is None
            else SimpleNamespace(
                optimizer_name=optimizer_name,
                function_evaluation_count=optimizer_step_count,
                converged=converged,
                hit_lower_boundary=hit_lower_parameter_boundary,
                hit_upper_boundary=hit_upper_parameter_boundary,
            )
        ),
        boundary_assessment=(
            None
            if stable_conclusion_supported is None
            else SimpleNamespace(
                stable_conclusion_supported=stable_conclusion_supported,
            )
        ),
        identifiability_warnings=[],
    )


def optional_float(value: object) -> float | None:
    if isinstance(value, (float, int)):
        return float(value)
    return None


def optional_int(value: object) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return None


def optional_string(value: object) -> str | None:
    if isinstance(value, str):
        return value
    return None
