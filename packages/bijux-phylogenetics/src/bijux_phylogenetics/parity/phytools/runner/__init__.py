from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import subprocess  # nosec B404 - parity helpers invoke repository-owned reference commands
import tempfile

from ..registry import PhytoolsParityCase, _selected_cases
from .comparison import (
    load_json as _load_json,
    load_rows_table as _load_rows_table,
    mismatch_reason as _mismatch_reason,
    reference_rows_filename as _reference_rows_filename,
    row_mismatch_reason as _row_mismatch_reason,
)
from .comparative_payloads import build_comparative_case_payload
from .continuous_payloads import build_continuous_case_payload
from .discrete_payloads import build_discrete_case_payload
from .models import (
    PhytoolsParityObservation,
    PhytoolsParityReport,
    PhytoolsParitySummaryRow as PhytoolsParitySummaryRow,
)
from .runtime import (
    bijux_commit as _bijux_commit,
    bijux_version as _bijux_version,
    failure_root as _failure_root,
    phytools_runner_path as _phytools_runner_path,
    reference_environment as _reference_environment,
    repository_root as _repository_root,
)
from .reporting import (
    build_phytools_parity_report,
    persist_failure_bundle as _persist_failure_bundle,
    write_phytools_parity_observation_table as write_phytools_parity_observation_table,
    write_phytools_parity_summary_table as write_phytools_parity_summary_table,
)
from .signal_payloads import build_signal_case_payload


def _optional_payload_string(payload: dict[str, object], key: str) -> str | None:
    value = payload.get(key)
    return value if isinstance(value, str) else None
def _write_case_file(path: Path, case: PhytoolsParityCase) -> Path:
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


def _build_bijux_case_payload(
    case: PhytoolsParityCase,
) -> tuple[dict[str, object], list[dict[str, object]] | None]:
    tree_path = case.input_fixtures[0]
    traits_path = case.input_fixtures[1] if len(case.input_fixtures) > 1 else None
    comparative_payload = build_comparative_case_payload(
        case,
        tree_path=tree_path,
        traits_path=traits_path,
    )
    if comparative_payload is not None:
        return comparative_payload
    discrete_payload = build_discrete_case_payload(
        case,
        tree_path=tree_path,
        traits_path=traits_path,
    )
    if discrete_payload is not None:
        return discrete_payload
    continuous_payload = build_continuous_case_payload(
        case,
        tree_path=tree_path,
        traits_path=traits_path,
    )
    if continuous_payload is not None:
        return continuous_payload
    signal_payload = build_signal_case_payload(
        case,
        tree_path=tree_path,
        traits_path=traits_path,
    )
    if signal_payload is not None:
        return signal_payload
    raise ValueError(f"unsupported phytools parity operation: {case.operation}")
def run_phytools_parity_cases(
    *,
    case_ids: list[str] | None = None,
    rscript_executable: str = "Rscript",
    failure_root: Path | None = None,
) -> PhytoolsParityReport:
    """Run governed live `phytools` parity cases through the checked-in R runner."""
    selected = _selected_cases(case_ids)
    observations: list[PhytoolsParityObservation] = []
    active_failure_root = _failure_root() if failure_root is None else failure_root
    bijux_version = _bijux_version()
    bijux_commit = _bijux_commit()
    for case in selected:
        with tempfile.TemporaryDirectory(
            prefix=f"bijux-phytools-parity-{case.case_id}-"
        ) as tmpdir:
            working_root = Path(tmpdir)
            case_file = _write_case_file(working_root / "case.json", case)
            execution_root = working_root / "reference"
            execution_root.mkdir(parents=True, exist_ok=True)
            bijux_summary: dict[str, object] | None = None
            bijux_rows: list[dict[str, object]] | None = None
            bijux_error: dict[str, object] | None = None
            try:
                bijux_summary, bijux_rows = _build_bijux_case_payload(case)
            except Exception as error:
                bijux_error = {
                    "error_type": type(error).__name__,
                    "message": str(error),
                }
            execution_payload: dict[str, object] | None = None
            reference_summary: dict[str, object] | None = None
            reference_rows: list[dict[str, object]] | None = None
            reference_error: dict[str, object] | None = None
            status = "failed"
            mismatch_reason: str | None = None
            artifact_root: Path | None = None
            r_version: str | None = None
            phytools_version: str | None = None
            process_stdout = ""
            process_stderr = ""
            try:
                # Repository-owned R parity runner.
                process = subprocess.run(  # nosec
                    [
                        rscript_executable,
                        str(_phytools_runner_path()),
                        str(case_file),
                        str(execution_root),
                    ],
                    capture_output=True,
                    check=False,
                    cwd=_repository_root(),
                    env=_reference_environment(),
                    text=True,
                )
                process_stdout = process.stdout
                process_stderr = process.stderr
            except FileNotFoundError:
                process = None
                status = "skipped"
                mismatch_reason = "rscript_unavailable"
            if process is not None and process.returncode == 0:
                execution_path = execution_root / "reference-execution.json"
                summary_path = execution_root / "reference-summary.json"
                if not execution_path.exists():
                    mismatch_reason = "reference_execution_failed"
                else:
                    execution_payload = _load_json(execution_path)
                    r_version = _optional_payload_string(execution_payload, "r_version")
                    phytools_version = _optional_payload_string(
                        execution_payload, "phytools_version"
                    )
                    execution_status = execution_payload.get("status")
                    if execution_status == "unavailable":
                        status = "skipped"
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason",
                                "phytools_package_unavailable",
                            )
                        )
                    elif execution_status != "ok":
                        reference_error = {
                            "error_type": str(
                                execution_payload.get(
                                    "error_type",
                                    execution_payload.get(
                                        "mismatch_reason",
                                        "reference_execution_failed",
                                    ),
                                )
                            ),
                            "message": str(execution_payload.get("message", "")),
                        }
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason",
                                "reference_execution_failed",
                            )
                        )
                    elif not summary_path.exists():
                        mismatch_reason = "reference_summary_missing"
                    else:
                        reference_summary = _load_json(summary_path)
                        mismatch_reason = _mismatch_reason(
                            case,
                            reference_summary=reference_summary,
                            bijux_summary=bijux_summary,
                        )
                        rows_filename = _reference_rows_filename(case)
                        if mismatch_reason is None and rows_filename is not None:
                            rows_path = execution_root / rows_filename
                            if not rows_path.exists():
                                mismatch_reason = "reference_rows_missing"
                            else:
                                reference_rows = _load_rows_table(rows_path)
                                mismatch_reason = _row_mismatch_reason(
                                    case,
                                    reference_rows=reference_rows,
                                    bijux_rows=bijux_rows,
                                )
                        if mismatch_reason is None:
                            status = "passed"
            elif process is not None and process.returncode != 0:
                mismatch_reason = "reference_execution_failed"
            if status != "passed":
                artifact_root = _persist_failure_bundle(
                    failure_root=active_failure_root,
                    case=case,
                    case_file=case_file,
                    execution_root=execution_root,
                    execution_payload=execution_payload,
                    reference_summary=reference_summary,
                    bijux_summary=bijux_summary,
                    reference_rows=reference_rows,
                    bijux_rows=bijux_rows,
                    reference_error=reference_error,
                    bijux_error=bijux_error,
                    mismatch_reason=mismatch_reason or "reference_execution_failed",
                )
                if process_stdout:
                    (artifact_root / "reference-stdout.txt").write_text(
                        process_stdout,
                        encoding="utf-8",
                    )
                if process_stderr:
                    (artifact_root / "reference-stderr.txt").write_text(
                        process_stderr,
                        encoding="utf-8",
                    )
            observations.append(
                PhytoolsParityObservation(
                    case_id=case.case_id,
                    fixture_id=case.fixture_id,
                    function_name=case.function_name,
                    python_function_name=case.python_function_name,
                    input_fixtures=case.input_fixtures,
                    tolerance=case.tolerance,
                    r_version=r_version,
                    phytools_version=phytools_version,
                    bijux_version=bijux_version,
                    bijux_commit=bijux_commit,
                    status=status,
                    passed=status == "passed",
                    mismatch_reason=mismatch_reason,
                    reproducible_artifact_root=artifact_root,
                    reference_summary=reference_summary,
                    bijux_summary=bijux_summary,
                    reference_rows=reference_rows,
                    bijux_rows=bijux_rows,
                    reference_error=reference_error,
                    bijux_error=bijux_error,
                )
            )
    return build_phytools_parity_report(observations)
