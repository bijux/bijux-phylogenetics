from __future__ import annotations

from pathlib import Path

# Parity helpers invoke repository-owned reference commands under governed paths.
import subprocess  # nosec B404
import tempfile

from ..boundary_warning_registry import build_geiger_boundary_warning_rows
from ..likelihood_policy import build_geiger_likelihood_policy_rows
from ..model_confidence import build_geiger_model_confidence_rows
from ..optimizer_triage import build_geiger_optimizer_triage_rows
from ..parameterization_registry import build_geiger_parameterization_registry_rows
from ..registry import (
    GeigerParityCase,
)
from ..registry import (
    list_geiger_parity_cases as list_geiger_parity_cases,
)
from .comparison import (
    load_json as _load_json,
)
from .comparison import (
    load_rows_table as _load_rows_table,
)
from .comparison import (
    mismatch_reason as _mismatch_reason,
)
from .comparison import (
    optional_payload_string as _optional_payload_string,
)
from .comparison import (
    row_mismatch_reason as _row_mismatch_reason,
)
from .continuous_payloads import (
    build_bijux_continuous_case_payload as _build_bijux_continuous_case_payload,
)
from .continuous_payloads import (
    build_bijux_model_comparison_payload as _build_bijux_model_comparison_payload,
)
from .discrete_payloads import (
    build_bijux_discrete_case_payload as _build_bijux_discrete_case_payload,
)
from .models import (
    GeigerParityObservation,
    GeigerParityReport,
)
from .models import (
    GeigerParitySummaryRow as GeigerParitySummaryRow,
)
from .reporting import (
    persist_failure_bundle as _persist_failure_bundle,
)
from .reporting import (
    summary_rows as _summary_rows,
)
from .reporting import (
    write_geiger_parity_observation_table as write_geiger_parity_observation_table,
)
from .reporting import (
    write_geiger_parity_summary_table as write_geiger_parity_summary_table,
)
from .runtime import (
    bijux_commit as _bijux_commit,
)
from .runtime import (
    bijux_version as _bijux_version,
)
from .runtime import (
    failure_root as _failure_root,
)
from .runtime import (
    geiger_runner_path as _geiger_runner_path,
)
from .runtime import (
    reference_environment as _reference_environment,
)
from .runtime import (
    repository_root as _repository_root,
)
from .runtime import (
    selected_cases as _selected_cases,
)
from .runtime import (
    write_case_file as _write_case_file,
)


def _build_bijux_case_payload(
    case: GeigerParityCase,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    if case.operation == "fit-discrete-mk":
        return _build_bijux_discrete_case_payload(case)
    if case.operation == "compare-fitcontinuous-models":
        return _build_bijux_model_comparison_payload(case)
    return _build_bijux_continuous_case_payload(case)


def run_geiger_parity_cases(
    *,
    case_ids: list[str] | None = None,
    rscript_executable: str = "Rscript",
    failure_root: Path | None = None,
) -> GeigerParityReport:
    """Run governed live `geiger` parity cases through the checked-in R runner."""
    selected = _selected_cases(case_ids)
    observations: list[GeigerParityObservation] = []
    active_failure_root = _failure_root() if failure_root is None else failure_root
    bijux_version = _bijux_version()
    bijux_commit = _bijux_commit()
    for case in selected:
        with tempfile.TemporaryDirectory(
            prefix=f"bijux-geiger-parity-{case.case_id}-"
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
            geiger_version: str | None = None
            process_stdout = ""
            process_stderr = ""
            try:
                process = subprocess.run(  # nosec B603
                    [
                        rscript_executable,
                        str(_geiger_runner_path()),
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
                rows_path = execution_root / "reference-parameters.tsv"
                if not execution_path.exists():
                    mismatch_reason = "reference_execution_failed"
                else:
                    execution_payload = _load_json(execution_path)
                    r_version = _optional_payload_string(execution_payload, "r_version")
                    geiger_version = _optional_payload_string(
                        execution_payload,
                        "geiger_version",
                    )
                    execution_status = execution_payload.get("status")
                    if execution_status == "unavailable":
                        status = "skipped"
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason",
                                "geiger_package_unavailable",
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
                        if rows_path.exists():
                            reference_rows = _load_rows_table(rows_path)
                        mismatch_reason = _mismatch_reason(
                            case,
                            reference_summary=reference_summary,
                            bijux_summary=bijux_summary,
                        )
                        if mismatch_reason is None:
                            if reference_rows is None:
                                mismatch_reason = "reference_rows_missing"
                            else:
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
                GeigerParityObservation(
                    case_id=case.case_id,
                    fixture_id=case.fixture_id,
                    function_name=case.function_name,
                    python_function_name=case.python_function_name,
                    input_fixtures=case.input_fixtures,
                    model_name=case.model_name,
                    optimizer_settings=case.optimizer_settings,
                    tolerance=case.tolerance,
                    r_version=r_version,
                    geiger_version=geiger_version,
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
    case_count = len(observations)
    passed_case_count = sum(1 for item in observations if item.status == "passed")
    failed_case_count = sum(1 for item in observations if item.status == "failed")
    skipped_case_count = sum(1 for item in observations if item.status == "skipped")
    return GeigerParityReport(
        observations=observations,
        optimizer_triage_rows=build_geiger_optimizer_triage_rows(observations),
        boundary_warning_rows=build_geiger_boundary_warning_rows(observations),
        likelihood_policy_rows=build_geiger_likelihood_policy_rows(observations),
        model_confidence_rows=build_geiger_model_confidence_rows(observations),
        parameterization_registry_rows=build_geiger_parameterization_registry_rows(
            observations
        ),
        summary_rows=_summary_rows(observations),
        case_count=case_count,
        passed_case_count=passed_case_count,
        failed_case_count=failed_case_count,
        skipped_case_count=skipped_case_count,
        all_passed=case_count > 0
        and passed_case_count == case_count
        and failed_case_count == 0
        and skipped_case_count == 0,
        limitations=[
            "The governed live `geiger` parity registry stays narrow in this harness tranche and only covers the currently owned continuous-mode fit surfaces.",
            "This harness requires Rscript plus the `geiger`, `ape`, and `jsonlite` R packages for live reference execution.",
        ],
    )
