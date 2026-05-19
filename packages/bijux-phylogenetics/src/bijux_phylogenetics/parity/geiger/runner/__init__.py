from __future__ import annotations

import csv
from dataclasses import asdict
import json
from pathlib import Path
import shutil
import subprocess  # nosec B404 - parity helpers invoke repository-owned reference commands
import tempfile

from ..optimizer_triage import build_geiger_optimizer_triage_rows
from ..boundary_warning_registry import build_geiger_boundary_warning_rows
from ..likelihood_policy import build_geiger_likelihood_policy_rows
from ..model_confidence import build_geiger_model_confidence_rows
from ..parameterization_registry import build_geiger_parameterization_registry_rows
from ..registry import (
    GeigerParityCase,
    list_geiger_parity_cases as list_geiger_parity_cases,
)
from .models import (
    GeigerParityObservation,
    GeigerParityReport,
    GeigerParitySummaryRow,
)
from .comparison import (
    load_json as _load_json,
    load_rows_table as _load_rows_table,
    mismatch_reason as _mismatch_reason,
    optional_payload_string as _optional_payload_string,
    row_mismatch_reason as _row_mismatch_reason,
)
from .continuous_payloads import (
    build_bijux_continuous_case_payload as _build_bijux_continuous_case_payload,
    build_bijux_model_comparison_payload as _build_bijux_model_comparison_payload,
)
from .discrete_payloads import (
    build_bijux_discrete_case_payload as _build_bijux_discrete_case_payload,
)
from .runtime import (
    bijux_commit as _bijux_commit,
    bijux_version as _bijux_version,
    failure_root as _failure_root,
    geiger_runner_path as _geiger_runner_path,
    reference_environment as _reference_environment,
    repository_root as _repository_root,
    selected_cases as _selected_cases,
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




def _write_rows_table(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = list(rows[0].keys()) if rows else ["parameter", "value"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _persist_failure_bundle(
    *,
    failure_root: Path,
    case: GeigerParityCase,
    case_file: Path,
    execution_root: Path,
    execution_payload: dict[str, object] | None,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
    reference_rows: list[dict[str, object]] | None,
    bijux_rows: list[dict[str, object]] | None,
    reference_error: dict[str, object] | None,
    bijux_error: dict[str, object] | None,
    mismatch_reason: str,
) -> Path:
    artifact_root = failure_root / case.case_id
    if artifact_root.exists():
        shutil.rmtree(artifact_root)
    artifact_root.mkdir(parents=True, exist_ok=True)
    shutil.copy2(case_file, artifact_root / "case.json")
    if execution_root.exists():
        shutil.copytree(
            execution_root,
            artifact_root / "reference",
            dirs_exist_ok=True,
        )
    if execution_payload is not None:
        (artifact_root / "reference-execution.json").write_text(
            json.dumps(execution_payload, indent=2),
            encoding="utf-8",
        )
    if reference_summary is not None:
        (artifact_root / "reference-summary.json").write_text(
            json.dumps(reference_summary, indent=2),
            encoding="utf-8",
        )
    if bijux_summary is not None:
        (artifact_root / "bijux-summary.json").write_text(
            json.dumps(bijux_summary, indent=2),
            encoding="utf-8",
        )
    if reference_rows is not None:
        _write_rows_table(artifact_root / "reference-parameters.tsv", reference_rows)
    if bijux_rows is not None:
        _write_rows_table(artifact_root / "bijux-parameters.tsv", bijux_rows)
    if reference_error is not None:
        (artifact_root / "reference-error.json").write_text(
            json.dumps(reference_error, indent=2),
            encoding="utf-8",
        )
    if bijux_error is not None:
        (artifact_root / "bijux-error.json").write_text(
            json.dumps(bijux_error, indent=2),
            encoding="utf-8",
        )
    (artifact_root / "mismatch-reason.txt").write_text(
        mismatch_reason,
        encoding="utf-8",
    )
    return artifact_root


def _summary_rows(
    observations: list[GeigerParityObservation],
) -> list[GeigerParitySummaryRow]:
    by_function: dict[str, list[GeigerParityObservation]] = {}
    for observation in observations:
        by_function.setdefault(observation.function_name, []).append(observation)
    return [
        GeigerParitySummaryRow(
            function_name=function_name,
            case_count=len(items),
            passed_case_count=sum(1 for item in items if item.status == "passed"),
            failed_case_count=sum(1 for item in items if item.status == "failed"),
            skipped_case_count=sum(1 for item in items if item.status == "skipped"),
        )
        for function_name, items in sorted(by_function.items())
    ]


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
                process = subprocess.run(  # nosec
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


def write_geiger_parity_summary_table(
    path: Path,
    report: GeigerParityReport,
) -> Path:
    """Write one row per governed `geiger` function summary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "function_name",
                "case_count",
                "passed_case_count",
                "failed_case_count",
                "skipped_case_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.summary_rows:
            writer.writerow(asdict(row))
    return path


def write_geiger_parity_observation_table(
    path: Path,
    report: GeigerParityReport,
) -> Path:
    """Write one row per governed `geiger` parity observation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_id",
                "fixture_id",
                "function_name",
                "python_function_name",
                "input_fixtures",
                "model_name",
                "optimizer_settings",
                "tolerance",
                "r_version",
                "geiger_version",
                "bijux_version",
                "bijux_commit",
                "status",
                "passed",
                "mismatch_reason",
                "reproducible_artifact_root",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for observation in report.observations:
            writer.writerow(
                {
                    "case_id": observation.case_id,
                    "fixture_id": observation.fixture_id,
                    "function_name": observation.function_name,
                    "python_function_name": observation.python_function_name,
                    "input_fixtures": json.dumps(
                        [str(path) for path in observation.input_fixtures]
                    ),
                    "model_name": observation.model_name,
                    "optimizer_settings": json.dumps(
                        observation.optimizer_settings,
                        sort_keys=True,
                    ),
                    "tolerance": format(observation.tolerance, ".12g"),
                    "r_version": observation.r_version or "",
                    "geiger_version": observation.geiger_version or "",
                    "bijux_version": observation.bijux_version,
                    "bijux_commit": observation.bijux_commit or "",
                    "status": observation.status,
                    "passed": str(observation.passed).lower(),
                    "mismatch_reason": observation.mismatch_reason or "",
                    "reproducible_artifact_root": ""
                    if observation.reproducible_artifact_root is None
                    else str(observation.reproducible_artifact_root),
                }
            )
    return path
