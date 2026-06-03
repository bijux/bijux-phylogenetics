from __future__ import annotations

from dataclasses import replace
from pathlib import Path

# Parity helpers invoke repository-owned reference commands under governed paths.
import subprocess  # nosec B404
import tempfile

from ..registry import _selected_cases, _write_case_file
from .comparison import (
    _apply_expected_status_contract,
    _determine_reference_mismatch_reason,
)
from .dispatch import _build_bijux_case_payload
from .failure_artifacts import _persist_failure_bundle
from .models import ApeParityObservation, ApeParityReport
from .normalization import _load_json, _optional_payload_string
from .reference_payloads import _load_reference_case_payload
from .reporting import build_ape_parity_report
from .runtime import (
    ape_runner_path as _ape_runner_path,
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
    reference_environment as _reference_environment,
)
from .runtime import (
    repository_root as _repository_root,
)
from .tree_payloads import _materialize_reference_input


def _write_reference_process_streams(
    artifact_root: Path,
    *,
    process_stdout: str,
    process_stderr: str,
) -> None:
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


def run_ape_parity_cases(
    *,
    case_ids: list[str] | None = None,
    rscript_executable: str = "Rscript",
    failure_root: Path | None = None,
    fixtures_root: Path | None = None,
) -> ApeParityReport:
    """Run governed live `ape` parity cases through the checked-in R runner."""
    selected = _selected_cases(case_ids=case_ids, fixtures_root=fixtures_root)
    observations: list[ApeParityObservation] = []
    active_failure_root = _failure_root() if failure_root is None else failure_root
    bijux_version = _bijux_version()
    bijux_commit = _bijux_commit()
    for case in selected:
        with tempfile.TemporaryDirectory(
            prefix=f"bijux-ape-parity-{case.case_id}-"
        ) as tmpdir:
            working_root = Path(tmpdir)
            reference_input_path = _materialize_reference_input(case, working_root)
            reference_case = replace(case, input_fixture=reference_input_path)
            case_file = _write_case_file(working_root / "case.json", reference_case)
            execution_root = working_root / "reference"
            execution_root.mkdir(parents=True, exist_ok=True)
            bijux_summary: dict[str, object] | None = None
            bijux_rows: list[dict[str, object]] | None = None
            bijux_normalized_text: str | None = None
            bijux_error: dict[str, object] | None = None
            try:
                (
                    bijux_summary,
                    bijux_rows,
                    bijux_normalized_text,
                ) = _build_bijux_case_payload(case)
            except Exception as error:
                bijux_error = {
                    "error_type": type(error).__name__,
                    "message": str(error),
                }
            execution_payload: dict[str, object] | None = None
            reference_summary: dict[str, object] | None = None
            reference_error: dict[str, object] | None = None
            reference_rows: list[dict[str, object]] | None = None
            reference_normalized_text: str | None = None
            status = "failed"
            mismatch_reason: str | None = None
            artifact_root: Path | None = None
            r_version: str | None = None
            ape_version: str | None = None
            process_stdout = ""
            process_stderr = ""
            try:
                # Repository-owned R parity runner.
                process = subprocess.run(  # nosec B603
                    [
                        rscript_executable,
                        str(_ape_runner_path()),
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
            if process is None:
                pass
            elif process.returncode != 0:
                mismatch_reason = "reference_execution_failed"
            else:
                execution_path = execution_root / "reference-execution.json"
                if not execution_path.exists():
                    mismatch_reason = "reference_execution_failed"
                else:
                    execution_payload = _load_json(execution_path)
                    r_version = _optional_payload_string(execution_payload, "r_version")
                    ape_version = _optional_payload_string(
                        execution_payload,
                        "ape_version",
                    )
                    execution_status = execution_payload.get("status")
                    if execution_status == "unavailable":
                        status = "skipped"
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason",
                                "ape_package_unavailable",
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
                    else:
                        (
                            reference_summary,
                            reference_rows,
                            reference_normalized_text,
                        ) = _load_reference_case_payload(case, execution_root)
                        mismatch_reason = _determine_reference_mismatch_reason(
                            case=case,
                            execution_root=execution_root,
                            reference_summary=reference_summary,
                            bijux_summary=bijux_summary,
                            reference_rows=reference_rows,
                            bijux_rows=bijux_rows,
                            reference_normalized_text=reference_normalized_text,
                            bijux_normalized_text=bijux_normalized_text,
                        )
                        if mismatch_reason is None:
                            status = "passed"
            status, mismatch_reason = _apply_expected_status_contract(
                case=case,
                bijux_error=bijux_error,
                reference_error=reference_error,
                status=status,
                mismatch_reason=mismatch_reason,
            )
            if status != "passed":
                artifact_root = _persist_failure_bundle(
                    failure_root=active_failure_root,
                    case=case,
                    case_file=case_file,
                    execution_root=execution_root,
                    execution_payload=execution_payload,
                    reference_summary=reference_summary,
                    bijux_summary=bijux_summary,
                    reference_error=reference_error,
                    bijux_error=bijux_error,
                    reference_rows=reference_rows,
                    bijux_rows=bijux_rows,
                    bijux_normalized_text=bijux_normalized_text,
                    mismatch_reason=mismatch_reason or "reference_execution_failed",
                )
                _write_reference_process_streams(
                    artifact_root,
                    process_stdout=process_stdout,
                    process_stderr=process_stderr,
                )
            observations.append(
                ApeParityObservation(
                    case_id=case.case_id,
                    fixture_kind=case.fixture_kind,
                    fixture_id=case.fixture_id,
                    function_name=case.function_name,
                    python_function_name=case.python_function_name,
                    input_fixture=case.input_fixture,
                    tolerance=case.tolerance,
                    r_version=r_version,
                    ape_version=ape_version,
                    bijux_version=bijux_version,
                    bijux_commit=bijux_commit,
                    status=status,
                    passed=status == "passed",
                    mismatch_reason=mismatch_reason,
                    reproducible_artifact_root=artifact_root,
                    reference_summary=reference_summary,
                    bijux_summary=bijux_summary,
                    reference_error=reference_error,
                    bijux_error=bijux_error,
                )
            )
    return build_ape_parity_report(observations)
