from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess  # nosec B404
import time
from typing import Any

from bijux_phylogenetics.runtime.errors import (
    EngineUnavailableError,
    EngineWorkflowError,
    InvalidAlignmentError,
)

_WARNING_PATTERN = re.compile(r"\b(warn(?:ing)?|caution)\b", re.IGNORECASE)


@dataclass(slots=True)
class EngineVersionInfo:
    engine_name: str
    executable: str
    command: list[str]
    text: str


@dataclass(slots=True)
class EngineRunReport:
    engine_name: str
    workflow: str
    executable: str
    working_directory: Path
    version: EngineVersionInfo
    command: list[str]
    exit_code: int
    stdout_path: Path
    stderr_path: Path
    output_paths: dict[str, Path]
    warning_lines: list[str]
    started_at_utc: str
    ended_at_utc: str
    runtime_seconds: float
    timeout_seconds: float | None
    timed_out: bool


@dataclass(slots=True)
class EngineResumeCheck:
    resume_allowed: bool
    reason: str


@dataclass(slots=True)
class EngineOutputObservation:
    output_name: str
    path: Path
    exists: bool
    path_kind: str
    size_bytes: int | None
    sha256: str | None


@dataclass(slots=True)
class EngineIncompleteRunRecord:
    engine_name: str
    workflow: str
    executable: str
    working_directory: Path
    manifest_path: Path
    command: list[str]
    stdout_path: Path
    stderr_path: Path
    output_paths: dict[str, Path]
    started_at_utc: str
    ended_at_utc: str | None
    timeout_seconds: float | None
    timed_out: bool
    exit_code: int | None
    failure_reason: str | None
    failure_message: str | None
    missing_output_names: list[str] = field(default_factory=list)
    observed_outputs: list[EngineOutputObservation] = field(default_factory=list)


@dataclass(slots=True)
class EngineActiveRunRecord:
    engine_name: str
    workflow: str
    executable: str
    working_directory: Path
    manifest_path: Path
    command: list[str]
    process_id: int
    started_at_utc: str


def build_engine_output_error(
    message: str,
    *,
    code: str,
    engine_name: str,
    workflow: str,
    path: Path | None = None,
    output_name: str | None = None,
    artifact_kind: str | None = None,
    details: dict[str, Any] | None = None,
) -> EngineWorkflowError:
    """Build one structured engine-output validation error."""
    payload: dict[str, Any] = {
        "engine_name": engine_name,
        "workflow": workflow,
    }
    if path is not None:
        payload["path"] = str(path)
    if output_name is not None:
        payload["output_name"] = output_name
    if artifact_kind is not None:
        payload["artifact_kind"] = artifact_kind
    if details is not None:
        payload.update(details)
    return EngineWorkflowError(message, code=code, details=payload)


def resolve_engine_executable(executable: str | Path) -> str:
    """Resolve an engine executable from PATH or an explicit filesystem path."""
    candidate = Path(executable)
    if candidate.is_absolute() or candidate.parent != Path("."):
        if not candidate.exists():
            raise EngineUnavailableError(
                f"engine executable was not found: {candidate}"
            )
        if not candidate.is_file():
            raise EngineUnavailableError(
                f"engine executable is not a file: {candidate}"
            )
        if not shutil.which(str(candidate)):
            raise EngineUnavailableError(
                f"engine executable is not runnable: {candidate}"
            )
        return str(candidate)

    resolved = shutil.which(str(executable))
    if resolved is None:
        raise EngineUnavailableError(
            f"engine executable is not available on PATH: {executable}"
        )
    return resolved


def summarize_engine_warnings(*, stdout_text: str, stderr_text: str) -> list[str]:
    """Extract warning-like log lines from engine stdout and stderr."""
    warnings: list[str] = []
    for text in (stderr_text, stdout_text):
        for line in text.splitlines():
            message = line.strip()
            if message and _WARNING_PATTERN.search(message) and message not in warnings:
                warnings.append(message)
    return warnings


def validate_timeout_seconds(timeout_seconds: float | None) -> float | None:
    """Validate one optional engine timeout setting."""
    if timeout_seconds is not None and timeout_seconds <= 0:
        raise ValueError(
            f"timeout_seconds must be positive when provided, got {timeout_seconds}"
        )
    return timeout_seconds


def utc_now_text() -> str:
    """Return one stable UTC timestamp for engine manifests and ledgers."""
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def read_engine_version(
    engine_name: str,
    executable: str | Path,
    *,
    version_args: tuple[str, ...],
    timeout_seconds: float | None = None,
) -> EngineVersionInfo:
    """Run the engine-specific version command and capture its text."""
    validate_timeout_seconds(timeout_seconds)
    resolved = resolve_engine_executable(executable)
    command = [resolved, *version_args]
    try:
        result = subprocess.run(  # nosec B603
            command,
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as error:
        budget = "unspecified" if timeout_seconds is None else f"{timeout_seconds:.3f}"
        raise EngineWorkflowError(
            f"{engine_name} version command timed out after {budget} seconds: {' '.join(command)}"
        ) from error
    fragments = [result.stdout.strip(), result.stderr.strip()]
    version_text = "\n".join(fragment for fragment in fragments if fragment).strip()
    if not version_text:
        raise EngineWorkflowError(
            f"{engine_name} did not produce version text for command: {' '.join(command)}"
        )
    return EngineVersionInfo(
        engine_name=engine_name,
        executable=resolved,
        command=command,
        text=version_text,
    )


def engine_incomplete_marker_path(manifest_path: Path) -> Path:
    """Return the sidecar path used to mark one incomplete engine run."""
    return manifest_path.with_suffix(".incomplete.json")


def engine_active_marker_path(manifest_path: Path) -> Path:
    """Return the sidecar path used to mark one actively running engine workflow."""
    return manifest_path.with_suffix(".running.json")


def _observe_engine_output(
    output_name: str,
    path: Path,
) -> EngineOutputObservation:
    exists = path.exists()
    if not exists:
        return EngineOutputObservation(
            output_name=output_name,
            path=path,
            exists=False,
            path_kind="missing",
            size_bytes=None,
            sha256=None,
        )
    if path.is_file():
        return EngineOutputObservation(
            output_name=output_name,
            path=path,
            exists=True,
            path_kind="file",
            size_bytes=path.stat().st_size,
            sha256=file_sha256(path),
        )
    if path.is_dir():
        return EngineOutputObservation(
            output_name=output_name,
            path=path,
            exists=True,
            path_kind="directory",
            size_bytes=None,
            sha256=None,
        )
    return EngineOutputObservation(
        output_name=output_name,
        path=path,
        exists=True,
        path_kind="other",
        size_bytes=None,
        sha256=None,
    )


def observe_engine_outputs(
    output_paths: dict[str, Path],
) -> list[EngineOutputObservation]:
    """Capture the current on-disk state of declared engine outputs."""
    return [
        _observe_engine_output(output_name, output_paths[output_name])
        for output_name in sorted(output_paths)
    ]


def restore_active_engine_run(
    payload: dict[str, Any],
) -> EngineActiveRunRecord:
    """Restore one active-run marker payload into a typed record."""
    return EngineActiveRunRecord(
        engine_name=str(payload["engine_name"]),
        workflow=str(payload["workflow"]),
        executable=str(payload["executable"]),
        working_directory=Path(payload["working_directory"]),
        manifest_path=Path(payload["manifest_path"]),
        command=[str(item) for item in payload["command"]],
        process_id=int(payload["process_id"]),
        started_at_utc=str(payload["started_at_utc"]),
    )


def load_active_engine_run(
    manifest_path: Path,
) -> EngineActiveRunRecord | None:
    """Load one active-run marker when it exists."""
    marker_path = engine_active_marker_path(manifest_path)
    if not marker_path.exists():
        return None
    return restore_active_engine_run(load_engine_manifest(marker_path))


def write_active_engine_run(record: EngineActiveRunRecord) -> Path:
    """Persist one active-run marker."""
    return write_engine_manifest(
        engine_active_marker_path(record.manifest_path),
        record,
    )


def clear_active_engine_run(manifest_path: Path) -> None:
    """Remove one active-run marker when a workflow is no longer executing."""
    engine_active_marker_path(manifest_path).unlink(missing_ok=True)


def _process_is_alive(process_id: int) -> bool:
    if process_id <= 0:
        return False
    try:
        os.kill(process_id, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False
    return True


def active_engine_run_is_live(record: EngineActiveRunRecord) -> bool:
    """Report whether an active-run marker still belongs to a live process."""
    return _process_is_alive(record.process_id)


def acquire_active_engine_run(record: EngineActiveRunRecord) -> Path:
    """Claim exclusive execution ownership for one engine-workflow manifest."""
    marker_path = engine_active_marker_path(record.manifest_path)
    if marker_path.exists():
        try:
            active = load_active_engine_run(record.manifest_path)
        except Exception as error:  # pragma: no cover - defensive malformed-marker path
            raise EngineWorkflowError(
                "engine workflow could not verify an existing active-run marker",
                code="engine_workflow_running_marker_invalid",
                details={
                    "manifest_path": str(record.manifest_path),
                    "marker_path": str(marker_path),
                },
            ) from error
        if active is not None and active_engine_run_is_live(active):
            raise EngineWorkflowError(
                "engine workflow is already running for the requested output manifest",
                code="engine_workflow_already_running",
                details={
                    "manifest_path": str(record.manifest_path),
                    "marker_path": str(marker_path),
                    "running_process_id": active.process_id,
                    "running_workflow": active.workflow,
                    "running_engine_name": active.engine_name,
                },
            )
        marker_path.unlink(missing_ok=True)
    payload_text = (
        json.dumps(asdict(record), default=str, indent=2, sort_keys=True) + "\n"
    )
    try:
        marker_fd = os.open(marker_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError as error:
        raise EngineWorkflowError(
            "engine workflow is already running for the requested output manifest",
            code="engine_workflow_already_running",
            details={
                "manifest_path": str(record.manifest_path),
                "marker_path": str(marker_path),
            },
        ) from error
    try:
        with os.fdopen(marker_fd, "w", encoding="utf-8") as handle:
            handle.write(payload_text)
    except Exception:
        marker_path.unlink(missing_ok=True)
        raise
    return marker_path


def restore_incomplete_engine_run(
    payload: dict[str, Any],
) -> EngineIncompleteRunRecord:
    """Restore one incomplete-run marker payload into a typed record."""
    return EngineIncompleteRunRecord(
        engine_name=str(payload["engine_name"]),
        workflow=str(payload["workflow"]),
        executable=str(payload["executable"]),
        working_directory=Path(payload["working_directory"]),
        manifest_path=Path(payload["manifest_path"]),
        command=[str(item) for item in payload["command"]],
        stdout_path=Path(payload["stdout_path"]),
        stderr_path=Path(payload["stderr_path"]),
        output_paths={
            str(key): Path(value)
            for key, value in dict(payload["output_paths"]).items()
        },
        started_at_utc=str(payload["started_at_utc"]),
        ended_at_utc=(
            None
            if payload.get("ended_at_utc") is None
            else str(payload["ended_at_utc"])
        ),
        timeout_seconds=(
            None
            if payload.get("timeout_seconds") is None
            else float(payload["timeout_seconds"])
        ),
        timed_out=bool(payload.get("timed_out", False)),
        exit_code=(
            None if payload.get("exit_code") is None else int(payload["exit_code"])
        ),
        failure_reason=(
            None
            if payload.get("failure_reason") is None
            else str(payload["failure_reason"])
        ),
        failure_message=(
            None
            if payload.get("failure_message") is None
            else str(payload["failure_message"])
        ),
        missing_output_names=[
            str(item) for item in payload.get("missing_output_names", [])
        ],
        observed_outputs=[
            EngineOutputObservation(
                output_name=str(dict(item)["output_name"]),
                path=Path(dict(item)["path"]),
                exists=bool(dict(item)["exists"]),
                path_kind=str(dict(item)["path_kind"]),
                size_bytes=(
                    None
                    if dict(item).get("size_bytes") is None
                    else int(dict(item)["size_bytes"])
                ),
                sha256=(
                    None
                    if dict(item).get("sha256") is None
                    else str(dict(item)["sha256"])
                ),
            )
            for item in payload.get("observed_outputs", [])
        ],
    )


def load_incomplete_engine_run(
    manifest_path: Path,
) -> EngineIncompleteRunRecord | None:
    """Load one incomplete-run marker when it exists."""
    marker_path = engine_incomplete_marker_path(manifest_path)
    if not marker_path.exists():
        return None
    return restore_incomplete_engine_run(load_engine_manifest(marker_path))


def write_incomplete_engine_run(record: EngineIncompleteRunRecord) -> Path:
    """Persist one incomplete-run marker."""
    return write_engine_manifest(
        engine_incomplete_marker_path(record.manifest_path),
        record,
    )


def clear_incomplete_engine_run(manifest_path: Path) -> None:
    """Remove one incomplete-run marker when a workflow completes cleanly."""
    engine_incomplete_marker_path(manifest_path).unlink(missing_ok=True)


def cleanup_incomplete_engine_run(
    manifest_path: Path,
) -> EngineIncompleteRunRecord | None:
    """Delete one incomplete-run marker and the outputs it recorded."""
    record = load_incomplete_engine_run(manifest_path)
    if record is None:
        return None
    for path in {
        record.manifest_path,
        record.stdout_path,
        record.stderr_path,
        *record.output_paths.values(),
    }:
        if path.exists():
            path.unlink()
    clear_incomplete_engine_run(manifest_path)
    return record


def update_incomplete_engine_run(
    manifest_path: Path,
    *,
    ended_at_utc: str | None = None,
    timed_out: bool | None = None,
    exit_code: int | None = None,
    failure_reason: str | None = None,
    failure_message: str | None = None,
    missing_output_names: list[str] | None = None,
    observed_outputs: list[EngineOutputObservation] | None = None,
) -> EngineIncompleteRunRecord | None:
    """Refresh one incomplete-run marker with the latest failure details."""
    record = load_incomplete_engine_run(manifest_path)
    if record is None:
        return None
    if ended_at_utc is not None:
        record.ended_at_utc = ended_at_utc
    if timed_out is not None:
        record.timed_out = timed_out
    if exit_code is not None:
        record.exit_code = exit_code
    if failure_reason is not None:
        record.failure_reason = failure_reason
    if failure_message is not None:
        record.failure_message = failure_message
    if missing_output_names is not None:
        record.missing_output_names = list(missing_output_names)
    if observed_outputs is not None:
        record.observed_outputs = list(observed_outputs)
    write_incomplete_engine_run(record)
    return record


def execute_engine_command(
    *,
    engine_name: str,
    workflow: str,
    executable: str,
    version: EngineVersionInfo,
    command_args: list[str],
    work_dir: Path,
    stdout_path: Path,
    stderr_path: Path,
    output_paths: dict[str, Path],
    manifest_path: Path,
    timeout_seconds: float | None = None,
) -> EngineRunReport:
    """Execute an engine command, capture logs, and validate expected outputs."""
    validate_timeout_seconds(timeout_seconds)
    work_dir.mkdir(parents=True, exist_ok=True)
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    command = [executable, *command_args]
    started_at_utc = utc_now_text()
    started = time.perf_counter()
    incomplete_record = EngineIncompleteRunRecord(
        engine_name=engine_name,
        workflow=workflow,
        executable=executable,
        working_directory=work_dir,
        manifest_path=manifest_path,
        command=command,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        output_paths=output_paths,
        started_at_utc=started_at_utc,
        ended_at_utc=None,
        timeout_seconds=timeout_seconds,
        timed_out=False,
        exit_code=None,
        failure_reason=None,
        failure_message=None,
        missing_output_names=[],
        observed_outputs=[],
    )
    active_record = EngineActiveRunRecord(
        engine_name=engine_name,
        workflow=workflow,
        executable=executable,
        working_directory=work_dir,
        manifest_path=manifest_path,
        command=command,
        process_id=os.getpid(),
        started_at_utc=started_at_utc,
    )
    acquire_active_engine_run(active_record)
    write_incomplete_engine_run(incomplete_record)
    try:
        with (
            stdout_path.open("w", encoding="utf-8") as stdout_handle,
            stderr_path.open("w", encoding="utf-8") as stderr_handle,
        ):
            try:
                result = subprocess.run(  # nosec B603
                    command,
                    cwd=work_dir,
                    stdout=stdout_handle,
                    stderr=stderr_handle,
                    text=True,
                    check=False,
                    timeout=timeout_seconds,
                )
            except subprocess.TimeoutExpired as error:
                stdout_text = (
                    stdout_path.read_text(encoding="utf-8")
                    if stdout_path.exists()
                    else ""
                )
                stderr_text = (
                    stderr_path.read_text(encoding="utf-8")
                    if stderr_path.exists()
                    else ""
                )
                budget = (
                    "unspecified"
                    if timeout_seconds is None
                    else f"{timeout_seconds:.3f}"
                )
                incomplete_record.ended_at_utc = utc_now_text()
                incomplete_record.timed_out = True
                incomplete_record.failure_reason = "engine_command_timeout"
                incomplete_record.failure_message = (
                    f"{engine_name} {workflow} timed out after {budget} seconds"
                )
                incomplete_record.observed_outputs = observe_engine_outputs(
                    output_paths
                )
                incomplete_record.missing_output_names = [
                    observation.output_name
                    for observation in incomplete_record.observed_outputs
                    if not observation.exists
                ]
                write_incomplete_engine_run(incomplete_record)
                warning_lines = summarize_engine_warnings(
                    stdout_text=stdout_text,
                    stderr_text=stderr_text,
                )
                raise EngineWorkflowError(
                    f"{engine_name} {workflow} timed out after {budget} seconds; "
                    f"stderr log: {stderr_path}; warnings: {len(warning_lines)}"
                ) from error
        stdout_text = (
            stdout_path.read_text(encoding="utf-8") if stdout_path.exists() else ""
        )
        stderr_text = (
            stderr_path.read_text(encoding="utf-8") if stderr_path.exists() else ""
        )
        ended_at_utc = utc_now_text()
        warning_lines = summarize_engine_warnings(
            stdout_text=stdout_text, stderr_text=stderr_text
        )
        if result.returncode != 0:
            incomplete_record.ended_at_utc = ended_at_utc
            incomplete_record.exit_code = result.returncode
            incomplete_record.failure_reason = "engine_command_exit_nonzero"
            incomplete_record.failure_message = (
                f"{engine_name} {workflow} failed with exit code {result.returncode}"
            )
            incomplete_record.observed_outputs = observe_engine_outputs(output_paths)
            incomplete_record.missing_output_names = [
                observation.output_name
                for observation in incomplete_record.observed_outputs
                if not observation.exists
            ]
            write_incomplete_engine_run(incomplete_record)
            raise EngineWorkflowError(
                f"{engine_name} {workflow} failed with exit code {result.returncode}; stderr log: {stderr_path}"
            )
        missing_outputs = [
            {"output_name": output_name, "path": str(path)}
            for output_name, path in output_paths.items()
            if not path.exists()
        ]
        if missing_outputs:
            incomplete_record.ended_at_utc = ended_at_utc
            incomplete_record.exit_code = result.returncode
            incomplete_record.failure_reason = "engine_required_output_missing"
            incomplete_record.failure_message = (
                f"{engine_name} {workflow} did not produce expected outputs"
            )
            incomplete_record.observed_outputs = observe_engine_outputs(output_paths)
            incomplete_record.missing_output_names = [
                str(entry["output_name"]) for entry in missing_outputs
            ]
            write_incomplete_engine_run(incomplete_record)
            missing_paths = ", ".join(
                entry["path"]
                for entry in missing_outputs
                if isinstance(entry["path"], str)
            )
            raise build_engine_output_error(
                f"{engine_name} {workflow} did not produce expected outputs: {missing_paths}",
                code="engine_required_output_missing",
                engine_name=engine_name,
                workflow=workflow,
                details={
                    "missing_outputs": missing_outputs,
                    "declared_output_count": len(output_paths),
                },
            )
        return EngineRunReport(
            engine_name=engine_name,
            workflow=workflow,
            executable=executable,
            working_directory=work_dir,
            version=version,
            command=command,
            exit_code=result.returncode,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            output_paths=output_paths,
            warning_lines=warning_lines,
            started_at_utc=started_at_utc,
            ended_at_utc=ended_at_utc,
            runtime_seconds=max(0.0, round(time.perf_counter() - started, 6)),
            timeout_seconds=timeout_seconds,
            timed_out=False,
        )
    finally:
        clear_active_engine_run(manifest_path)


def load_unaligned_fasta(path: Path) -> list[tuple[str, str]]:
    """Load a FASTA file while allowing unequal sequence lengths."""
    records: list[tuple[str, str]] = []
    current_identifier: str | None = None
    current_sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_identifier is not None:
                    records.append((current_identifier, "".join(current_sequence)))
                current_identifier = line[1:].strip()
                current_sequence = []
                continue
            if current_identifier is None:
                raise InvalidAlignmentError(
                    f"alignment sequence appears before any FASTA header in {path}"
                )
            current_sequence.append(line)
    if current_identifier is not None:
        records.append((current_identifier, "".join(current_sequence)))
    if not records:
        raise InvalidAlignmentError(f"alignment contains no FASTA records: {path}")

    ids = [identifier for identifier, _sequence in records]
    duplicate_ids = sorted(
        identifier for identifier in set(ids) if ids.count(identifier) > 1
    )
    if duplicate_ids:
        raise InvalidAlignmentError(
            f"alignment contains duplicate sequence ids: {', '.join(duplicate_ids)}"
        )
    return records


def file_sha256(path: Path) -> str:
    """Compute a deterministic SHA-256 checksum for one file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_file_checksums(paths: list[Path]) -> dict[str, str]:
    """Hash every existing file path in a stable key order."""
    return {
        str(path): file_sha256(path)
        for path in sorted(paths, key=str)
        if path.exists() and path.is_file()
    }


def write_engine_manifest(path: Path, payload: Any) -> Path:
    """Write a deterministic engine-workflow manifest."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(payload), default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def load_engine_manifest(path: Path) -> dict[str, Any]:
    """Load a previously written engine-workflow manifest."""
    return json.loads(path.read_text(encoding="utf-8"))
