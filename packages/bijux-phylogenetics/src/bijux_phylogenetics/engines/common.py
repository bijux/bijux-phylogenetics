from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import UTC, datetime
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess  # nosec B404
import time
from typing import Any

from bijux_phylogenetics.errors import (
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
    failure_message: str | None


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
        failure_message=(
            None
            if payload.get("failure_message") is None
            else str(payload["failure_message"])
        ),
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
    failure_message: str | None = None,
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
    if failure_message is not None:
        record.failure_message = failure_message
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
        failure_message=None,
    )
    write_incomplete_engine_run(incomplete_record)
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
                stdout_path.read_text(encoding="utf-8") if stdout_path.exists() else ""
            )
            stderr_text = (
                stderr_path.read_text(encoding="utf-8") if stderr_path.exists() else ""
            )
            budget = (
                "unspecified" if timeout_seconds is None else f"{timeout_seconds:.3f}"
            )
            incomplete_record.ended_at_utc = utc_now_text()
            incomplete_record.timed_out = True
            incomplete_record.failure_message = (
                f"{engine_name} {workflow} timed out after {budget} seconds"
            )
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
        incomplete_record.failure_message = (
            f"{engine_name} {workflow} failed with exit code {result.returncode}"
        )
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
        incomplete_record.failure_message = (
            f"{engine_name} {workflow} did not produce expected outputs"
        )
        write_incomplete_engine_run(incomplete_record)
        missing_paths = ", ".join(
            entry["path"] for entry in missing_outputs if isinstance(entry["path"], str)
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
