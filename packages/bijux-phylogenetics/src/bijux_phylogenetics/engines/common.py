from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess  # nosec B404
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


@dataclass(slots=True)
class EngineResumeCheck:
    resume_allowed: bool
    reason: str


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


def read_engine_version(
    engine_name: str,
    executable: str | Path,
    *,
    version_args: tuple[str, ...],
) -> EngineVersionInfo:
    """Run the engine-specific version command and capture its text."""
    resolved = resolve_engine_executable(executable)
    command = [resolved, *version_args]
    result = subprocess.run(command, capture_output=True, text=True, check=False)  # nosec B603
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
) -> EngineRunReport:
    """Execute an engine command, capture logs, and validate expected outputs."""
    work_dir.mkdir(parents=True, exist_ok=True)
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    command = [executable, *command_args]
    with (
        stdout_path.open("w", encoding="utf-8") as stdout_handle,
        stderr_path.open("w", encoding="utf-8") as stderr_handle,
    ):
        result = subprocess.run(  # nosec B603
            command,
            cwd=work_dir,
            stdout=stdout_handle,
            stderr=stderr_handle,
            text=True,
            check=False,
        )
    stdout_text = (
        stdout_path.read_text(encoding="utf-8") if stdout_path.exists() else ""
    )
    stderr_text = (
        stderr_path.read_text(encoding="utf-8") if stderr_path.exists() else ""
    )
    warning_lines = summarize_engine_warnings(
        stdout_text=stdout_text, stderr_text=stderr_text
    )
    if result.returncode != 0:
        raise EngineWorkflowError(
            f"{engine_name} {workflow} failed with exit code {result.returncode}; stderr log: {stderr_path}"
        )
    missing_outputs = [str(path) for path in output_paths.values() if not path.exists()]
    if missing_outputs:
        raise EngineWorkflowError(
            f"{engine_name} {workflow} did not produce expected outputs: {', '.join(missing_outputs)}"
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
