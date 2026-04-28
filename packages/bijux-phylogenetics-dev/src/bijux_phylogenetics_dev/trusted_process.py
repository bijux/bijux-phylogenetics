"""Run repository-owned commands through a small validated subprocess wrapper."""

from __future__ import annotations

from collections.abc import Sequence
import os
from pathlib import Path
import subprocess  # nosec B404

CommandArg = str | os.PathLike[str]


class TrustedCommandError(RuntimeError):
    """Raised when a trusted command exits with a non-zero status."""

    def __init__(
        self, command: Sequence[CommandArg], error: subprocess.CalledProcessError
    ):
        """Capture the failed trusted command with its subprocess details."""
        self.command = [os.fspath(part) for part in command]
        self.returncode = error.returncode
        self.stdout = error.stdout
        self.stderr = error.stderr
        super().__init__(f"trusted command failed: {' '.join(self.command)}")


def _normalize_command(command: Sequence[CommandArg]) -> list[str]:
    """Normalize command."""
    if not command:
        raise ValueError("trusted command is empty")
    executable = Path(os.fspath(command[0])).expanduser()
    if not executable.is_absolute():
        raise ValueError(
            f"trusted command requires an absolute executable path: {command[0]!r}"
        )
    return [os.fspath(executable), *(os.fspath(part) for part in command[1:])]


def run_text(
    command: Sequence[CommandArg],
    *,
    cwd: Path | None = None,
    check: bool,
    capture_output: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run a trusted command and return its completed process record."""
    normalized = _normalize_command(command)
    try:
        return subprocess.run(  # nosec B603
            normalized,
            cwd=cwd,
            check=check,
            capture_output=capture_output,
            text=True,
        )
    except subprocess.CalledProcessError as error:
        raise TrustedCommandError(normalized, error) from error
