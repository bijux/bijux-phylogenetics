from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from bijux_phylogenetics import __version__
from bijux_phylogenetics.runtime.error_explanations import (
    explain_phylogenetics_error,
    explanation_payload,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError
from bijux_phylogenetics.runtime.identity import CLI_NAME, IMPORT_NAME, PACKAGE_NAME


@dataclass(slots=True)
class CommandResult:
    """Structured envelope for CLI command responses."""

    status: str
    command: str
    inputs: list[str]
    outputs: list[str]
    warnings: list[str] = field(default_factory=list)
    errors: list[dict[str, Any]] = field(default_factory=list)
    metrics: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)
    data: Any | None = None


def build_command_result(
    *,
    command: str,
    inputs: list[Path | str],
    outputs: list[Path | str] | None = None,
    warnings: list[str] | None = None,
    metrics: dict[str, Any] | None = None,
    data: Any | None = None,
) -> CommandResult:
    """Build a successful CLI result envelope."""
    return CommandResult(
        status="ok",
        command=command,
        inputs=[str(item) for item in inputs],
        outputs=[str(item) for item in outputs or []],
        warnings=list(warnings or []),
        metrics=dict(metrics or {}),
        provenance={
            "cli": CLI_NAME,
            "package": PACKAGE_NAME,
            "import_name": IMPORT_NAME,
            "version": __version__,
        },
        data=data,
    )


def build_error_result(
    *,
    command: str,
    inputs: list[Path | str],
    error: PhylogeneticsError,
) -> CommandResult:
    """Build a failed CLI result envelope."""
    explanation = explain_phylogenetics_error(
        error,
        inputs=[Path(item) if not isinstance(item, Path) else item for item in inputs],
    )
    details = dict(error.details)
    details.update(
        {
            key: value
            for key, value in explanation_payload(explanation).items()
            if key not in details
        }
    )
    return CommandResult(
        status="error",
        command=command,
        inputs=[str(item) for item in inputs],
        outputs=[],
        warnings=[],
        errors=[
            {
                "code": error.code,
                "message": error.message,
                **({"details": details} if details else {}),
            }
        ],
        metrics={},
        provenance={
            "cli": CLI_NAME,
            "package": PACKAGE_NAME,
            "import_name": IMPORT_NAME,
            "version": __version__,
        },
        data=None,
    )
