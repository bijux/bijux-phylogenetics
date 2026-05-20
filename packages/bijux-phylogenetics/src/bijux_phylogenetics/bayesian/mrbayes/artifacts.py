from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.runtime.errors import EngineWorkflowError


def _mrbayes_artifact_error(
    message: str,
    *,
    code: str,
    path: Path,
    artifact_kind: str,
    details: dict[str, object] | None = None,
) -> EngineWorkflowError:
    payload: dict[str, object] = {
        "path": str(path),
        "artifact_kind": artifact_kind,
    }
    if details is not None:
        payload.update(details)
    return EngineWorkflowError(message, code=code, details=payload)
