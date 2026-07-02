from __future__ import annotations

import os
from pathlib import Path

from bijux_phylogenetics.runtime.errors import (
    ExternalEngineAccessError,
    PhylogeneticsError,
)

RUN_LIVE_PARITY_ENV_VAR = "RUN_LIVE_PARITY"
_TRUTHY_VALUES = frozenset({"1", "true", "yes", "on"})
_FALSEY_VALUES = frozenset({"0", "false", "no", "off", ""})


def _read_live_parity_env() -> bool:
    raw_value = os.environ.get(RUN_LIVE_PARITY_ENV_VAR)
    if raw_value is None:
        return False
    normalized_value = raw_value.strip().lower()
    if normalized_value in _TRUTHY_VALUES:
        return True
    if normalized_value in _FALSEY_VALUES:
        return False
    raise PhylogeneticsError(
        f"{RUN_LIVE_PARITY_ENV_VAR} must be one of "
        f"{', '.join(sorted(_TRUTHY_VALUES | _FALSEY_VALUES))}, got {raw_value!r}",
        code="live_parity_opt_in_invalid",
        details={
            "environment_variable": RUN_LIVE_PARITY_ENV_VAR,
            "received_value": raw_value,
            "allowed_values": sorted(_TRUTHY_VALUES | _FALSEY_VALUES),
        },
    )


def is_live_parity_enabled() -> bool:
    """Return whether live reference parity is explicitly enabled."""
    return _read_live_parity_env()


def require_live_parity_access(
    *,
    surface: str,
    reference_source: str | None = None,
    executable: str | Path | None = None,
) -> None:
    """Require explicit opt-in before launching live reference parity."""
    if is_live_parity_enabled():
        return
    details: dict[str, object] = {
        "surface": surface,
        "environment_variable": RUN_LIVE_PARITY_ENV_VAR,
        "policy": "live_parity_opt_in",
    }
    if reference_source is not None:
        details["reference_source"] = reference_source
    if executable is not None:
        details["executable"] = str(executable)
    reference_fragment = (
        ""
        if reference_source is None
        else f" for {reference_source}"
    )
    raise ExternalEngineAccessError(
        "live reference parity requires explicit opt-in"
        f"{reference_fragment}: set {RUN_LIVE_PARITY_ENV_VAR}=1 to continue",
        details=details,
    )
