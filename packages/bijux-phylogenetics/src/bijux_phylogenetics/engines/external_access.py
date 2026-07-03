from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from contextvars import ContextVar
import os
from pathlib import Path

from bijux_phylogenetics.runtime.errors import (
    ExternalEngineAccessError,
    PhylogeneticsError,
)

NATIVE_MAXIMUM_LIKELIHOOD_ONLY_ENV_VAR = (
    "BIJUX_PHYLOGENETICS_NATIVE_MAXIMUM_LIKELIHOOD_ONLY"
)
_NATIVE_MAXIMUM_LIKELIHOOD_ONLY_OVERRIDE: ContextVar[bool | None] = ContextVar(
    "native_maximum_likelihood_only_override",
    default=None,
)
_TRUTHY_VALUES = frozenset({"1", "true", "yes", "on"})
_FALSEY_VALUES = frozenset({"0", "false", "no", "off", ""})

__all__ = [
    "NATIVE_MAXIMUM_LIKELIHOOD_ONLY_ENV_VAR",
    "is_native_maximum_likelihood_only_mode_enabled",
    "native_maximum_likelihood_only_mode",
    "require_external_engine_access",
]


def _read_native_maximum_likelihood_only_env() -> bool:
    raw_value = os.environ.get(NATIVE_MAXIMUM_LIKELIHOOD_ONLY_ENV_VAR)
    if raw_value is None:
        return False
    normalized_value = raw_value.strip().lower()
    if normalized_value in _TRUTHY_VALUES:
        return True
    if normalized_value in _FALSEY_VALUES:
        return False
    raise PhylogeneticsError(
        f"{NATIVE_MAXIMUM_LIKELIHOOD_ONLY_ENV_VAR} must be one of "
        f"{', '.join(sorted(_TRUTHY_VALUES | _FALSEY_VALUES))}, got {raw_value!r}",
        code="native_maximum_likelihood_only_invalid",
        details={
            "environment_variable": NATIVE_MAXIMUM_LIKELIHOOD_ONLY_ENV_VAR,
            "received_value": raw_value,
            "allowed_values": sorted(_TRUTHY_VALUES | _FALSEY_VALUES),
        },
    )


def is_native_maximum_likelihood_only_mode_enabled() -> bool:
    """Return whether native-only mode currently forbids external-engine access."""
    override = _NATIVE_MAXIMUM_LIKELIHOOD_ONLY_OVERRIDE.get()
    if override is not None:
        return override
    return _read_native_maximum_likelihood_only_env()


@contextmanager
def native_maximum_likelihood_only_mode() -> Iterator[None]:
    """Temporarily forbid external-engine access within the active execution context."""
    token = _NATIVE_MAXIMUM_LIKELIHOOD_ONLY_OVERRIDE.set(True)
    try:
        yield
    finally:
        _NATIVE_MAXIMUM_LIKELIHOOD_ONLY_OVERRIDE.reset(token)


def require_external_engine_access(
    *,
    surface: str,
    engine_name: str | None = None,
    executable: str | Path | None = None,
) -> None:
    """Require the current runtime policy to permit one external-engine access attempt."""
    if not is_native_maximum_likelihood_only_mode_enabled():
        return
    details: dict[str, object] = {
        "surface": surface,
        "environment_variable": NATIVE_MAXIMUM_LIKELIHOOD_ONLY_ENV_VAR,
        "policy": "native_maximum_likelihood_only",
    }
    if engine_name is not None:
        details["engine_name"] = engine_name
    if executable is not None:
        details["executable"] = str(executable)
    engine_fragment = "" if engine_name is None else f" for {engine_name}"
    raise ExternalEngineAccessError(
        "native maximum-likelihood only mode forbids external engine access"
        f"{engine_fragment}: {surface}",
        details=details,
    )
