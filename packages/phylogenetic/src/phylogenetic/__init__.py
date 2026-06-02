"""Compatibility alias module for bijux-phylogenetics."""

from __future__ import annotations

from importlib import import_module, metadata
from typing import Any

from bijux_phylogenetics.comparative.common import (
    ComparativeDataset as ComparativeDataset,
)

from .runtime_alias import install_runtime_aliases

_ALIAS_PACKAGE = "phylogenetic"
_RUNTIME_PACKAGE = "bijux_phylogenetics"
_LOCAL_SUBMODULES = frozenset({"__main__", "cli", "runtime_alias"})
_runtime_module = import_module(_RUNTIME_PACKAGE)

install_runtime_aliases(
    alias_package=_ALIAS_PACKAGE,
    runtime_package=_RUNTIME_PACKAGE,
    local_submodules=_LOCAL_SUBMODULES,
)

for _name in getattr(_runtime_module, "__all__", ()):
    if _name == "__version__":
        continue
    globals()[_name] = getattr(_runtime_module, _name)

try:
    __version__ = metadata.version(_ALIAS_PACKAGE)
except metadata.PackageNotFoundError:
    __version__ = "0.1.0"

__all__ = list(getattr(_runtime_module, "__all__", ()))
if "ComparativeDataset" not in __all__:
    __all__.append("ComparativeDataset")


def __getattr__(name: str) -> Any:
    """Forward top-level compatibility lookups to the canonical runtime package."""
    return getattr(_runtime_module, name)


def __dir__() -> list[str]:
    """Expose the canonical runtime attributes in interactive discovery."""
    return sorted(set(globals()) | set(dir(_runtime_module)))
