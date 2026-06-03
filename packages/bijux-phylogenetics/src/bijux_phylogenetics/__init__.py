"""Curated package gateway for Bijux phylogenetics."""

from __future__ import annotations

from importlib import import_module, metadata
from typing import TYPE_CHECKING, Any

_PUBLIC_MODULES = (
    "ancestral",
    "api",
    "bayesian",
    "biogeography",
    "comparative",
    "datasets",
    "distance",
    "evidence",
    "parsimony",
    "parity",
    "phylo",
    "trees",
)

try:
    __version__ = metadata.version("bijux-phylogenetics")
except metadata.PackageNotFoundError:
    __version__ = "0.1.0"

__all__ = [
    "__version__",
    "ancestral",
    "api",
    "bayesian",
    "biogeography",
    "comparative",
    "datasets",
    "distance",
    "evidence",
    "parsimony",
    "parity",
    "phylo",
    "trees",
]

if TYPE_CHECKING:
    from . import (
        ancestral,
        api,
        bayesian,
        biogeography,
        comparative,
        datasets,
        distance,
        evidence,
        parity,
        parsimony,
        phylo,
        trees,
    )


def __getattr__(name: str) -> Any:
    """Load curated top-level modules on demand."""
    if name in _PUBLIC_MODULES:
        module = import_module(f".{name}", __name__)
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """Expose the curated package gateway in interactive discovery."""
    return sorted(set(globals()) | set(__all__))
