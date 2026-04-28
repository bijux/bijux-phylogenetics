from __future__ import annotations

from dataclasses import dataclass
from importlib import metadata
import platform
import sys


@dataclass(slots=True)
class DependencyStatus:
    name: str
    version: str
    available: bool


@dataclass(slots=True)
class EnvironmentInventory:
    python_version: str
    host_platform: str
    dependencies: list[DependencyStatus]


_DEPENDENCIES = (
    "biopython",
    "dendropy",
    "scikit-bio",
    "ete3",
    "ete4",
    "numpy",
    "pandas",
    "bijux-phylogenetics",
)


def _dependency_status(name: str) -> DependencyStatus:
    try:
        return DependencyStatus(name=name, version=metadata.version(name), available=True)
    except metadata.PackageNotFoundError:
        return DependencyStatus(name=name, version="unavailable", available=False)


def inspect_environment() -> EnvironmentInventory:
    """Inspect runtime and optional backend dependency availability."""
    return EnvironmentInventory(
        python_version=sys.version.split()[0],
        host_platform=platform.platform(),
        dependencies=[_dependency_status(name) for name in _DEPENDENCIES],
    )
