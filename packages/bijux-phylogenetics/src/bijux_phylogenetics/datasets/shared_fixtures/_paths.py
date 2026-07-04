from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.runtime.repository_layout import (
    governed_fixture_root as locate_governed_fixture_root,
    package_root as locate_package_root,
)


def package_root() -> Path:
    return locate_package_root(start=Path(__file__))


def governed_fixture_root() -> Path:
    return locate_governed_fixture_root(start=Path(__file__))
