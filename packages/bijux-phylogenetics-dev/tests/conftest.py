"""Test-path bootstrap for repository-owned maintainer package tests."""

from __future__ import annotations

from configparser import ConfigParser
from pathlib import Path
import sys

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
PYTEST_CONFIG_PATH = REPO_ROOT / "configs" / "pytest.ini"

for source_root in (
    REPO_ROOT / "packages" / "bijux-phylogenetics-dev" / "src",
    REPO_ROOT / "packages" / "bijux-phylogenetics" / "src",
    REPO_ROOT / "packages" / "phylogenetic" / "src",
):
    source_root_text = str(source_root)
    if source_root_text not in sys.path:
        sys.path.insert(0, source_root_text)


def pytest_configure(config: pytest.Config) -> None:
    """Register repository-owned shared pytest markers for maintainer tests."""
    parser = ConfigParser()
    parser.read(PYTEST_CONFIG_PATH, encoding="utf-8")
    configured_markers = {
        marker for marker in parser["pytest"]["markers"].splitlines() if marker.strip()
    }

    existing_markers = set(config.getini("markers"))
    for marker in sorted(configured_markers - existing_markers):
        config.addinivalue_line("markers", marker)
