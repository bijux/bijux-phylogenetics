from __future__ import annotations

from configparser import ConfigParser
from pathlib import Path

import pytest

from bijux_phylogenetics.parity.live_access import RUN_LIVE_PARITY_ENV_VAR

PACKAGE_ROOT = Path(__file__).resolve().parent.parent
REPO_ROOT = PACKAGE_ROOT.parents[1]
PYTEST_CONFIG_PATH = REPO_ROOT / "configs" / "pytest.ini"


def pytest_configure(config: object) -> None:
    parser = ConfigParser()
    parser.read(PYTEST_CONFIG_PATH, encoding="utf-8")
    configured_markers = {
        marker for marker in parser["pytest"]["markers"].splitlines() if marker.strip()
    }

    existing_markers = set(config.getini("markers"))
    for marker in sorted(configured_markers - existing_markers):
        config.addinivalue_line("markers", marker)


@pytest.fixture
def enable_live_parity_access(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(RUN_LIVE_PARITY_ENV_VAR, "1")
