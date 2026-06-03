from __future__ import annotations

from pathlib import Path


def package_root() -> Path:
    return Path(__file__).resolve().parents[4]


def governed_fixture_root() -> Path:
    return package_root() / "tests" / "fixtures"
