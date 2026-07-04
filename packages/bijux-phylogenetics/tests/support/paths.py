from __future__ import annotations

from pathlib import Path

TESTS_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = TESTS_ROOT.parent
PACKAGES_ROOT = PACKAGE_ROOT.parent
REPOSITORY_ROOT = PACKAGES_ROOT.parent
FIXTURES_ROOT = TESTS_ROOT / "fixtures"
