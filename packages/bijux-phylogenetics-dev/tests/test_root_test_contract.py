from __future__ import annotations

from configparser import ConfigParser
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]


def _pytest_config() -> ConfigParser:
    parser = ConfigParser()
    parser.read(REPO_ROOT / "configs" / "pytest.ini", encoding="utf-8")
    return parser


def _pytest_testpaths() -> set[str]:
    raw_value = _pytest_config()["pytest"]["testpaths"]
    return {line.strip() for line in raw_value.splitlines() if line.strip()}


def test_root_pytest_configuration_covers_all_package_test_directories() -> None:
    testpaths = _pytest_testpaths()
    expected_paths = {
        test_dir.relative_to(REPO_ROOT).as_posix()
        for test_dir in (REPO_ROOT / "packages").glob("*/tests")
    }
    assert testpaths == expected_paths
