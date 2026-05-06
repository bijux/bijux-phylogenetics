from __future__ import annotations

from configparser import ConfigParser
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]


def _tox_config() -> ConfigParser:
    parser = ConfigParser()
    parser.read(REPO_ROOT / "tox.ini", encoding="utf-8")
    return parser


def _envlist() -> set[str]:
    envlist = _tox_config()["tox"]["envlist"]
    return {line.strip() for line in envlist.splitlines() if line.strip()}


def test_root_tox_keeps_the_shared_env_families_and_drops_proteomics_only_ones() -> None:
    envlist = _envlist()

    assert "security" in envlist
    assert "docs" in envlist
    assert "fmt-{dev,core}" not in envlist
    assert "api-freeze-core" not in envlist
    assert "openapi-drift-core" not in envlist


def test_root_make_declares_shared_maintainer_commands() -> None:
    root_make = (REPO_ROOT / "makes" / "root.mk").read_text(encoding="utf-8")

    assert "check:" in root_make
    assert "sync-badges:" in root_make
    assert "check-badges:" in root_make
