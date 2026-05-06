from __future__ import annotations

from pathlib import Path
import tomllib

REPO_ROOT = Path(__file__).resolve().parents[3]
PROTEOMICS_ONLY_EXTENSION_COMMANDS = {
    "ensure-venv:",
    "nlenv:",
    "manage_examples:",
    "manage_models:",
    "api-freeze:",
    "openapi-drift:",
    "architecture-check:",
}


def _root_pyproject() -> dict[str, object]:
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)


def test_root_pyproject_declares_shared_quality_tooling() -> None:
    tool_section = _root_pyproject()["tool"]
    interrogate = tool_section["interrogate"]
    bandit = tool_section["bandit"]

    assert interrogate == {"fail-under": 32, "color": True}
    assert bandit == {
        "skips": ["B404", "B311"],
        "exclude_dirs": [
            ".venv",
            "tests",
            "artifacts",
            ".pytest_cache",
            ".ruff_cache",
        ],
    }


def test_root_pyproject_uses_only_the_shared_dev_group() -> None:
    assert set(_root_pyproject()["dependency-groups"]) == {"dev"}


def test_root_make_does_not_declare_proteomics_only_extensions() -> None:
    root_make = (REPO_ROOT / "makes" / "root.mk").read_text(encoding="utf-8")

    assert not any(command in root_make for command in PROTEOMICS_ONLY_EXTENSION_COMMANDS)
