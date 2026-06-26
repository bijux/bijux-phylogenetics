from __future__ import annotations

from pathlib import Path
import tomllib
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[3]
WORKSPACE_TOOL = "bijux_phylogenetics"


def _root_pyproject() -> dict[str, Any]:
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)


def test_root_pyproject_uses_shared_workspace_build_contract() -> None:
    pyproject = _root_pyproject()

    assert pyproject["build-system"] == {
        "requires": ["hatchling>=1.27.0,<1.31", "hatch-vcs>=0.4.0,<1.0"],
        "build-backend": "hatchling.build",
    }

    project = pyproject["project"]
    assert project["dynamic"] == ["version"]
    assert "version" not in project

    hatch_version = pyproject["tool"]["hatch"]["version"]
    assert hatch_version["source"] == "vcs"
    assert hatch_version["tag-pattern"] == "^v(?P<version>.*)$"

    assert pyproject["tool"]["uv"]["workspace"]["members"] == ["packages/*"]
    assert pyproject["tool"]["hatch"]["build"]["targets"]["wheel"] == {
        "bypass-selection": True
    }


def test_root_pyproject_exposes_all_workspace_packages_to_root_dev_installs() -> None:
    pyproject = _root_pyproject()
    workspace_packages = set(pyproject["tool"][WORKSPACE_TOOL]["packages"])
    dev_group = pyproject["dependency-groups"]["dev"]
    dev_group_entries = {entry.split("[", 1)[0] for entry in dev_group}

    assert dev_group_entries.issuperset(workspace_packages)


def test_root_pyproject_uses_workspace_sources_for_every_workspace_package() -> None:
    pyproject = _root_pyproject()
    workspace_packages = set(pyproject["tool"][WORKSPACE_TOOL]["packages"])
    uv_sources = pyproject["tool"]["uv"]["sources"]

    assert set(uv_sources) == workspace_packages
    assert {
        name for name, config in uv_sources.items() if config == {"workspace": True}
    } == (workspace_packages)


def test_root_pyproject_declares_workspace_package_dirs_for_shared_make_tooling() -> (
    None
):
    pyproject = _root_pyproject()
    workspace_packages = set(pyproject["tool"][WORKSPACE_TOOL]["packages"])
    package_dirs = pyproject["tool"]["bijux_canon"]["package_dirs"]

    assert set(package_dirs) == workspace_packages
    assert all(
        (REPO_ROOT / relative_path).is_dir() for relative_path in package_dirs.values()
    )
