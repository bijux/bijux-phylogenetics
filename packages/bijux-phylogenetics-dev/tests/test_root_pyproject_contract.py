from __future__ import annotations

from pathlib import Path
import tomllib
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[3]
WORKSPACE_TOOL = "bijux_phylogenetics"


def _root_pyproject() -> dict[str, Any]:
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)


def test_root_pyproject_keeps_distribution_builds_in_workspace_packages() -> None:
    pyproject = _root_pyproject()
    assert "project" not in pyproject
    assert "build-system" not in pyproject
    assert "hatch" not in pyproject["tool"]
    assert pyproject["tool"]["uv"]["workspace"]["members"] == ["packages/*"]
    names = []
    for path in sorted((REPO_ROOT / "packages").glob("*/pyproject.toml")):
        package = tomllib.loads(path.read_text())
        names.append(package["project"]["name"])
        assert package["build-system"]["build-backend"] == "hatchling.build"
    assert set(names) == set(pyproject["tool"][WORKSPACE_TOOL]["packages"])
    assert len(names) == len(set(names))


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
    package_dirs = pyproject["tool"][WORKSPACE_TOOL]["package_dirs"]

    assert set(package_dirs) == workspace_packages
    assert all(
        (REPO_ROOT / relative_path).is_dir() for relative_path in package_dirs.values()
    )
