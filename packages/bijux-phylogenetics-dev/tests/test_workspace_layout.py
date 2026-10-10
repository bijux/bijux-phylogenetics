"""Workspace package layout coverage."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tomllib
from typing import Any, cast

REPO_ROOT = Path(__file__).resolve().parents[3]
ARTIFACT_ALIAS_SCRIPT = (
    REPO_ROOT
    / ".bijux"
    / "shared"
    / "bijux-makes-py"
    / "repository"
    / "artifact_aliases.py"
)
ROOT_ARTIFACT_ALIAS_WRAPPER = REPO_ROOT / "makes" / "repository_artifact_aliases.py"
PACKAGE_ARTIFACT_LINKS = {
    "artifacts": "",
    ".venv": "venv",
    ".hypothesis": "hypothesis",
    ".benchmarks": "benchmarks",
}
ROOT_ARTIFACT_LINKS = {
    ".venv": "artifacts/root/check-venv",
    ".hypothesis": "artifacts/root/hypothesis",
    ".benchmarks": "artifacts/root/benchmarks",
    ".tox": "artifacts/root/tox",
}
ROOT_FORBIDDEN_CACHE_PATHS = (
    ".pytest_cache",
    ".ruff_cache",
)


def _workspace_metadata() -> dict[str, Any]:
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        data = tomllib.load(handle)
    return cast(dict[str, Any], data["tool"]["bijux_phylogenetics"])


def _package_path(package_name: str) -> Path:
    return REPO_ROOT / "packages" / package_name


def _expected_package_link_target(package_name: str, link_name: str) -> str:
    suffix = PACKAGE_ARTIFACT_LINKS[link_name]
    base = Path("..") / ".." / "artifacts" / package_name
    return str(base / suffix) if suffix else str(base)


def _artifact_alias_paths() -> list[str]:
    workspace = _workspace_metadata()
    package_names = sorted(cast(list[str], workspace["packages"]))
    paths = sorted(ROOT_ARTIFACT_LINKS)
    for package_name in package_names:
        for link_name in PACKAGE_ARTIFACT_LINKS:
            paths.append(f"packages/{package_name}/{link_name}")
    return paths


def _assert_symlink(*, link_path: Path, expected_target: str) -> None:
    assert link_path.is_symlink()
    assert link_path.readlink().as_posix() == expected_target


def test_setup_materializes_governed_artifact_aliases(tmp_path: Path) -> None:
    workspace = _workspace_metadata()
    repo_root = tmp_path / "repo"
    packages_dir = repo_root / "packages"
    packages_dir.mkdir(parents=True)

    for package_name in sorted(cast(list[str], workspace["packages"])):
        package_root = packages_dir / package_name
        package_root.mkdir()
        (package_root / "pyproject.toml").write_text("[project]\nname='test'\n")

    subprocess.run(
        [
            sys.executable,
            str(ARTIFACT_ALIAS_SCRIPT),
            "root",
            "--repo-root",
            str(repo_root),
            "--packages-dir",
            str(packages_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    for link_name, expected_target in ROOT_ARTIFACT_LINKS.items():
        _assert_symlink(
            link_path=repo_root / link_name,
            expected_target=expected_target,
        )

    for package_name in sorted(cast(list[str], workspace["packages"])):
        package_root = packages_dir / package_name
        for link_name in PACKAGE_ARTIFACT_LINKS:
            _assert_symlink(
                link_path=package_root / link_name,
                expected_target=_expected_package_link_target(package_name, link_name),
            )


def test_root_setup_preserves_existing_root_environment(tmp_path: Path) -> None:
    workspace = _workspace_metadata()
    repo_root = tmp_path / "repo"
    packages_dir = repo_root / "packages"
    packages_dir.mkdir(parents=True)
    conflicting_root_venv = repo_root / ".venv"
    conflicting_root_venv.mkdir()
    (conflicting_root_venv / "pyvenv.cfg").write_text("home = /tmp/python\n")

    for package_name in sorted(cast(list[str], workspace["packages"])):
        package_root = packages_dir / package_name
        package_root.mkdir()
        (package_root / "pyproject.toml").write_text("[project]\nname='test'\n")

    completed = subprocess.run(
        [
            sys.executable,
            str(ROOT_ARTIFACT_ALIAS_WRAPPER),
            "root",
            "--repo-root",
            str(repo_root),
            "--packages-dir",
            str(packages_dir),
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0
    assert completed.stderr == ""
    assert conflicting_root_venv.is_dir()
    assert not conflicting_root_venv.is_symlink()
    assert (conflicting_root_venv / "pyvenv.cfg").read_text() == "home = /tmp/python\n"


def test_artifact_alias_paths_are_ignored_by_git() -> None:
    paths = _artifact_alias_paths()
    result = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "check-ignore", "--no-index", "--stdin"],
        input="\n".join(paths),
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert set(result.stdout.splitlines()) == set(paths)


def test_artifact_alias_paths_stay_untracked_by_git() -> None:
    result = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "ls-files", "-z", "--", *_artifact_alias_paths()],
        check=True,
        capture_output=True,
    )

    assert result.stdout == b""


def test_repository_root_keeps_non_artifact_caches_out_of_tree() -> None:
    stray_configs_artifacts = REPO_ROOT / "configs" / "artifacts"
    assert not stray_configs_artifacts.exists()
    for path_name in ROOT_FORBIDDEN_CACHE_PATHS:
        assert not (REPO_ROOT / path_name).exists()
