"""Workspace package layout coverage."""

from __future__ import annotations

import os
from pathlib import Path
import tomllib
from typing import Any, cast

REPO_ROOT = Path(__file__).resolve().parents[3]
PACKAGE_ARTIFACT_LINKS = {
    "artifacts": "",
    ".venv": "venv",
    ".hypothesis": "hypothesis",
    ".benchmarks": "benchmarks",
}


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


def test_package_roots_use_repository_artifact_symlinks() -> None:
    workspace = _workspace_metadata()
    failures: list[str] = []

    for package_name in sorted(cast(list[str], workspace["packages"])):
        package_root = _package_path(package_name)
        for link_name in PACKAGE_ARTIFACT_LINKS:
            link_path = package_root / link_name
            if not link_path.is_symlink():
                failures.append(f"{package_name}: missing symlink {link_name}")
                continue
            target = os.readlink(link_path)
            expected = _expected_package_link_target(package_name, link_name)
            if target != expected:
                failures.append(
                    f"{package_name}: {link_name} -> {target!r}, expected {expected!r}"
                )

    assert not failures, "package artifact symlink contract failed:\n" + "\n".join(
        failures
    )
