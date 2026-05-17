from __future__ import annotations

from pathlib import Path
import tomllib

from bijux_phylogenetics_dev.release.license_assets import (
    managed_assets,
    synchronize_license_assets,
)


def test_managed_assets_cover_every_workspace_package() -> None:
    assets = managed_assets()
    targets = {asset.target.parent.name for asset in assets}
    assert targets == {
        "bijux-phylogenetics",
        "phylogenetic",
        "bijux-phylogenetics-dev",
    }


def test_license_assets_match_repository_root_contents() -> None:
    failures: list[str] = []

    for asset in managed_assets():
        if asset.target.is_symlink() or not asset.target.is_file():
            failures.append(f"{asset.target}: expected synchronized regular file")
            continue
        if asset.target.read_bytes() != asset.source.read_bytes():
            failures.append(f"{asset.target}: content drift from {asset.source}")

    assert not failures, "managed legal asset linkage failed:\n" + "\n".join(failures)


def test_license_assets_are_synchronized() -> None:
    assert synchronize_license_assets(check=True) == []


def test_python_packages_declare_package_local_legal_files() -> None:
    package_projects = (
        Path("packages/bijux-phylogenetics/pyproject.toml"),
        Path("packages/bijux-phylogenetics-dev/pyproject.toml"),
        Path("packages/phylogenetic/pyproject.toml"),
    )

    for project_path in package_projects:
        project = tomllib.loads(project_path.read_text(encoding="utf-8"))
        assert project["project"]["license-files"] == ["LICENSE", "NOTICE"]
        sdist_config = project["tool"]["hatch"]["build"]["targets"]["sdist"]
        assert sdist_config["exclude"] == ["tests"]
