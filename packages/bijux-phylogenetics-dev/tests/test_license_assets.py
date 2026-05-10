from __future__ import annotations

from bijux_phylogenetics_dev.release.license_assets import (
    ROOT_LEGAL_ARTIFACTS,
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


def test_license_assets_link_back_to_repository_root() -> None:
    failures: list[str] = []

    for asset in managed_assets():
        if not asset.target.is_symlink():
            failures.append(f"{asset.target}: expected symlink")
            continue
        expected = ROOT_LEGAL_ARTIFACTS[asset.target.name]
        target = asset.target.readlink()
        if target != expected:
            failures.append(f"{asset.target}: {target!s} != {expected!s}")

    assert not failures, "managed legal asset linkage failed:\n" + "\n".join(failures)


def test_license_assets_are_synchronized() -> None:
    assert synchronize_license_assets(check=True) == []
