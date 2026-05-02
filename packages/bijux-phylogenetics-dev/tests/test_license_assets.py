from __future__ import annotations

from bijux_phylogenetics_dev.release.license_assets import managed_assets


def test_managed_assets_cover_every_workspace_package() -> None:
    assets = managed_assets()
    targets = {asset.target.parent.name for asset in assets}
    assert targets == {
        "bijux-phylogenetics",
        "bijux-phylogenetics-dev",
    }
