from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]


class MkDocsLoader(yaml.SafeLoader):
    """SafeLoader that tolerates MkDocs python/name tags."""


def _construct_unknown(
    loader: MkDocsLoader,
    tag_suffix: str,
    node: yaml.Node,
) -> object | None:
    del tag_suffix
    if isinstance(node, yaml.ScalarNode):
        return loader.construct_scalar(node)
    if isinstance(node, yaml.SequenceNode):
        return loader.construct_sequence(node)
    if isinstance(node, yaml.MappingNode):
        return loader.construct_mapping(node)
    return None


cast(Any, MkDocsLoader).add_multi_constructor("", _construct_unknown)


def _load_yaml(path: Path) -> dict[str, object]:
    payload = yaml.load(path.read_text(encoding="utf-8"), Loader=MkDocsLoader) or {}
    return cast(dict[str, object], payload)


def test_mkdocs_extra_javascript_assets_exist() -> None:
    shared_config = _load_yaml(REPO_ROOT / "mkdocs.shared.yml")
    javascript_assets = shared_config["extra_javascript"]

    assert isinstance(javascript_assets, list)

    missing_assets = [
        asset_path
        for asset_path in javascript_assets
        if isinstance(asset_path, str)
        and asset_path.startswith("assets/")
        and not (REPO_ROOT / "docs" / asset_path).exists()
    ]

    assert not missing_assets, (
        "mkdocs.shared.yml references missing docs assets:\n"
        + "\n".join(sorted(missing_assets))
    )
