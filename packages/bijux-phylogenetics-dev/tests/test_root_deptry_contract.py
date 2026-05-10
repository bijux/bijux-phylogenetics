from __future__ import annotations

from pathlib import Path
import tomllib
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[3]


def _deptry_config() -> dict[str, Any]:
    with (REPO_ROOT / "configs" / "deptry.toml").open("rb") as handle:
        return tomllib.load(handle)


def test_root_deptry_configuration_uses_supported_dev_group_contract() -> None:
    deptry_config = _deptry_config()["tool"]["deptry"]

    assert deptry_config["optional_dependencies_dev_groups"] == ["dev"]
    assert "pep621_dev_dependency_groups" not in deptry_config
    assert deptry_config["ignore"] == ["DEP002", "DEP003"]
    assert deptry_config["extend_exclude"] == ["docs", "apis", "artifacts", "site"]
