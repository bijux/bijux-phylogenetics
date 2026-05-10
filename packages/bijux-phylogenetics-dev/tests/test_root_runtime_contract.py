from __future__ import annotations

import ast
from configparser import ConfigParser
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]


def _tox_config() -> ConfigParser:
    parser = ConfigParser()
    parser.read(REPO_ROOT / "tox.ini", encoding="utf-8")
    return parser


def _envlist() -> set[str]:
    envlist = _tox_config()["tox"]["envlist"]
    return {line.strip() for line in envlist.splitlines() if line.strip()}


def test_root_tox_keeps_the_shared_env_families_and_drops_proteomics_only_ones() -> (
    None
):
    envlist = _envlist()

    assert _tox_config()["tox"]["toxworkdir"] == "{tox_root}/artifacts/root/tox"
    assert "security" in envlist
    assert "docs" in envlist
    assert "fmt-{dev,core}" not in envlist
    assert "api-freeze-core" not in envlist
    assert "openapi-drift-core" not in envlist


def test_root_make_declares_shared_maintainer_commands() -> None:
    root_make = (REPO_ROOT / "makes" / "root.mk").read_text(encoding="utf-8")

    assert "check:" in root_make
    assert "sync-badges:" in root_make
    assert "check-badges:" in root_make


def test_phylogenetic_alias_security_audits_the_installed_environment() -> None:
    package_make = (REPO_ROOT / "makes" / "packages" / "phylogenetic.mk").read_text(
        encoding="utf-8"
    )

    assert "SECURITY_AUDIT_PREPARE_MODE = pyproject" not in package_make
    assert 'PIP_AUDIT_INPUTS = -r "$(SECURITY_REQS)"' not in package_make


def test_top_level_runtime_exports_cover_every_relative_import() -> None:
    package_init = (
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "src"
        / "bijux_phylogenetics"
        / "__init__.py"
    )
    module = ast.parse(package_init.read_text(encoding="utf-8"))

    imported_names = [
        alias.asname or alias.name
        for node in module.body
        if isinstance(node, ast.ImportFrom) and node.level > 0
        for alias in node.names
    ]
    exported_names: list[str] = []
    for node in module.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(
            isinstance(target, ast.Name) and target.id == "__all__"
            for target in node.targets
        ):
            continue
        assert isinstance(node.value, ast.List)
        exported_names = [
            element.value
            for element in node.value.elts
            if isinstance(element, ast.Constant) and isinstance(element.value, str)
        ]
        break

    missing_exports = [name for name in imported_names if name not in exported_names]
    assert missing_exports == []
