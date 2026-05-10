from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics_dev.release.publication_guard import (
    assert_publishable_repository,
)


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _minimal_repo(tmp_path: Path) -> Path:
    repo_root = tmp_path / "repo"
    _write(
        repo_root / "pyproject.toml",
        """
[tool.bijux_phylogenetics]
config_dir = "configs"
make_dir = "makes"
""".strip()
        + "\n",
    )
    _write(
        repo_root / "configs" / "config_ssot.toml",
        """
[tool.bijux_phylogenetics.config_ssot]
required_root_files = ["configs/mypy.ini", "configs/pytest.ini"]
forbidden_package_config_filenames = ["mypy.ini", "pytest.ini"]
allowed_package_config_paths = []
audit_paths = ["tox.ini", "makes/packages/runtime.mk"]
expected_root_config_dir = "configs"
expected_root_make_dir = "makes"
expected_mypy_config_path = "configs/mypy.ini"
""".strip()
        + "\n",
    )
    _write(repo_root / "configs" / "mypy.ini", "[mypy]\nstrict = true\n")
    _write(repo_root / "configs" / "pytest.ini", "[pytest]\naddopts = -ra\n")
    _write(repo_root / "tox.ini", "[tox]\nminversion = 4.11\n")
    _write(
        repo_root / "makes" / "packages" / "runtime.mk",
        'MYPY_CONFIG = $(MONOREPO_ROOT)/configs/mypy.ini\n',
    )
    _write(repo_root / "packages" / "runtime" / "pyproject.toml", "[project]\nname='runtime'\n")
    return repo_root


def test_assert_publishable_repository_allows_clean_config_ssot_repo(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)

    assert_publishable_repository(repo_root=repo_root, require_config_ssot=True)


def test_assert_publishable_repository_rejects_config_ssot_drift(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    _write(repo_root / "packages" / "runtime" / "mypy.ini", "[mypy]\nstrict = false\n")

    with pytest.raises(SystemExit, match="config SSOT audit failed"):
        assert_publishable_repository(repo_root=repo_root, require_config_ssot=True)


def test_assert_publishable_repository_requires_repo_root_for_config_ssot() -> None:
    with pytest.raises(
        ValueError, match="repo_root is required when require_config_ssot is enabled"
    ):
        assert_publishable_repository(require_config_ssot=True)
