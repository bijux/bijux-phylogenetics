from __future__ import annotations

from pathlib import Path
import tomllib

from bijux_phylogenetics_dev.quality.config_ssot import load_config_ssot_policy

REPO_ROOT = Path(__file__).resolve().parents[3]


def _root_pyproject() -> dict[str, object]:
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)


def test_root_pyproject_points_to_repository_owned_config_and_make_dirs() -> None:
    workspace = _root_pyproject()["tool"]["bijux_phylogenetics"]

    assert workspace["config_dir"] == "configs"
    assert workspace["make_dir"] == "makes"


def test_config_ssot_policy_required_root_files_exist() -> None:
    policy = load_config_ssot_policy(REPO_ROOT)

    assert policy.required_root_files
    assert all((REPO_ROOT / relative_path).is_file() for relative_path in policy.required_root_files)


def test_package_makefiles_keep_mypy_on_the_root_config_surface() -> None:
    policy = load_config_ssot_policy(REPO_ROOT)
    expected = policy.expected_mypy_config_path

    for relative_path in policy.audit_paths:
        text = (REPO_ROOT / relative_path).read_text(encoding="utf-8")
        assert "$(PROJECT_DIR)/mypy.ini" not in text
        if "MYPY_CONFIG" in text or "QUALITY_MYPY_CONFIG" in text:
            assert expected in text


def test_forbidden_package_local_config_files_do_not_exist() -> None:
    policy = load_config_ssot_policy(REPO_ROOT)

    forbidden_paths = [
        package_root / filename
        for package_root in (REPO_ROOT / "packages").glob("*")
        if package_root.is_dir()
        for filename in policy.forbidden_package_config_filenames
    ]
    allowlist = {
        (REPO_ROOT / relative_path).resolve()
        for relative_path in policy.allowed_package_config_paths
    }

    offenders = [
        path.relative_to(REPO_ROOT).as_posix()
        for path in forbidden_paths
        if path.exists() and path.resolve() not in allowlist
    ]
    assert offenders == []
