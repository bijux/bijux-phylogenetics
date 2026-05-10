from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics_dev.quality.config_ssot import (
    build_config_ssot_report,
    check_config_ssot,
    load_config_ssot_policy,
)

REPO_ROOT = Path(__file__).resolve().parents[3]


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


def test_load_config_ssot_policy_reads_repo_owned_policy() -> None:
    policy = load_config_ssot_policy(REPO_ROOT)

    assert policy.expected_root_config_dir == "configs"
    assert policy.expected_root_make_dir == "makes"
    assert policy.expected_mypy_config_path == "configs/mypy.ini"
    assert "mypy.ini" in policy.forbidden_package_config_filenames
    assert "makes/packages/bijux-phylogenetics.mk" in policy.audit_paths


def test_build_config_ssot_report_is_clean_for_the_repository() -> None:
    report = build_config_ssot_report(REPO_ROOT)

    assert report.issue_count == 0, report.to_dict()
    assert report.package_config_allowlist == []
    assert report.package_local_config_files == []


def test_build_config_ssot_report_flags_forbidden_package_local_config(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    _write(repo_root / "packages" / "runtime" / "mypy.ini", "[mypy]\nstrict = false\n")

    report = build_config_ssot_report(repo_root)

    assert report.issue_count == 1
    assert report.issues[0].code == "forbidden-package-config"
    assert report.issues[0].path == "packages/runtime/mypy.ini"


def test_build_config_ssot_report_flags_package_local_mypy_reference(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    _write(
        repo_root / "makes" / "packages" / "runtime.mk",
        "MYPY_CONFIG = $(PROJECT_DIR)/mypy.ini\n",
    )

    report = build_config_ssot_report(repo_root)

    assert report.issue_count == 2
    issue_codes = {issue.code for issue in report.issues}
    assert issue_codes == {"local-mypy-reference", "mypy-config-path-drift"}


def test_check_config_ssot_writes_json_report(tmp_path: Path) -> None:
    repo_root = _minimal_repo(tmp_path)
    report_path = repo_root / "artifacts" / "root" / "config-ssot-audit.json"

    report = check_config_ssot(repo_root, json_out=report_path)

    assert report.issue_count == 0
    assert report_path.is_file()
