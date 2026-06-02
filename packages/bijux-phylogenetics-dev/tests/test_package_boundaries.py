from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics_dev.quality.package_boundaries import (
    build_package_boundary_report,
    load_package_boundary_policy,
)
from bijux_phylogenetics_dev.quality.policies import PACKAGE_BOUNDARIES_POLICY_PATH

REPO_ROOT = Path(__file__).resolve().parents[3]


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _minimal_repo(tmp_path: Path) -> Path:
    repo_root = tmp_path / "repo"
    _write(
        repo_root / PACKAGE_BOUNDARIES_POLICY_PATH,
        """
        [tool.bijux_phylogenetics.package_boundaries]
known_repo_module_roots = ["demo_runtime", "demo_alias", "demo_dev"]
forbidden_runtime_top_level_exports = ["EvidenceBundleReport"]
alias_allowed_local_files = ["__init__.py", "__main__.py", "cli.py", "py.typed", "runtime_alias.py"]

        [tool.bijux_phylogenetics.package_boundaries.runtime_evidence_compatibility]
runtime_version_spec = ">=0.1.0,<1.0"
supported_api_modules = ["demo_runtime.api"]
supported_api_locators = ["demo_runtime.api:run_pgls"]
notes = "demo"

        [tool.bijux_phylogenetics.package_boundaries.package_roles."demo-runtime"]
role = "runtime"
package_dir = "packages/demo-runtime"
module_root = "demo_runtime"
allowed_repo_import_roots = ["demo_runtime"]
owned_module_prefixes = ["demo_runtime"]
required_install_dependencies = ["numpy>=1.0"]

        [tool.bijux_phylogenetics.package_boundaries.package_roles."demo-alias"]
role = "compatibility-alias"
package_dir = "packages/demo-alias"
module_root = "demo_alias"
allowed_repo_import_roots = ["demo_alias", "demo_runtime"]
owned_module_prefixes = ["demo_alias"]
required_install_dependencies = ["demo-runtime>=0.1.0,<1.0"]

        [tool.bijux_phylogenetics.package_boundaries.package_roles."demo-dev"]
role = "maintainer"
package_dir = "packages/demo-dev"
module_root = "demo_dev"
allowed_repo_import_roots = ["demo_dev", "demo_runtime"]
owned_module_prefixes = ["demo_dev"]
required_install_dependencies = ["PyYAML>=6.0"]

        [tool.bijux_phylogenetics.package_boundaries.target_package_roles."demo-evidence"]
role = "evidence-consumer"
target_module_root = "demo_evidence"
required_runtime_dependency = "demo-runtime>=0.1.0,<1.0"
""".strip()
        + "\n",
    )
    _write(
        repo_root / "packages" / "demo-runtime" / "pyproject.toml",
        """
[project]
name = "demo-runtime"
version = "0.1.0"
dependencies = ["numpy>=1.0"]
""".strip()
        + "\n",
    )
    _write(
        repo_root / "packages" / "demo-alias" / "pyproject.toml",
        """
[project]
name = "demo-alias"
version = "0.1.0"
dependencies = ["demo-runtime>=0.1.0,<1.0"]
""".strip()
        + "\n",
    )
    _write(
        repo_root / "packages" / "demo-dev" / "pyproject.toml",
        """
[project]
name = "demo-dev"
version = "0.1.0"
dependencies = ["PyYAML>=6.0"]
""".strip()
        + "\n",
    )
    _write(
        repo_root
        / "packages"
        / "demo-runtime"
        / "src"
        / "demo_runtime"
        / "__init__.py",
        """
from .api import run_pgls

__all__ = ["run_pgls"]
""".strip()
        + "\n",
    )
    _write(
        repo_root / "packages" / "demo-runtime" / "src" / "demo_runtime" / "api.py",
        "def run_pgls() -> str:\n    return 'ok'\n",
    )
    _write(
        repo_root
        / "packages"
        / "demo-runtime"
        / "src"
        / "demo_runtime"
        / "comparative"
        / "evidence_contract.py",
        """
SUPPORTED_EVIDENCE_API_MODULES = ("demo_runtime.api",)
SUPPORTED_EVIDENCE_API_LOCATORS = ("demo_runtime.api:run_pgls",)
""".strip()
        + "\n",
    )
    _write(
        repo_root / "packages" / "demo-alias" / "src" / "demo_alias" / "__init__.py",
        "from demo_runtime import run_pgls\n",
    )
    _write(
        repo_root / "packages" / "demo-alias" / "src" / "demo_alias" / "__main__.py",
        "from .cli import main\n",
    )
    _write(
        repo_root / "packages" / "demo-alias" / "src" / "demo_alias" / "cli.py",
        "def main() -> int:\n    return 0\n",
    )
    _write(
        repo_root
        / "packages"
        / "demo-alias"
        / "src"
        / "demo_alias"
        / "runtime_alias.py",
        "def install_runtime_aliases() -> None:\n    return None\n",
    )
    _write(
        repo_root / "packages" / "demo-alias" / "src" / "demo_alias" / "py.typed",
        "",
    )
    _write(
        repo_root / "packages" / "demo-dev" / "src" / "demo_dev" / "__init__.py",
        "from demo_runtime import run_pgls\n",
    )
    return repo_root


def test_load_package_boundary_policy_reads_repo_owned_policy() -> None:
    policy = load_package_boundary_policy(REPO_ROOT)

    assert "bijux-phylogenetics" in policy.package_roles
    assert "phylogenetic" in policy.package_roles
    assert policy.target_package_roles == {}
    assert policy.runtime_evidence_compatibility.runtime_version_spec == ">=0.1.0,<1.0"


def test_runtime_package_boundary_policy_allows_secured_xml_dependency() -> None:
    policy = load_package_boundary_policy(REPO_ROOT)

    assert policy.package_roles[
        "bijux-phylogenetics"
    ].required_install_dependencies == (
        "biopython>=1.87,<2.0",
        "cairosvg>=2.9.0,<3.0",
        "defusedxml>=0.7.1,<1.0",
        "PyYAML>=6.0,<7.0",
    )


def test_build_package_boundary_report_accepts_governed_repo(tmp_path: Path) -> None:
    repo_root = _minimal_repo(tmp_path)

    report = build_package_boundary_report(repo_root)

    assert report["issue_count"] == 0
    assert report["runtime_public_api"]["forbidden_top_level_exports_present"] == []
    assert report["runtime_evidence_compatibility"]["contract_locators"] == [
        "demo_runtime.api:run_pgls"
    ]


def test_build_package_boundary_report_flags_forbidden_runtime_exports(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    _write(
        repo_root
        / "packages"
        / "demo-runtime"
        / "src"
        / "demo_runtime"
        / "__init__.py",
        "__all__ = ['run_pgls', 'EvidenceBundleReport']\n",
    )

    report = build_package_boundary_report(repo_root)

    issue_codes = {issue["code"] for issue in report["issues"]}
    assert "forbidden-runtime-top-level-export" in issue_codes


def test_build_package_boundary_report_flags_cross_package_import_drift(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    _write(
        repo_root / "packages" / "demo-runtime" / "src" / "demo_runtime" / "bad.py",
        "from demo_dev import helper\n",
    )

    report = build_package_boundary_report(repo_root)

    issue_codes = {issue["code"] for issue in report["issues"]}
    assert "repo-import-boundary-drift" in issue_codes


def test_build_package_boundary_report_flags_runtime_evidence_contract_drift(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    _write(
        repo_root
        / "packages"
        / "demo-runtime"
        / "src"
        / "demo_runtime"
        / "comparative"
        / "evidence_contract.py",
        """
SUPPORTED_EVIDENCE_API_MODULES = ("demo_runtime.api",)
SUPPORTED_EVIDENCE_API_LOCATORS = ("demo_runtime.api:missing_export",)
""".strip()
        + "\n",
    )

    report = build_package_boundary_report(repo_root)

    issue_codes = {issue["code"] for issue in report["issues"]}
    assert "runtime-evidence-locator-contract-drift" in issue_codes


def test_build_package_boundary_report_allows_runtime_contract_supersets(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    _write(
        repo_root
        / "packages"
        / "demo-runtime"
        / "src"
        / "demo_runtime"
        / "comparative"
        / "evidence_contract.py",
        """
SUPPORTED_EVIDENCE_API_MODULES = ("demo_runtime.api", "demo_runtime.extra")
SUPPORTED_EVIDENCE_API_LOCATORS = ("demo_runtime.api:run_pgls", "demo_runtime.extra:run_extra")
""".strip()
        + "\n",
    )
    _write(
        repo_root / "packages" / "demo-runtime" / "src" / "demo_runtime" / "extra.py",
        "def run_extra() -> str:\n    return 'extra'\n",
    )

    report = build_package_boundary_report(repo_root)

    issue_codes = {issue["code"] for issue in report["issues"]}
    assert "runtime-evidence-module-contract-drift" not in issue_codes
    assert "runtime-evidence-locator-contract-drift" not in issue_codes


@pytest.mark.slow
def test_repository_package_boundary_report_is_clean() -> None:
    report = build_package_boundary_report(REPO_ROOT)

    assert report["issue_count"] == 0
