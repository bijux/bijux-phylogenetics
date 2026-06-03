from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import pytest

from bijux_phylogenetics_dev.quality.evidence_artifacts import sync_evidence_artifacts
from bijux_phylogenetics_dev.quality.evidence_inputs import build_inputs_manifest
from bijux_phylogenetics_dev.quality.policies import (
    CONFIG_SSOT_POLICY_PATH,
    EXECUTION_SURFACES_POLICY_PATH,
    PACKAGE_BOUNDARIES_POLICY_PATH,
    PUBLICATION_READINESS_POLICY_PATH,
)
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
        repo_root / CONFIG_SSOT_POLICY_PATH,
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
    _write(
        repo_root / PACKAGE_BOUNDARIES_POLICY_PATH,
        """
[tool.bijux_phylogenetics.package_boundaries]
known_repo_module_roots = ["bijux_phylogenetics", "phylogenetic", "bijux_phylogenetics_dev"]
forbidden_runtime_top_level_exports = ["EvidenceBundleReport", "bundle_directory"]
alias_allowed_local_files = ["__init__.py", "__main__.py", "cli.py", "py.typed", "runtime_alias.py"]

[tool.bijux_phylogenetics.package_boundaries.runtime_evidence_compatibility]
runtime_version_spec = ">=0.1.0,<1.0"
supported_api_modules = ["bijux_phylogenetics.comparative"]
supported_api_locators = ["bijux_phylogenetics.comparative:inspect_pgls_inputs", "bijux_phylogenetics.comparative:run_pgls"]
notes = "Minimal runtime consumer contract."

[tool.bijux_phylogenetics.package_boundaries.package_roles."bijux-phylogenetics"]
role = "runtime"
package_dir = "packages/bijux-phylogenetics"
module_root = "bijux_phylogenetics"
allowed_repo_import_roots = ["bijux_phylogenetics"]
owned_module_prefixes = ["bijux_phylogenetics"]
required_install_dependencies = ["biopython>=1.0"]

[tool.bijux_phylogenetics.package_boundaries.package_roles."phylogenetic"]
role = "compatibility-alias"
package_dir = "packages/phylogenetic"
module_root = "phylogenetic"
allowed_repo_import_roots = ["phylogenetic", "bijux_phylogenetics"]
owned_module_prefixes = ["phylogenetic"]
required_install_dependencies = ["bijux-phylogenetics>=0.1.0,<1.0"]

[tool.bijux_phylogenetics.package_boundaries.package_roles."bijux-phylogenetics-dev"]
role = "maintainer"
package_dir = "packages/bijux-phylogenetics-dev"
module_root = "bijux_phylogenetics_dev"
allowed_repo_import_roots = ["bijux_phylogenetics_dev", "bijux_phylogenetics"]
owned_module_prefixes = ["bijux_phylogenetics_dev"]
required_install_dependencies = ["PyYAML>=6.0"]

[tool.bijux_phylogenetics.package_boundaries.target_package_roles."bijux-phylogenetics-evidence"]
role = "evidence-consumer"
target_module_root = "bijux_phylogenetics_evidence"
required_runtime_dependency = "bijux-phylogenetics>=0.1.0,<1.0"
""".strip()
        + "\n",
    )
    _write(
        repo_root / EXECUTION_SURFACES_POLICY_PATH,
        """
[tool.bijux_phylogenetics.execution_surfaces]
required_root_make_targets = [
  "check-evidence-governance:",
  "report-evidence-governance:",
  "check-evidence-completeness:",
  "report-evidence-completeness:",
  "check-execution-surfaces:",
  "report-execution-surfaces:",
]
required_tox_envs = [
  "repository-contracts",
  "config-ssot",
  "evidence-governance",
  "evidence-completeness",
  "publish-readiness",
  "release-readiness-gate",
]

[tool.bijux_phylogenetics.execution_surfaces.tox_commands]
repository-contracts = ["make check-shared-bijux-py check-config-layout check-make-layout help"]
config-ssot = ["make check-config-ssot"]
evidence-governance = ["make check-evidence-governance"]
evidence-completeness = ["make check-evidence-completeness"]
publish-readiness = ["make report-release-readiness"]
release-readiness-gate = ["make check-release-readiness"]
""".strip()
        + "\n",
    )
    _write(
        repo_root / PUBLICATION_READINESS_POLICY_PATH,
        """
[tool.bijux_phylogenetics.publication_readiness]
required_evidence_input_manifest = "inputs.manifest.json"
required_evidence_bundle_code_files = ["reference.R", "analysis.py"]
required_evidence_bundle_artifacts = ["reference.R", "analysis.py", "checks.json", "report.md", "provenance.json"]
expected_publishable_packages = [
    "bijux-phylogenetics",
    "phylogenetic",
    "bijux-phylogenetics-dev",
]
target_shape_packages = [
    "bijux-phylogenetics",
    "phylogenetic",
    "bijux-phylogenetics-dev",
]
forbidden_runtime_subpackages = ["evidence"]
required_root_make_targets = [
    "validate-evidence-book:",
    "report-evidence-completeness:",
    "check-evidence-completeness:",
    "report-evidence-governance:",
    "check-evidence-governance:",
    "sync-evidence-artifacts:",
    "check-evidence-artifacts:",
    "check-artifact-governance:",
    "report-execution-surfaces:",
    "check-execution-surfaces:",
    "report-package-boundaries:",
    "check-package-boundaries:",
    "report-package-bundles:",
    "check-package-bundles:",
    "report-publish-readiness:",
    "check-publish-readiness:",
    "report-release-readiness:",
    "check-release-readiness:",
]

[tool.bijux_phylogenetics.publication_readiness.package_policy."bijux-phylogenetics"]
package_dir = "packages/bijux-phylogenetics"
wheel_module_root = "bijux_phylogenetics"
allowed_dependencies = ["biopython>=1.0"]
required_sdist_entries = ["src/bijux_phylogenetics/__init__.py"]
required_wheel_entries = ["bijux_phylogenetics/__init__.py"]
forbidden_archive_prefixes = ["tests/", "docs/"]

[tool.bijux_phylogenetics.publication_readiness.package_policy."phylogenetic"]
package_dir = "packages/phylogenetic"
wheel_module_root = "phylogenetic"
allowed_dependencies = ["bijux-phylogenetics>=0.1.0,<1.0"]
required_sdist_entries = ["src/phylogenetic/__init__.py"]
required_wheel_entries = ["phylogenetic/__init__.py"]
forbidden_archive_prefixes = ["tests/", "docs/"]

[tool.bijux_phylogenetics.publication_readiness.package_policy."bijux-phylogenetics-dev"]
package_dir = "packages/bijux-phylogenetics-dev"
wheel_module_root = "bijux_phylogenetics_dev"
allowed_dependencies = ["PyYAML>=6.0"]
required_sdist_entries = ["src/bijux_phylogenetics_dev/__init__.py"]
required_wheel_entries = ["bijux_phylogenetics_dev/__init__.py"]
forbidden_archive_prefixes = ["tests/", "docs/"]
""".strip()
        + "\n",
    )
    _write(repo_root / "configs" / "mypy.ini", "[mypy]\nstrict = true\n")
    _write(repo_root / "configs" / "pytest.ini", "[pytest]\naddopts = -ra\n")
    _write(
        repo_root / "tox.ini",
        """
[tox]
minversion = 4.11
toxworkdir = {tox_root}/artifacts/root/tox
envlist =
    repository-contracts
    config-ssot
    evidence-governance
    evidence-completeness
    publish-readiness
    release-readiness-gate

[testenv:repository-contracts]
commands =
    make check-shared-bijux-py check-config-layout check-make-layout help

[testenv:config-ssot]
commands =
    make check-config-ssot

[testenv:evidence-governance]
commands =
    make check-evidence-governance

[testenv:evidence-completeness]
commands =
    make check-evidence-completeness

[testenv:publish-readiness]
commands =
    make report-release-readiness

[testenv:release-readiness-gate]
commands =
    make check-release-readiness
""".strip()
        + "\n",
    )
    _write(
        repo_root / "makes" / "root.mk",
        "\n".join(
            [
                "sync-evidence-artifacts:",
                "check-evidence-artifacts:",
                "validate-evidence-book:",
                "report-evidence-completeness:",
                "check-evidence-completeness:",
                "report-evidence-governance:",
                "check-evidence-governance:",
                "check-artifact-governance:",
                "report-execution-surfaces:",
                "check-execution-surfaces:",
                "report-package-boundaries:",
                "check-package-boundaries:",
                "report-package-bundles:",
                "check-package-bundles:",
                "report-publish-readiness:",
                "check-publish-readiness:",
                "report-release-readiness:",
                "check-release-readiness:",
            ]
        )
        + "\n",
    )
    _write(
        repo_root / "makes" / "packages" / "runtime.mk",
        "MYPY_CONFIG = $(MONOREPO_ROOT)/configs/mypy.ini\n",
    )
    _write(
        repo_root / "packages" / "bijux-phylogenetics" / "pyproject.toml",
        """
[project]
name = "bijux-phylogenetics"
dependencies = ["biopython>=1.0"]

[tool.hatch.build.targets.sdist]
include = ["src/bijux_phylogenetics/**"]

[tool.hatch.build.targets.wheel]
packages = ["src/bijux_phylogenetics"]
""".strip()
        + "\n",
    )
    _write(
        repo_root
        / "packages"
        / "bijux-phylogenetics"
        / "src"
        / "bijux_phylogenetics"
        / "__init__.py",
        """
from .comparative import inspect_pgls_inputs, run_pgls

__version__ = "0.1.0"
__all__ = ["inspect_pgls_inputs", "run_pgls"]
""".strip()
        + "\n",
    )
    _write(
        repo_root
        / "packages"
        / "bijux-phylogenetics"
        / "src"
        / "bijux_phylogenetics"
        / "comparative"
        / "evidence_contract.py",
        """
SUPPORTED_EVIDENCE_API_MODULES = ("bijux_phylogenetics.comparative",)
SUPPORTED_EVIDENCE_API_LOCATORS = (
    "bijux_phylogenetics.comparative:inspect_pgls_inputs",
    "bijux_phylogenetics.comparative:run_pgls",
)
""".strip()
        + "\n",
    )
    _write(
        repo_root
        / "packages"
        / "bijux-phylogenetics"
        / "src"
        / "bijux_phylogenetics"
        / "comparative"
        / "__init__.py",
        """
def inspect_pgls_inputs() -> dict[str, str]:
    return {"status": "ok"}


def run_pgls() -> dict[str, str]:
    return {"status": "ok"}


__all__ = ["inspect_pgls_inputs", "run_pgls"]
""".strip()
        + "\n",
    )
    _write(
        repo_root / "packages" / "phylogenetic" / "pyproject.toml",
        """
[project]
name = "phylogenetic"
dependencies = ["bijux-phylogenetics>=0.1.0,<1.0"]

[tool.hatch.build.targets.sdist]
include = ["src/phylogenetic/**"]

[tool.hatch.build.targets.wheel]
packages = ["src/phylogenetic"]
""".strip()
        + "\n",
    )
    _write(
        repo_root
        / "packages"
        / "phylogenetic"
        / "src"
        / "phylogenetic"
        / "__init__.py",
        'from bijux_phylogenetics import inspect_pgls_inputs, run_pgls\n\n__version__ = "0.1.0"\n',
    )
    _write(
        repo_root
        / "packages"
        / "phylogenetic"
        / "src"
        / "phylogenetic"
        / "__main__.py",
        "from .cli import main\n",
    )
    _write(
        repo_root / "packages" / "phylogenetic" / "src" / "phylogenetic" / "cli.py",
        "def main() -> int:\n    return 0\n",
    )
    _write(
        repo_root
        / "packages"
        / "phylogenetic"
        / "src"
        / "phylogenetic"
        / "runtime_alias.py",
        "def install_runtime_aliases() -> None:\n    return None\n",
    )
    _write(
        repo_root / "packages" / "phylogenetic" / "src" / "phylogenetic" / "py.typed",
        "",
    )
    _write(
        repo_root / "packages" / "bijux-phylogenetics-dev" / "pyproject.toml",
        """
[project]
name = "bijux-phylogenetics-dev"
dependencies = ["PyYAML>=6.0"]

[tool.hatch.build.targets.sdist]
include = ["src/bijux_phylogenetics_dev/**"]

[tool.hatch.build.targets.wheel]
packages = ["src/bijux_phylogenetics_dev"]
""".strip()
        + "\n",
    )
    _write(
        repo_root
        / "packages"
        / "bijux-phylogenetics-dev"
        / "src"
        / "bijux_phylogenetics_dev"
        / "__init__.py",
        "__version__ = '0.1.0'\n",
    )
    _write(
        repo_root / "evidence-book" / "studies" / "demo-study" / "README.md",
        "# Demo Study\n\nMinimal governed study.\n",
    )
    _write(
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "provenance"
        / "sources.json",
        """
{
  "intake_policy": "repository-owned-source",
  "schema_version": 1,
  "source_count": 1,
  "sources": [
    {
      "source_id": "demo-fixture",
      "kind": "repository-fixture",
      "label": "Demo fixture",
      "locator": "packages/bijux-phylogenetics/tests/fixtures/demo.tsv",
      "read_only": true
    }
  ],
  "study_id": "demo-study"
}
""".strip()
        + "\n",
    )
    _write(
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "datasets"
        / "registry.json",
        """
{
  "dataset_count": 1,
  "datasets": [
    {
      "dataset_id": "dataset-001",
      "kind": "repository-fixture",
      "label": "Demo fixture",
      "locator": "packages/bijux-phylogenetics/tests/fixtures/demo.tsv",
      "schema_summary": "Small governed fixture."
    }
  ],
  "schema_version": 1,
  "study_id": "demo-study"
}
""".strip()
        + "\n",
    )
    _write(
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "evidence-001"
        / "manifest.json",
        """
{
  "schema_version": 1,
  "evidence_id": "evidence-001",
  "evidence_title": "Demo evidence bundle",
  "study_id": "demo-study",
  "summary": "Minimal governed evidence bundle.",
  "owner_package": "bijux-phylogenetics",
  "claim_ids": ["demo-claim"],
  "freshness": {
    "governed_code_paths": ["packages/bijux-phylogenetics/src/bijux_phylogenetics"],
    "last_generated_on": "2026-05-10",
    "source_basis_locators": [
      "packages/bijux-phylogenetics/tests/fixtures/demo.tsv",
      "evidence-book/studies/demo-study/evidence-001/parity.json"
    ]
  },
  "ownership": {
    "analytical_surfaces": ["demo-analysis"],
    "owner_package": "bijux-phylogenetics"
  },
  "claim_tags": ["demo"],
  "comparison_mode": "direct_r_parity",
  "verdict": {
    "status": "matched",
    "summary": "Demo parity matched."
  },
  "limitations": [],
  "source_basis": [
    {
      "kind": "repository-fixture",
      "label": "Demo fixture",
      "locator": "packages/bijux-phylogenetics/tests/fixtures/demo.tsv"
    },
    {
      "kind": "repository-reference",
      "label": "Local parity payload",
      "locator": "evidence-book/studies/demo-study/evidence-001/parity.json"
    }
  ]
}
""".strip()
        + "\n",
    )
    _write(
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "evidence-001"
        / "analysis.py",
        "def run() -> dict[str, str]:\n    return {'status': 'ok'}\n",
    )
    _write(
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "evidence-001"
        / "reference.R",
        "status <- 'ok'\n",
    )
    _write(
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "evidence-001"
        / "parity.json",
        '{"status": "ok"}\n',
    )
    _write(
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "evidence-001"
        / "report.md",
        "# Demo report\n",
    )
    _write(
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "evidence-001"
        / "results.json",
        '{"status": "ok"}\n',
    )
    _write(
        repo_root
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "fixtures"
        / "demo.tsv",
        "species\tvalue\nA\t1\n",
    )
    sync_evidence_artifacts(repo_root)
    bundle_root = (
        repo_root / "evidence-book" / "studies" / "demo-study" / "evidence-001"
    )
    inputs_manifest = build_inputs_manifest(repo_root, bundle_root)
    _write(
        bundle_root / "inputs.manifest.json",
        json.dumps(inputs_manifest, indent=2, sort_keys=True) + "\n",
    )
    return repo_root


def test_assert_publishable_repository_allows_clean_config_ssot_repo(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)

    assert_publishable_repository(
        repo_root=repo_root,
        require_config_ssot=True,
        require_package_boundaries=True,
        require_publish_readiness=True,
    )


@pytest.mark.slow
def test_assert_publishable_repository_allows_clean_package_bundle_repo() -> None:
    assert_publishable_repository(
        repo_root=REPO_ROOT,
        require_config_ssot=True,
        require_package_bundles=True,
    )


@pytest.mark.slow
def test_assert_publishable_repository_allows_clean_publish_readiness_repo() -> None:
    assert_publishable_repository(
        repo_root=REPO_ROOT,
        require_package_boundaries=True,
        require_publish_readiness=True,
    )


def test_assert_publishable_repository_rejects_config_ssot_drift(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    _write(
        repo_root / "packages" / "bijux-phylogenetics" / "mypy.ini",
        "[mypy]\nstrict = false\n",
    )

    with pytest.raises(SystemExit, match="config SSOT audit failed"):
        assert_publishable_repository(
            repo_root=repo_root,
            require_config_ssot=True,
            require_publish_readiness=True,
        )


def test_assert_publishable_repository_requires_repo_root_for_repo_level_guards() -> (
    None
):
    with pytest.raises(
        ValueError,
        match="repo_root is required when repository-level publish guards are enabled",
    ):
        assert_publishable_repository(require_publish_readiness=True)


def test_assert_publishable_repository_rejects_package_boundary_drift(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    _write(
        repo_root / "packages" / "phylogenetic" / "src" / "phylogenetic" / "notes.py",
        "ALIASED = True\n",
    )

    with pytest.raises(SystemExit, match="package boundary audit failed"):
        assert_publishable_repository(
            repo_root=repo_root,
            require_package_boundaries=True,
        )


def test_assert_publishable_repository_rejects_publish_readiness_drift(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    _write(repo_root / "docs" / ".DS_Store", "")

    with pytest.raises(SystemExit, match="publish-readiness report failed"):
        assert_publishable_repository(
            repo_root=repo_root,
            require_publish_readiness=True,
        )


@pytest.mark.slow
def test_publication_guard_module_runs_without_runpy_warning() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "-Werror",
            "-m",
            "bijux_phylogenetics_dev.release.publication_guard",
            "--pyproject",
            str(REPO_ROOT / "packages" / "bijux-phylogenetics" / "pyproject.toml"),
            "--package-name",
            "bijux-phylogenetics",
            "--repo-root",
            str(REPO_ROOT),
            "--require-config-ssot",
            "--require-package-bundles",
            "--allow-local-version",
            "--allow-prerelease",
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr


REPO_ROOT = Path(__file__).resolve().parents[3]
