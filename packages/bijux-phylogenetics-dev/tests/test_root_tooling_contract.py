from __future__ import annotations

from pathlib import Path
import tomllib
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[3]
PROTEOMICS_ONLY_EXTENSION_COMMANDS = {
    "ensure-venv:",
    "nlenv:",
    "manage_examples:",
    "manage_models:",
    "api-freeze:",
    "openapi-drift:",
    "architecture-check:",
}
REQUIRED_ROOT_TARGET_SNIPPETS = {
    "list-evidence-studies:",
    "build-evidence-book:",
    "build-evidence-study:",
    "build-evidence-unit:",
    "validate-evidence-book:",
    "sync-evidence-artifacts:",
    "sync-evidence-unit-artifacts:",
    "sync-evidence-unit-inputs:",
    "check-evidence-artifacts:",
    "check-evidence-unit-artifacts:",
    "check-evidence-unit-inputs:",
    "report-evidence-completeness:",
    "check-evidence-completeness:",
    "report-evidence-governance:",
    "check-evidence-governance:",
    "report-artifact-governance:",
    "check-artifact-governance:",
    "report-execution-surfaces:",
    "check-execution-surfaces:",
    "report-package-boundaries:",
    "check-package-boundaries:",
    "report-package-bundles:",
    "check-package-bundles:",
    "check-config-ssot:",
    "report-publish-readiness:",
    "check-publish-readiness:",
    "report-release-readiness:",
    "check-release-readiness:",
    "test-all:",
    "test-all-plus-run-time:",
    "check: sync-license-assets lock-check check-config-ssot check-evidence-governance check-execution-surfaces check-package-boundaries lint test quality security docs build sbom",
}


def _root_pyproject() -> dict[str, Any]:
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)


def test_root_pyproject_declares_shared_quality_tooling() -> None:
    tool_section = _root_pyproject()["tool"]
    interrogate = tool_section["interrogate"]
    bandit = tool_section["bandit"]

    assert interrogate == {"fail-under": 32, "color": True}
    assert bandit == {
        "skips": ["B404", "B311"],
        "exclude_dirs": [
            ".venv",
            "tests",
            "artifacts",
            ".pytest_cache",
            ".ruff_cache",
        ],
    }


def test_root_pyproject_uses_only_the_shared_dev_group() -> None:
    assert set(_root_pyproject()["dependency-groups"]) == {"dev"}


def test_root_make_does_not_declare_proteomics_only_extensions() -> None:
    root_make = (REPO_ROOT / "makes" / "root.mk").read_text(encoding="utf-8")

    assert not any(
        command in root_make for command in PROTEOMICS_ONLY_EXTENSION_COMMANDS
    )


def test_root_make_wires_config_ssot_into_repository_checks() -> None:
    root_make = (REPO_ROOT / "makes" / "root.mk").read_text(encoding="utf-8")

    assert all(snippet in root_make for snippet in REQUIRED_ROOT_TARGET_SNIPPETS)


def test_root_make_routes_test_all_across_repository_packages() -> None:
    root_make = (REPO_ROOT / "makes" / "root.mk").read_text(encoding="utf-8")

    assert "ROOT_PACKAGE_TARGETS += test-all test-all-plus-run-time" in root_make
    assert "ROOT_TARGET_GROUPS_test-all ?= check" in root_make
    assert "ROOT_TARGET_SHARED_ENV_test-all ?= 1" in root_make
    assert "ROOT_TARGET_PACKAGES_test-all := $(CHECK_PACKAGES)" in root_make
    assert (
        "test-all: ## Run every repository test surface, including slow, evaluation, and real-local tests"
        in root_make
    )


def test_root_make_routes_test_all_plus_run_time_across_repository_packages() -> None:
    root_make = (REPO_ROOT / "makes" / "root.mk").read_text(encoding="utf-8")

    assert "ROOT_PACKAGE_TARGETS += test-all test-all-plus-run-time" in root_make
    assert "ROOT_TARGET_GROUPS_test-all-plus-run-time ?= check" in root_make
    assert "ROOT_TARGET_SHARED_ENV_test-all-plus-run-time ?= 1" in root_make
    assert (
        "ROOT_TARGET_PACKAGES_test-all-plus-run-time := $(CHECK_PACKAGES)" in root_make
    )
    assert (
        "test-all-plus-run-time: ## Run every repository test surface and report per-test durations"
        in root_make
    )


def test_package_makefiles_defer_monorepo_root_dependent_paths() -> None:
    package_makefiles = [
        REPO_ROOT / "makes" / "packages" / "bijux-phylogenetics.mk",
        REPO_ROOT / "makes" / "packages" / "bijux-phylogenetics-dev.mk",
        REPO_ROOT / "makes" / "packages" / "phylogenetic.mk",
    ]

    for makefile in package_makefiles:
        text = makefile.read_text(encoding="utf-8")
        assert ":= $(MONOREPO_ROOT)" not in text
        assert ':= "$(MONOREPO_ROOT)' not in text


def test_root_apis_surface_has_no_placeholder_readme() -> None:
    apis_root = REPO_ROOT / "apis"

    assert not (apis_root / "README.md").exists()
    if not apis_root.exists():
        return

    assert all(path.is_dir() for path in apis_root.iterdir())
