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
    assert "repository-contracts" in envlist
    assert "config-ssot" in envlist
    assert "evidence-governance" in envlist
    assert "evidence-completeness" in envlist
    assert "publish-readiness" in envlist
    assert "release-readiness-gate" in envlist
    assert "security" in envlist
    assert "docs" in envlist
    assert "fmt-{dev,core}" not in envlist
    assert "api-freeze-core" not in envlist
    assert "openapi-drift-core" not in envlist


def test_root_make_declares_shared_maintainer_commands() -> None:
    root_make = (REPO_ROOT / "makes" / "root.mk").read_text(encoding="utf-8")

    assert "check:" in root_make
    assert "test-all:" in root_make
    assert "sync-badges:" in root_make
    assert "check-badges:" in root_make
    assert "validate-evidence-book:" in root_make
    assert "report-release-readiness:" in root_make


def test_root_tox_isolates_repository_evidence_and_publish_surfaces() -> None:
    config = _tox_config()

    assert config["testenv:repository-contracts"]["change_dir"] == "{tox_root}"
    assert (
        config["testenv:repository-contracts"]["commands"].strip()
        == "make check-shared-bijux-py check-config-layout check-make-layout help"
    )
    assert config["testenv:evidence-governance"]["change_dir"] == "{tox_root}"
    assert (
        config["testenv:evidence-governance"]["commands"].strip()
        == "make check-evidence-governance\nmake rerun-governed-evidence-cleanroom"
    )
    assert config["testenv:evidence-completeness"]["change_dir"] == "{tox_root}"
    assert (
        config["testenv:evidence-completeness"]["commands"].strip()
        == "make check-evidence-completeness"
    )
    assert config["testenv:config-ssot"]["change_dir"] == "{tox_root}"
    assert config["testenv:config-ssot"]["commands"].strip() == "make check-config-ssot"
    assert config["testenv:publish-readiness"]["change_dir"] == "{tox_root}"
    assert (
        config["testenv:publish-readiness"]["commands"].strip()
        == "make report-release-readiness"
    )
    assert config["testenv:release-readiness-gate"]["change_dir"] == "{tox_root}"
    assert (
        config["testenv:release-readiness-gate"]["commands"].strip()
        == "make check-release-readiness"
    )


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


def test_top_level_runtime_exports_do_not_leak_evidence_book_helpers() -> None:
    package_init = (
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "src"
        / "bijux_phylogenetics"
        / "__init__.py"
    )
    module = ast.parse(package_init.read_text(encoding="utf-8"))

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

    assert "EvidenceBundleReport" not in exported_names
    assert "bundle_directory" not in exported_names


def test_runtime_workflows_use_provenance_bundle_contracts_instead_of_evidence_modules() -> (
    None
):
    runtime_paths = [
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "src"
        / "bijux_phylogenetics"
        / "bayesian"
        / "evidence.py",
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "src"
        / "bijux_phylogenetics"
        / "engines"
        / "evidence.py",
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "src"
        / "bijux_phylogenetics"
        / "core"
        / "demo.py",
    ]
    for path in runtime_paths:
        text = path.read_text(encoding="utf-8")
        assert "bijux_phylogenetics.evidence.bundles" not in text


def test_runtime_package_make_exposes_unfiltered_test_all_surface() -> None:
    runtime_make = (
        REPO_ROOT / "makes" / "packages" / "bijux-phylogenetics.mk"
    ).read_text(encoding="utf-8")

    assert (
        'TEST_MAIN_ARGS = -m "not slow and not real_local and not evaluation"'
        in runtime_make
    )
    assert "test-all: TEST_MAIN_ARGS =" in runtime_make
    assert "test-all: PYTEST_ADDOPTS_EXTRA = -o timeout=0" in runtime_make
    assert "test-all: test" in runtime_make


def test_avian_dataset_export_regression_surfaces_stay_slow_marked() -> None:
    module_path = (
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_avian_reproductive_trait_dataset.py"
    )
    module = ast.parse(module_path.read_text(encoding="utf-8"))

    expected_slow_functions = {
        "test_write_avian_reproductive_trait_workflow_bundle_matches_packaged_expected_outputs",
        "test_export_avian_reproductive_trait_dataset_copies_expected_outputs",
    }
    slow_functions: set[str] = set()
    for node in module.body:
        if not isinstance(node, ast.FunctionDef):
            continue
        for decorator in node.decorator_list:
            if not isinstance(decorator, ast.Attribute):
                continue
            if (
                isinstance(decorator.value, ast.Attribute)
                and isinstance(decorator.value.value, ast.Name)
                and decorator.value.value.id == "pytest"
                and decorator.value.attr == "mark"
                and decorator.attr == "slow"
            ):
                slow_functions.add(node.name)

    assert expected_slow_functions <= slow_functions


def test_repository_test_all_surface_disables_pytest_timeout_in_all_packages() -> None:
    root_make = (REPO_ROOT / "makes" / "root.mk").read_text(encoding="utf-8")
    dev_make = (
        REPO_ROOT / "makes" / "packages" / "bijux-phylogenetics-dev.mk"
    ).read_text(encoding="utf-8")
    alias_make = (REPO_ROOT / "makes" / "packages" / "phylogenetic.mk").read_text(
        encoding="utf-8"
    )

    assert "PYTEST_ADDOPTS_EXTRA='-o timeout=0'" in root_make
    assert "test-all: PYTEST_ADDOPTS_EXTRA = -o timeout=0" in dev_make
    assert "test-all: PYTEST_ADDOPTS_EXTRA = -o timeout=0" in alias_make


def test_root_conftest_registers_markers_from_repository_pytest_config() -> None:
    conftest_path = REPO_ROOT / "packages" / "bijux-phylogenetics" / "conftest.py"
    conftest_text = conftest_path.read_text(encoding="utf-8")

    assert 'PYTEST_CONFIG_PATH = REPO_ROOT / "configs" / "pytest.ini"' in conftest_text
    assert 'config.getini("markers")' in conftest_text
    assert 'config.addinivalue_line("markers", marker)' in conftest_text
