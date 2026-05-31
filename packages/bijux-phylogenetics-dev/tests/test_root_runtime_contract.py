from __future__ import annotations

import ast
from configparser import ConfigParser
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
PACKAGE_INIT = (
    REPO_ROOT
    / "packages"
    / "bijux-phylogenetics"
    / "src"
    / "bijux_phylogenetics"
    / "__init__.py"
)
EXPECTED_ROOT_GATEWAYS = [
    "__version__",
    "ancestral",
    "api",
    "bayesian",
    "biogeography",
    "comparative",
    "datasets",
    "distance",
    "evidence",
    "parsimony",
    "parity",
    "phylo",
    "trees",
]


def _tox_config() -> ConfigParser:
    parser = ConfigParser()
    parser.read(REPO_ROOT / "tox.ini", encoding="utf-8")
    return parser


def _envlist() -> set[str]:
    envlist = _tox_config()["tox"]["envlist"]
    return {line.strip() for line in envlist.splitlines() if line.strip()}


def _is_pytest_marker(node: ast.AST, marker_name: str) -> bool:
    return (
        isinstance(node, ast.Attribute)
        and node.attr == marker_name
        and isinstance(node.value, ast.Attribute)
        and node.value.attr == "mark"
        and isinstance(node.value.value, ast.Name)
        and node.value.value.id == "pytest"
    )


def _contains_pytest_slow_marker(node: ast.AST) -> bool:
    if _is_pytest_marker(node, "slow"):
        return True
    if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
        return any(_contains_pytest_slow_marker(element) for element in node.elts)
    return False


def _slow_marked_functions(module_path: Path) -> set[str]:
    module = ast.parse(module_path.read_text(encoding="utf-8"))
    module_has_slow_mark = False
    for node in module.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(
            isinstance(target, ast.Name) and target.id == "pytestmark"
            for target in node.targets
        ):
            continue
        module_has_slow_mark = _contains_pytest_slow_marker(node.value)
        if module_has_slow_mark:
            break

    function_names = {
        node.name for node in module.body if isinstance(node, ast.FunctionDef)
    }
    if module_has_slow_mark:
        return function_names

    slow_functions: set[str] = set()
    for node in module.body:
        if not isinstance(node, ast.FunctionDef):
            continue
        for decorator in node.decorator_list:
            if _is_pytest_marker(decorator, "slow"):
                slow_functions.add(node.name)
                break
    return slow_functions


def _stress_marked_functions(module_path: Path) -> set[str]:
    module = ast.parse(module_path.read_text(encoding="utf-8"))
    stress_functions: set[str] = set()
    for node in module.body:
        if not isinstance(node, ast.FunctionDef):
            continue
        for decorator in node.decorator_list:
            if _is_pytest_marker(decorator, "stress_small") or _is_pytest_marker(
                decorator, "stress_heavy"
            ):
                stress_functions.add(node.name)
                break
    return stress_functions


def _package_init_module() -> ast.Module:
    return ast.parse(PACKAGE_INIT.read_text(encoding="utf-8"))


def _package_init_exports() -> list[str]:
    module = _package_init_module()
    for node in module.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(
            isinstance(target, ast.Name) and target.id == "__all__"
            for target in node.targets
        ):
            continue
        assert isinstance(node.value, ast.List)
        return [
            element.value
            for element in node.value.elts
            if isinstance(element, ast.Constant) and isinstance(element.value, str)
        ]
    raise AssertionError("__all__ definition not found in package root")


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
    assert "ROOT_PACKAGE_TARGETS += test-all test-all-plus-run-time" in root_make
    assert "ROOT_TARGET_GROUPS_test-all ?= check" in root_make
    assert "ROOT_TARGET_GROUPS_test-all-plus-run-time ?= check" in root_make
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


def test_top_level_runtime_exports_define_curated_domain_gateways() -> None:
    assert _package_init_exports() == EXPECTED_ROOT_GATEWAYS


def test_top_level_runtime_exports_use_lazy_module_gateway_pattern() -> None:
    module = _package_init_module()

    import_from_names = [
        alias.asname or alias.name
        for node in module.body
        if isinstance(node, ast.ImportFrom) and node.level > 0
        for alias in node.names
    ]

    assert import_from_names == []


def test_top_level_runtime_exports_do_not_leak_leaf_level_symbols() -> None:
    exported_names = _package_init_exports()

    assert "EvidenceBundleReport" not in exported_names
    assert "bundle_directory" not in exported_names
    assert "trim_alignment" not in exported_names
    assert "run_pgls" not in exported_names


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


def test_runtime_package_make_exposes_unfiltered_test_all_plus_run_time_surface() -> (
    None
):
    runtime_make = (
        REPO_ROOT / "makes" / "packages" / "bijux-phylogenetics.mk"
    ).read_text(encoding="utf-8")

    assert (
        'TEST_MAIN_ARGS = -m "not slow and not real_local and not evaluation"'
        in runtime_make
    )
    assert "test-all-plus-run-time: TEST_MAIN_ARGS =" in runtime_make
    assert (
        "test-all-plus-run-time: PYTEST_ADDOPTS_EXTRA = -o timeout=0 --durations=0 --durations-min=0"
        in runtime_make
    )
    assert "test-all-plus-run-time: test" in runtime_make


def test_dev_package_make_excludes_slow_tests_by_default() -> None:
    dev_make = (
        REPO_ROOT / "makes" / "packages" / "bijux-phylogenetics-dev.mk"
    ).read_text(encoding="utf-8")

    assert 'TEST_MAIN_ARGS = -m "not slow"' in dev_make


def test_dev_package_make_exposes_unfiltered_test_all_surface() -> None:
    dev_make = (
        REPO_ROOT / "makes" / "packages" / "bijux-phylogenetics-dev.mk"
    ).read_text(encoding="utf-8")

    assert 'TEST_MAIN_ARGS = -m "not slow"' in dev_make
    assert "test-all: TEST_MAIN_ARGS =" in dev_make
    assert "test-all: PYTEST_ADDOPTS_EXTRA = -o timeout=0" in dev_make
    assert "test-all: test" in dev_make


def test_dev_package_make_exposes_unfiltered_test_all_plus_run_time_surface() -> None:
    dev_make = (
        REPO_ROOT / "makes" / "packages" / "bijux-phylogenetics-dev.mk"
    ).read_text(encoding="utf-8")

    assert 'TEST_MAIN_ARGS = -m "not slow"' in dev_make
    assert "test-all-plus-run-time: TEST_MAIN_ARGS =" in dev_make
    assert (
        "test-all-plus-run-time: PYTEST_ADDOPTS_EXTRA = -o timeout=0 --durations=0 --durations-min=0"
        in dev_make
    )
    assert "test-all-plus-run-time: test" in dev_make


def test_long_running_runtime_workflows_stay_slow_marked() -> None:
    expected_slow_functions_by_module = {
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_avian_reproductive_trait_dataset.py": {
            "test_write_avian_reproductive_trait_workflow_bundle_matches_packaged_expected_outputs",
            "test_run_avian_reproductive_trait_demo_materializes_dataset_and_workflow",
            "test_cli_demo_avian_reproductive_traits_json_output_reports_dataset_and_workflow",
        },
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_fasta_to_tree_real_workflows.py": {
            "test_run_fasta_to_tree_workflow_matches_real_output_golden",
            "test_adapter_fasta_to_tree_cli_matches_real_output_golden",
        },
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_large_dataset_stress.py": {
            "test_benchmark_large_dataset_stress_suite_small_tier_reports_all_workloads",
            "test_cli_benchmark_stress_suite_reports_tier_and_observation_count",
            "test_benchmark_large_dataset_stress_suite_heavy_tier_meets_large_input_thresholds",
        },
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_validation_corpus_iteration.py": {
            "test_build_method_accuracy_dashboard_summarizes_fixture_and_corpus_pass_rates",
            "test_build_runtime_and_memory_dashboards_cover_sites_and_posterior_samples",
            "test_build_scientific_validation_report_separates_claim_statuses",
        },
        REPO_ROOT / "packages" / "bijux-phylogenetics" / "tests" / "test_runtime.py": {
            "test_cli_benchmark_tree_validation_reports_observations",
            "test_cli_evidence_book_validate_json_output",
            "test_cli_evidence_book_build_json_output",
            "test_cli_evidence_book_build_selected_evidence_json_output",
            "test_cli_evidence_book_rerun_json_output",
        },
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_primate_comparative_dataset.py": {
            "test_write_primate_comparative_workflow_bundle_matches_packaged_expected_outputs",
            "test_run_primate_comparative_demo_materializes_dataset_and_workflow",
            "test_cli_demo_primate_comparative_json_output_reports_dataset_and_workflow",
        },
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_pgls_lambda_fit.py": {
            "test_run_pgls_estimated_lambda_matches_primate_reference_bundle",
        },
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_catarrhine_mitogenome_five_locus_panel_dataset.py": {
            "test_write_catarrhine_mitogenome_five_locus_panel_workflow_bundle_matches_packaged_expected_outputs",
            "test_run_catarrhine_mitogenome_five_locus_panel_demo_materializes_dataset_and_workflow",
            "test_cli_demo_catarrhine_mitogenome_five_locus_panel_json_output_reports_multilocus_review",
        },
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_pleistocene_bear_cytb_fragment_dataset.py": {
            "test_write_pleistocene_bear_cytb_fragment_workflow_bundle_matches_packaged_expected_outputs",
            "test_run_pleistocene_bear_cytb_fragment_demo_materializes_dataset_and_workflow",
            "test_cli_demo_pleistocene_bear_cytb_fragments_json_output_reports_missingness_review",
        },
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_central_european_seashore_flora_dataset.py": {
            "test_write_central_european_seashore_flora_workflow_bundle_matches_packaged_expected_outputs",
            "test_run_central_european_seashore_flora_demo_materializes_dataset_and_workflow",
            "test_cli_demo_central_european_seashore_flora_json_output_reports_dataset_and_workflow",
        },
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_rabies_cross_host_geography_panel_dataset.py": {
            "test_write_rabies_cross_host_geography_panel_workflow_bundle_matches_packaged_expected_outputs",
            "test_run_rabies_cross_host_geography_panel_demo_materializes_dataset_and_workflow",
            "test_cli_demo_rabies_cross_host_geography_panel_json_output_reports_integrated_workflow",
        },
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_influenza_a_ha_reference_dataset.py": {
            "test_write_influenza_a_ha_reference_workflow_bundle_matches_packaged_expected_outputs",
            "test_run_influenza_a_ha_reference_demo_materializes_dataset_and_workflow",
            "test_cli_demo_influenza_a_ha_reference_panel_json_output_reports_dataset_and_workflow",
        },
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_known_answer_reference_dataset.py": {
            "test_write_known_answer_reference_workflow_bundle_matches_packaged_expected_outputs",
            "test_run_known_answer_reference_demo_materializes_dataset_and_workflow",
            "test_public_runtime_exports_include_known_answer_reference_surface",
            "test_cli_demo_known_answer_reference_panel_json_output_reports_recovery_metrics",
        },
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_fasta_to_tree_workflow.py": {
            "test_run_fasta_to_tree_workflow_materializes_expected_outputs_for_three_datasets",
        },
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_evidence_workbench.py": {
            "test_refresh_evidence_book_writes_repo_index_outputs",
        },
    }

    for (
        module_path,
        expected_slow_functions,
    ) in expected_slow_functions_by_module.items():
        assert expected_slow_functions <= _slow_marked_functions(module_path)


def test_governed_stress_tests_stay_slow_marked() -> None:
    stress_modules = [
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_large_dataset_stress.py",
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "test_large_fasta_streaming.py",
    ]

    for module_path in stress_modules:
        assert _stress_marked_functions(module_path) <= _slow_marked_functions(
            module_path
        )


def test_long_running_maintainer_governance_surfaces_stay_slow_marked() -> None:
    expected_slow_functions_by_module = {
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics-dev"
        / "tests"
        / "test_evidence_cleanroom.py": {
            "test_repository_cleanroom_rerun_keeps_primate_longevity_selection_clean",
            "test_repository_selected_cleanroom_reruns_keep_governed_selections_clean",
        },
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics-dev"
        / "tests"
        / "test_publication_guard.py": {
            "test_assert_publishable_repository_allows_clean_package_bundle_repo",
            "test_assert_publishable_repository_allows_clean_publish_readiness_repo",
            "test_publication_guard_module_runs_without_runpy_warning",
        },
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics-dev"
        / "tests"
        / "test_publish_readiness.py": {
            "test_build_publish_readiness_report_exposes_repository_blockers",
            "test_build_publish_readiness_report_accepts_secured_runtime_dependency_contract",
        },
        REPO_ROOT
        / "packages"
        / "bijux-phylogenetics-dev"
        / "tests"
        / "test_package_bundles.py": {
            "test_check_package_bundles_builds_and_audits_publishable_packages",
            "test_check_package_bundles_rebuilds_cleanly_into_staged_output_directories",
            "test_check_package_bundles_keeps_published_bundle_outputs_unchanged",
        },
    }

    for (
        module_path,
        expected_slow_functions,
    ) in expected_slow_functions_by_module.items():
        assert expected_slow_functions <= _slow_marked_functions(module_path)


def test_repository_test_all_surface_disables_pytest_timeout_in_all_packages() -> None:
    root_make = (REPO_ROOT / "makes" / "root.mk").read_text(encoding="utf-8")
    dev_make = (
        REPO_ROOT / "makes" / "packages" / "bijux-phylogenetics-dev.mk"
    ).read_text(encoding="utf-8")
    alias_make = (REPO_ROOT / "makes" / "packages" / "phylogenetic.mk").read_text(
        encoding="utf-8"
    )

    assert "ROOT_PACKAGE_TARGETS += test-all test-all-plus-run-time" in root_make
    assert "ROOT_TARGET_PACKAGES_test-all := $(CHECK_PACKAGES)" in root_make
    assert "test-all: TEST_MAIN_ARGS =" in dev_make
    assert "test-all: PYTEST_ADDOPTS_EXTRA = -o timeout=0" in dev_make
    assert "test-all: PYTEST_ADDOPTS_EXTRA = -o timeout=0" in alias_make


def test_repository_test_all_plus_run_time_surface_disables_timeout_and_reports_durations_in_all_packages() -> (
    None
):
    root_make = (REPO_ROOT / "makes" / "root.mk").read_text(encoding="utf-8")
    dev_make = (
        REPO_ROOT / "makes" / "packages" / "bijux-phylogenetics-dev.mk"
    ).read_text(encoding="utf-8")
    alias_make = (REPO_ROOT / "makes" / "packages" / "phylogenetic.mk").read_text(
        encoding="utf-8"
    )

    assert "ROOT_PACKAGE_TARGETS += test-all test-all-plus-run-time" in root_make
    assert (
        "ROOT_TARGET_PACKAGES_test-all-plus-run-time := $(CHECK_PACKAGES)" in root_make
    )
    assert "test-all-plus-run-time: TEST_MAIN_ARGS =" in dev_make
    assert (
        "test-all-plus-run-time: PYTEST_ADDOPTS_EXTRA = -o timeout=0 --durations=0 --durations-min=0"
        in dev_make
    )
    assert (
        "test-all-plus-run-time: PYTEST_ADDOPTS_EXTRA = -o timeout=0 --durations=0 --durations-min=0"
        in alias_make
    )


def test_root_conftest_registers_markers_from_repository_pytest_config() -> None:
    conftest_path = (
        REPO_ROOT / "packages" / "bijux-phylogenetics" / "tests" / "conftest.py"
    )
    conftest_text = conftest_path.read_text(encoding="utf-8")

    assert 'PYTEST_CONFIG_PATH = REPO_ROOT / "configs" / "pytest.ini"' in conftest_text
    assert 'config.getini("markers")' in conftest_text
    assert 'config.addinivalue_line("markers", marker)' in conftest_text
