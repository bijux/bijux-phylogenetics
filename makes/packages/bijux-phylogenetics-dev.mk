PACKAGE_KIND := python
PACKAGE_IMPORT_NAME := bijux_phylogenetics_dev
PACKAGE_INSTALL_SPEC := .[dev]
PACKAGE_INSTALL_PYTHON_PACKAGES = $(MONOREPO_ROOT)/packages/bijux-phylogenetics
RUFF_CONFIG = $(MONOREPO_ROOT)/configs/ruff.toml
BUILD_PRE_TARGETS := sync-license-assets-package
TEST_PATHS := tests
TEST_PATHS_UNIT := tests
TEST_SOURCE_PATHS := src ../bijux-phylogenetics/src
TEST_MAIN_ARGS = -m "not slow"
TEST_CLEAN_PATHS = "$(MONOREPO_ROOT)/.pytest_cache" "$(MONOREPO_ROOT)/.ruff_cache"
INTERROGATE_PATHS := src
QUALITY_PATHS := src tests
QUALITY_MYPY_TARGETS := src ../bijux-phylogenetics/src
SECURITY_AUDIT_PREPARE_MODE = pyproject
PIP_AUDIT_INPUTS = -r "$(SECURITY_REQS)"
TEST_PRE_TARGETS := sync-license-assets-package
LINT_PRE_TARGETS := sync-license-assets-package
QUALITY_PRE_TARGETS := sync-license-assets-package
SECURITY_EXTRA_TARGETS := sync-license-assets-package
ENABLE_PYDOCSTYLE := 1
SKIP_MYPY := 0
PACKAGE_ALL_TARGETS := clean install test lint quality security build sbom

sync-license-assets-package:
	@"$(VENV_PYTHON)" -m bijux_phylogenetics_dev.release.license_assets sync
.PHONY: sync-license-assets-package

test-all: TEST_MAIN_ARGS =
test-all: PYTEST_ADDOPTS_EXTRA = -o timeout=0
test-all: test
.PHONY: test-all

test-all-plus-run-time: TEST_MAIN_ARGS =
test-all-plus-run-time: PYTEST_ADDOPTS_EXTRA = -o timeout=0 --durations=0 --durations-min=0
test-all-plus-run-time: test
.PHONY: test-all-plus-run-time

include $(abspath $(dir $(firstword $(MAKEFILE_LIST))))/../bijux-py/package.mk
