PACKAGE_KIND := repository-python
PACKAGE_IMPORT_NAME := bijux_phylogenetics
PACKAGE_INSTALL_PYTHON_PACKAGES = "$(MONOREPO_ROOT)/packages/bijux-phylogenetics-dev[dev]"
LINT_DIRS = $(MONOREPO_ROOT)/packages/bijux-phylogenetics/src $(MONOREPO_ROOT)/packages/bijux-phylogenetics/tests
MYPY_TARGETS = src
ENABLE_MYPY := 1
ENABLE_CODESPELL := 0
ENABLE_RADON := 0
ENABLE_PYDOCSTYLE := 0
TEST_PATHS := tests
TEST_PATHS_UNIT := tests
TEST_PATHS_EVALUATION = $(MONOREPO_ROOT)/packages/bijux-phylogenetics/tests
TEST_REAL_LOCAL_PATH = $(MONOREPO_ROOT)/packages/bijux-phylogenetics/tests/real_local
TEST_MAIN_ARGS = -m "not slow and not real_local and not evaluation" --maxfail=1 -q
TEST_UNIT_DIR_ARGS = -m "not slow and not real_local and not evaluation" --maxfail=1 -q
TEST_UNIT_FALLBACK_ARGS = -k "not e2e and not integration and not functional" -m "not slow and not real_local and not evaluation" --maxfail=1 -q
TEST_EVALUATION_ARGS = -m "evaluation and scientific_validation and not slow" -s -p no:cov
TEST_REAL_LOCAL_ARGS = -m "real_local and engine_real and not scientific_validation and not slow" -s -p no:cov
TEST_CLEAN_PATHS = "$(MONOREPO_ROOT)/.pytest_cache" "$(MONOREPO_ROOT)/.ruff_cache"
QUALITY_PATHS = src tests
MYPY_CONFIG = $(MONOREPO_ROOT)/configs/mypy.ini
QUALITY_MYPY_CONFIG = $(MONOREPO_ROOT)/configs/mypy.ini
QUALITY_MYPY_TARGETS = src
SECURITY_AUDIT_PREPARE_MODE = pyproject
PIP_AUDIT_INPUTS = -r "$(SECURITY_REQS)"
TEST_PRE_TARGETS := sync-license-assets-package
LINT_PRE_TARGETS := sync-license-assets-package
QUALITY_PRE_TARGETS := sync-license-assets-package
SECURITY_EXTRA_TARGETS := sync-license-assets-package
QUALITY_POST_TARGETS := quality-compileall
BUILD_PACKAGE_NAME := bijux-phylogenetics
BUILD_PRE_TARGETS := sync-license-assets-package
BUILD_TEMP_CLEAN_PATHS := build dist *.egg-info
BUILD_POST_TARGETS := build-install-smoke
PACKAGE_NAME := bijux-phylogenetics
SBOM_REQUIREMENTS_WRITER := -m bijux_phylogenetics_dev.quality.requirements_writer

quality-compileall:
	@"$(VENV_PYTHON)" -m compileall src | tee "$(PROJECT_ARTIFACTS_DIR)/quality/compileall.log"
.PHONY: quality-compileall

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

build-install-smoke:
	@rm -rf "$(PROJECT_ARTIFACTS_DIR)/tmp/build-install-smoke"
	@PYTHONPATH="$(MONOREPO_ROOT)/packages/bijux-phylogenetics-dev/src$${PYTHONPATH:+:$${PYTHONPATH}}" \
	"$(VENV_PYTHON)" -m bijux_phylogenetics_dev.quality.package_install_smoke \
		--repo-root "$(MONOREPO_ROOT)" \
		--dist-dir "$(BUILD_DIR_ABS)" \
		--artifacts-root "$(PROJECT_ARTIFACTS_DIR)/tmp/build-install-smoke" \
		--build-python "$(VENV_PYTHON)" \
		--artifact-kind both
.PHONY: build-install-smoke

external-engine-lane:
	@$(SELF_MAKE) real-local
	@$(SELF_MAKE) scientific-validation-lane
	@echo "✔ external engine execution lanes completed"
.PHONY: external-engine-lane

scientific-validation-lane:
	@echo "→ Running scientific validation tests (manual only)"
	@$(PYTEST) $(PYTEST_INFO_FLAGS) --version
	@mkdir -p "$(TEST_ARTIFACTS_DIR)" "$(HYPOTHESIS_DB_DIR)" "$(TMP_DIR)"
	@( cd "$(PYTEST_ROOTDIR_ABS)" && \
	  PYTHONPATH="$(TEST_SOURCE_PATH_ABS)$${PYTHONPATH:+:$${PYTHONPATH}}" \
	  PYTHONDONTWRITEBYTECODE=1 \
	  HYPOTHESIS_DATABASE_DIRECTORY="$(HYPOTHESIS_DB_ABS)" \
	  $(TEST_PYCACHE_ENV) \
	  sh -c '$(PYTEST) --rootdir "$(PYTEST_ROOTDIR_ABS)" -c "$(PYTEST_INI_ABS)" "packages/bijux-phylogenetics/tests" -m "evaluation and scientific_validation and not slow" -o addopts= -s -p no:cov' )
	@echo "✔ scientific validation lane completed"
.PHONY: scientific-validation-lane

scientific-validation-slow:
	@echo "✘ slow tests are reserved for 'make test-all' or 'make test-all-plus-run-time'"
	@exit 2
.PHONY: scientific-validation-slow

stress-small-lane:
	@echo "✘ slow tests are reserved for 'make test-all' or 'make test-all-plus-run-time'"
	@exit 2
.PHONY: stress-small-lane

stress-heavy-lane:
	@echo "✘ slow tests are reserved for 'make test-all' or 'make test-all-plus-run-time'"
	@exit 2
.PHONY: stress-heavy-lane

include $(abspath $(dir $(firstword $(MAKEFILE_LIST))))/../bijux-py/package.mk
