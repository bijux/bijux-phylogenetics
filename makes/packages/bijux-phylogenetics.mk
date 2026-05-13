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
TEST_PATHS_EVALUATION := $(MONOREPO_ROOT)/packages/bijux-phylogenetics/tests
TEST_REAL_LOCAL_PATH := $(MONOREPO_ROOT)/packages/bijux-phylogenetics/tests/real_local
TEST_MAIN_ARGS = -m "not slow and not real_local and not evaluation"
TEST_UNIT_DIR_ARGS = -m "not slow and not real_local and not evaluation" --maxfail=1 -q
TEST_UNIT_FALLBACK_ARGS = -k "not e2e and not integration and not functional" -m "not slow and not real_local and not evaluation" --maxfail=1 -q
TEST_EVALUATION_ARGS = -m "evaluation and scientific_validation" -s -p no:cov
TEST_REAL_LOCAL_ARGS = -m "real_local and engine_real and not scientific_validation" -s -p no:cov
TEST_CLEAN_PATHS := "$(MONOREPO_ROOT)/.pytest_cache" "$(MONOREPO_ROOT)/.ruff_cache"
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

quality-compileall:
	@"$(VENV_PYTHON)" -m compileall src | tee "$(PROJECT_ARTIFACTS_DIR)/quality/compileall.log"
.PHONY: quality-compileall

sync-license-assets-package:
	@"$(VENV_PYTHON)" -m bijux_phylogenetics_dev.release.license_assets sync
.PHONY: sync-license-assets-package

build-install-smoke:
	@tmp_root="$(PROJECT_ARTIFACTS_DIR)/tmp/build-install-smoke"; \
	dist_name="$$(printf '%s' "$(BUILD_PACKAGE_NAME)" | tr '-' '_')"; \
	wheel_path="$$(ls -1t "$(BUILD_DIR_ABS)/$${dist_name}"-*.whl | head -n 1)"; \
	sdist_path="$$(ls -1t "$(BUILD_DIR_ABS)/$${dist_name}"-*.tar.gz | head -n 1)"; \
	export PIP_DISABLE_PIP_VERSION_CHECK=1; \
	if [ -z "$$wheel_path" ] || [ -z "$$sdist_path" ]; then \
	  echo "✘ Missing build artifacts for $(BUILD_PACKAGE_NAME) in $(BUILD_DIR_ABS)"; \
	  exit 1; \
	fi; \
	rm -rf "$$tmp_root"; \
	"$(BUILD_PYTHON)" -m venv "$$tmp_root/smoke"; \
	"$$tmp_root/smoke/bin/python" -m pip install "$$wheel_path"; \
	"$$tmp_root/smoke/bin/bijux-phylogenetics" --version; \
	"$$tmp_root/smoke/bin/bijux-phylogenetics" --help >/dev/null; \
	"$$tmp_root/smoke/bin/python" -m pip uninstall -y "$(BUILD_PACKAGE_NAME)" >/dev/null; \
	"$$tmp_root/smoke/bin/python" -m pip install "$$sdist_path"; \
	"$$tmp_root/smoke/bin/bijux-phylogenetics" --version; \
	"$$tmp_root/smoke/bin/bijux-phylogenetics" --help >/dev/null
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
	@echo "→ Running slow scientific validation tests (manual only)"
	@$(PYTEST) $(PYTEST_INFO_FLAGS) --version
	@mkdir -p "$(TEST_ARTIFACTS_DIR)" "$(HYPOTHESIS_DB_DIR)" "$(TMP_DIR)"
	@( cd "$(PYTEST_ROOTDIR_ABS)" && \
	  PYTHONPATH="$(TEST_SOURCE_PATH_ABS)$${PYTHONPATH:+:$${PYTHONPATH}}" \
	  PYTHONDONTWRITEBYTECODE=1 \
	  HYPOTHESIS_DATABASE_DIRECTORY="$(HYPOTHESIS_DB_ABS)" \
	  $(TEST_PYCACHE_ENV) \
	  sh -c '$(PYTEST) --rootdir "$(PYTEST_ROOTDIR_ABS)" -c "$(PYTEST_INI_ABS)" "packages/bijux-phylogenetics/tests" -m "evaluation and scientific_validation and slow" -o addopts= -o timeout=600 -s -p no:cov' )
	@echo "✔ slow scientific validation lane completed"
.PHONY: scientific-validation-slow

stress-small-lane:
	@echo "→ Running governed small stress tier"
	@$(PYTEST) $(PYTEST_INFO_FLAGS) --version
	@mkdir -p "$(TEST_ARTIFACTS_DIR)" "$(HYPOTHESIS_DB_DIR)" "$(TMP_DIR)"
	@( cd "$(PYTEST_ROOTDIR_ABS)" && \
	  PYTHONPATH="$(TEST_SOURCE_PATH_ABS)$${PYTHONPATH:+:$${PYTHONPATH}}" \
	  PYTHONDONTWRITEBYTECODE=1 \
	  HYPOTHESIS_DATABASE_DIRECTORY="$(HYPOTHESIS_DB_ABS)" \
	  $(TEST_PYCACHE_ENV) \
	  sh -c '$(PYTEST) --rootdir "$(PYTEST_ROOTDIR_ABS)" -c "$(PYTEST_INI_ABS)" "packages/bijux-phylogenetics/tests/test_large_dataset_stress.py" -m "evaluation and stress_small" -o addopts= -s -p no:cov' )
	@echo "✔ small stress tier completed"
.PHONY: stress-small-lane

stress-heavy-lane:
	@echo "→ Running governed heavy stress tier"
	@$(PYTEST) $(PYTEST_INFO_FLAGS) --version
	@mkdir -p "$(TEST_ARTIFACTS_DIR)" "$(HYPOTHESIS_DB_DIR)" "$(TMP_DIR)"
	@( cd "$(PYTEST_ROOTDIR_ABS)" && \
	  PYTHONPATH="$(TEST_SOURCE_PATH_ABS)$${PYTHONPATH:+:$${PYTHONPATH}}" \
	  PYTHONDONTWRITEBYTECODE=1 \
	  HYPOTHESIS_DATABASE_DIRECTORY="$(HYPOTHESIS_DB_ABS)" \
	  $(TEST_PYCACHE_ENV) \
	  sh -c '$(PYTEST) --rootdir "$(PYTEST_ROOTDIR_ABS)" -c "$(PYTEST_INI_ABS)" "packages/bijux-phylogenetics/tests/test_large_dataset_stress.py" -m "evaluation and stress_heavy" -o addopts= -o timeout=600 -s -p no:cov' )
	@echo "✔ heavy stress tier completed"
.PHONY: stress-heavy-lane

include $(abspath $(dir $(firstword $(MAKEFILE_LIST))))/../bijux-py/package.mk
