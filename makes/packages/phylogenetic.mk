PACKAGE_KIND := python
PACKAGE_IMPORT_NAME := phylogenetic
PACKAGE_INSTALL_SPEC := .[dev]
PACKAGE_INSTALL_PYTHON_PACKAGES = "$(MONOREPO_ROOT)/packages/bijux-phylogenetics[dev]"
RUFF_CONFIG = $(MONOREPO_ROOT)/configs/ruff.toml
BUILD_PRE_TARGETS := sync-license-assets-package
TEST_PATHS := tests
TEST_PATHS_UNIT := tests
TEST_SOURCE_PATHS := src
TEST_CLEAN_PATHS = "$(MONOREPO_ROOT)/.pytest_cache" "$(MONOREPO_ROOT)/.ruff_cache"
INTERROGATE_PATHS := src
QUALITY_PATHS := src tests
# Audit the installed environment so the local workspace dependency on
# bijux-phylogenetics is not misclassified as an external PyPI requirement.
ENABLE_PYDOCSTYLE := 1
SKIP_MYPY := 0
PACKAGE_ALL_TARGETS := clean install test lint quality security build sbom
FROZEN_TEST_REPORT_PYTHON = $(MONOREPO_ROOT)/artifacts/root/check-venv/bin/python
FULL_SUITE_PYTEST_WORKERS ?= 8
FULL_SUITE_PYTEST_DISTRIBUTION ?= worksteal

sync-license-assets-package:
	@for file_name in LICENSE NOTICE; do \
	  target_path="$(PROJECT_DIR)/$$file_name"; \
	  expected_target="../../$$file_name"; \
	  if [ -L "$$target_path" ] && [ "$$(readlink "$$target_path")" = "$$expected_target" ]; then \
	    continue; \
	  fi; \
	  if [ -L "$$target_path" ] || [ -e "$$target_path" ]; then \
	    rm -f "$$target_path"; \
	  fi; \
	  ln -s "$$expected_target" "$$target_path"; \
	done
.PHONY: sync-license-assets-package

test-all: PYTEST_ADDOPTS_EXTRA = -o timeout=0 -n $(FULL_SUITE_PYTEST_WORKERS) --dist $(FULL_SUITE_PYTEST_DISTRIBUTION)
test-all: test
.PHONY: test-all

test-all-plus-run-time: PYTEST_ADDOPTS_EXTRA = -o timeout=0 -n $(FULL_SUITE_PYTEST_WORKERS) --dist $(FULL_SUITE_PYTEST_DISTRIBUTION) --durations=0 --durations-min=0
test-all-plus-run-time: test
.PHONY: test-all-plus-run-time

test-all-frozen: ## Compare the full-suite output against TEST_ALL_FROZEN_REF
	@test -n "$(TEST_ALL_FROZEN_REF)" || { echo "TEST_ALL_FROZEN_REF is required"; exit 2; }
	@$(MAKE) -C "$(MONOREPO_ROOT)" root-check-env >/dev/null
	@"$(FROZEN_TEST_REPORT_PYTHON)" -m bijux_phylogenetics_dev.quality.frozen_test_report check --repo-root "$(MONOREPO_ROOT)" --package "$(PROJECT_SLUG)" --frozen-ref "$(TEST_ALL_FROZEN_REF)" --artifacts-root "$(PROJECT_ARTIFACTS_DIR)/test-frozen" --json-out "$(PROJECT_ARTIFACTS_DIR)/test-frozen/report.json"
.PHONY: test-all-frozen

include $(abspath $(dir $(firstword $(MAKEFILE_LIST))))/../bijux-py/package.mk
