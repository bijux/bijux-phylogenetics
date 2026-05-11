PACKAGE_KIND := python
PACKAGE_IMPORT_NAME := phylogenetic
PACKAGE_INSTALL_SPEC := .[dev]
PACKAGE_INSTALL_PYTHON_PACKAGES = "$(MONOREPO_ROOT)/packages/bijux-phylogenetics[dev]"
RUFF_CONFIG = $(MONOREPO_ROOT)/configs/ruff.toml
BUILD_PRE_TARGETS := sync-license-assets-package
TEST_PATHS := tests
TEST_PATHS_UNIT := tests
TEST_SOURCE_PATHS := src
TEST_CLEAN_PATHS := "$(MONOREPO_ROOT)/.pytest_cache" "$(MONOREPO_ROOT)/.ruff_cache"
INTERROGATE_PATHS := src
QUALITY_PATHS := src tests
# Audit the installed environment so the local workspace dependency on
# bijux-phylogenetics is not misclassified as an external PyPI requirement.
ENABLE_PYDOCSTYLE := 1
SKIP_MYPY := 0
PACKAGE_ALL_TARGETS := clean install test lint quality security build sbom

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

test-all: test
.PHONY: test-all

include $(abspath $(dir $(firstword $(MAKEFILE_LIST))))/../bijux-py/package.mk
