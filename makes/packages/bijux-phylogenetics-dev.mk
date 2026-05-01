PACKAGE_KIND := python
PACKAGE_IMPORT_NAME := bijux_phylogenetics_dev
PACKAGE_INSTALL_SPEC := .[dev]
RUFF_CONFIG = $(MONOREPO_ROOT)/configs/ruff.toml
BUILD_PRE_TARGETS := sync-license-assets-package
TEST_PATHS := tests
TEST_PATHS_UNIT := tests
TEST_SOURCE_PATHS := src
INTERROGATE_PATHS := src
QUALITY_PATHS := src tests
SECURITY_AUDIT_PREPARE_MODE = pyproject
PIP_AUDIT_INPUTS = -r "$(SECURITY_REQS)"
ENABLE_PYDOCSTYLE := 1
SKIP_MYPY := 0
PACKAGE_ALL_TARGETS := clean install test lint quality security build sbom

sync-license-assets-package:
	@for file_name in LICENSE NOTICE; do \
	  source_path="$(MONOREPO_ROOT)/$$file_name"; \
	  target_path="$(PROJECT_DIR)/$$file_name"; \
	  if [ -L "$$target_path" ] || [ ! -f "$$target_path" ] || ! cmp -s "$$source_path" "$$target_path"; then \
	    rm -f "$$target_path"; \
	    cp "$$source_path" "$$target_path"; \
	  fi; \
	done
.PHONY: sync-license-assets-package

include $(abspath $(dir $(firstword $(MAKEFILE_LIST))))/../bijux-py/package.mk
