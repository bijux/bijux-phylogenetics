ROOT_MAKEFILE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))

include $(ROOT_MAKEFILE_DIR)/bijux-py/root/env.mk
include $(ROOT_MAKEFILE_DIR)/env.mk
include $(ROOT_MAKEFILE_DIR)/packages.mk

ROOT_DEV_PYTHONPATH := $(CURDIR)/packages/bijux-phylogenetics-dev/src
BIJUX_PY_SYSTEM_REL ?= .bijux/shared/bijux-makes-py
ROOT_CHECK_VENV := $(ROOT_ARTIFACTS_DIR)/check-venv
ROOT_DOCS_DEV_ADDR ?= 127.0.0.1:8000
UV_SYNC := UV_PROJECT_ENVIRONMENT="$(ROOT_CHECK_VENV)" $(UV) sync --frozen --group dev --python "$(PYTHON)"
CLI := $(ROOT_CHECK_VENV)/bin/bijux-phylogenetics
DEV_RUN = PYTHONPATH="$(CURDIR)/packages/bijux-phylogenetics-dev/src$${PYTHONPATH:+:$$PYTHONPATH}" "$(ROOT_CHECK_PYTHON)"
DOCS_RENDER_SERVE_CONFIG := 0

include $(ROOT_MAKEFILE_DIR)/bijux-py/repository/root.mk

include $(ROOT_MAKEFILE_DIR)/bijux-py/root/package-dispatch.mk
include $(ROOT_MAKEFILE_DIR)/bijux-py/root/docs.mk
include $(ROOT_MAKEFILE_DIR)/bijux-docs.mk
include $(ROOT_MAKEFILE_DIR)/bijux-std.mk
include $(ROOT_MAKEFILE_DIR)/bijux-py/repository/config-layout.mk
include $(ROOT_MAKEFILE_DIR)/bijux-py/repository/make-layout.mk
include $(ROOT_MAKEFILE_DIR)/bijux-py/bijux.mk

DOCS_BUILD_PREPARE_TARGETS := bijux-docs-sync docs-prepare-source
DOCS_CHECK_PREPARE_TARGETS := bijux-docs-sync docs-prepare-source
DOCS_SERVE_PREPARE_TARGETS := bijux-docs-sync docs-render-serve-config

.PHONY: \
	help list list-all install lock lock-check lint quality security test docs docs-check docs-serve api build sbom clean all \
	check package-check package-smoke package-source-smoke package-verify sync-badges sync-license-assets test-goldens demo \
	clean-root-artifacts root-check-env check-shared-bijux-py

check: sync-license-assets lock-check lint test quality security docs build sbom ## Run the full repository verification flow

sync-badges: root-check-env ## Render shared badge blocks into managed README surfaces
	@$(DEV_RUN) -m bijux_phylogenetics_dev.docs.badge_sync sync

check-badges: root-check-env ## Verify README badge blocks match docs/badges.md
	@$(DEV_RUN) -m bijux_phylogenetics_dev.docs.badge_sync check

sync-license-assets: root-check-env ## Sync package LICENSE and NOTICE links from root sources
	@$(DEV_RUN) -m bijux_phylogenetics_dev.release.license_assets sync

package-check: build ## Validate built distributions with twine
	@"$(ROOT_CHECK_PYTHON)" -m twine check "$(CURDIR)/artifacts/bijux-phylogenetics/build"/*.whl "$(CURDIR)/artifacts/bijux-phylogenetics/build"/*.tar.gz
.PHONY: package-check

package-smoke: build ## Install the wheel into a temp venv and run the CLI
	@rm -rf "$(CURDIR)/artifacts/tmp/package-smoke"
	@"$(UV)" venv --python "$(PYTHON)" "$(CURDIR)/artifacts/tmp/package-smoke"
	@"$(UV)" pip install --python "$(CURDIR)/artifacts/tmp/package-smoke/bin/python" "$(CURDIR)/artifacts/bijux-phylogenetics/build"/*.whl
	@"$(CURDIR)/artifacts/tmp/package-smoke/bin/bijux-phylogenetics" --version
	@"$(CURDIR)/artifacts/tmp/package-smoke/bin/bijux-phylogenetics" --help >/dev/null
.PHONY: package-smoke

package-source-smoke: build ## Install the sdist into a temp venv and run the CLI
	@rm -rf "$(CURDIR)/artifacts/tmp/package-source-smoke"
	@"$(UV)" venv --python "$(PYTHON)" "$(CURDIR)/artifacts/tmp/package-source-smoke"
	@"$(UV)" pip install --python "$(CURDIR)/artifacts/tmp/package-source-smoke/bin/python" "$(CURDIR)/artifacts/bijux-phylogenetics/build"/*.tar.gz
	@"$(CURDIR)/artifacts/tmp/package-source-smoke/bin/bijux-phylogenetics" --version
	@"$(CURDIR)/artifacts/tmp/package-source-smoke/bin/bijux-phylogenetics" --help >/dev/null
.PHONY: package-source-smoke

package-verify: package-check package-smoke package-source-smoke ## Run the full packaging proof surface
.PHONY: package-verify

test-goldens: root-check-env ## Compare checked-in golden outputs
	@"$(ROOT_CHECK_VENV)/bin/python" -m pytest packages/bijux-phylogenetics/tests/test_goldens.py
.PHONY: test-goldens

demo: root-check-env ## Run the public capability demo into repository artifacts
	@"$(CLI)" demo run --out "$(CURDIR)/artifacts/bijux-phylogenetics/demo"
.PHONY: demo

HELP_WIDTH := 22
include $(ROOT_MAKEFILE_DIR)/bijux-py/ci/help.mk

##@ Repository
help: ## Show generated repository commands from included make modules
check-shared-bijux-py: ## Verify shared bijux-py make modules match across sibling repositories
check-config-layout: ## Validate the repository config tree shape and required tool configs
check-make-layout: ## Validate the repository make tree shape and required entrypoints
