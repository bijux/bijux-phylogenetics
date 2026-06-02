ROOT_MAKEFILE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))

include $(ROOT_MAKEFILE_DIR)/bijux-py/root/env.mk
include $(ROOT_MAKEFILE_DIR)/env.mk
include $(ROOT_MAKEFILE_DIR)/packages.mk

ROOT_DEV_PYTHONPATH := $(CURDIR)/packages/bijux-phylogenetics-dev/src:$(CURDIR)/packages/bijux-phylogenetics/src
BIJUX_PY_SYSTEM_REL ?= .bijux/shared/bijux-makes-py
ROOT_CHECK_VENV := $(ROOT_ARTIFACTS_DIR)/check-venv
ROOT_DOCS_DEV_ADDR ?= 127.0.0.1:8000
UV_SYNC := UV_PROJECT_ENVIRONMENT="$(ROOT_CHECK_VENV)" $(UV) sync --frozen --group dev --python "$(PYTHON)"
CLI := $(ROOT_CHECK_VENV)/bin/bijux-phylogenetics
DEV_RUN = PYTHONPATH="$(ROOT_DEV_PYTHONPATH)$${PYTHONPATH:+:$$PYTHONPATH}" "$(ROOT_CHECK_PYTHON)"
DOCS_RENDER_SERVE_CONFIG := 0
ROOT_PACKAGE_TARGETS += test-all test-all-plus-run-time
ROOT_TARGET_GROUPS_test-all ?= check
ROOT_TARGET_GROUPS_test-all-plus-run-time ?= check
ROOT_TARGET_SHARED_ENV_test-all ?= 1
ROOT_TARGET_SHARED_ENV_test-all-plus-run-time ?= 1

include $(ROOT_MAKEFILE_DIR)/bijux-py/repository/root.mk

ROOT_FORBIDDEN_ARTIFACTS := $(filter-out \
	"$(CURDIR)/.hypothesis" \
	"$(CURDIR)/.benchmarks",$(ROOT_FORBIDDEN_ARTIFACTS))

include $(ROOT_MAKEFILE_DIR)/bijux-py/root/package-dispatch.mk
ROOT_TARGET_PACKAGES_test-all := $(CHECK_PACKAGES)
ROOT_TARGET_PACKAGES_test-all-plus-run-time := $(CHECK_PACKAGES)
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
	help list list-all install lock lock-check lint quality security test test-all test-all-plus-run-time docs docs-check docs-serve api build sbom clean all \
	check package-check package-smoke package-source-smoke package-verify sync-badges sync-license-assets test-goldens demo \
	install-external-engine-runtime test-external-engines \
	clean-root-artifacts root-check-env check-shared-bijux-py check-config-ssot \
	list-evidence-studies build-evidence-book build-evidence-study build-evidence-unit validate-evidence-book \
	sync-evidence-artifacts sync-evidence-unit-artifacts check-evidence-artifacts check-evidence-unit-artifacts sync-evidence-unit-inputs check-evidence-unit-inputs report-evidence-completeness check-evidence-completeness report-evidence-governance check-evidence-governance \
	report-artifact-governance check-artifact-governance report-execution-surfaces check-execution-surfaces \
	report-package-boundaries check-package-boundaries report-package-bundles check-package-bundles report-publish-readiness check-publish-readiness \
	report-release-readiness check-release-readiness

EVIDENCE_STUDY_ID ?=
EVIDENCE_ID ?=
EVIDENCE_IDS ?=

check: sync-license-assets lock-check check-config-ssot check-evidence-governance check-execution-surfaces check-package-boundaries lint test quality security docs build sbom ## Run the full repository verification flow

test-all: ## Run every repository test surface, including slow, evaluation, and real-local tests
test-all-plus-run-time: ## Run every repository test surface and report per-test durations

install-external-engine-runtime: ## Install the local external engine lane dependencies on macOS with Homebrew
	@command -v brew >/dev/null || { echo "Homebrew is required for install-external-engine-runtime"; exit 2; }
	@brew install mafft trimal mrbayes brewsci/bio/fasttree brewsci/bio/iqtree2
	@brew install --cask beast2
.PHONY: install-external-engine-runtime

test-external-engines: root-check-env ## Run real external engine integration and scientific validation lanes
	@$(MAKE) -f "$(CURDIR)/makes/packages/bijux-phylogenetics.mk" -C "$(CURDIR)/packages/bijux-phylogenetics" external-engine-lane
.PHONY: test-external-engines

test-scientific-validation-slow: root-check-env ## Run the slow governed external-engine validation lane
	@echo "✘ slow tests are reserved for 'make test-all' or 'make test-all-plus-run-time'"
	@exit 2
.PHONY: test-scientific-validation-slow

test-stress-small: root-check-env ## Run the governed small large-dataset stress tier
	@echo "✘ slow tests are reserved for 'make test-all' or 'make test-all-plus-run-time'"
	@exit 2
.PHONY: test-stress-small

test-stress-heavy: root-check-env ## Run the governed heavy large-dataset stress tier
	@echo "✘ slow tests are reserved for 'make test-all' or 'make test-all-plus-run-time'"
	@exit 2
.PHONY: test-stress-heavy

sync-badges: root-check-env ## Render shared badge blocks into managed README surfaces
	@$(DEV_RUN) -m bijux_phylogenetics_dev.docs.badge_sync sync

check-badges: root-check-env ## Verify README badge blocks match docs/badges.md
	@$(DEV_RUN) -m bijux_phylogenetics_dev.docs.badge_sync check

sync-license-assets: root-check-env ## Sync package LICENSE and NOTICE links from root sources
	@$(DEV_RUN) -m bijux_phylogenetics_dev.release.license_assets sync

check-config-ssot: root-check-env ## Validate repository-owned config source-of-truth rules
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.config_ssot check --repo-root "$(CURDIR)" --json-out "$(ROOT_ARTIFACTS_DIR)/config-ssot-audit.json"
.PHONY: check-config-ssot

list-evidence-studies: root-check-env ## Report registered evidence studies into repository artifacts
	@"$(CLI)" evidence book studies --json > "$(ROOT_ARTIFACTS_DIR)/evidence-studies.json"
.PHONY: list-evidence-studies

build-evidence-book: root-check-env ## Refresh the full evidence-book and store the governed build report
	@"$(CLI)" evidence book build --json > "$(ROOT_ARTIFACTS_DIR)/evidence-book-build.json"
.PHONY: build-evidence-book

build-evidence-study: root-check-env ## Refresh one governed evidence study selected by EVIDENCE_STUDY_ID
	@test -n "$(EVIDENCE_STUDY_ID)" || { echo "EVIDENCE_STUDY_ID is required"; exit 2; }
	@"$(CLI)" evidence book build "$(EVIDENCE_STUDY_ID)" --json > "$(ROOT_ARTIFACTS_DIR)/evidence-book-build.json"
.PHONY: build-evidence-study

build-evidence-unit: root-check-env ## Refresh one governed Evidence ID selected by EVIDENCE_STUDY_ID and EVIDENCE_ID
	@test -n "$(EVIDENCE_STUDY_ID)" || { echo "EVIDENCE_STUDY_ID is required"; exit 2; }
	@test -n "$(EVIDENCE_ID)" || { echo "EVIDENCE_ID is required"; exit 2; }
	@"$(CLI)" evidence book build "$(EVIDENCE_STUDY_ID)" --evidence-id "$(EVIDENCE_ID)" --json > "$(ROOT_ARTIFACTS_DIR)/evidence-unit-build.json"
.PHONY: build-evidence-unit

validate-evidence-book: root-check-env ## Validate the governed evidence-book structure and index surfaces
	@"$(CLI)" evidence book validate --json > "$(ROOT_ARTIFACTS_DIR)/evidence-book-validation.json"
.PHONY: validate-evidence-book

sync-evidence-artifacts: root-check-env ## Render governed local artifact surfaces for every evidence bundle
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.evidence_artifacts sync --repo-root "$(CURDIR)"
.PHONY: sync-evidence-artifacts

sync-evidence-unit-artifacts: root-check-env ## Render governed local artifact surfaces for one Evidence ID
	@test -n "$(EVIDENCE_STUDY_ID)" || { echo "EVIDENCE_STUDY_ID is required"; exit 2; }
	@test -n "$(EVIDENCE_ID)" || { echo "EVIDENCE_ID is required"; exit 2; }
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.evidence_artifacts sync --repo-root "$(CURDIR)" --study-id "$(EVIDENCE_STUDY_ID)" --evidence-id "$(EVIDENCE_ID)" > "$(ROOT_ARTIFACTS_DIR)/evidence-unit-artifacts.json"
.PHONY: sync-evidence-unit-artifacts

check-evidence-artifacts: root-check-env ## Validate governed local artifact surfaces for every evidence bundle
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.evidence_artifacts check --repo-root "$(CURDIR)"
.PHONY: check-evidence-artifacts

check-evidence-unit-artifacts: root-check-env ## Validate governed local artifact surfaces for one Evidence ID
	@test -n "$(EVIDENCE_STUDY_ID)" || { echo "EVIDENCE_STUDY_ID is required"; exit 2; }
	@test -n "$(EVIDENCE_ID)" || { echo "EVIDENCE_ID is required"; exit 2; }
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.evidence_artifacts check --repo-root "$(CURDIR)" --study-id "$(EVIDENCE_STUDY_ID)" --evidence-id "$(EVIDENCE_ID)" > "$(ROOT_ARTIFACTS_DIR)/evidence-unit-artifacts.json"
.PHONY: check-evidence-unit-artifacts

sync-evidence-unit-inputs: root-check-env ## Render governed input manifests for one Evidence ID
	@test -n "$(EVIDENCE_STUDY_ID)" || { echo "EVIDENCE_STUDY_ID is required"; exit 2; }
	@test -n "$(EVIDENCE_ID)" || { echo "EVIDENCE_ID is required"; exit 2; }
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.evidence_inputs sync --repo-root "$(CURDIR)" --study-id "$(EVIDENCE_STUDY_ID)" --evidence-id "$(EVIDENCE_ID)" > "$(ROOT_ARTIFACTS_DIR)/evidence-unit-inputs.json"
.PHONY: sync-evidence-unit-inputs

check-evidence-unit-inputs: root-check-env ## Validate governed input manifests for one Evidence ID
	@test -n "$(EVIDENCE_STUDY_ID)" || { echo "EVIDENCE_STUDY_ID is required"; exit 2; }
	@test -n "$(EVIDENCE_ID)" || { echo "EVIDENCE_ID is required"; exit 2; }
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.evidence_inputs check --repo-root "$(CURDIR)" --study-id "$(EVIDENCE_STUDY_ID)" --evidence-id "$(EVIDENCE_ID)" > "$(ROOT_ARTIFACTS_DIR)/evidence-unit-inputs.json"
.PHONY: check-evidence-unit-inputs

report-evidence-completeness: root-check-env ## Audit evidence bundle completeness into repository artifacts
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.evidence_completeness report --repo-root "$(CURDIR)" --json-out "$(ROOT_ARTIFACTS_DIR)/evidence-completeness.json"
.PHONY: report-evidence-completeness

check-evidence-completeness: root-check-env ## Fail when any evidence bundle is structurally incomplete
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.evidence_completeness check --repo-root "$(CURDIR)" --json-out "$(ROOT_ARTIFACTS_DIR)/evidence-completeness.json"
.PHONY: check-evidence-completeness

report-evidence-governance: root-check-env ## Build evidence governance report surfaces without package-quality checks
	@$(MAKE) validate-evidence-book
	@$(MAKE) check-evidence-artifacts
	@$(MAKE) report-evidence-completeness
	@$(MAKE) report-artifact-governance
.PHONY: report-evidence-governance

check-evidence-governance: root-check-env ## Enforce evidence-only governance checks
	@$(MAKE) validate-evidence-book
	@$(MAKE) check-evidence-artifacts
	@$(MAKE) check-evidence-completeness
	@$(MAKE) check-artifact-governance
.PHONY: check-evidence-governance

report-artifact-governance: root-check-env ## Audit tox, make, and workflow artifact output discipline
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.artifact_governance report --repo-root "$(CURDIR)" --json-out "$(ROOT_ARTIFACTS_DIR)/artifact-governance.json"
.PHONY: report-artifact-governance

check-artifact-governance: root-check-env ## Fail when repo execution surfaces drift from governed artifact output discipline
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.artifact_governance check --repo-root "$(CURDIR)" --json-out "$(ROOT_ARTIFACTS_DIR)/artifact-governance.json"
.PHONY: check-artifact-governance

report-execution-surfaces: root-check-env ## Audit make, tox, and workflow execution-surface ownership
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.execution_surfaces report --repo-root "$(CURDIR)" --json-out "$(ROOT_ARTIFACTS_DIR)/execution-surfaces.json"
.PHONY: report-execution-surfaces

check-execution-surfaces: root-check-env ## Fail when repository execution surfaces mix ownership or governance concerns
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.execution_surfaces check --repo-root "$(CURDIR)" --json-out "$(ROOT_ARTIFACTS_DIR)/execution-surfaces.json"
.PHONY: check-execution-surfaces

report-package-boundaries: root-check-env ## Audit package ownership, exports, and cross-package boundary contracts
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.package_boundaries report --repo-root "$(CURDIR)" --json-out "$(ROOT_ARTIFACTS_DIR)/package-boundaries.json"
.PHONY: report-package-boundaries

check-package-boundaries: root-check-env ## Fail when runtime, alias, or maintainer package boundaries drift from owned contracts
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.package_boundaries check --repo-root "$(CURDIR)" --json-out "$(ROOT_ARTIFACTS_DIR)/package-boundaries.json"
.PHONY: check-package-boundaries

report-package-bundles: root-check-env ## Build package bundle audit reports for all publishable packages
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.package_bundles report --repo-root "$(CURDIR)" --artifacts-root "$(ROOT_ARTIFACTS_DIR)/package-bundles" --json-out "$(ROOT_ARTIFACTS_DIR)/package-bundles.json"
.PHONY: report-package-bundles

check-package-bundles: root-check-env ## Fail when publishable package bundles drift from governed policy
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.package_bundles check --repo-root "$(CURDIR)" --artifacts-root "$(ROOT_ARTIFACTS_DIR)/package-bundles" --json-out "$(ROOT_ARTIFACTS_DIR)/package-bundles.json"
.PHONY: check-package-bundles

report-publish-readiness: root-check-env ## Build the repository publish-readiness scorecard
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.publish_readiness --repo-root "$(CURDIR)" --json-out "$(ROOT_ARTIFACTS_DIR)/publish-readiness.json"
.PHONY: report-publish-readiness

check-publish-readiness: root-check-env ## Fail when the repository is not publish-ready
	@$(DEV_RUN) -m bijux_phylogenetics_dev.quality.publish_readiness --check --repo-root "$(CURDIR)" --json-out "$(ROOT_ARTIFACTS_DIR)/publish-readiness.json"
.PHONY: check-publish-readiness

report-release-readiness: root-check-env ## Build the full release-readiness evidence and governance report surface
	@$(MAKE) check-config-ssot
	@$(MAKE) report-evidence-governance
	@$(MAKE) report-execution-surfaces
	@$(MAKE) report-package-boundaries
	@$(MAKE) report-package-bundles
	@$(MAKE) report-publish-readiness
.PHONY: report-release-readiness

check-release-readiness: root-check-env ## Enforce the full release-readiness gate
	@$(MAKE) check-config-ssot
	@$(MAKE) check-evidence-governance
	@$(MAKE) check-execution-surfaces
	@$(MAKE) check-package-boundaries
	@$(MAKE) check-package-bundles
	@$(MAKE) check-publish-readiness
.PHONY: check-release-readiness

package-check: build ## Validate built distributions with twine
	@"$(ROOT_CHECK_PYTHON)" -m twine check "$(CURDIR)/artifacts/bijux-phylogenetics/build"/*.whl "$(CURDIR)/artifacts/bijux-phylogenetics/build"/*.tar.gz
.PHONY: package-check

package-smoke: build ## Install the wheel into a temp venv and run the governed CLI smoke proof
	@rm -rf "$(CURDIR)/artifacts/tmp/package-smoke"
	@PYTHONPATH="$(CURDIR)/packages/bijux-phylogenetics-dev/src$${PYTHONPATH:+:$${PYTHONPATH}}" \
	"$(ROOT_CHECK_PYTHON)" -m bijux_phylogenetics_dev.quality.package_install_smoke \
		--repo-root "$(CURDIR)" \
		--dist-dir "$(CURDIR)/artifacts/bijux-phylogenetics/build" \
		--artifacts-root "$(CURDIR)/artifacts/tmp/package-smoke" \
		--build-python "$(ROOT_CHECK_PYTHON)" \
		--artifact-kind wheel
.PHONY: package-smoke

package-source-smoke: build ## Install the sdist into a temp venv and run the governed CLI smoke proof
	@rm -rf "$(CURDIR)/artifacts/tmp/package-source-smoke"
	@PYTHONPATH="$(CURDIR)/packages/bijux-phylogenetics-dev/src$${PYTHONPATH:+:$${PYTHONPATH}}" \
	"$(ROOT_CHECK_PYTHON)" -m bijux_phylogenetics_dev.quality.package_install_smoke \
		--repo-root "$(CURDIR)" \
		--dist-dir "$(CURDIR)/artifacts/bijux-phylogenetics/build" \
		--artifacts-root "$(CURDIR)/artifacts/tmp/package-source-smoke" \
		--build-python "$(ROOT_CHECK_PYTHON)" \
		--artifact-kind sdist
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
