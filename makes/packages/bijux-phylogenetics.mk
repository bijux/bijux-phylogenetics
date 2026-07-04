PACKAGE_KIND := repository-python
PACKAGE_IMPORT_NAME := bijux_phylogenetics
PACKAGE_INSTALL_PYTHON_PACKAGES = "$(MONOREPO_ROOT)/packages/bijux-phylogenetics-dev[dev]"
LINT_DIRS = $(MONOREPO_ROOT)/packages/bijux-phylogenetics/src $(MONOREPO_ROOT)/packages/bijux-phylogenetics/tests
MYPY_TARGETS = src
ENABLE_MYPY := 1
ENABLE_CODESPELL := 0
ENABLE_RADON := 0
ENABLE_PYDOCSTYLE := 0
TEST_SMOKE_PATHS := tests/smoke
TEST_FULL_PATHS := tests
TEST_PATHS := $(TEST_SMOKE_PATHS)
TEST_PATHS_UNIT := $(TEST_FULL_PATHS)
TEST_PATHS_EVALUATION = $(MONOREPO_ROOT)/packages/bijux-phylogenetics/tests
TEST_REAL_LOCAL_PATH = $(MONOREPO_ROOT)/packages/bijux-phylogenetics/tests/real_local
TEST_MAIN_ARGS = -m "not slow and not real_local and not evaluation" --maxfail=1 -q
TEST_UNIT_DIR_ARGS = -m "not slow and not real_local and not evaluation" --maxfail=1 -q
TEST_UNIT_FALLBACK_ARGS = -k "not e2e and not integration and not functional" -m "not slow and not real_local and not evaluation" --maxfail=1 -q
TEST_LABEL = development smoke suite
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
FULL_SUITE_TEST_PRE_TARGETS := sync-license-assets-package ensure-live-r-test-runtime
BUILD_PACKAGE_NAME := bijux-phylogenetics
BUILD_PRE_TARGETS := sync-license-assets-package
BUILD_TEMP_CLEAN_PATHS := build dist *.egg-info
BUILD_POST_TARGETS := build-install-smoke
PACKAGE_NAME := bijux-phylogenetics
SBOM_REQUIREMENTS_WRITER := -m bijux_phylogenetics_dev.quality.requirements_writer
FROZEN_TEST_REPORT_PYTHON = $(MONOREPO_ROOT)/artifacts/root/check-venv/bin/python
FULL_SUITE_PYTEST_WORKERS ?= 8
FULL_SUITE_PYTEST_DISTRIBUTION ?= worksteal
CORE_TEST_FILTER_ARGS = -m "not slow and not real_local and not evaluation" --maxfail=1 -q

define run_unit_partition
	@$(SELF_MAKE) test \
		TEST_LABEL="$(1)" \
		TEST_PATHS="$(2)" \
		TEST_MAIN_ARGS='$(CORE_TEST_FILTER_ARGS)' \
		TEST_ARTIFACTS_DIR="$(PROJECT_ARTIFACTS_DIR)/test/$(3)"
endef

quality-compileall:
	@"$(VENV_PYTHON)" -m compileall src | tee "$(PROJECT_ARTIFACTS_DIR)/quality/compileall.log"
.PHONY: quality-compileall

sync-license-assets-package:
	@"$(VENV_PYTHON)" -m bijux_phylogenetics_dev.release.license_assets sync
.PHONY: sync-license-assets-package

ensure-live-r-test-runtime:
	@"$(VENV_PYTHON)" -m bijux_phylogenetics_dev.quality.live_r_runtime --repo-root "$(MONOREPO_ROOT)" --json-out "$(PROJECT_ARTIFACTS_DIR)/test/live-r-runtime.json"
.PHONY: ensure-live-r-test-runtime

test-all: TEST_LABEL = full test suite
test-all: TEST_PATHS = $(TEST_FULL_PATHS)
test-all: TEST_PRE_TARGETS = $(FULL_SUITE_TEST_PRE_TARGETS)
test-all: TEST_MAIN_ARGS =
test-all: PYTEST_ADDOPTS_EXTRA = -o timeout=0 -n $(FULL_SUITE_PYTEST_WORKERS) --dist $(FULL_SUITE_PYTEST_DISTRIBUTION)
test-all: test
.PHONY: test-all

test-all-plus-run-time: TEST_LABEL = full test suite
test-all-plus-run-time: TEST_PATHS = $(TEST_FULL_PATHS)
test-all-plus-run-time: TEST_PRE_TARGETS = $(FULL_SUITE_TEST_PRE_TARGETS)
test-all-plus-run-time: TEST_MAIN_ARGS =
test-all-plus-run-time: PYTEST_ADDOPTS_EXTRA = -o timeout=0 -n $(FULL_SUITE_PYTEST_WORKERS) --dist $(FULL_SUITE_PYTEST_DISTRIBUTION) --durations=0 --durations-min=0
test-all-plus-run-time: test
.PHONY: test-all-plus-run-time

test-all-frozen: ## Compare the full-suite output against TEST_ALL_FROZEN_REF
	@test -n "$(TEST_ALL_FROZEN_REF)" || { echo "TEST_ALL_FROZEN_REF is required"; exit 2; }
	@$(MAKE) -C "$(MONOREPO_ROOT)" root-check-env >/dev/null
	@"$(FROZEN_TEST_REPORT_PYTHON)" -m bijux_phylogenetics_dev.quality.frozen_test_report check --repo-root "$(MONOREPO_ROOT)" --package "$(PROJECT_SLUG)" --frozen-ref "$(TEST_ALL_FROZEN_REF)" --artifacts-root "$(PROJECT_ARTIFACTS_DIR)/test-frozen" --json-out "$(PROJECT_ARTIFACTS_DIR)/test-frozen/report.json" --full-suite-workers "$(FULL_SUITE_PYTEST_WORKERS)"
.PHONY: test-all-frozen

test-unit-sequence-data-foundations:
	$(call run_unit_partition,unit sequence data foundations,tests/alignment tests/datasets tests/distance,sequence-data-foundations)
.PHONY: test-unit-sequence-data-foundations

test-unit-ancestral-dating-evidence:
	$(call run_unit_partition,unit ancestral dating and evidence,tests/ancestral tests/dating tests/evidence,ancestral-dating-evidence)
.PHONY: test-unit-ancestral-dating-evidence

test-unit-bayesian-model-definition:
	$(call run_unit_partition,unit bayesian model definition,tests/bayesian/models tests/bayesian/priors,bayesian-model-definition)
.PHONY: test-unit-bayesian-model-definition

test-unit-bayesian-inference-runtime:
	$(call run_unit_partition,unit bayesian inference runtime,tests/bayesian/engines tests/bayesian/posterior tests/bayesian/proposals tests/bayesian/sampling,bayesian-inference-runtime)
.PHONY: test-unit-bayesian-inference-runtime

test-unit-comparative-regression-and-signal:
	$(call run_unit_partition,unit comparative regression and signal,tests/comparative/assessment tests/comparative/contracts tests/comparative/foundations tests/comparative/optimization tests/comparative/regression tests/comparative/trait_evolution,comparative-regression-and-signal)
.PHONY: test-unit-comparative-regression-and-signal

test-unit-comparative-biogeography-and-publication:
	$(call run_unit_partition,unit comparative biogeography and publication,tests/comparative/biogeography tests/comparative/macroevolution tests/comparative/publication tests/publication,comparative-biogeography-and-publication)
.PHONY: test-unit-comparative-biogeography-and-publication

test-unit-likelihood-search-and-fit:
	$(call run_unit_partition,unit likelihood search and fit,tests/likelihood/assessment tests/likelihood/optimization tests/likelihood/search,likelihood-search-and-fit)
.PHONY: test-unit-likelihood-search-and-fit

test-unit-likelihood-placements-and-topology:
	$(call run_unit_partition,unit likelihood placements and topology,tests/likelihood/inference tests/likelihood/placement tests/likelihood/starting_trees tests/likelihood/support tests/likelihood/wrappers tests/parsimony tests/topology,likelihood-placements-and-topology)
.PHONY: test-unit-likelihood-placements-and-topology

test-unit-state-model-foundations:
	$(call run_unit_partition,unit state model foundations,tests/state_models/compiled tests/state_models/ctmc tests/state_models/discrete_traits tests/state_models/finite_state tests/state_models/parameters,state-model-foundations)
.PHONY: test-unit-state-model-foundations

test-unit-state-model-mixtures-and-tree-workflows:
	$(call run_unit_partition,unit state model mixtures and tree workflows,tests/state_models/families tests/state_models/mixtures tests/state_models/optimization tests/trees tests/workflows tests/parity,state-model-mixtures-and-tree-workflows)
.PHONY: test-unit-state-model-mixtures-and-tree-workflows

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
