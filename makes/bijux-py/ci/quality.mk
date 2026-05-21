INTERROGATE_PATHS              ?= src
QUALITY_PATHS                  ?= $(INTERROGATE_PATHS)
QUALITY_ARTIFACTS_DIR          ?= $(PROJECT_ARTIFACTS_DIR)/quality
QUALITY_OK_MARKER              ?= $(QUALITY_ARTIFACTS_DIR)/_passed
QUALITY_MYPY_CONFIG            ?=
QUALITY_MYPY_FLAGS             ?= --strict
QUALITY_MYPY_TARGETS           ?= $(QUALITY_PATHS)
QUALITY_VULTURE_MIN_CONFIDENCE ?= 80
QUALITY_PRE_TARGETS            ?=
QUALITY_POST_TARGETS           ?=
QUALITY_RUN_MKDOCS             ?= 0
QUALITY_MKDOCS_CONFIG          ?= $(MKDOCS_CFG)
QUALITY_MYPY_CACHE_DIR         ?= $(QUALITY_ARTIFACTS_DIR)/.mypy_cache
QUALITY_MKDOCS_SITE_DIR        ?= $(QUALITY_ARTIFACTS_DIR)/docs/site
QUALITY_MKDOCS_CACHE_DIR       ?= $(QUALITY_ARTIFACTS_DIR)/docs/.cache
QUALITY_MKDOCS_PYCACHE_DIR     ?= $(QUALITY_ARTIFACTS_DIR)/docs/pycache
QUALITY_MKDOCS_PYTHON          ?= $(VENV_PYTHON)
QUALITY_MKDOCS_BUILD_FLAGS     ?= --strict --config-file "$(QUALITY_MKDOCS_CONFIG)" --site-dir "$(QUALITY_MKDOCS_SITE_DIR)"
QUALITY_DEPTRY_TARGET          ?= $(PROJECT_DIR)
QUALITY_DEPTRY_COMMAND         ?= $(DEPTRY) --config "$(DEPTRY_CONFIG)" "$(QUALITY_DEPTRY_TARGET)"
QUALITY_DEPTRY_VERSION_COMMAND ?= $(DEPTRY) --version
QUALITY_INTERROGATE_FLAGS      ?=
QUALITY_SELF_MAKE              ?= $(SELF_MAKE)
DOCSTRING_COVERAGE            ?= $(VENV_PYTHON) -m bijux_phylogenetics_dev.quality.docstring_coverage

PYTHON      ?= $(shell command -v python3 || command -v python)
VULTURE     ?= $(VENV_PYTHON) -m vulture
DEPTRY      ?= $(VENV_PYTHON) -m deptry
MYPY        ?= $(VENV_PYTHON) -m mypy

include $(abspath $(dir $(lastword $(MAKEFILE_LIST))))/util.mk

SKIP_DEPTRY      ?= 0
SKIP_INTERROGATE ?= 0
SKIP_MYPY        ?= 0

ifeq ($(shell uname -s),Darwin)
  BREW_PREFIX  := $(shell command -v brew >/dev/null 2>&1 && brew --prefix)
  CAIRO_PREFIX := $(shell test -n "$(BREW_PREFIX)" && brew --prefix cairo)
  QUALITY_ENV  := DYLD_FALLBACK_LIBRARY_PATH="$(BREW_PREFIX)/lib:$(CAIRO_PREFIX)/lib:$$DYLD_FALLBACK_LIBRARY_PATH"
else
  QUALITY_ENV  :=
endif

define run_docstring_coverage_report
	echo "→ Generating docstring coverage report (<100%)"; \
	mkdir -p "$(QUALITY_ARTIFACTS_DIR)"; \
	set +e; \
	  OUT="$$( $(QUALITY_ENV) $(DOCSTRING_COVERAGE) $(QUALITY_INTERROGATE_FLAGS) --verbose $(INTERROGATE_PATHS) )"; \
	  rc=$$?; \
	  printf '%s\n' "$$OUT" >"$(QUALITY_ARTIFACTS_DIR)/docstring_coverage.full.txt"; \
	  OFF="$$(printf '%s\n' "$$OUT" | awk -F' \\| ' '/^FILE: / { \
	    name=$$1; sub(/^FILE: /, "", name); cov=$$2; \
	    gsub(/^[ \t]+|[ \t]+$$/, "", name); gsub(/^[ \t]+|[ \t]+$$/, "", cov); \
	    if (cov != "100.0%" && cov != "100%") printf("  - %s (%s)\n", name, cov); \
	  }')"; \
	  printf '%s\n' "$$OFF" >"$(QUALITY_ARTIFACTS_DIR)/docstring_coverage.offenders.txt"; \
	  if [ $$rc -eq 0 ]; then \
	    RESULT_LINE="$$(printf '%s\n' "$$OUT" | awk '/^RESULT: / {print; exit}')"; \
	    if [ -n "$$RESULT_LINE" ]; then printf '%s\n' "$$RESULT_LINE"; else echo "✔ Docstring coverage passed"; fi; \
	  elif [ -n "$$OFF" ]; then \
	    printf '%s\n' "$$OFF"; \
	  else \
	    printf '%s\n' "$$OUT"; \
	  fi; \
	  exit $$rc
endef

.PHONY: quality interrogate-report docs-links quality-clean

quality:
	@echo "→ Running quality checks..."
	@mkdir -p "$(QUALITY_ARTIFACTS_DIR)" "$(QUALITY_MYPY_CACHE_DIR)"
	$(call run_make_targets,$(QUALITY_PRE_TARGETS),$(QUALITY_SELF_MAKE))
	@echo "   - Dead code analysis (Vulture)"
	@set -eu; \
	  { $(VULTURE) --version 2>/dev/null || echo vulture; } >"$(QUALITY_ARTIFACTS_DIR)/vulture.log"; \
	  OUT="$$( $(VULTURE) $(QUALITY_PATHS) --min-confidence $(QUALITY_VULTURE_MIN_CONFIDENCE) 2>&1 || true )"; \
	  printf '%s\n' "$$OUT" >>"$(QUALITY_ARTIFACTS_DIR)/vulture.log"; \
	  if [ -z "$$OUT" ]; then echo "✔ Vulture: no dead code found." >>"$(QUALITY_ARTIFACTS_DIR)/vulture.log"; fi
	@echo "   - Dependency hygiene (Deptry)"
	@if [ "$(SKIP_DEPTRY)" = "1" ]; then \
	  echo "✖ Deptry must remain enabled for $(PROJECT_SLUG)" | tee "$(QUALITY_ARTIFACTS_DIR)/deptry.log"; \
	  exit 1; \
	fi
	@set -euo pipefail; \
	  if [ -n "$(strip $(QUALITY_DEPTRY_VERSION_COMMAND))" ]; then \
	    { $(QUALITY_DEPTRY_VERSION_COMMAND) 2>/dev/null || true; } >"$(QUALITY_ARTIFACTS_DIR)/deptry.log"; \
	  else \
	    : >"$(QUALITY_ARTIFACTS_DIR)/deptry.log"; \
	  fi; \
	  $(QUALITY_DEPTRY_COMMAND) 2>&1 | tee -a "$(QUALITY_ARTIFACTS_DIR)/deptry.log"
	@echo "   - Static typing (Mypy)"
	@if [ "$(SKIP_MYPY)" = "1" ]; then \
	  echo "✖ Mypy must remain enabled for $(PROJECT_SLUG)" | tee "$(QUALITY_ARTIFACTS_DIR)/mypy.log"; \
	  exit 1; \
	fi
	@if [ -z "$(QUALITY_MYPY_CONFIG)" ]; then \
	  echo "✖ QUALITY_MYPY_CONFIG is required for $(PROJECT_SLUG)" | tee "$(QUALITY_ARTIFACTS_DIR)/mypy.log"; \
	  exit 1; \
	fi
	@set -euo pipefail; \
	  $(MYPY) --config-file "$(QUALITY_MYPY_CONFIG)" $(QUALITY_MYPY_FLAGS) --cache-dir "$(QUALITY_MYPY_CACHE_DIR)" $(QUALITY_MYPY_TARGETS) 2>&1 | tee "$(QUALITY_ARTIFACTS_DIR)/mypy.log"
	@echo "   - Documentation coverage"
	@if [ "$(SKIP_INTERROGATE)" = "1" ]; then \
	  echo "✖ Docstring coverage must remain enabled for $(PROJECT_SLUG)" | tee "$(QUALITY_ARTIFACTS_DIR)/docstring_coverage.full.txt"; \
	  exit 1; \
	fi
	@$(call run_docstring_coverage_report)
	$(call run_make_targets,$(QUALITY_POST_TARGETS),$(QUALITY_SELF_MAKE))
	@if [ "$(QUALITY_RUN_MKDOCS)" = "1" ]; then \
	  echo "   - MkDocs build"; \
	  mkdir -p "$(QUALITY_MKDOCS_SITE_DIR)" "$(QUALITY_MKDOCS_CACHE_DIR)"; \
	  XDG_CACHE_HOME="$(QUALITY_MKDOCS_CACHE_DIR)" PYTHONPYCACHEPREFIX="$(QUALITY_MKDOCS_PYCACHE_DIR)" "$(QUALITY_MKDOCS_PYTHON)" -m mkdocs build $(QUALITY_MKDOCS_BUILD_FLAGS); \
	fi
	@echo "✔ Quality checks passed"
	@printf "OK\n" >"$(QUALITY_OK_MARKER)"

interrogate-report:
	$(call run_docstring_coverage_report)

docs-links:
	@echo "→ docs-links is not configured for $(PROJECT_SLUG)"

quality-clean:
	@echo "→ Cleaning quality artifacts"
	@rm -rf "$(QUALITY_ARTIFACTS_DIR)"

##@ Quality
quality: ## Run repository quality checks and write logs to $(QUALITY_ARTIFACTS_DIR)
interrogate-report: ## Save the full docstring coverage report and offenders list
docs-links: ## Run the package docs link checker when configured
quality-clean: ## Remove quality artifacts
