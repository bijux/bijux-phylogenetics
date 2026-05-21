BIJUX_REPOSITORY_ENV_OVERLAY_INCLUDED := 1

MONOREPO_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/..)
ROOT_MAKE_DIR := $(MONOREPO_ROOT)/makes
QUALITY_GATE_PYTHON ?= $(if $(wildcard $(VENV_PYTHON)),$(abspath $(VENV_PYTHON)),$(if $(wildcard $(VENV)/bin/python),$(abspath $(VENV)/bin/python),$(PYTHON)))
DEPTRY_SCAN_SCRIPT ?= PYTHONPATH="$(MONOREPO_ROOT)/packages/bijux-phylogenetics-dev/src$${PYTHONPATH:+:$$PYTHONPATH}" "$(QUALITY_GATE_PYTHON)" -m bijux_phylogenetics_dev.quality.deptry_scan
DEPTRY_CONFIG ?= $(MONOREPO_ROOT)/configs/deptry.toml
QUALITY_DEPTRY_COMMAND ?= $(DEPTRY_SCAN_SCRIPT) --config "$(DEPTRY_CONFIG)" --project-dir . $(QUALITY_PATHS)
QUALITY_DEPTRY_VERSION_COMMAND ?=
CODESPELL ?= $(VENV_PYTHON) -m codespell_lib

include $(ROOT_MAKE_DIR)/bijux-py/repository/env.mk

# Prefer the active package virtualenv's uv executable when available so the
# shared make layers do not depend on a separately installed global uv binary.
UV := $(if $(wildcard $(ACT)/uv),$(ACT)/uv,$(if $(wildcard $(VENV)/bin/uv),$(VENV)/bin/uv,uv))

# Package roots expose tracked symlink aliases for repository-owned artifact
# locations. Keep the aliases stable and let clean targets operate on the
# canonical repository artifact tree instead of deleting the links themselves.
COMMON_PYTHON_CLEAN_PATHS := $(filter-out .hypothesis .benchmarks,$(COMMON_PYTHON_CLEAN_PATHS))
PROJECT_ARTIFACT_PRESERVE_DIRS ?= venv hypothesis benchmarks
PROJECT_ARTIFACT_CHILD_CLEAN_PATHS := $(shell if [ -d "$(PROJECT_ARTIFACTS_DIR)" ]; then find "$(PROJECT_ARTIFACTS_DIR)" -mindepth 1 -maxdepth 1 $(foreach dir,$(PROJECT_ARTIFACT_PRESERVE_DIRS),! -name "$(dir)") -print; fi)
COMMON_ARTIFACT_CLEAN_PATHS := $(PROJECT_ARTIFACT_CHILD_CLEAN_PATHS)
