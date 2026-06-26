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
