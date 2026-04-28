API_FREEZE_COMMAND ?= $(VENV_PYTHON) -m bijux_phylogenetics_dev.api.freeze_contracts
API_OPENAPI_DRIFT_COMMAND ?= $(VENV_PYTHON) -m bijux_phylogenetics_dev.api.openapi_drift
API_INSTALL_PYTHON_PACKAGES ?= --reinstall click prance openapi-spec-validator

include $(ROOT_MAKE_DIR)/bijux-py/api-freeze.mk
