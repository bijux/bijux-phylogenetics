PACKAGE_KIND := repository-python
PACKAGE_IMPORT_NAME := bijux_phylogenetics
PACKAGE_INSTALL_PYTHON_PACKAGES = "$(MONOREPO_ROOT)/packages/bijux-phylogenetics-dev[dev]"
LINT_DIRS = src tests
MYPY_TARGETS = src
ENABLE_MYPY := 1
ENABLE_CODESPELL := 0
ENABLE_RADON := 0
ENABLE_PYDOCSTYLE := 0
TEST_PATHS := tests
TEST_PATHS_UNIT := tests
TEST_MAIN_ARGS = -m "not slow"
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
BUILD_PACKAGE_NAME := bijux-phylogenetics
BUILD_PRE_TARGETS := sync-license-assets-package
BUILD_TEMP_CLEAN_PATHS := build dist *.egg-info
BUILD_POST_TARGETS := build-install-smoke
PACKAGE_NAME := bijux-phylogenetics

quality-compileall:
	@"$(VENV_PYTHON)" -m compileall src | tee "$(PROJECT_ARTIFACTS_DIR)/quality/compileall.log"
.PHONY: quality-compileall

sync-license-assets-package:
	@"$(VENV_PYTHON)" -m bijux_phylogenetics_dev.release.license_assets sync
.PHONY: sync-license-assets-package

build-install-smoke:
	@tmp_root="$(PROJECT_ARTIFACTS_DIR)/tmp/build-install-smoke"; \
	dist_name="$$(printf '%s' "$(BUILD_PACKAGE_NAME)" | tr '-' '_')"; \
	wheel_path="$$(ls -1t "$(BUILD_DIR_ABS)/$${dist_name}"-*.whl | head -n 1)"; \
	sdist_path="$$(ls -1t "$(BUILD_DIR_ABS)/$${dist_name}"-*.tar.gz | head -n 1)"; \
	export PIP_DISABLE_PIP_VERSION_CHECK=1; \
	if [ -z "$$wheel_path" ] || [ -z "$$sdist_path" ]; then \
	  echo "✘ Missing build artifacts for $(BUILD_PACKAGE_NAME) in $(BUILD_DIR_ABS)"; \
	  exit 1; \
	fi; \
	rm -rf "$$tmp_root"; \
	"$(BUILD_PYTHON)" -m venv "$$tmp_root/smoke"; \
	"$$tmp_root/smoke/bin/python" -m pip install "$$wheel_path"; \
	"$$tmp_root/smoke/bin/bijux-phylogenetics" --version; \
	"$$tmp_root/smoke/bin/bijux-phylogenetics" --help >/dev/null; \
	"$$tmp_root/smoke/bin/python" -m pip uninstall -y "$(BUILD_PACKAGE_NAME)" >/dev/null; \
	"$$tmp_root/smoke/bin/python" -m pip install "$$sdist_path"; \
	"$$tmp_root/smoke/bin/bijux-phylogenetics" --version; \
	"$$tmp_root/smoke/bin/bijux-phylogenetics" --help >/dev/null
.PHONY: build-install-smoke

include $(abspath $(dir $(firstword $(MAKEFILE_LIST))))/../bijux-py/package.mk
