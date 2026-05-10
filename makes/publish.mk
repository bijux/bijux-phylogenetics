PUBLISH_VERSION_RESOLVER ?= -m bijux_phylogenetics_dev.release.version_resolver
PUBLISH_VERSION_GUARD ?= -m bijux_phylogenetics_dev.release.publication_guard
PUBLISH_VERSION_GUARD_FLAGS += --repo-root "$(MONOREPO_ROOT)" --require-config-ssot --require-publish-readiness
PUBLISH_UPLOAD_ENABLED ?= 1
PUBLISH_ALLOW_PRERELEASE ?= 0
PUBLISH_ALLOW_LOCAL_VERSION ?= 0

include $(ROOT_MAKE_DIR)/bijux-py/repository/publish.mk
