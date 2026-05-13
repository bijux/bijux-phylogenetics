ROOT_PACKAGE_PROFILE_DIR ?= $(ROOT_MAKEFILE_DIR)/packages

PACKAGE_RECORDS := \
	bijux-phylogenetics|primary,check,buildable,sbom|bijux-phylogenetics.mk \
	phylogenetic|compat,check,buildable,sbom|phylogenetic.mk \
	bijux-phylogenetics-dev|check|bijux-phylogenetics-dev.mk

include $(ROOT_MAKEFILE_DIR)/bijux-py/package-catalog.mk
