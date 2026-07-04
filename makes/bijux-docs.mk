# Shared Bijux docs shell synchronization and contract enforcement.

# Prefer the repository root virtualenv when present, then fall back to the
# quality-gate interpreter selection used elsewhere in this repository.
PYTHON_BIN ?= $(if $(wildcard $(MONOREPO_ROOT)/.venv/bin/python),$(abspath $(MONOREPO_ROOT)/.venv/bin/python),$(QUALITY_GATE_PYTHON))
BIJUX_DOCS_SYNC_SCRIPT ?= .bijux/shared/bijux-docs/tooling/scripts/sync_bijux_docs.sh
BIJUX_DOCS_SOT_GUARD ?= .bijux/shared/bijux-docs/tooling/scripts/verify_bijux_docs_source_of_truth.sh
BIJUX_DOCS_CONTRACT_GUARD ?= .bijux/shared/bijux-docs/tooling/quality/validate_bijux_docs_contract.py

.PHONY: bijux-docs-sync bijux-docs-check shell-sync shell-check

bijux-docs-sync: ## Synchronize shared Bijux docs shell into docs assets
	@bash "$(BIJUX_DOCS_SYNC_SCRIPT)"

bijux-docs-check: ## Validate Bijux docs shell contract and drift checks
	@"$(PYTHON_BIN)" "$(BIJUX_DOCS_CONTRACT_GUARD)" .
	@bash "$(BIJUX_DOCS_SOT_GUARD)"

# Backward-compatible aliases.
shell-sync: bijux-docs-sync ## Alias for bijux-docs-sync
shell-check: bijux-docs-check ## Alias for bijux-docs-check
