# bijux-phylogenetics

<!-- bijux-phylogenetics-badges:generated:start -->
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-3776AB?logo=python&logoColor=white)](https://pypi.org/project/bijux-phylogenetics/)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-0F766E)](https://github.com/bijux/bijux-phylogenetics/blob/main/LICENSE)
[![Verify](https://github.com/bijux/bijux-phylogenetics/workflows/repo%20/%20verify/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml?query=branch%3Amain)
[![Release PyPI](https://github.com/bijux/bijux-phylogenetics/workflows/release-pypi/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-pypi.yml)
[![Release GHCR](https://github.com/bijux/bijux-phylogenetics/workflows/release-ghcr/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-ghcr.yml)
[![Release GitHub](https://github.com/bijux/bijux-phylogenetics/workflows/release-github/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-github.yml)
[![Docs](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml)
[![Release](https://img.shields.io/github/v/release/bijux/bijux-phylogenetics?display_name=tag&label=release)](https://github.com/bijux/bijux-phylogenetics/releases)
[![GHCR packages](https://img.shields.io/badge/ghcr-2%20packages-181717?logo=github)](https://github.com/bijux?tab=packages&repo_name=bijux-phylogenetics)
[![Published packages](https://img.shields.io/badge/published%20packages-2-2563EB)](https://github.com/bijux/bijux-phylogenetics/tree/main/packages)

[![bijux-phylogenetics](https://img.shields.io/pypi/v/bijux-phylogenetics?label=bijux--phylogenetics&logo=pypi)](https://pypi.org/project/bijux-phylogenetics/)
[![phylogenetics](https://img.shields.io/pypi/v/phylogenetics?label=phylogenetics&logo=pypi)](https://pypi.org/project/phylogenetics/)

[![bijux-phylogenetics](https://img.shields.io/badge/bijux--phylogenetics-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fbijux-phylogenetics)
[![phylogenetics](https://img.shields.io/badge/phylogenetics-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fphylogenetics)

[![bijux-phylogenetics docs](https://img.shields.io/badge/docs-bijux--phylogenetics-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
[![phylogenetics docs](https://img.shields.io/badge/docs-phylogenetics-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
<!-- bijux-phylogenetics-badges:generated:end -->

`bijux-phylogenetics` is a reproducible phylogenetics workbench for tree
inspection, validation, comparison, metadata linkage, evidence capture, and
publishable reporting.

## Repository Layout

The repository keeps these durable top-level surfaces:

- `packages/` for published runtime, alias, and maintainer packages
- `docs/` for documentation source
- `examples/` for tracked workflow examples
- `datasets/` for durable repository-owned reference inputs
- `reports/` for tracked proof outputs and reviewer-facing artifacts
- `tests/` for cross-package test assets as the repository grows

The repository follows the same Bijux Python workspace pattern used by
`bijux-canon`, `bijux-proteomics`, and `bijux-pollenomics`: shared `.bijux`
assets, standardized `makes/` and `configs/`, a root `uv` workspace, a primary
runtime package, a compatibility alias package, and a repository-owned
maintainer package.




This repository publishes `2` public packages. The current runtime focuses on
reproducible tree hygiene, taxon normalization, tree and trait pruning,
metadata linkage, shared-clade comparison, evidence manifests, and HTML report
generation rather than tree inference.

## Start Here

- read the docs home: [Documentation home](https://bijux.io/bijux-phylogenetics/)
- inspect the runtime package source: [packages/bijux-phylogenetics](packages/bijux-phylogenetics)
- inspect the alias package source: [packages/phylogenetics](packages/phylogenetics)
- inspect maintainer tooling: [packages/bijux-phylogenetics-dev](packages/bijux-phylogenetics-dev)

## What This Repository Produces

Today, the checked-in repository produces these durable outcomes:

- a Python runtime package for phylogenetic validation, inspection, comparison, and reporting
- a compatibility alias distribution for the shorter `phylogenetics` command
- a repository-owned maintainer package for docs, release, and quality automation
- a MkDocs documentation site that builds into `artifacts/root/docs/site/`

## Common Workflows

- `make install` syncs the editable environment from the tracked `uv.lock`
- `make check` runs the main repository verification pass: lock check, lint, tests, docs, and distribution verification
- `make docs-serve` serves the docs locally at `http://127.0.0.1:8000/`
- `make package-verify` runs wheel, sdist, and smoke-install package proof targets
- `make sync-badges` renders the shared badge catalog into managed README surfaces

## Working Rules


- treat `.bijux/shared/`, `makes/`, and `configs/` as managed standard surfaces
- keep generated outputs under `artifacts/` unless the task explicitly governs another tracked destination
- update README and docs claims only when the current repository state actually supports them
