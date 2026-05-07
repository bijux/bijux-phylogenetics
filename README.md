# bijux-phylogenetics

<!-- bijux-phylogenetics-badges:generated:start -->
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-3776AB?logo=python&logoColor=white)](https://pypi.org/project/bijux-phylogenetics/)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-0F766E)](https://github.com/bijux/bijux-phylogenetics/blob/main/LICENSE)
[![Verify](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml/badge.svg?branch=main)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml?query=branch%3Amain)
[![Release PyPI](https://img.shields.io/badge/release-pypi%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-pypi.yml)
[![Release GHCR](https://img.shields.io/badge/release-ghcr%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-ghcr.yml)
[![Release GitHub](https://img.shields.io/badge/release-github%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-github.yml)
[![Docs](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml)
[![Release](https://img.shields.io/github/v/release/bijux/bijux-phylogenetics?display_name=tag&label=release)](https://github.com/bijux/bijux-phylogenetics/releases)
[![GHCR packages](https://img.shields.io/badge/ghcr-1%20packages-181717?logo=github)](https://github.com/bijux?tab=packages&repo_name=bijux-phylogenetics)
[![Published packages](https://img.shields.io/badge/published%20packages-1-2563EB)](https://github.com/bijux/bijux-phylogenetics/tree/main/packages)

[![bijux-phylogenetics](https://img.shields.io/pypi/v/bijux-phylogenetics?label=bijux--phylogenetics&logo=pypi)](https://pypi.org/project/bijux-phylogenetics/)

[![bijux-phylogenetics](https://img.shields.io/badge/bijux--phylogenetics-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fbijux-phylogenetics)

[![bijux-phylogenetics docs](https://img.shields.io/badge/docs-bijux--phylogenetics-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
<!-- bijux-phylogenetics-badges:generated:end -->

`bijux-phylogenetics` is a contract-first Python repository for reproducible
phylogenetic validation, alignment diagnostics, comparative analysis, and
reviewable reporting around trees, alignments, and trait tables.

The repository follows the same managed workspace pattern as
`bijux-canon`, `bijux-proteomics`, and `bijux-pollenomics`: shared standards
under `.bijux/`, repository-owned `makes/` and `configs/`, a root `uv`
workspace, one primary runtime package, and one repository-owned maintainer
package.

This repository publishes `1` package. Each release tag builds one staged
bundle, uploads the distribution to PyPI, publishes the release bundle to its
exact GHCR package page under the `bijux` account, and attaches the same staged
assets to the GitHub Release.

## Start Here

- read the public docs home: [Documentation home](https://bijux.io/bijux-phylogenetics/)
- inspect the runtime handbook: [Repository overview](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
- inspect the maintainer handbook: [Maintainer overview](https://bijux.io/bijux-phylogenetics/03-bijux-phylogenetics-maintain/)
- inspect the runtime package source: [`packages/bijux-phylogenetics`](packages/bijux-phylogenetics)

## What This Repository Produces

Today, the checked-in repository produces three durable outcomes:

- a Python runtime package for tree validation, comparison, analysis, and reporting
- a repository-owned maintainer package for docs, release, and quality automation
- a MkDocs documentation site that builds into `artifacts/root/docs/site/`

## Package Map

The `1` publishable package in this repository is:

| Package | Role | Links |
| --- | --- | --- |
| `bijux-phylogenetics` | Canonical runtime package for phylogenetic diagnostics, comparative workflows, and reporting surfaces | <a href="https://pypi.org/project/bijux-phylogenetics/"><img alt="PyPI" src="https://img.shields.io/badge/pypi-3775A9?logo=pypi&logoColor=white" height="18"></a> <a href="https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/"><img alt="Docs" src="https://img.shields.io/badge/docs-2563EB?logo=materialformkdocs&logoColor=white" height="18"></a> <a href="https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fbijux-phylogenetics"><img alt="GHCR" src="https://img.shields.io/badge/ghcr-181717?logo=github&logoColor=white" height="18"></a> <a href="https://github.com/bijux/bijux-phylogenetics/tree/main/packages/bijux-phylogenetics"><img alt="Source" src="https://img.shields.io/badge/source-181717?logo=github&logoColor=white" height="18"></a> |

Repository-owned developer tooling also lives here in
[`packages/bijux-phylogenetics-dev`](packages/bijux-phylogenetics-dev), but it
is for maintaining the workspace rather than for end-user installation.

## Current Scope

The current repository scope is deliberate. `bijux-phylogenetics` is meant to
govern evidence, diagnostics, and reproducible workflow surfaces around
phylogenetic work. It is not trying to replace full inference engines.

What exists today:

- tree validation, inspection, comparison, rooting, and rendering workflows
- alignment diagnostics, filtering, trimming, and translation workflows
- comparative trait analysis, ancestral-state workflows, and distance-based analysis
- governed adapter surfaces for external tools such as alignment and inference engines
- checked-in API contracts, package handbooks, and release automation

What does not exist today:

- a native maximum-likelihood or Bayesian inference engine implemented inside this repository
- a broad checked-in benchmark corpus that replaces external method validation; targeted comparative trust reports may live under `reports/`
- a claim that one passing workflow is sufficient scientific evidence on its own

## Common Commands

- `make help` to list repository automation targets
- `make install` to sync the editable environment from the tracked `uv.lock`
- `make check` to run the main verification pass: lock check, lint, tests, docs, and distribution verification
- `make package-verify` to run wheel, sdist, and smoke-install package proofs
- `make docs-serve` to serve the docs locally at `http://127.0.0.1:8000/`
- `make sync-badges` to render the shared badge catalog into managed README surfaces

## Local Artifact Contract

- transient local outputs belong under `artifacts/`, not as ad hoc root-level
  cache or build directories
- curated checked-in comparative trust reports belong under `reports/`; they are
  repository evidence bundles rather than transient execution products
- the shared root environment lives at `artifacts/root/check-venv/`
- the MkDocs site builds to `artifacts/root/docs/site/`

## Repository Layout

The root keeps repository-owned concerns explicit:

- `apis/` for checked-in API contracts, pinned canonical JSON, and schema digests
- `configs/` for shared tool configuration
- `docs/` for the repository handbook and package handbook index
- `makes/` for automation and orchestration
- `reports/` for checked-in comparative validation studies and trust reports
- `.github/workflows/` for CI, release, and docs deployment pipelines
- `packages/` for publishable package boundaries and maintainer tooling
- `artifacts/` for transient local outputs such as built docs, package bundles, and test products

That split is intentional: runtime code stays local to packages, and repository
governance stays visible and reviewable at the root.

## Documentation

The canonical project documentation lives in `docs/` and is built with MkDocs.

Useful entry points:

- docs home: [`docs/index.md`](docs/index.md)
- runtime handbook: [`docs/01-bijux-phylogenetics/index.md`](docs/01-bijux-phylogenetics/index.md)
- maintainer handbook: [`docs/03-bijux-phylogenetics-maintain/index.md`](docs/03-bijux-phylogenetics-maintain/index.md)
- package source: [`packages/bijux-phylogenetics/src`](packages/bijux-phylogenetics/src)
- package tests: [`packages/bijux-phylogenetics/tests`](packages/bijux-phylogenetics/tests)

## License

This repository is licensed under the Apache License 2.0. See
[`LICENSE`](LICENSE) and [`NOTICE`](NOTICE).
