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
[![GHCR packages](https://img.shields.io/badge/ghcr-2%20packages-181717?logo=github)](https://github.com/bijux?tab=packages&repo_name=bijux-phylogenetics)
[![Published packages](https://img.shields.io/badge/published%20packages-2-2563EB)](https://github.com/bijux/bijux-phylogenetics/tree/main/packages)

[![bijux-phylogenetics](https://img.shields.io/pypi/v/bijux-phylogenetics?label=bijux--phylogenetics&logo=pypi)](https://pypi.org/project/bijux-phylogenetics/)
[![phylogenetic](https://img.shields.io/pypi/v/phylogenetic?label=phylogenetic&logo=pypi)](https://pypi.org/project/phylogenetic/)

[![bijux-phylogenetics](https://img.shields.io/badge/bijux--phylogenetics-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fbijux-phylogenetics)
[![phylogenetic](https://img.shields.io/badge/phylogenetic-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fphylogenetic)

[![bijux-phylogenetics docs](https://img.shields.io/badge/docs-bijux--phylogenetics-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/public/phylogenetics/)
[![phylogenetic docs](https://img.shields.io/badge/docs-phylogenetic-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/public/phylogenetics/)
<!-- bijux-phylogenetics-badges:generated:end -->

`bijux-phylogenetics` is a contract-first Python repository for reproducible
phylogenetic validation, alignment diagnostics, comparative analysis, and
reviewable reporting around trees, alignments, and trait tables.

The repository follows the same managed workspace pattern as
`bijux-canon`, `bijux-proteomics`, and `bijux-pollenomics`: shared standards
under `.bijux/`, repository-owned `makes/` and `configs/`, a root `uv`
workspace, one canonical runtime package, one compatibility alias package, and
one repository-owned maintainer package.

This repository publishes `2` packages. Each release tag builds two staged
bundles, uploads both distributions to PyPI, publishes each release bundle to
its exact GHCR package page under the `bijux` account, and attaches the same
staged assets to the GitHub Release.

## Start Here

- read the public docs home: [Documentation home](https://bijux.io/bijux-phylogenetics/)
- inspect the public runtime guide: [Phylogenetics product guide](https://bijux.io/bijux-phylogenetics/public/phylogenetics/)
- inspect the public evidence guide: [Evidence Book](https://bijux.io/bijux-phylogenetics/public/phylogenetics-evidence-book/)
- inspect the internal maintainer guide: [Internal documentation](https://bijux.io/bijux-phylogenetics/internal/)
- inspect the runtime package source: [`packages/bijux-phylogenetics`](packages/bijux-phylogenetics)

## What This Repository Produces

Today, the checked-in repository produces four durable outcomes:

- a Python runtime package for tree validation, comparison, analysis, and reporting
- a shorter compatibility alias package that installs the same runtime under the `phylogenetic` distribution name
- a repository-owned maintainer package for docs, release, and quality automation
- a MkDocs documentation site that builds into `artifacts/root/docs/site/`

## Package Map

The `2` publishable packages in this repository are:

| Package | Role | Links |
| --- | --- | --- |
| `bijux-phylogenetics` | Canonical runtime package for phylogenetic diagnostics, comparative workflows, and reporting surfaces | <a href="https://pypi.org/project/bijux-phylogenetics/"><img alt="PyPI" src="https://img.shields.io/badge/pypi-3775A9?logo=pypi&logoColor=white" height="18"></a> <a href="https://bijux.io/bijux-phylogenetics/public/phylogenetics/"><img alt="Docs" src="https://img.shields.io/badge/docs-2563EB?logo=materialformkdocs&logoColor=white" height="18"></a> <a href="https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fbijux-phylogenetics"><img alt="GHCR" src="https://img.shields.io/badge/ghcr-181717?logo=github&logoColor=white" height="18"></a> <a href="https://github.com/bijux/bijux-phylogenetics/tree/main/packages/bijux-phylogenetics"><img alt="Source" src="https://img.shields.io/badge/source-181717?logo=github&logoColor=white" height="18"></a> |
| `phylogenetic` | Compatibility alias package that installs the canonical runtime with a shorter distribution and CLI name | <a href="https://pypi.org/project/phylogenetic/"><img alt="PyPI" src="https://img.shields.io/badge/pypi-3775A9?logo=pypi&logoColor=white" height="18"></a> <a href="https://bijux.io/bijux-phylogenetics/public/phylogenetics/"><img alt="Docs" src="https://img.shields.io/badge/docs-2563EB?logo=materialformkdocs&logoColor=white" height="18"></a> <a href="https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fphylogenetic"><img alt="GHCR" src="https://img.shields.io/badge/ghcr-181717?logo=github&logoColor=white" height="18"></a> <a href="https://github.com/bijux/bijux-phylogenetics/tree/main/packages/phylogenetic"><img alt="Source" src="https://img.shields.io/badge/source-181717?logo=github&logoColor=white" height="18"></a> |

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
- partition parsing, validation, summary tables, and partitioned adapter entrypoints for multi-locus inference
- comparative trait analysis, ancestral-state workflows, and distance-based analysis
- explicit supported-versus-excluded distance-tree method policy, including a governed `bionj` exclusion for this round
- governed adapter surfaces for external tools such as alignment and inference engines
- package handbooks and release automation

Current evidence closure is narrower than raw runtime breadth. The governed
status lives in the evidence-book:

- analytical surface coverage: [`evidence-book/index/analytical-surface-coverage.md`](evidence-book/index/analytical-surface-coverage.md)
- claim re-audit: [`evidence-book/index/claim-reaudit.md`](evidence-book/index/claim-reaudit.md)
- closure criteria: [`evidence-book/index/closure-criteria.md`](evidence-book/index/closure-criteria.md)
- maturity scorecard: [`evidence-book/index/evidence-maturity-scorecard.md`](evidence-book/index/evidence-maturity-scorecard.md)
- completion gates: [`evidence-book/index/completion-gates.md`](evidence-book/index/completion-gates.md)

What does not exist today:

- a native maximum-likelihood or Bayesian inference engine implemented inside this repository
- a broad checked-in benchmark corpus that replaces external method validation; targeted comparative trust studies may live under `evidence-book/`
- a claim that one passing workflow is sufficient scientific evidence on its own

## Common Commands

- `make help` to list repository automation targets
- `make install` to sync the editable environment from the tracked `uv.lock`
- `make check` to run the main verification pass: lock check, lint, tests, docs, and distribution verification
- `make package-verify` to build the wheel and sdist, validate them with Twine,
  install each into a clean virtual environment, copy the packaged example
  inputs through the runtime API, run CLI help, validate the copied example
  FASTA, render a packaged tree report bundle, and fit a comparative PGLS
  model on the packaged primate dataset
- `make docs-serve` to serve the docs locally at `http://127.0.0.1:8000/`
- `make sync-badges` to render the shared badge catalog into managed README surfaces

## Local Artifact Contract

- transient local outputs belong under `artifacts/`, not as ad hoc root-level
  cache or build directories
- curated checked-in comparative trust studies belong under `evidence-book/`; they are
  repository evidence bundles rather than transient execution products
- the shared root environment lives at `artifacts/root/check-venv/`
- the MkDocs site builds to `artifacts/root/docs/site/`

## Repository Layout

The root keeps repository-owned concerns explicit:

- `apis/` for versioned checked-in API contract bundles when this repository owns them
- `configs/` for shared tool configuration
- `docs/` for the repository handbook and package handbook index
- `makes/` for automation and orchestration
- `evidence-book/` for checked-in comparative parity studies and trust evidence
- `.github/workflows/` for CI, release, and docs deployment pipelines
- `packages/` for publishable package boundaries and maintainer tooling
- `artifacts/` for transient local outputs such as built docs, package bundles, and test products

That split is intentional: runtime code stays local to packages, and repository
governance stays visible and reviewable at the root.

## Documentation

The canonical project documentation lives in `docs/` and is built with MkDocs.

Useful entry points:

- docs home: [`docs/index.md`](docs/index.md)
- public runtime guide: [`docs/public/phylogenetics/index.md`](docs/public/phylogenetics/index.md)
- public evidence guide: [`docs/public/phylogenetics-evidence-book/index.md`](docs/public/phylogenetics-evidence-book/index.md)
- internal maintainer guide: [`docs/internal/maintain/index.md`](docs/internal/maintain/index.md)
- package source: [`packages/bijux-phylogenetics/src`](packages/bijux-phylogenetics/src)
- package tests: [`packages/bijux-phylogenetics/tests`](packages/bijux-phylogenetics/tests)

## License

This repository is licensed under the Apache License 2.0. See
[`LICENSE`](LICENSE) and [`NOTICE`](NOTICE).
