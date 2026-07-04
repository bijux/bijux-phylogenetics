# phylogenetic

<!-- bijux-phylogenetics-badges:generated:start -->
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-3776AB?logo=python&logoColor=white)](https://pypi.org/project/phylogenetic/)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-0F766E)](https://github.com/bijux/bijux-phylogenetics/blob/main/LICENSE)
[![Verify](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml/badge.svg?branch=main)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml?query=branch%3Amain)
[![Release PyPI](https://img.shields.io/badge/release-pypi%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-pypi.yml)
[![Release GHCR](https://img.shields.io/badge/release-ghcr%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-ghcr.yml)
[![Release GitHub](https://img.shields.io/badge/release-github%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-github.yml)
[![Docs](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml)

[![phylogenetic](https://img.shields.io/pypi/v/phylogenetic?label=phylogenetic&logo=pypi)](https://pypi.org/project/phylogenetic/)
[![bijux-phylogenetics](https://img.shields.io/pypi/v/bijux-phylogenetics?label=bijux--phylogenetics&logo=pypi)](https://pypi.org/project/bijux-phylogenetics/)

[![phylogenetic](https://img.shields.io/badge/phylogenetic-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fphylogenetic)
[![bijux-phylogenetics](https://img.shields.io/badge/bijux--phylogenetics-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fbijux-phylogenetics)

[![phylogenetic docs](https://img.shields.io/badge/docs-phylogenetic-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
[![bijux-phylogenetics docs](https://img.shields.io/badge/docs-bijux--phylogenetics-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
<!-- bijux-phylogenetics-badges:generated:end -->

Compatibility alias package for `bijux-phylogenetics`.

This distribution installs the canonical phylogenetics runtime under the
shorter `phylogenetic` package and CLI names. It is not a reduced runtime and
it is not a fork. It points at the same benchmark, workflow, and reporting
surfaces that the canonical package publishes.

The full handbook includes diagrams and architecture maps on the documentation
site. This PyPI page stays text-first so it renders cleanly in the package
index.

## Install

`phylogenetic` supports Python 3.11 and newer.

```bash
python3.11 -m pip install phylogenetic
phylogenetic --help
```

## What You Still Get

| Capability family | What carries through the alias |
| --- | --- |
| Trees and alignments | inspection, validation, comparison, and preparation workflows |
| Native runtime | native likelihood and native inference surfaces |
| Comparative analysis | comparative and ancestral workflows |
| Reporting | benchmark, workflow, and reporting surfaces |
| External engines | governed wrapper-backed external-engine workflows |

## What The Alias Does Not Change

- it does not remove native inference surfaces
- it does not remove comparative, ancestral, or parsimony families
- it does not change the evidence-book or quality reading model
- it does not publish a smaller scientific contract than the canonical package

Choose this package when you want the shorter install and command names while
keeping the same public runtime contract.

Choose the canonical `bijux-phylogenetics` package instead when you want the
most explicit package name for notebooks, papers, release notes, or downstream
documentation.

## Trust And Evidence

The alias changes naming, not scientific scope. Runtime breadth, benchmark
surfaces, and evidence-book boundaries all come from the canonical runtime and
documentation set.

## Read Next

- [Runtime product guide](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
- [Runtime architecture](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/architecture/)
- [Surface selection guide](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/surface-selection/)
- [Artifact consumption guide](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/artifact-consumption-guide/)
- [Native Inference And Benchmark Surfaces](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/native-inference-and-benchmarks/)
- [Native maximum-likelihood workflows](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/operations/native-maximum-likelihood-workflows/)
- [Native Bayesian workflows](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/operations/native-bayesian-workflows/)
- [Native benchmark review](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/operations/native-benchmark-review/)
- [Release review workflow](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/quality/release-review-workflow/)
- [Evidence-book guide](https://bijux.io/bijux-phylogenetics/02-bijux-phylogenetics-evidence-book/)
- [Canonical runtime package on GitHub](https://github.com/bijux/bijux-phylogenetics/tree/main/packages/bijux-phylogenetics)
