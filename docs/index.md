---
title: Bijux Phylogenetics
audience: public
type: overview
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-03
---

# Bijux Phylogenetics

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

`bijux-phylogenetics` is a reproducible phylogenetics workbench for tree
inspection, validation, comparison, alignment trimming, coding-sequence
translation, native likelihood and inference work, publication figure
packaging, evidence capture, and publishable HTML reports.

The repository keeps the same managed Python workspace shape as the other
Bijux scientific repos: shared standards under `.bijux/`, repository-owned
`makes/` and `configs/`, one canonical runtime package, one compatibility
alias package, and one maintainer package.

The public distributions are `bijux-phylogenetics` for the canonical runtime
surface and `phylogenetic` for the shorter compatibility alias. The same
repository also carries a maintainer package for documentation, release, and
quality automation.

## Read this next

- public runtime guide: [Phylogenetics product guide](public/phylogenetics/index.md)
- runtime architecture guide: [Runtime architecture](public/phylogenetics/architecture/index.md)
- native maximum-likelihood workflows: [Native maximum-likelihood workflows](public/phylogenetics/operations/native-maximum-likelihood-workflows.md)
- native Bayesian workflows: [Native Bayesian workflows](public/phylogenetics/operations/native-bayesian-workflows.md)
- native benchmark review: [Native benchmark review](public/phylogenetics/operations/native-benchmark-review.md)
- release review workflow: [Release review workflow](public/phylogenetics/quality/release-review-workflow.md)
- public evidence guide: [Evidence-book guide](public/phylogenetics-evidence-book/index.md)
- internal maintainer guide: [Internal documentation](internal/index.md)


## Notable Runtime Surfaces

- tree, split, and tree-set review on one owned `PhyloTree` runtime
- native finite-state likelihood evaluation across nucleotide, protein, codon, and discrete Mk families
- native maximum-likelihood inference result objects, support evidence, and wrapper-comparison benchmark reports
- native Bayesian entry points for supported DNA fixed-topology and topology-sampling workflows
- comparative and ancestral analysis surfaces that link directly to governed artifacts and reviewer-facing reports
- external-engine workflow orchestration where the repository intentionally treats wrappers as product surfaces rather than pretending they are native implementations

## Evidence Closure

The runtime surface is broader than the currently closed evidence surface.
The governed evidence-book records that difference explicitly, and the public
handbook now distinguishes native ownership, wrapper orchestration, and
evidence closure rather than flattening them into one claim:

- start from the public evidence guide: [Evidence Book](public/phylogenetics-evidence-book/index.md)
- review analytical surface coverage and closure status from the linked evidence-book index outputs there
