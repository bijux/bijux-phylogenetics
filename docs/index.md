---
title: Bijux Phylogenetics
audience: public
type: overview
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-04
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

[![bijux-phylogenetics docs](https://img.shields.io/badge/docs-bijux--phylogenetics-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/phylogenetics/)
[![phylogenetic docs](https://img.shields.io/badge/docs-phylogenetic-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/phylogenetics/)
<!-- bijux-phylogenetics-badges:generated:end -->

`bijux-phylogenetics` is a reproducible phylogenetics toolkit for tree
inspection, alignment review, native inference, comparative analysis,
ancestral reconstruction, evidence-linked reporting, and governed workflow
bundles.

This documentation is split into three reader-facing surfaces so claims stay
clear:

- the product guide explains what the runtime does and how to use it
- the evidence book explains what has been checked against governed study
  records
- the maintainer guide explains how documentation, release, and repository
  health stay aligned

```mermaid
flowchart LR
    A[Runtime package] --> B[Product guide]
    A --> C[Workflow outputs and reports]
    C --> D[Evidence book]
    B --> E[Public users]
    D --> F[Scientific reviewers]
    G[Maintainer tooling] --> H[Maintainer guide]
    H --> I[Release and docs stewardship]
```

## Start By Goal

- I want to understand the runtime surface and choose the right interface.
  Start with the [Phylogenetics product guide](phylogenetics/index.md).
- I want to understand what has been checked on course-derived comparative
  studies.
  Start with the [Evidence Book](phylogenetics-evidence-book/index.md).
- I want to maintain this repository without letting docs and release claims
  drift.
  Start with the [Maintainer guide](phylogenetics-maintainer/index.md).

## What This Repository Publishes

- `bijux-phylogenetics`
  The canonical runtime distribution and the documentation center of gravity.
- `phylogenetic`
  A shorter compatibility distribution that maps to the same runtime story.
- `bijux-phylogenetics-dev`
  Repository-owned maintainer tooling for docs, release, and quality
  automation.

## Product Highlights

- owned tree, split, and tree-set review on the `PhyloTree` runtime
- native finite-state likelihood across nucleotide, protein, codon, and
  discrete Mk families
- native maximum-likelihood result objects with support and benchmark evidence
- native Bayesian DNA workflows for supported fixed-topology and
  joint-topology runs
- comparative and ancestral workflows that write reviewable artifacts instead
  of opaque notebook state
- governed wrapper orchestration for external scientific engines when the
  wrapper is itself the honest public surface

## Read In Order

1. [Product guide](phylogenetics/index.md)
2. [Surface selection guide](phylogenetics/interfaces/surface-selection.md)
3. [Installation and setup](phylogenetics/operations/installation-and-setup.md)
4. [Claim reading guide](phylogenetics/foundation/claim-reading-guide.md)
5. [Evidence Book](phylogenetics-evidence-book/index.md)
6. [Release review workflow](phylogenetics/quality/release-review-workflow.md)

## Important Boundary

The runtime surface is broader than the currently closed evidence surface. The
documentation keeps that distinction explicit on purpose. A feature can be
published, documented, and usable without implying that every scientific claim
around it is already closed by study evidence.
