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

[![bijux-phylogenetics docs](https://img.shields.io/badge/docs-bijux--phylogenetics-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
[![phylogenetic docs](https://img.shields.io/badge/docs-phylogenetic-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
<!-- bijux-phylogenetics-badges:generated:end -->

`bijux-phylogenetics` is a large phylogenetics runtime and evidence-bearing
repository. The checked-in tree already covers owned tree semantics, alignment
processing, native likelihood and inference, comparative analysis, ancestral
reconstruction, explicit parsimony, report generation, benchmark review, and
governed study evidence.

This site exists to make that depth visible without overstating what has been
scientifically closed.

```mermaid
flowchart TD
    A[Canonical runtime] --> B[Product handbook]
    A --> C[CLI, Python, and native contracts]
    C --> D[Governed artifacts and reports]
    D --> E[Evidence book]
    F[Maintainer package] --> G[Maintainer handbook]
    B --> H[Users and downstream developers]
    E --> I[Scientific reviewers]
    G --> J[Release and trust stewardship]
```

## What Readers Should Understand Immediately

| Surface | What it is for |
| --- | --- |
| [Product handbook](01-bijux-phylogenetics/index.md) | runtime ownership, interfaces, workflows, outputs, and quality reading |
| [Evidence book](02-bijux-phylogenetics-evidence-book/index.md) | governed comparative-study dossiers, bundle IDs, provenance, and bounded study claims |
| [Maintainer handbook](03-bijux-phylogenetics-maintain/index.md) | documentation truth, release alignment, artifact hygiene, and public trust stewardship |

## What The Checked-In Runtime Already Covers

| Capability family | Publicly visible depth |
| --- | --- |
| Trees | `PhyloTree`, rootedness, support normalization, topology comparison, clades, MRCA review, tree sets, and rendering |
| Alignments | FASTA validation, trimming, translation, coding checks, partition-aware handling, and packaged example inputs |
| Native inference | finite-state likelihood, native maximum-likelihood result objects, supported native Bayesian DNA workflows, and benchmark surfaces |
| Comparative biology | PGLS, signal, regression interpretation, discrete evolution, diversification, and model-selection reading |
| Ancestral and parsimony analysis | continuous and discrete reconstruction, uncertainty-ledger outputs, and explicit parsimony families |
| Datasets and reporting | packaged datasets, study inputs, benchmark corpora, HTML and TSV outputs, manifests, and evidence-linked bundles |

## Start By Goal

- I want to understand what the runtime owns and which interface I should use.
  Start with the [Phylogenetics product handbook](01-bijux-phylogenetics/index.md).
- I want to understand what the governed comparative-study record actually
  proves.
  Start with the [Evidence book](02-bijux-phylogenetics-evidence-book/index.md).
- I want to understand how the repository keeps documentation and release
  language honest.
  Start with the [Maintainer handbook](03-bijux-phylogenetics-maintain/index.md).

## Reader Routes

```mermaid
flowchart LR
    A[What can I run?] --> B[Product handbook]
    C[What has been checked?] --> D[Evidence book]
    E[Why should I trust the release story?] --> F[Maintainer handbook]
```

- Read the product handbook when the question is about runtime behavior,
  package roles, interfaces, workflows, outputs, or validation lanes.
- Read the evidence book when the question is about provenance, bundle IDs,
  study dossiers, or bounded parity claims.
- Read the maintainer handbook when the question is about narrative
  alignment, artifact hygiene, metadata truth, or release readiness.

## Flagship Starting Points

- [Surface selection guide](01-bijux-phylogenetics/interfaces/surface-selection.md)
- [Native maximum-likelihood workflows](01-bijux-phylogenetics/operations/native-maximum-likelihood-workflows.md)
- [Native Bayesian workflows](01-bijux-phylogenetics/operations/native-bayesian-workflows.md)
- [Validation lanes](01-bijux-phylogenetics/quality/validation-lanes.md)
- [Evidence book overview](02-bijux-phylogenetics-evidence-book/index.md)
- [Maintainer handbook overview](03-bijux-phylogenetics-maintain/index.md)

## Important Boundary

The code surface is broader than the currently closed evidence surface, and
the evidence surface is narrower than the full runtime surface. These docs
keep that distinction explicit on purpose. A feature can be documented and
usable without claiming that every surrounding scientific interpretation has
already been closed by governed study evidence.
