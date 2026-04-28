---
title: Bijux Phylogenetics
audience: public
type: overview
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-04-28
---

# Bijux Phylogenetics

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
inspection, validation, comparison, alignment trimming, coding-sequence
translation, explicit rooting transforms, evidence capture, and publishable
HTML reports.

The repository keeps the same managed Python workspace shape as the other
Bijux scientific repos: shared standards under `.bijux/`, repository-owned
`makes/` and `configs/`, one primary runtime package, one alias package, and
one maintainer package.

## Read this next

- runtime package handbook: [Repository overview](01-bijux-phylogenetics/index.md)
- maintainer handbook: [Maintainer overview](03-bijux-phylogenetics-maintain/index.md)


## Notable Workflows

- trim all-gap or all-missing alignment sites and remove high-missingness sequences
- inspect coding alignments for stop codons and frameshift-like sequence lengths
- translate nucleotide coding alignments into amino-acid alignments
- export pairwise sequence identity matrices for downstream review
- root trees on explicit outgroups or reroot trees by midpoint
