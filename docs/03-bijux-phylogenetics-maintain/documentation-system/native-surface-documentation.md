---
title: Native Surface Documentation
audience: maintainer
type: handbook
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-04
---

# Native Surface Documentation

When a native surface becomes part of the public release story, documentation
is not finished when one reference page mentions the API name.

## Minimum Public Coverage For One Flagship Native Surface

Each flagship native surface should have all of the following:

- one reference page that names the stable contract boundary
- one practical workflow page that shows how a serious user runs it
- one quality or evidence page that explains how to interpret trust,
  benchmark, or validation status
- one route from the repository entry surfaces to those pages

## Current Example Pattern

The native maximum-likelihood and native Bayesian families follow that shape:

- reference:
  `docs/phylogenetics/interfaces/native-inference-and-benchmarks.md`
- workflows:
  the native maximum-likelihood and native Bayesian pages under
  `docs/phylogenetics/operations/`
- benchmark interpretation:
  `docs/phylogenetics/quality/native-benchmark-evidence.md`

## Entry Surfaces To Check

- root `README.md`
- `packages/bijux-phylogenetics/README.md`
- `packages/phylogenetic/README.md` when the alias package should advertise
  the same runtime story
- `docs/index.md`
- `docs/phylogenetics/index.md`
- the matching reference, operations, and quality pages

## Weak Pattern To Prevent

- the API exists
- one reference page names it
- no workflow page shows serious use
- no quality page explains how to trust it

That state is better than silence, but it is still not release-grade
documentation for a flagship scientific surface.
