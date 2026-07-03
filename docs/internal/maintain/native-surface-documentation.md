---
title: Native Surface Documentation
audience: maintainer
type: handbook
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-03
---

# Native Surface Documentation

When one native surface becomes part of the public release story, documentation
work is not finished when one reference page mentions the API name.

## Minimum Public Coverage For One Flagship Native Surface

Each flagship native surface should have all of the following:

- one reference page that names the stable contract boundary
- one practical workflow page that shows how a serious user actually runs it
- one quality or evidence page that explains how to interpret trust,
  benchmark, or validation status
- one route from the root README or package README surfaces to those pages

The current native maximum-likelihood and native Bayesian surfaces now follow
that shape:

- reference: `docs/public/phylogenetics/interfaces/native-inference-and-benchmarks.md`
- workflows: the native maximum-likelihood and native Bayesian pages under
  `docs/public/phylogenetics/operations/`
- benchmark interpretation: `docs/public/phylogenetics/quality/native-benchmark-evidence.md`

## Entry Surfaces That Usually Need Updates

When one native surface becomes public-facing, check all of these:

- root `README.md`
- `packages/bijux-phylogenetics/README.md`
- `packages/phylogenetic/README.md` when the alias package should advertise
  the same runtime story
- `docs/index.md`
- `docs/public/phylogenetics/index.md`
- the matching public reference, operations, and quality pages

## Failure Mode To Prevent

The weak pattern is:

- one API exists
- one reference page names it
- no public workflow page shows how to use it
- no public quality page explains how to trust it

That state is better than silence, but it is still not release-grade
documentation for a flagship scientific surface.
