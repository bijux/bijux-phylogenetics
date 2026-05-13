---
title: Runtime Purpose and Scope
audience: public
type: explanation
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-10
---

# What This Runtime Is For

`bijux-phylogenetics` is a reproducible phylogenetics runtime for inspection,
comparison, annotation, comparative analysis, ancestral reconstruction,
adapter-driven engine workflows, and evidence-linked reporting.

It is not trying to replace every upstream inference engine. Its value comes
from making important phylogenetic operations reviewable, repeatable, and
packaged behind one governed Python surface.

## The Core Job

The runtime owns the operational layer between tracked inputs and reviewable
outputs:

- tree validation, rooting, comparison, and rendering
- alignment diagnostics, trimming, coding checks, and translation
- comparative trait analysis and ancestral-state reconstruction
- adapter workflows for external alignment and inference tools
- report, bundle, and artifact generation for downstream review

## What The Runtime Does Not Claim

- it does not claim to be a native maximum-likelihood or Bayesian inference
  engine
- it does not claim that the evidence-book already closes every runtime surface
- it does not flatten serious scientific questions into one generic success
  message

## Use This Section Next

- [repository scope and limits](repository-scope-and-limits.md) for the hard
  boundary
- [product surface and ownership](product-surface-and-ownership.md) for the
  runtime/package split
- [evidence boundary](evidence-boundary.md) for the difference between runtime
  capability and evidence closure
