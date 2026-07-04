---
title: Public Trust Surface Maintenance
audience: maintainer
type: handbook
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-04
---

# Public Trust Surface Maintenance

This repository has several public trust surfaces that must move together.

## Surfaces To Review Together

- runtime architecture pages
- native workflow guides
- release review guides
- evidence-book navigation guides
- package metadata and package READMEs

## Why This Matters

If those surfaces drift, the repository can start looking broader, cleaner, or
more mature than its governed evidence actually supports.

## Questions To Ask

- Does the architecture section still describe the real package and runtime
  ownership split?
- Do workflow guides still match the supported native and wrapper-backed
  surfaces?
- Does the release review guide still point to the real governed report and
  gate surfaces?
- Does the evidence-book overview avoid stale or contradictory summary claims?
- Do package metadata and package README links reinforce the same public story?

## Common Failure Pattern

The common failure is not one missing page. It is partial alignment:

- workflow pages get updated
- architecture pages stay stale
- evidence pages still summarize old status
- release pages still describe the previous review flow

That is when the repository starts looking polished in fragments instead of
coherent as one release surface.
