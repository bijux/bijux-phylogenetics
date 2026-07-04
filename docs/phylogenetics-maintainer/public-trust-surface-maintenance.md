---
title: Public Trust Surface Maintenance
audience: maintainer
type: handbook
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-03
---

# Public Trust Surface Maintenance

This repository now has multiple public trust surfaces that must move together:

- runtime architecture pages
- native workflow guides
- release review guides
- evidence-book navigation guides
- package metadata and package READMEs

If those surfaces drift, the repository starts looking broader or cleaner than
its governed evidence actually supports.

## Surfaces To Review Together

When one major public capability changes, review all of these together:

- `README.md`
- package `README.md` files
- `docs/phylogenetics/architecture/`
- `docs/phylogenetics/operations/`
- `docs/phylogenetics/quality/`
- `docs/phylogenetics-evidence-book/`
- package `pyproject.toml` metadata

## Questions To Ask

- does the architecture section still describe the real package and runtime
  ownership split
- do workflow guides still match the supported native and wrapper-backed
  surfaces
- does the release review guide still point to the real governed report and
  gate surfaces
- does the evidence-book overview avoid stale or contradictory summary claims
- do package metadata and package README links still reinforce the same public
  story

## Common Failure Pattern

The common failure is not one missing page. It is partial alignment:

- workflow pages get updated
- architecture pages stay stale
- evidence pages still summarize old status
- release pages still describe the previous review flow

That is the point where the repository starts looking polished in fragments
instead of coherent as one release surface.
