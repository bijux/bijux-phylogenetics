---
title: Product Surface and Ownership
audience: public
type: explanation
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-10
---

# Product Surface and Ownership

The public runtime surface is intentionally narrower than the full repository.

## Public Product Surface

- the canonical runtime package: `bijux-phylogenetics`
- the compatibility alias package: `phylogenetic`
- the documented CLI and Python entrypoints
- reviewable artifacts emitted from governed workflows

## Non-Public Ownership

- maintainer automation in `bijux-phylogenetics-dev`
- repository-only make, standards, and release wiring
- evidence-book generation and validation machinery

That split keeps the runtime understandable. Users should not need maintainer
tooling to understand what the runtime does.
