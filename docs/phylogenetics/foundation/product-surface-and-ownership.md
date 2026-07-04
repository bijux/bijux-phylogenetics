---
title: Product Surface and Ownership
audience: public
type: explanation
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Product Surface and Ownership

The public runtime surface is intentionally narrower than the full repository.

## Public Product Surface

- the canonical runtime package: `bijux-phylogenetics`
- the compatibility alias package: `phylogenetic`
- the documented CLI and Python entrypoints
- the documented native likelihood, native inference, and benchmark contracts
- reviewable artifacts emitted from governed workflows

## Public Ownership Boundaries

- owned native tree, alignment, likelihood, comparative, and ancestral
  runtime surfaces are part of the product when they are named in the public
  docs and import surfaces
- owned native maximum-likelihood and Bayesian surfaces are public when the
  docs name their result contracts, limits, and evidence boundaries directly
- wrapper orchestration is also part of the serious public product surface when
  the repository intentionally exposes it as governed workflow execution rather
  than pretending it is an implementation detail

## What Public Ownership Does Not Mean

- it does not mean every internal helper or report row is a stable public API
- it does not mean every public workflow is native-first
- it does not mean the evidence-book has already closed every owned surface

## Non-Public Ownership

- maintainer automation in `bijux-phylogenetics-dev`
- repository-only make, standards, and release wiring
- evidence-book generation and validation machinery

That split keeps the runtime understandable. Users should not need maintainer
tooling to understand what the runtime does.
