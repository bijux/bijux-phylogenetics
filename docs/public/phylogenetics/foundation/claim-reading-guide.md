---
title: Claim Reading Guide
audience: public
type: explanation
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Claim Reading Guide

This repository publishes multiple kinds of claim, and they do not all mean the
same thing.

## Main Claim Types

- runtime claim
  The package owns and documents a surface that users can run.
- native ownership claim
  The repository implements the named surface locally rather than only wrapping
  an external engine.
- wrapper claim
  The repository intentionally exposes an external-engine workflow as a governed
  product surface without pretending the implementation is native.
- evidence claim
  The repository has bounded trust material showing how one surface behaved on
  governed studies, parity cases, or benchmark lanes.
- release claim
  The repository currently looks publishable across docs, metadata, readiness
  reports, and evidence boundaries.

## Reading Rule

Do not convert one claim type into another automatically.

- a documented runtime surface is not automatically evidence-closed
- a native implementation is not automatically broader than the documented
  public contract
- a wrapper-backed workflow is not automatically a weaker product surface, but
  it also is not a native claim
- one evidence study does not automatically close the whole runtime family
- one green release gate does not erase bounded or open evidence limits

## Practical Examples

- native maximum-likelihood result objects are public runtime and native
  ownership claims
- external FASTA-to-tree orchestration is a public workflow claim with wrapper
  ownership where that remains the honest boundary
- the evidence-book records where those surfaces have bounded study, parity, or
  benchmark support
- the release review workflow records whether the repository tells that story
  coherently at one release point

## Read This With

- [runtime and evidence boundary](evidence-boundary.md)
- [product surface and ownership](product-surface-and-ownership.md)
- [validation lanes](../quality/validation-lanes.md)
- [evidence-book overview](../../phylogenetics-evidence-book/index.md)
