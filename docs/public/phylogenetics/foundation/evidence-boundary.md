---
title: Runtime and Evidence Boundary
audience: public
type: explanation
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Runtime And Evidence Boundary

The runtime handbook and the evidence-book answer different questions.

## The Runtime Handbook Answers

- what the product surface is
- how to run it
- which commands and packages matter
- which native-versus-wrapper ownership boundary applies
- what operational limits apply

## The Evidence-Book Answers

- which Lund-derived study slices have governed parity evidence
- where trust is bounded, partial, or still open
- how R reference behavior and Bijux behavior compare on explicit cases

## What This Boundary Protects

- the public docs can name a native surface without implying that every model,
  benchmark lane, or study claim is already closed
- the public docs can describe wrapper-backed workflows honestly without
  downgrading the native runtime that now exists
- the evidence-book can stay narrower than the runtime without making the
  runtime docs false, as long as the narrower evidence boundary stays explicit

## Reading Rule

If you need a more explicit reading model for those different claim types, use
the [claim reading guide](claim-reading-guide.md). This page states the core
boundary; the reading guide explains how to interpret runtime, native, wrapper,
evidence, and release claims without collapsing them together.

The runtime can be broader than the currently closed evidence surface. That is
not a contradiction as long as the documentation says so plainly.
