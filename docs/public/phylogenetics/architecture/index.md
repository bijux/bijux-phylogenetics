---
title: Runtime Architecture
audience: public
type: explanation
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Runtime Architecture

The runtime architecture is designed to keep public behavior understandable and
package ownership explicit.

## Architectural Priorities

- one canonical runtime package
- one thin compatibility alias package
- explicit boundaries between runtime, maintainer tooling, and evidence
- public-facing capability families grouped by owned responsibility rather than
  by module history
- deterministic artifacts instead of ad hoc local outputs

## What This Section Answers

Use the architecture section when the question is:

- which package should I install or depend on
- where native runtime ownership ends and wrapper orchestration begins
- how major capability families are grouped inside the public runtime story
- why the evidence-book is linked to the repository but not treated as the
  runtime itself

## Read This Section Next

- [runtime system map](runtime-system-map.md) for the main capability families
- [package map](package-map.md) for package ownership
- [evidence consumer boundary](evidence-consumer-boundary.md) for the runtime
  versus evidence split
