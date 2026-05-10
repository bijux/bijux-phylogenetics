---
title: Runtime Architecture
audience: public
type: explanation
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-10
---

# Runtime Architecture

The runtime architecture is designed to keep public behavior understandable and
package ownership explicit.

## Architectural Priorities

- one canonical runtime package
- one thin compatibility alias package
- explicit boundaries between runtime, maintainer tooling, and evidence
- deterministic artifacts instead of ad hoc local outputs

## Read This Section Next

- [runtime system map](runtime-system-map.md) for the main capability families
- [package map](package-map.md) for package ownership
- [evidence consumer boundary](evidence-consumer-boundary.md) for the runtime
  versus evidence split
